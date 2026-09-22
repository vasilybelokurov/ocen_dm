"""Joint photometric/kinematic fitting for the positive stellar DF branch.

Photometric uncertainties must be supplied explicitly. The legacy composite
profile provides relative weights, not calibrated errors. Its optional error
scale below defines a diagnostic objective, not a calibrated evidence model.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import copy

import numpy as np

from .likelihood import ARCSEC_PER_RAD, KinematicData, ProfileLikelihood
from .positive_df import DFModelConfig, PositiveDFModel


def model_config_from_dict(payload):
    """Dispatch explicit new families; absent family preserves legacy schema."""
    if payload.get("family") == "regularized_exponential":
        from .regularized_df import RegularizedDFConfig
        return RegularizedDFConfig.from_dict(payload)
    if payload.get("family") == "lowered_isothermal":
        from .lowered_isothermal import LoweredIsothermalConfig
        return LoweredIsothermalConfig.from_dict(payload)
    if "family" in payload:
        raise ValueError(f"unsupported DF family {payload['family']}")
    return DFModelConfig.from_dict(payload)


def build_df_model(config):
    if isinstance(config, DFModelConfig):
        return PositiveDFModel(config)
    from .regularized_df import RegularizedDFConfig, RegularizedDFModel
    if isinstance(config, RegularizedDFConfig):
        return RegularizedDFModel(config)
    from .lowered_isothermal import LoweredIsothermalConfig, LoweredIsothermalModel
    if isinstance(config, LoweredIsothermalConfig):
        return LoweredIsothermalModel(config)
    raise ValueError(f"unsupported DF configuration {type(config).__name__}")


@dataclass(frozen=True)
class PhotometricData:
    r_arcsec: np.ndarray
    mu: np.ndarray
    sigma_mag: np.ndarray
    source: str
    error_model: str

    def __post_init__(self):
        arrays = (self.r_arcsec, self.mu, self.sigma_mag)
        if any(np.asarray(a).ndim != 1 for a in arrays) or len(self.mu) < 2:
            raise ValueError("photometry needs at least two points in 1D arrays")
        if len({len(a) for a in arrays}) != 1 or any(np.any(~np.isfinite(a)) for a in arrays):
            raise ValueError("invalid photometric arrays")
        if np.any(self.r_arcsec <= 0) or np.any(self.sigma_mag <= 0):
            raise ValueError("photometric radii/errors must be positive")

    @classmethod
    def from_relative_weights(cls, profile, sigma_mag):
        if not np.isfinite(sigma_mag) or sigma_mag <= 0:
            raise ValueError("an explicit positive photometric error scale is required")
        keep = profile.weight > 0
        return cls(profile.r_arcsec[keep], profile.mu[keep],
                   sigma_mag/np.sqrt(profile.weight[keep]), profile.source,
                   f"ADOPTED diagnostic errors: {sigma_mag:g} mag / sqrt(relative weight); "
                   "not calibrated measurement errors; splice covariance omitted")

    def to_dict(self):
        return dict(schema_version=1, r_arcsec=self.r_arcsec.tolist(), mu=self.mu.tolist(),
                    sigma_mag=self.sigma_mag.tolist(), source=self.source, error_model=self.error_model)

    @classmethod
    def from_dict(cls, row):
        row = dict(row)
        if row.pop("schema_version") != 1:
            raise ValueError("unsupported photometric snapshot")
        return cls(**{k: np.array(v, float) if k in ("r_arcsec", "mu", "sigma_mag") else v
                      for k, v in row.items()})

    def compare(self, sigma):
        sigma = np.asarray(sigma)
        if sigma.shape != self.mu.shape or np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
            raise ValueError("invalid projected DF surface density")
        shape = -2.5*np.log10(sigma)
        # Profile out one common magnitude zero-point (constant stellar M/L).
        offset = np.average(self.mu-shape, weights=1/self.sigma_mag**2)
        predicted = shape+offset
        residual = (predicted-self.mu)/self.sigma_mag
        chi2 = float(residual @ residual)
        loglike = -0.5*chi2-float(np.log(np.sqrt(2*np.pi)*self.sigma_mag).sum())
        return dict(chi2=chi2, loglike=loglike, offset=float(offset), prediction=predicted)


class DFJointProblem:
    """One evaluation shares the converged model and its cached DF moments.

    Uses the saved selection nodes, units and split-normal errors of the Jeans
    analysis. Published streaming subtraction is an explicit moment-level
    approximation: the current DF is the nonrotating even part, not a fitted
    rotating DF. Infeasible streaming subtraction is rejected, never floored.
    """
    def __init__(self, data: KinematicData, photometry: PhotometricData):
        self.data, self.photometry = data, photometry
        self.likelihood = ProfileLikelihood(KinematicData(tuple(
            replace(p, streaming2=None) for p in data.profiles)))

    def evaluate(self, config: DFModelConfig, model=None):
        model = model or build_df_model(config)
        if model.config != config:
            raise ValueError("supplied model/config mismatch")
        pred = self.likelihood.predict(model, config.distance_kpc)
        terms = {}
        for p in self.data.profiles:
            if p.streaming2 is not None:
                variance = pred[p.name]**2-p.streaming2
                if np.any(variance <= 0):
                    raise ValueError(f"{p.name}: adopted streaming exceeds DF second moment")
                pred[p.name] = np.sqrt(variance)
            error = np.where(pred[p.name] > p.value, p.err_hi, p.err_lo)
            terms[p.name] = dict(
                chi2=float(np.sum(((pred[p.name]-p.value)/error)**2)), n=p.n,
                loglike=float(ProfileLikelihood._split_normal_lnlike(pred[p.name], p).sum()))
        R = self.photometry.r_arcsec*config.distance_kpc*1000/ARCSEC_PER_RAD
        photo = self.photometry.compare(model.projected_moments(R)["Sigma"])
        chi2_kin = sum(t["chi2"] for t in terms.values())
        lnlike_kin = sum(t["loglike"] for t in terms.values())
        return dict(model=model, predictions=pred, photometry=photo, terms=terms,
                    chi2_kinematic=chi2_kin, loglike_kinematic=lnlike_kin,
                    objective=chi2_kin+photo["chi2"],
                    loglike_profiled=lnlike_kin+photo["loglike"])


@dataclass(frozen=True)
class FitCoordinate:
    """Bounded optimizer coordinate (not an inference prior). Dotted config path."""
    path: str
    lower: float
    upper: float
    log: bool = False

    def __post_init__(self):
        if not np.isfinite(self.lower+self.upper) or self.lower >= self.upper:
            raise ValueError("invalid fit bounds")
        if self.log and self.lower <= 0:
            raise ValueError("log coordinate needs positive bounds")
        if self.path.startswith("numerics.") or self.path.endswith("fraction"):
            raise ValueError("numerics are fixed; use mixture_log_ratio.i to fit positive fractions")
        if self.path.startswith("mixture_log_ratio."):
            index = self.path.removeprefix("mixture_log_ratio.")
            if not index.isdigit() or int(index) < 1 or self.log:
                raise ValueError("mixture_log_ratio.i needs i >= 1 and log=False")

    def encode(self, value):
        return float(np.log(value) if self.log else value)

    def decode(self, value):
        return float(np.exp(value) if self.log else value)

    @property
    def bounds(self):
        return (self.encode(self.lower), self.encode(self.upper))

    def get(self, config):
        if self.path.startswith("mixture_log_ratio."):
            if not hasattr(config, "components"):
                raise ValueError("mixture coordinates require a mixture model")
            i = int(self.path.split(".")[1])
            if i >= len(config.components):
                raise ValueError("mixture index exceeds the number of components")
            node = np.log(config.components[i].fraction/config.components[0].fraction)
        else:
            node = config.to_dict()
            for key in self.path.split("."):
                node = node[int(key)] if isinstance(node, (tuple, list)) else node[key]
        if not self.lower <= node <= self.upper:
            raise ValueError(f"initial {self.path} lies outside fit bounds")
        return self.encode(node)


def config_at(config, coordinates, x):
    """Decode into a new validated config, leaving the starting model untouched."""
    if len(coordinates) != len(x) or len({c.path for c in coordinates}) != len(coordinates):
        raise ValueError("invalid coordinate vector")
    row = copy.deepcopy(config.to_dict())
    logits = np.log([c.fraction for c in getattr(config, "components", ())])
    if len(logits):
        logits -= logits[0]
    change_weights = False
    for c, value in zip(coordinates, x):
        if not c.bounds[0] <= value <= c.bounds[1]:
            raise ValueError(f"{c.path} outside bounds")
        if c.path.startswith("mixture_log_ratio."):
            i = int(c.path.split(".")[1])
            if i >= len(logits):
                raise ValueError("mixture index exceeds the number of components")
            logits[i] = value
            change_weights = True
            continue
        keys = c.path.split(".")
        node = row
        for key in keys[:-1]:
            node = node[int(key)] if isinstance(node, (tuple, list)) else node[key]
        if keys[-1] not in node:
            raise ValueError(f"unknown parameter {c.path}")
        node[keys[-1]] = c.decode(value)
    if change_weights:
        fractions = np.exp(logits-max(logits))
        fractions /= sum(fractions)
        for component, fraction in zip(row["components"], fractions):
            component["fraction"] = float(fraction)
    return type(config).from_dict(row)
