"""Model families, priors and the nested-sampling driver for the Jeans fits.

Experiments follow the specification (section 6):

* **K1** -- no dark matter: stars (Trager MGE scaled by ``M_star``), a
  mass-segregated remnant Plummer sphere, a central point mass (``M_bh``, whose
  log-uniform prior reaches down to a floor that is dynamically zero), the
  three-parameter anisotropy and one multiplicative nuisance per Gaia instrument.
* **K2** -- K1 plus a truncated generalized-NFW halo, cored (``gamma = 0``) or
  cuspy (``gamma = 1``), parametrised by the halo mass within 100 pc and the
  scale radius so that the prior is flat in the quantity the data constrain.

Both are Bayesian models with the same likelihood; the sampler returns the
evidence so the pair can be compared. All parameters are transformed from the
unit hypercube, so the same families serve ultranest and dynesty.

Provenance: every run writes ``run.yaml`` with the git commit, the sha256 of each
processed table used, the prior table and the sampler settings.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
from astropy.table import Table
from scipy.special import erfinv

from ..light_model import MGEFit, build_stellar_mge, fit_mge_projected, load_tracer_profile
from ..mass_models import CompositeMassModel, PointMass, RemnantPlummer, TruncatedGNFW
from ..paths import processed_dir, results_dir
from .anisotropy import Anisotropy
from .jeans import SphericalJeans
from .likelihood import KinematicData, ProfileLikelihood

__all__ = ["Prior", "Parameter", "NoDarkMatterModel", "DarkMatterModel", "FitProblem",
           "run_nested", "maximum_likelihood", "OCEN_DISTANCE_KPC"]

OCEN_DISTANCE_KPC = 5.43          # Baumgardt & Vasiliev 2021; oMEGACat kinematic 5.494 differs by 1.2 per cent
OCEN_DISTANCE_ERR_KPC = 0.05
#: default distance prior: the literature value. Distance is degenerate with the
#: LOS/PM dispersion ratio, which is exactly what the MUSE and HST samples disagree
#: on at the 4 per cent level (JOURNAL 2026-09-17), so it is never fixed in a fit.
DEFAULT_DISTANCE_PRIOR = None     # set below once Prior is defined
#: smallest Gaussian allowed in the light model: the innermost Trager (1995) datum.
#: Unconstrained NNLS put 0.24 per cent of the light in 5-arcsec Gaussians inside
#: the data; such a compact tracer in a harmonic core is dynamically cold and the
#: fit compensated with a spurious 3e4 Msun point mass (JOURNAL 2026-09-17).
MGE_SIGMA_RANGE_ARCSEC = (10.5, 3000.0)
#: multiplicative nuisance per instrument relative to the HST PM sample. MUSE
#: targets bright giants (median F625W 17.4) whereas the HST PM sample is 3 mag
#: fainter (20.7); energy equipartition makes the giants ~3-4 per cent colder at
#: the same radius (measured on the oMEGACat catalogue, 2026-09-17). Gaia scales
#: absorb the crowding systematics seen in the overlap with HST.
NUISANCE_INSTRUMENTS = ("MUSE", "GaiaDR2", "GaiaEDR3")
NUISANCE_RANGE = {"MUSE": (0.85, 1.15), "GaiaDR2": (0.7, 1.3), "GaiaEDR3": (0.7, 1.3)}
PROFILE_GRID_PC = np.geomspace(0.1, 500.0, 60)


# ------------------------------------------------------------------ priors ---
@dataclass(frozen=True)
class Prior:
    """A one-dimensional prior with a unit-cube transform.

    ``kind``: ``uniform(lo, hi)``, ``loguniform(lo, hi)`` (uniform in log10 of the
    value; ``lo``/``hi`` are the values themselves), ``normal(mu, sigma)`` or
    ``truncnormal(mu, sigma, lo, hi)``.
    """

    kind: str
    lo: float = -np.inf
    hi: float = np.inf
    mu: float = 0.0
    sigma: float = 1.0

    def __post_init__(self) -> None:
        if self.kind not in ("uniform", "loguniform", "normal", "truncnormal"):
            raise ValueError(f"unknown prior kind {self.kind!r}")
        if self.kind in ("uniform", "loguniform", "truncnormal") and not self.lo < self.hi:
            raise ValueError("lo must be < hi")
        if self.kind == "loguniform" and self.lo <= 0:
            raise ValueError("loguniform needs lo > 0")
        if self.kind in ("normal", "truncnormal") and self.sigma <= 0:
            raise ValueError("sigma must be positive")

    def transform(self, u: np.ndarray) -> np.ndarray:
        """Map ``u`` in (0, 1) to the parameter value."""
        u = np.asarray(u, float)
        if self.kind == "uniform":
            return self.lo + (self.hi - self.lo) * u
        if self.kind == "loguniform":
            return 10 ** (np.log10(self.lo) + (np.log10(self.hi) - np.log10(self.lo)) * u)
        if self.kind == "normal":
            return self.mu + self.sigma * np.sqrt(2) * erfinv(2 * u - 1)
        # truncated normal: map u onto the CDF interval [F(lo), F(hi)]
        from scipy.special import erf
        F = lambda x: 0.5 * (1 + erf((x - self.mu) / (self.sigma * np.sqrt(2))))
        return self.mu + self.sigma * np.sqrt(2) * erfinv(2 * (F(self.lo) + (F(self.hi) - F(self.lo)) * u) - 1)

    def describe(self) -> str:
        if self.kind in ("uniform", "loguniform"):
            return f"{self.kind}({self.lo:g}, {self.hi:g})"
        if self.kind == "normal":
            return f"normal({self.mu:g}, {self.sigma:g})"
        return f"truncnormal({self.mu:g}, {self.sigma:g}; [{self.lo:g}, {self.hi:g}])"


DEFAULT_DISTANCE_PRIOR = Prior("normal", mu=OCEN_DISTANCE_KPC, sigma=OCEN_DISTANCE_ERR_KPC)


@dataclass(frozen=True)
class Parameter:
    name: str
    prior: Prior
    unit: str = ""
    latex: str = ""


# ---------------------------------------------------------------- families ---
class NoDarkMatterModel:
    """Experiment K1. See the module docstring for the physical content."""

    label = "K1_noDM"

    def __init__(self, mge_fit: MGEFit | None = None, instruments: Sequence[str] = NUISANCE_INSTRUMENTS,
                 distance_prior: Prior | None = "default", fix_distance: bool = False,
                 tracer: str = "trager", backend: str = "jeans", fixed: dict[str, float] | None = None,
                 constant_beta: bool = False, beta0_max: float = 0.0, distance_kpc: float | None = None) -> None:
        """
        Parameters
        ----------
        fixed : dict, optional
            Parameters held at given values and removed from the sampled vector
            (e.g. ``{"M_rem": 1e4, "M_bh": 1e2}`` to switch remnants and BH off).
        constant_beta : bool
            Use a single anisotropy ``beta_0`` (``beta_inf = beta_0``); literature
            comparisons with constant-anisotropy models need this.
        beta0_max : float
            Upper edge of the ``beta_0`` prior (default 0, the An & Evans bound for a
            cored tracer; literature presets that fitted radial anisotropy raise it).
        distance_kpc : float, optional
            Fixed distance to use instead of :data:`OCEN_DISTANCE_KPC` when
            ``fix_distance`` is true.
        """
        if backend not in ("jeans", "jam", "agama"):
            raise ValueError("backend must be 'jeans', 'jam' or 'agama'")
        self.backend = backend
        self.tracer = tracer
        self.fixed = dict(fixed or {})
        self.constant_beta = constant_beta
        self.beta0_max = float(beta0_max)
        self.fixed_distance_kpc = float(distance_kpc) if distance_kpc is not None else OCEN_DISTANCE_KPC
        if mge_fit is None:
            profile = load_tracer_profile(tracer)
            # smallest Gaussian = innermost datum of whichever profile is used
            sigma_range = (max(MGE_SIGMA_RANGE_ARCSEC[0] if tracer == "trager" else float(profile.r_arcsec.min()), 0.5),
                           MGE_SIGMA_RANGE_ARCSEC[1])
            mge_fit = fit_mge_projected(profile, sigma_range_arcsec=sigma_range)
        self.mge_fit = mge_fit
        if tracer != "trager" and not hasattr(self, "gamma"):
            self.label = f"{self.label}_{tracer}"
        if backend != "jeans" and not hasattr(self, "gamma"):
            self.label = f"{self.label}_{backend}"
        self.nuisance_instruments = tuple(instruments)
        if fix_distance:
            distance_prior = None                      # fixed at OCEN_DISTANCE_KPC (tests, comparisons)
        elif distance_prior == "default":
            distance_prior = DEFAULT_DISTANCE_PRIOR
        self.distance_prior = distance_prior
        self.parameters: tuple[Parameter, ...] = self._parameters()
        self.names = tuple(p.name for p in self.parameters)

    def _physical_parameters(self) -> list[Parameter]:
        params = [
            Parameter("M_star", Prior("loguniform", 1e6, 1e7), "Msun", r"M_\star"),
            Parameter("M_rem", Prior("loguniform", 1e4, 3e6), "Msun", r"M_{\rm rem}"),
            Parameter("a_rem", Prior("loguniform", 0.3, 20.0), "pc", r"a_{\rm rem}"),
            Parameter("M_bh", Prior("loguniform", 1e2, 3e5), "Msun", r"M_\bullet"),
        ]
        b0max = getattr(self, "beta0_max", 0.0)
        if getattr(self, "backend", "jeans") == "agama":
            # Cuddeford-Osipkov-Merritt family of the positive DF: beta -> 1 beyond r_a
            params += [Parameter("beta_0", Prior("uniform", -1.0, b0max), "", r"\beta_0"),
                       Parameter("r_a", Prior("loguniform", 1.0, 1000.0), "pc", r"r_a")]
        elif getattr(self, "constant_beta", False):
            params += [Parameter("beta_0", Prior("uniform", -1.0, b0max), "", r"\beta")]
        else:
            params += [Parameter("beta_0", Prior("uniform", -1.0, b0max), "", r"\beta_0"),      # An & Evans: cored tracer => beta_0 <= 0
                       Parameter("beta_inf", Prior("uniform", -1.0, 1.0), "", r"\beta_\infty"),
                       Parameter("r_beta", Prior("loguniform", 0.5, 100.0), "pc", r"r_\beta")]
        return [p for p in params if p.name not in getattr(self, "fixed", {})]

    def _parameters(self) -> tuple[Parameter, ...]:
        params = self._physical_parameters()
        for inst in self.nuisance_instruments:
            lo, hi = NUISANCE_RANGE.get(inst, (0.7, 1.3))
            params.append(Parameter(f"s_{inst}", Prior("uniform", lo, hi), "", rf"s_{{\rm {inst}}}"))
        if self.distance_prior is not None:
            params.append(Parameter("distance", self.distance_prior, "kpc", "D"))
        return tuple(params)

    # -- construction ---------------------------------------------------------
    def distance(self, theta: dict[str, float]) -> float:
        return float(theta.get("distance", self.fixed_distance_kpc))

    def complete(self, theta: dict[str, float]) -> dict[str, float]:
        """Sampled parameters plus the fixed ones and the constant-beta aliases."""
        full = {**self.fixed, **theta}
        if self.constant_beta and "beta_inf" not in full:
            full["beta_inf"] = full["beta_0"]
            full.setdefault("r_beta", 1.0)
        return full

    def scales(self, theta: dict[str, float]) -> dict[str, float]:
        return {inst: float(theta[f"s_{inst}"]) for inst in self.nuisance_instruments}

    def extra_components(self, theta: dict[str, float], distance_kpc: float) -> list:
        return []

    def build(self, theta: dict[str, float]) -> tuple[Any, float, dict[str, float]]:
        """Return ``(model, distance_kpc, instrument_scales)`` for one parameter vector.

        ``model`` is a :class:`SphericalJeans` (backend ``'jeans'``), a
        :class:`~ocen_dm.kinematics.backends.JamBackend` (``'jam'``) or an
        :class:`~ocen_dm.kinematics.backends.AgamaDFBackend` (``'agama'``); all
        expose ``projected_moments`` and ``mass``.
        """
        theta = self.complete(theta)
        D = self.distance(theta)
        stars = build_stellar_mge(self.mge_fit, D, total_mass=theta["M_star"])
        comps = [stars, RemnantPlummer(theta["M_rem"], theta["a_rem"]), PointMass(theta["M_bh"])]
        comps += self.extra_components(theta, D)
        mass = CompositeMassModel(comps)
        if self.backend == "agama":
            from .backends import AgamaDFBackend
            return AgamaDFBackend(mass, stars, beta0=theta["beta_0"], r_a=theta["r_a"], distance_kpc=D), D, self.scales(theta)
        anis = Anisotropy(theta["beta_0"], theta["beta_inf"], theta["r_beta"])
        if self.backend == "jam":
            from .backends import JamBackend
            return JamBackend(mass, stars, anis, D), D, self.scales(theta)
        return SphericalJeans(mass, stars, anis), D, self.scales(theta)

    def transform(self, u: np.ndarray) -> np.ndarray:
        u = np.asarray(u, float)
        return np.array([p.prior.transform(u[..., i]) for i, p in enumerate(self.parameters)]).T

    def to_dict(self, x: np.ndarray) -> dict[str, float]:
        return dict(zip(self.names, np.asarray(x, float)))


class DarkMatterModel(NoDarkMatterModel):
    """Experiment K2: K1 plus a truncated gNFW halo.

    Parameters
    ----------
    gamma : float
        Inner slope: 0 (cored) or 1 (NFW).
    r_t : float
        Truncation radius in pc, fixed (the Jacobi radius is of order a few
        hundred pc; the data end well inside it).
    """

    def __init__(self, gamma: float = 0.0, r_t: float = 1000.0, **kwargs: Any) -> None:
        self.gamma = float(gamma)
        self.r_t = float(r_t)
        super().__init__(**kwargs)
        self.label = "K2_cored" if gamma == 0 else ("K2_nfw" if gamma == 1 else f"K2_gnfw{gamma:g}")
        if self.tracer != "trager":
            self.label += f"_{self.tracer}"
        if self.backend != "jeans":
            self.label += f"_{self.backend}"

    def _physical_parameters(self) -> list[Parameter]:
        return super()._physical_parameters() + [
            Parameter("M_dm_100", Prior("loguniform", 1e3, 3e7), "Msun", r"M_{\rm DM}(<100\,{\rm pc})"),
            Parameter("r_s", Prior("loguniform", 3.0, 1000.0), "pc", r"r_s"),
        ]

    def halo(self, theta: dict[str, float]) -> TruncatedGNFW:
        unit = TruncatedGNFW(rho_s=1.0, r_s=theta["r_s"], gamma=self.gamma, r_t=self.r_t)
        rho_s = theta["M_dm_100"] / float(np.asarray(unit.enclosed_mass(100.0)).reshape(-1)[0])
        return TruncatedGNFW(rho_s=rho_s, r_s=theta["r_s"], gamma=self.gamma, r_t=self.r_t,
                             name="dm_cored" if self.gamma == 0 else "dm_nfw")

    def extra_components(self, theta: dict[str, float], distance_kpc: float) -> list:
        return [self.halo(theta)]


# ---------------------------------------------------------------- problem ---
class FitProblem:
    """A model family bound to data: the callable the samplers need."""

    def __init__(self, family: NoDarkMatterModel, data: KinematicData) -> None:
        self.family = family
        self.data = data
        self.likelihood = ProfileLikelihood(data)
        self.n_calls = 0

    def loglike_vector(self, x: np.ndarray) -> float:
        self.n_calls += 1
        theta = self.family.to_dict(x)
        try:
            jeans, D, scales = self.family.build(theta)
            ln = self.likelihood.lnlike(jeans, D, scales)
        except (ValueError, FloatingPointError, ZeroDivisionError):
            return -1e100
        return ln if np.isfinite(ln) else -1e100

    def loglike_unit(self, u: np.ndarray) -> float:
        return self.loglike_vector(self.family.transform(u))

    def predict(self, x: np.ndarray) -> dict[str, np.ndarray]:
        jeans, D, scales = self.family.build(self.family.to_dict(x))
        return self.likelihood.predict(jeans, D, scales)

    def chi2(self, x: np.ndarray) -> dict[str, tuple[float, int]]:
        jeans, D, scales = self.family.build(self.family.to_dict(x))
        return self.likelihood.chi2_terms(jeans, D, scales)

    def radial_profile(self, x: np.ndarray, grid: np.ndarray = PROFILE_GRID_PC) -> dict[str, np.ndarray]:
        model, _, _ = self.family.build(self.family.to_dict(x))
        return model.mass.radial_profile(grid)

    def mock_data(self, x: np.ndarray, rng: np.random.Generator) -> KinematicData:
        """Data with the real bins and errors but values drawn from model ``x`` (injection tests)."""
        pred = self.predict(x)
        out = []
        for p in self.data.profiles:
            err = 0.5 * (p.err_lo + p.err_hi)
            value = pred[p.name] + rng.normal(0.0, err)
            out.append(BinnedProfile_replace(p, value=np.maximum(value, 0.05 * pred[p.name])))
        return KinematicData(tuple(out))


def BinnedProfile_replace(p, **changes):
    from dataclasses import replace
    return replace(p, **changes)


# -------------------------------------------------------------- optimiser ---
def maximum_likelihood(problem: FitProblem, n_starts: int = 8, seed: int = 0,
                       maxiter: int = 400) -> tuple[np.ndarray, float]:
    """Multi-start Nelder-Mead in the unit cube; returns ``(x_best, lnL_best)``.

    Used for smoke tests and injection recovery where the full posterior is not
    needed. Not a substitute for the nested-sampling run.
    """
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    n = len(problem.family.parameters)
    best_u, best = None, -np.inf
    starts = np.clip(rng.uniform(0.15, 0.85, size=(n_starts, n)), 1e-6, 1 - 1e-6)
    for u0 in starts:
        f = lambda u: -problem.loglike_unit(np.clip(u, 1e-9, 1 - 1e-9))
        res = minimize(f, u0, method="Nelder-Mead", options={"maxiter": maxiter, "xatol": 1e-4, "fatol": 1e-3})
        if -res.fun > best:
            best, best_u = -res.fun, np.clip(res.x, 1e-9, 1 - 1e-9)
    return problem.family.transform(best_u), float(best)


# ------------------------------------------------------------------ driver ---
def _git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True,
                              cwd=str(results_dir().parent)).stdout.strip()
    except Exception:                                  # noqa: BLE001 - provenance best effort
        return "unknown"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_nested(problem: FitProblem, out_dir: Path, *, n_live: int = 400, dlogz: float = 0.5,
               max_ncalls: int | None = None, seed: int = 42, n_profile_samples: int = 300,
               resume: str = "overwrite", verbose: bool = False, step_sampler: bool = False,
               data_provenance: dict[str, Any] | None = None) -> dict[str, Any]:
    """Run ultranest on ``problem`` and write posterior, summary and profiles to ``out_dir``.

    Files written: ``posterior.ecsv`` (equally weighted samples), ``summary.json``
    (evidence, quantiles, chi2 per dataset at the maximum-likelihood sample),
    ``profiles.npz`` (``radial_profile`` on :data:`PROFILE_GRID_PC` for a random
    subset of samples), ``run.yaml`` (provenance). ultranest's own output goes to
    ``out_dir/ultranest``.
    """
    import ultranest
    import yaml

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fam = problem.family
    t0 = time.time()
    np.random.seed(seed)
    sampler = ultranest.ReactiveNestedSampler(list(fam.names), problem.loglike_vector, fam.transform,
                                              log_dir=str(out_dir / "ultranest"), resume=resume)
    if step_sampler:
        # MLFriends rejection sampling needed ~150 calls per iteration on the 11-13-D
        # K1/K2 posteriors (2026-09-17); a slice sampler costs ~2 ndim steps per
        # iteration whatever the region shape.
        import ultranest.stepsampler as uss
        sampler.stepsampler = uss.SliceSampler(
            nsteps=2 * len(fam.names), generate_direction=uss.generate_mixture_random_direction)
    result = sampler.run(min_num_live_points=n_live, dlogz=dlogz, max_ncalls=max_ncalls,
                         viz_callback=None, show_status=verbose)
    elapsed = time.time() - t0

    samples = np.asarray(result["samples"])
    logl = np.asarray(result["weighted_samples"]["logl"])
    wsamp = np.asarray(result["weighted_samples"]["points"])
    x_ml = wsamp[np.argmax(logl)]

    post = Table(samples, names=list(fam.names))
    post.meta["family"] = fam.label
    post.write(out_dir / "posterior.ecsv", format="ascii.ecsv", overwrite=True)
    np.save(out_dir / "ml_x.npy", x_ml)                 # for --mock-from injection runs

    q = np.percentile(samples, [16, 50, 84], axis=0)
    chi2 = problem.chi2(x_ml)
    summary = {
        "family": fam.label, "logz": float(result["logz"]), "logzerr": float(result["logzerr"]),
        "n_calls": int(result["ncall"]), "n_samples": int(len(samples)), "elapsed_s": round(elapsed, 1),
        "parameters": {n: {"p16": float(q[0, i]), "p50": float(q[1, i]), "p84": float(q[2, i]),
                           "ml": float(x_ml[i]), "unit": fam.parameters[i].unit,
                           "prior": fam.parameters[i].prior.describe()} for i, n in enumerate(fam.names)},
        "lnL_max": float(logl.max()),
        "chi2_ml": {k: {"chi2": round(c, 2), "n": n} for k, (c, n) in chi2.items()},
        "chi2_ml_total": round(sum(c for c, _ in chi2.values()), 2),
        "n_points": problem.data.n_points,
        "datasets": {p.name: {"kind": p.kind, "n": p.n, "instrument": p.instrument, "note": p.note}
                     for p in problem.data.profiles},
        # what was actually fitted: real profiles, or a mock realisation. Without this a
        # report cannot know which data to draw the model against (JOURNAL 2026-09-18).
        "data": data_provenance or {"kind": "real"},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    rng = np.random.default_rng(seed)
    idx = rng.choice(len(samples), size=min(n_profile_samples, len(samples)), replace=False)
    profs = [problem.radial_profile(samples[i]) for i in idx]
    stack = {k: np.array([p[k] for p in profs]) for k in profs[0] if k != "r"}
    np.savez_compressed(out_dir / "profiles.npz", r=PROFILE_GRID_PC, sample_index=idx, **stack)

    kin = processed_dir() / "kinematics"
    run = {
        "family": fam.label, "git_commit": _git_commit(), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t0)),
        "elapsed_s": round(elapsed, 1), "sampler": {"name": "ultranest", "version": ultranest.__version__,
                                                    "n_live": n_live, "dlogz": dlogz, "seed": seed, "max_ncalls": max_ncalls},
        "data": data_provenance or {"kind": "real"},
        "priors": {p.name: p.prior.describe() for p in fam.parameters},
        "distance_kpc": None if fam.distance_prior else OCEN_DISTANCE_KPC,
        "inputs": {f.name: _sha256(f) for f in sorted(kin.glob("*.ecsv"))},
    }
    (out_dir / "run.yaml").write_text(yaml.safe_dump(run, sort_keys=False))
    return summary
