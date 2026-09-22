"""Positive stellar mass/light mixtures and free-potential mock recovery.

Every stellar basis DF contributes independently to gravity and observed light.
The common potential is solved from the mass DF before projecting the light DF.
Dark remnants/halo retain the explicit static-density model of the production
branch: spatial positivity, not a certified DF for those two components.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import copy

import numpy as np
from scipy.special import softmax
from scipy.interpolate import PchipInterpolator

from .df_capacity import CapacityObservations
from .df_fit import FitCoordinate, config_at
from .positive_df import DFModelConfig, PositiveDFModel, agama_pc


class DensityTemplate:
    """Positive compact density shape with independently free total mass and size.

Templates supply the known *profile family*, not its true amplitude or radius.
This matched-family control avoids confusing mass-family error with tracer-DF
error. Subsequent blind tests must free or change these source-density shapes.
    """

    def __init__(self, radii, density, total_mass, edge, name):
        self.r = np.asarray(radii, float)
        self.rho = np.asarray(density, float)
        if (self.r.ndim != 1 or self.r.shape != self.rho.shape or len(self.r) < 20 or
                np.any(np.diff(self.r) <= 0) or np.any(self.r <= 0) or
                np.any(~np.isfinite(self.rho)) or np.any(self.rho <= 0) or
                not np.isfinite(total_mass+edge) or total_mass <= 0 or edge <= self.r[-1]):
            raise ValueError("invalid positive density template")
        self.total_mass, self.edge, self.name = float(total_mass), float(edge), name
        self.spline = PchipInterpolator(np.log(self.r), np.log(self.rho))

    @classmethod
    def from_mock(cls, mock, population):
        index = next(i for i,p in enumerate(mock.populations) if p.name == population)
        r = np.geomspace(1e-5, mock.r_t*(1-1e-5), 1200)
        return cls(r, mock.intrinsic_components(r)[:, index, 0], mock.masses[index], mock.r_t, population)

    def density(self, radii, mass, scale=1.):
        if mass < 0 or scale <= 0 or not np.isfinite(mass+scale):
            raise ValueError("invalid template mass/scale")
        r = np.asarray(radii, float)/scale
        if np.any(~np.isfinite(r)) or np.any(r < 0):
            raise ValueError("invalid density radii")
        result = np.exp(self.spline(np.log(np.clip(r, self.r[0], self.r[-1]))))
        near_edge = r > self.r[-1]
        # The specified lowered-isothermal density vanishes as Psi^(5/2).
        result = np.where(near_edge, self.rho[-1]*(np.maximum(self.edge-r, 0)/(self.edge-self.r[-1]))**2.5, result)
        return result*(mass/self.total_mass)/scale**3

    def potential(self, mass, scale=1., nodes=160):
        if mass == 0:
            return None
        agama = agama_pc()
        return agama.Potential(type="Multipole", density=agama.Density(
            density=lambda xyz:self.density(np.linalg.norm(xyz,axis=1),mass,scale),symmetry="s"),
            lmax=0,gridSizeR=nodes,rmin=1e-5,rmax=2*self.edge*scale)

    def to_dict(self):
        return dict(name=self.name,radii_pc=self.r.tolist(),density=self.rho.tolist(),
                    total_mass=self.total_mass,edge_pc=self.edge)

    @classmethod
    def from_dict(cls,row):
        return cls(row["radii_pc"],row["density"],row["total_mass"],row["edge_pc"],row["name"])


@dataclass(frozen=True)
class TemplateMassConfig:
    remnant_mass: float = 5e4
    remnant_scale: float = 1.
    halo_mass: float = 3e5
    halo_scale: float = 1.

    def __post_init__(self):
        if (not all(np.isfinite([self.remnant_mass,self.remnant_scale,self.halo_mass,self.halo_scale]))
                or min(self.remnant_mass,self.halo_mass)<0 or min(self.remnant_scale,self.halo_scale)<=0):
            raise ValueError("invalid template mass parameters")


@dataclass(frozen=True)
class MassLightConfig:
    mass: DFModelConfig
    light_fractions: tuple[float, ...]
    dark: TemplateMassConfig | None = None

    def __post_init__(self):
        w = np.asarray(self.light_fractions)
        if (len(w) != len(self.mass.components) or np.any(~np.isfinite(w)) or
                np.any(w <= 0) or not np.isclose(w.sum(), 1., rtol=0, atol=1e-10)):
            raise ValueError("independent positive light fractions must sum to one")
        if self.dark is not None and (self.mass.M_rem or self.mass.M_dm_100):
            raise ValueError("choose template dark components or legacy densities, not both")

    def to_dict(self):
        return dict(schema_version=1, mass=self.mass.to_dict(), light_fractions=list(self.light_fractions),
                    dark=None if self.dark is None else dict(vars(self.dark)))

    @classmethod
    def from_dict(cls, payload):
        if payload["schema_version"] != 1:
            raise ValueError("unsupported mass/light configuration")
        return cls(DFModelConfig.from_dict(payload["mass"]), tuple(payload["light_fractions"]),
                   None if payload.get("dark") is None else TemplateMassConfig(**payload["dark"]))


class MassLightDFModel(PositiveDFModel):
    """Poisson uses the mass mixture; photometry and kinematics use the light mixture.

All component DFs are normalized to the same unit-integrated shape. The two
sets of weights specify mass and light fractions, not central-density fractions.
Their ratio gives each basis component's relative M/L. Absolute luminosity is
profiled out as in the existing photometric interface.
    """

    def __init__(self, config: MassLightConfig, *, initial_stellar_potential=None, templates=None):
        self.mass_light_config = config
        self.templates = templates
        base_config, extra = config.mass, []
        if templates is not None:
            if config.dark is None:
                raise ValueError("density templates need explicit template mass parameters")
            for template,mass,scale in [(templates["remnant"],config.dark.remnant_mass,config.dark.remnant_scale),
                                         (templates["halo"],config.dark.halo_mass,config.dark.halo_scale)]:
                potential = template.potential(mass,scale,config.mass.numerics.potential_nodes)
                if potential is not None:
                    extra.append(potential)
        elif config.dark is not None:
            raise ValueError("template masses require saved density templates")
        super().__init__(base_config, initial_stellar_potential=initial_stellar_potential,
                         extra_static_potentials=extra)
        self.config = config.mass
        self.mass_df = self.df
        self._set_light(config)

    def _set_light(self,config):
        agama = agama_pc()
        dfs = [replace(c, fraction=w).build(config.mass.M_star)
               for c, w in zip(config.mass.components, config.light_fractions)]
        self.df = dfs[0] if len(dfs) == 1 else agama.DistributionFunction(*dfs)
        self.galaxy = agama.GalaxyModel(self.potential, self.df, self.af)
        self.diagnostics["light_fractions"] = list(config.light_fractions)
        self.diagnostics["mass_fractions"] = [c.fraction for c in config.mass.components]
        self.diagnostics["relative_mass_to_light"] = [c.fraction/w for c, w in
                                                       zip(config.mass.components, config.light_fractions)]

    def reweight_light(self,config):
        if config.mass != self.config or config.dark != self.mass_light_config.dark:
            raise ValueError("light reweighting cannot change gravitating parameters")
        model=copy.copy(self)
        model.mass_light_config=config
        model.diagnostics=copy.deepcopy(self.diagnostics)
        for key in ("moment_table","_splines"):
            model.__dict__.pop(key,None)
        model._set_light(config)
        return model

    def mass_profile(self, radii):
        r = np.asarray(radii, float)
        xyz = np.column_stack((r, r*0, r*0))
        return -self.potential.force(xyz)[:, 0]*r*r/agama_pc().G


def recovery_coordinates(config, *, shapes="full", with_dm=True, template_mode=False):
    """Explicit optimizer limits; these do not define inference priors."""
    coords = [FitCoordinate("M_star", 5e5, 7e6, True)]
    for i in range(1, len(config.mass.components)):
        coords.append(FitCoordinate(f"mixture_log_ratio.{i}", -6., 6.))
    if shapes not in ("weights", "scales", "full"):
        raise ValueError("unknown shape stage")
    for i in range(len(config.mass.components)):
        if shapes in ("scales", "full"):
            coords.extend([FitCoordinate(f"components.{i}.J0", 5., 4000., True),
                           FitCoordinate(f"components.{i}.g_r", .15, 2.85)])
        if shapes == "full":
            coords.extend([FitCoordinate(f"components.{i}.slope_in", -3., 2.5),
                           FitCoordinate(f"components.{i}.slope_out", 3.6, 30.),
                           FitCoordinate(f"components.{i}.steepness", .3, 10., True),
                           FitCoordinate(f"components.{i}.h_r", .15, 2.85)])
    # Linear amplitudes include the exact zero-remnant and zero-DM boundaries.
    if template_mode:
        coords.extend([FitCoordinate("dark.remnant_mass",0.,8e5),FitCoordinate("dark.remnant_scale",.3,3.,True)])
        if with_dm:
            coords.extend([FitCoordinate("dark.halo_mass",0.,3e6),FitCoordinate("dark.halo_scale",.3,3.,True)])
    else:
        coords.extend([FitCoordinate("M_rem",0.,8e5),FitCoordinate("a_rem",.3,12.,True)])
        if with_dm:
            coords.extend([FitCoordinate("M_dm_100",0.,3e6),FitCoordinate("r_s",3.,150.,True)])
    return coords


class RecoveryProblem:
    def __init__(self, observations: CapacityObservations, config: MassLightConfig,
                 *, shapes="full", with_dm=True, noise=None, templates=None, shared_ml=False):
        if not with_dm and (config.mass.M_dm_100 != 0 or config.dark is not None and config.dark.halo_mass != 0):
            raise ValueError("no-DM recovery requires zero halo mass")
        self.observations, self.initial = observations, config
        self.templates=templates
        self.shared_ml=bool(shared_ml)
        self.coordinates = recovery_coordinates(config, shapes=shapes, with_dm=with_dm,template_mode=templates is not None)
        self.nmass = len(self.coordinates)
        nlight=0 if self.shared_ml else len(config.light_fractions)-1
        self.names = [c.path for c in self.coordinates]+[f"light_log_ratio.{i}" for i in range(1,nlight+1)]
        self.lower = np.r_[[c.bounds[0] for c in self.coordinates], np.full(nlight,-6.)]
        self.upper = np.r_[[c.bounds[1] for c in self.coordinates], np.full(nlight,6.)]
        self.noise = np.zeros(len(observations.truth)) if noise is None else np.asarray(noise, float).copy()
        if self.noise.shape != observations.truth.shape or np.any(~np.isfinite(self.noise)):
            raise ValueError("noise must have one finite entry per mock datum")
        # Noise is in physical observable units (including magnitudes).
        self.target = observations.truth+self.noise
        self.last_model = None

    def encode(self, config=None):
        c = config or self.initial
        w = np.asarray(c.light_fractions)
        return np.r_[[p.encode(getattr(c.dark,p.path.split(".")[1])) if p.path.startswith("dark.")
                      else p.get(c.mass) for p in self.coordinates], [] if self.shared_ml else np.log(w[1:]/w[0])]

    def decode(self, x):
        x = np.asarray(x)
        if x.shape != self.lower.shape or np.any(x < self.lower) or np.any(x > self.upper):
            raise ValueError("recovery coordinates outside bounds")
        pairs=[(p,v) for p,v in zip(self.coordinates,x[:self.nmass]) if not p.path.startswith("dark.")]
        mass = config_at(self.initial.mass, [p for p,v in pairs], [v for p,v in pairs])
        dark_changes={p.path.split(".")[1]:p.decode(v) for p,v in zip(self.coordinates,x[:self.nmass]) if p.path.startswith("dark.")}
        dark=replace(self.initial.dark,**dark_changes) if self.initial.dark is not None else None
        light = tuple(c.fraction for c in mass.components) if self.shared_ml else tuple(softmax(np.r_[0., x[self.nmass:]]))
        return MassLightConfig(mass, light, dark)

    def evaluate(self, x, *, warm=False):
        config = self.decode(x)
        seed = self.last_model.stellar_potential if warm and self.last_model is not None else None
        model = (self.last_model.reweight_light(config) if self.last_model is not None and
                 self.last_model.config == config.mass and self.last_model.mass_light_config.dark == config.dark else MassLightDFModel(
                 config, initial_stellar_potential=seed,templates=self.templates))
        moments = model.projected_moments(self.observations.radii)
        array = np.column_stack([moments[k] for k in ("Sigma", "los", "pmr", "pmt")])
        predicted = self.observations.reduce(array)
        residual = predicted-self.target
        nkin = self.observations.nkin
        residual[nkin:] -= np.average(residual[nkin:], weights=self.observations.photo_weights)
        residual /= self.observations.errors
        self.last_model = model
        return dict(model=model, config=config, projected=array, prediction=predicted,
                    residual=residual, score=float(residual @ residual))


def gaussian_mock_noise(observations, seed):
    """Independent Gaussian noise under the explicitly adopted mock error model."""
    return np.random.default_rng(seed).normal(size=len(observations.truth))*observations.errors
