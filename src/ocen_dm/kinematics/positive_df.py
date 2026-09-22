"""Positive spherical action DFs and self-consistent stellar gravity in AGAMA.

Separate from the legacy QuasiSpherical inversion backend: density and beta are
outputs of the same DF, without clipping. Units: pc, Msun, km/s, pc km/s.
Static remnant/halo densities are not themselves certified DFs.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from functools import cached_property
import importlib

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from ..mass_models import TruncatedGNFW


class DFConvergenceError(ValueError):
    """The model has not met its numerical accuracy requirements."""


def agama_pc():
    """Initialize once; never silently change units of existing AGAMA objects."""
    try:
        agama = importlib.import_module("agama")
    except ImportError as exc:
        raise ImportError("The positive-DF branch requires optional package AGAMA") from exc
    units = agama.getUnits()
    if not units:
        agama.setUnits(length=0.001, mass=1., velocity=1.)
        units = agama.getUnits()
    expected = {"length": 0.001, "mass": 1., "velocity": 1.}
    unit_names = {"length": "kpc", "mass": "Msun", "velocity": "km/s"}
    actual = {k: units[k].to_value(unit_names[k]) if hasattr(units[k], "to_value")
              else units[k] for k in expected}
    if any(not np.isclose(actual[k], v, rtol=1e-12, atol=0) for k, v in expected.items()):
        raise RuntimeError("Incompatible AGAMA units; run the DF branch in a fresh process")
    return agama


def _finite(**values):
    if not all(np.isfinite(v) for v in values.values()):
        raise ValueError(f"parameters must be finite: {values}")


@dataclass(frozen=True)
class ActionComponent:
    """Positive DoublePowerLaw DF; fractions are shared mass/light fractions.

    h = h_r Jr + (3-h_r)L/2 and g = g_r Jr + (3-g_r)L/2. These coefficients
    control orbital weights and are NOT beta parameters. Negative slope_in is
    allowed, suppressing the DF at small actions.
    """
    fraction: float = 1.
    J0: float = 300.
    slope_in: float = 0.
    slope_out: float = 7.
    steepness: float = 2.
    h_r: float = 1.5
    g_r: float = 1.

    def __post_init__(self):
        _finite(**asdict(self))
        if self.fraction <= 0 or self.J0 <= 0 or self.steepness <= 0:
            raise ValueError("fraction, J0 and steepness must be positive")
        if self.slope_in >= 3 or self.slope_out <= 3:
            raise ValueError("finite DF mass requires slope_in < 3 and slope_out > 3")
        if not (0 < self.h_r < 3 and 0 < self.g_r < 3):
            raise ValueError("radial action coefficients must lie strictly between 0 and 3")

    def build(self, stellar_mass):
        return agama_pc().DistributionFunction(
            type="DoublePowerLaw", mass=stellar_mass*self.fraction,
            J0=self.J0, slopeIn=self.slope_in, slopeOut=self.slope_out,
            steepness=self.steepness, coefJrIn=self.h_r, coefJzIn=(3-self.h_r)/2,
            coefJrOut=self.g_r, coefJzOut=(3-self.g_r)/2, rotFrac=0.)


@dataclass(frozen=True)
class DFNumerics:
    r_min: float = 1e-4
    r_max: float = 1e5
    potential_nodes: int = 160
    velocity_nodes: int = 64
    moment_nodes: int = 240
    projection_nodes: int = 128
    max_iterations: int = 50
    iteration_tolerance: float = 5e-4
    closure_tolerance: float = 0.005
    mass_tolerance: float = 0.005

    def __post_init__(self):
        _finite(**asdict(self))
        if not 0 < self.r_min < self.r_max:
            raise ValueError("invalid radial domain")
        for key in ("potential_nodes", "velocity_nodes", "moment_nodes", "projection_nodes", "max_iterations"):
            value = getattr(self, key)
            if not isinstance(value, int) or value < (2 if key == "max_iterations" else 32):
                raise ValueError(f"invalid {key}")
        if min(self.iteration_tolerance, self.closure_tolerance, self.mass_tolerance) <= 0:
            raise ValueError("tolerances must be positive")


@dataclass(frozen=True)
class DFModelConfig:
    M_star: float = 3e6
    M_bh: float = 1e4
    M_rem: float = 0.
    a_rem: float = 5.
    M_dm_100: float = 0.
    r_s: float = 30.
    gamma: float = 0.
    r_t: float = 1000.
    distance_kpc: float = 5.43
    seed_radius: float = 5.
    components: tuple[ActionComponent, ...] = field(default_factory=lambda: (ActionComponent(),))
    numerics: DFNumerics = field(default_factory=DFNumerics)

    def __post_init__(self):
        _finite(**{k: v for k, v in asdict(self).items() if k not in ("components", "numerics")})
        if min(self.M_star, self.a_rem, self.r_s, self.r_t,
               self.distance_kpc, self.seed_radius) <= 0:
            raise ValueError("stellar mass, scales and distance must be positive")
        if min(self.M_bh, self.M_rem, self.M_dm_100) < 0 or not 0 <= self.gamma < 3:
            raise ValueError("invalid static mass component")
        if not self.components or not np.isclose(sum(c.fraction for c in self.components), 1.,
                                                 rtol=0, atol=1e-10):
            raise ValueError("positive stellar component fractions must sum to one")

    def to_dict(self):
        return {"schema_version": 1, **asdict(self)}

    @classmethod
    def from_dict(cls, payload):
        row = dict(payload)
        if row.pop("schema_version") != 1:
            raise ValueError("unsupported DF model schema")
        row["components"] = tuple(ActionComponent(**c) for c in row["components"])
        row["numerics"] = DFNumerics(**row["numerics"])
        return cls(**row)


def _xyz(r):
    return np.column_stack((r, np.zeros((len(r), 2))))


def spherical_velocity_moments(potential, df, af, radii, order=64):
    """Positive 2D quadrature of a spherical even DF, using AGAMA actions.

    Integrate 4*pi*v^2 dv dmu with mu=|v_r|/v in [0,1], v in [0,v_escape].
    Tangential pressure is one component, hence v_t^2/2. No energy inversion,
    Jeans integration, or negative-value repair occurs here.
    """
    agama_pc()
    radii = np.asarray(radii, float)
    x, w = np.polynomial.legendre.leggauss(order)
    v, mu, w = (x+1)/2, (x+1)/2, w/2
    weight = v[:, None]**2*w[:, None]*w[None, :]
    outputs = []
    for start in range(0, len(radii), 32):
        r = radii[start:start+32]
        escape = np.sqrt(-2*potential.potential(_xyz(r)))
        posvel = np.zeros((len(r), order, order, 6))
        posvel[:, :, :, 0] = r[:, None, None]
        speed = escape[:, None, None]*v[None, :, None]
        posvel[:, :, :, 3] = speed*mu[None, None, :]
        posvel[:, :, :, 4] = speed*np.sqrt(1-mu[None, None, :]**2)
        f = df(af(posvel.reshape(-1, 6))).reshape(len(r), order, order)
        if np.any(~np.isfinite(f)) or np.any(f < 0):
            raise DFConvergenceError("non-finite/negative phase-space DF values")
        weighted = f*weight[None, :, :]*4*np.pi*escape[:, None, None]**3
        rho = weighted.sum(axis=(1, 2))
        pr = (weighted*posvel[:, :, :, 3]**2).sum(axis=(1, 2))
        pt = (weighted*posvel[:, :, :, 4]**2/2).sum(axis=(1, 2))
        outputs.append(np.column_stack((rho, pr, pt)))
    return np.concatenate(outputs)


def _density_interpolator(agama, radii, density):
    """Log cubic spline with explicit power-law continuation at either end."""
    logr, logrho = np.log(radii), np.log(density)
    spline = CubicSpline(logr, logrho)
    inner, outer = float(spline(logr[0], 1)), float(spline(logr[-1], 1))
    if inner <= -3 or outer >= -3:
        raise DFConvergenceError("density grid does not reach finite-mass asymptotes")

    def callback(xyz):
        r = np.linalg.norm(xyz, axis=1)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            lr = np.log(r)
            result = spline(np.clip(lr, logr[0], logr[-1]))
            lo, hi = lr < logr[0], lr > logr[-1]
            result[lo] = logrho[0]+inner*(lr[lo]-logr[0])
            result[hi] = logrho[-1]+outer*(lr[hi]-logr[-1])
            if inner == 0:
                result[r == 0] = logrho[0]
            return np.exp(result)
    return agama.Density(density=callback, symmetry="s")


class PositiveDFModel:
    """Converged equilibrium compatible with the binned moment likelihood.

    AGAMA integrates the analytic DF over velocities; projection integrates the
    resulting intrinsic density/pressures, not a new Jeans solution. Independent
    direct AGAMA projection is also available for accuracy checks.
    """
    def __init__(self, config: DFModelConfig, *, initial_stellar_potential=None,
                 extra_static_potentials=()):
        self.config = config
        agama = agama_pc()
        n = config.numerics
        dfs = [c.build(config.M_star) for c in config.components]
        self.df = dfs[0] if len(dfs) == 1 else agama.DistributionFunction(*dfs)
        self.static_potentials = list(extra_static_potentials)
        if config.M_bh:
            self.static_potentials.append(agama.Potential(
                type="Plummer", mass=config.M_bh, scaleRadius=0.))  # exact point mass
        if config.M_rem:
            self.static_potentials.append(agama.Potential(
                type="Plummer", mass=config.M_rem, scaleRadius=config.a_rem))
        self.halo = None
        if config.M_dm_100:
            unit = TruncatedGNFW(1., config.r_s, config.gamma, config.r_t)
            self.halo = TruncatedGNFW(config.M_dm_100/unit.enclosed_mass([100.])[0],
                                      config.r_s, config.gamma, config.r_t)
            density = agama.Density(
                density=lambda xyz: self.halo.density(np.linalg.norm(xyz, axis=1)), symmetry="s")
            self.static_potentials.append(agama.Potential(
                type="Multipole", density=density, lmax=0, gridSizeR=n.potential_nodes,
                rmin=n.r_min, rmax=n.r_max))
        stellar_potential = (initial_stellar_potential if initial_stellar_potential is not None else
                             agama.Potential(type="Plummer", mass=config.M_star,
                                             scaleRadius=config.seed_radius))
        self.potential = agama.Potential(stellar_potential, *self.static_potentials)
        radii = np.geomspace(n.r_min, n.r_max, n.potential_nodes)
        points = _xyz(radii)
        previous = None
        self.iterations = []
        for i in range(n.max_iterations):
            af = agama.ActionFinder(self.potential)
            density = spherical_velocity_moments(
                self.potential, self.df, af, radii, n.velocity_nodes)[:, 0]
            self._positive(density, "DF density")
            self.stellar_density = _density_interpolator(agama, radii, density)
            stellar_potential = agama.Potential(
                type="Multipole", density=self.stellar_density, lmax=0,
                gridSizeR=n.potential_nodes, rmin=n.r_min, rmax=n.r_max)
            self.potential = agama.Potential(stellar_potential, *self.static_potentials)
            force = -self.potential.force(points)[:, 0]
            # A fixed component can mask an unconverged DF in total force alone.
            current = np.concatenate((force, density))
            self._positive(current, "self-consistency force/density")
            if previous is not None:
                change = float(np.max(np.abs(current/previous - 1)))
                self.iterations.append(change)
                if change < n.iteration_tolerance:
                    break
            previous = current
        else:
            raise DFConvergenceError(f"DF equilibrium did not converge in {n.max_iterations} iterations")
        self.af = agama.ActionFinder(self.potential)
        self.stellar_potential = stellar_potential
        self.galaxy = agama.GalaxyModel(self.potential, self.df, self.af)
        # Double the velocity quadrature and sample BETWEEN density grid nodes.
        check_r = np.sqrt(radii[:-1]*radii[1:])
        direct = spherical_velocity_moments(
            self.potential, self.df, self.af, check_r, 2*n.velocity_nodes)[:, 0]
        closure = float(np.max(np.abs(direct/self.stellar_density.density(_xyz(check_r))-1)))
        mass = float(stellar_potential.totalMass())
        mass_error = abs(mass/config.M_star - 1)
        self.diagnostics = dict(iterations=i+1, last_change=self.iterations[-1],
                                density_closure=closure, stellar_mass=mass,
                                mass_error=mass_error, df_mass=float(self.df.totalMass()))
        if not np.isfinite(closure) or closure > n.closure_tolerance:
            raise DFConvergenceError(f"DF-density closure error {closure:g}")
        if not np.isfinite(mass_error) or mass_error > n.mass_tolerance:
            raise DFConvergenceError(f"stellar mass closure error {mass_error:g}; enlarge/refine grid")

    @staticmethod
    def _positive(a, label):
        if np.any(~np.isfinite(a)) or np.any(np.asarray(a) <= 0):
            raise DFConvergenceError(f"non-positive/non-finite {label}; no clipping is applied")

    def _check_inside_tracer(self, r):
        r = np.asarray(r)
        n = self.config.numerics
        if np.any(~np.isfinite(r)) or np.any(r < n.r_min) or np.any(r > n.r_max/100):
            raise ValueError("requested radii need a wider DF integration grid")

    @cached_property
    def moment_table(self):
        agama_pc()
        n = self.config.numerics
        r = np.geomspace(n.r_min, n.r_max, n.moment_nodes)
        values = spherical_velocity_moments(self.potential, self.df, self.af, r, n.velocity_nodes)
        self._positive(values, "DF density/pressures")
        rho, radial, tangential = values.T
        return dict(r=r, rho=rho, radial_pressure=radial,
                    tangential_pressure=tangential, beta=1-tangential/radial)

    @cached_property
    def _splines(self):
        t = self.moment_table
        return {k: PchipInterpolator(np.log(t["r"]), np.log(t[k]), extrapolate=False)
                for k in ("rho", "radial_pressure", "tangential_pressure")}

    def intrinsic_moments(self, r):
        agama_pc()
        r = np.atleast_1d(np.asarray(r, float))
        self._check_inside_tracer(r)
        values = {k: np.exp(s(np.log(r))) for k, s in self._splines.items()}
        values["beta"] = 1-values["tangential_pressure"]/values["radial_pressure"]
        return values

    def projected_moments(self, R):
        agama_pc()
        R = np.atleast_1d(np.asarray(R, float))
        self._check_inside_tracer(R)
        n = self.config.numerics
        x, w = np.polynomial.legendre.leggauss(n.projection_nodes)
        u_max = np.arccosh(n.r_max/R)
        u = u_max[:, None]*(x[None, :]+1)/2
        r = R[:, None]*np.cosh(u)
        weight = r*u_max[:, None]*w[None, :]
        rho, pr, pt = (np.exp(self._splines[k](np.log(r)))
                       for k in ("rho", "radial_pressure", "tangential_pressure"))
        q = (R[:, None]/r)**2
        out = {"Sigma": np.sum(weight*rho, axis=1),
               "los": np.sum(weight*((1-q)*pr+q*pt), axis=1),
               "pmr": np.sum(weight*(q*pr+(1-q)*pt), axis=1),
               "pmt": np.sum(weight*pt, axis=1)}
        for k, value in out.items():
            self._positive(value, k)
        return out

    def direct_projected_moments(self, R):
        """Independent AGAMA LOS/velocity integration for validation."""
        agama_pc()
        R = np.atleast_1d(np.asarray(R, float))
        self._check_inside_tracer(R)
        rho, tensor = self.galaxy.moments(np.column_stack((R, np.zeros(len(R)))))
        return {"Sigma": rho, "los": rho*tensor[:, 2],
                "pmr": rho*tensor[:, 0], "pmt": rho*tensor[:, 1]}

    def phase_space_density(self, posvel):
        """f(x,v) in mass/(pc^3 (km/s)^3), zero on unbound phase space.

        The central point is excluded because an unsoftened BH is singular.
        Input is an N x 6 array in project units, within the numerical aperture.
        """
        agama_pc()
        points = np.asarray(posvel, float)
        if points.ndim != 2 or points.shape[1] != 6 or np.any(~np.isfinite(points)):
            raise ValueError("phase-space coordinates must be a finite N x 6 array")
        self._check_inside_tracer(np.linalg.norm(points[:, :3], axis=1))
        energy = self.potential.potential(points[:, :3])+np.sum(points[:, 3:]**2, axis=1)/2
        bound = energy < 0
        values = np.zeros(len(points))
        if np.any(bound):
            values[bound] = self.df(self.af(points[bound]))
        if np.any(~np.isfinite(values)) or np.any(values < 0):
            raise DFConvergenceError("invalid phase-space DF evaluation")
        return values
