"""Independent generalized lowered-isothermal (Michie/LIMEPY) control.

f(E,L)=A exp[-L^2/(2 ra^2 s^2)] E_gamma(g,(Et-E)/s^2), E<Et.
This is the Gieles & Zocchi family, not the energy-shift/OM mock in df_mock.
We solve stellar Poisson gravity in the same prescribed halo/remnant potential;
an isolated King solution with a halo added afterwards is not used.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma, hyp1f1

from .positive_df import DFConvergenceError, DFNumerics, _xyz, agama_pc
from .regularized_df import PrescribedMatter


def lowered_exponential(g, x):
    """Stable exp(x)*P(g,x), zero for x<=0; includes Woolley and King limits."""
    x = np.asarray(x, float)
    positive = x > 0
    out = np.zeros_like(x)
    if g == 0:
        out[positive] = np.exp(x[positive])
    elif g == 1:
        out[positive] = np.expm1(x[positive])
    else:
        xp = x[positive]
        out[positive] = xp**g/gamma(g+1)*hyp1f1(1., g+1, xp)
    return out


@lru_cache(maxsize=16)
def _quadrature(order):
    z, w = np.polynomial.legendre.leggauss(order)
    return (z+1)/2, w/2


def lowered_moments(W, r, s, ra, g, order=64):
    """rho, radial pressure, one-component tangential pressure for A=1."""
    W, r = np.broadcast_arrays(np.asarray(W, float), np.asarray(r, float))
    shape = W.shape
    W, r = W.ravel(), r.ravel()
    v, w = _quadrature(order)
    output = np.zeros((len(W), 3))
    for start in range(0, len(W), 32):
        stop = start+32
        ws = np.maximum(W[start:stop], 0.)
        speed = s*np.sqrt(2*ws[:, None, None])*v[None, :, None]
        vt2 = speed**2*(1-v[None, None, :]**2)
        vr2 = speed**2*v[None, None, :]**2
        energy = ws[:, None, None]*(1-v[None, :, None]**2)
        orbital = np.exp(-r[start:stop, None, None]**2*vt2/(2*ra**2*s*s))
        f = lowered_exponential(g, energy)*orbital
        weight = 4*np.pi*(s*np.sqrt(2*ws))[:, None, None]**3* (
            v[None, :, None]**2*w[None, :, None]*w[None, None, :])
        output[start:stop] = np.column_stack(((f*weight).sum(axis=(1, 2)),
                                             (f*weight*vr2).sum(axis=(1, 2)),
                                             (f*weight*vt2/2).sum(axis=(1, 2))))
    return output.reshape(shape+(3,))


@dataclass(frozen=True)
class LoweredIsothermalConfig:
    # Central conditions specify A and Et through the embedded Poisson solve.
    W0: float = 7.
    velocity_scale: float = 15.
    rho0: float = 4000.
    anisotropy_radius: float = 40.
    truncation_index: float = 1.
    matter: PrescribedMatter = field(default_factory=PrescribedMatter)
    numerics: DFNumerics = field(default_factory=DFNumerics)
    distance_kpc: float = 5.43
    poisson_rtol: float = 2e-9

    def __post_init__(self):
        values = [self.W0, self.velocity_scale, self.rho0, self.anisotropy_radius,
                  self.distance_kpc, self.poisson_rtol]
        if not np.all(np.isfinite(values)) or min(values) <= 0:
            raise ValueError("lowered-isothermal scales must be finite and positive")
        if not np.isfinite(self.truncation_index) or not 0 <= self.truncation_index <= 3.5:
            raise ValueError("supported truncation index is 0 <= g <= 3.5")
        if self.W0 > 50 or self.poisson_rtol >= .01:
            raise ValueError("central depth or Poisson tolerance outside numerical support")

    def to_dict(self):
        return dict(schema_version=1, family="lowered_isothermal", **asdict(self))

    @classmethod
    def from_dict(cls, payload):
        row = dict(payload)
        if row.pop("schema_version") != 1 or row.pop("family") != "lowered_isothermal":
            raise ValueError("unsupported lowered-isothermal schema")
        row["matter"] = PrescribedMatter(**row["matter"])
        row["numerics"] = DFNumerics(**row["numerics"])
        return cls(**row)


class LoweredIsothermalModel:
    """Spherical embedded Poisson solution with a finite stellar energy boundary.

    W0, rho0, s, ra and g are the five stellar controls; Mstar, A and Et are
    derived. For an explicitly isotropic comparison fix ra very large and
    check convergence; it is not then a fitted coordinate.
    """
    def __init__(self, config: LoweredIsothermalConfig):
        self.config = config
        ag, n = agama_pc(), config.numerics
        self.G = ag.G
        self.static_potentials = config.matter.potentials(n)
        self.static_potential = ag.Potential(*self.static_potentials) if self.static_potentials else None
        s, ra, g = config.velocity_scale, config.anisotropy_radius, config.truncation_index
        self.radius_scale = np.sqrt(9*s*s/(4*np.pi*self.G*config.rho0))
        self.mass_scale = 4*np.pi*config.rho0*self.radius_scale**3
        central = lowered_moments(config.W0, 0., s, ra, g, n.velocity_nodes)[0]
        self.amplitude = config.rho0/central

        def dark_mass(r):
            return (-self.static_potential.force([r, 0., 0.])[0]*r*r/self.G
                    if self.static_potential is not None else 0.)

        def rhs(x, y):
            r = x*self.radius_scale
            rho = self.amplitude*lowered_moments(y[0], r, s, ra, g, n.velocity_nodes)[0]
            return [-9*(y[1]+dark_mass(r)/self.mass_scale)/x**2,
                    x*x*rho/config.rho0]

        def edge(x, y):
            return y[0]
        edge.terminal = True
        edge.direction = -1
        x0 = min(1e-6, n.r_min/self.radius_scale)
        r0 = x0*self.radius_scale
        dark_delta = (self.static_potential.potential([r0, 0, 0])-
                      self.static_potential.potential([0., 0., 0.])
                      if self.static_potential is not None else 0.)
        y0 = [config.W0-1.5*x0*x0-dark_delta/s**2, x0**3/3]
        self.solution = solve_ivp(rhs, (x0, n.r_max/self.radius_scale), y0,
                                  events=edge, dense_output=True, rtol=config.poisson_rtol,
                                  atol=config.poisson_rtol*.01, max_step=2.)
        if not self.solution.success or not len(self.solution.t_events[0]):
            raise DFConvergenceError("lowered-isothermal Poisson solve did not reach a finite edge")
        self.r_t = float(self.solution.t_events[0][0]*self.radius_scale)
        self.M_star = float(self.solution.y_events[0][0, 1]*self.mass_scale)
        phi_dark = (self.static_potential.potential([self.r_t, 0, 0])
                    if self.static_potential is not None else 0.)
        self.energy_cutoff = -self.G*self.M_star/self.r_t+phi_dark
        density = ag.Density(density=lambda xyz: self.intrinsic_moments(
            np.linalg.norm(xyz, axis=1))["rho"], symmetry="s")
        self.stellar_potential = ag.Potential(type="Multipole", density=density, lmax=0,
                                             gridSizeR=n.potential_nodes, rmin=n.r_min,
                                             rmax=2*self.r_t)
        self.potential = ag.Potential(self.stellar_potential, *self.static_potentials)
        r = np.geomspace(max(n.r_min, .001*self.radius_scale), .98*self.r_t, 80)
        reference = self.energy_cutoff-self._W(r)*s*s
        phi = self.potential.potential(_xyz(r))
        force = -self.potential.force(_xyz(r))[:, 0]
        expected_mass = self.enclosed_stellar_mass(r)+np.array([dark_mass(x) for x in r])
        force_error = float(np.max(abs(force*r*r/self.G/expected_mass-1)))
        mass_error = abs(self.stellar_potential.totalMass()/self.M_star-1)
        self.diagnostics = dict(stellar_mass=self.M_star, mass_error=mass_error,
                                force_closure=force_error,
                                potential_closure=float(np.max(abs(phi-reference))/s**2),
                                truncation_radius_pc=self.r_t, energy_cutoff=self.energy_cutoff,
                                poisson_function_evaluations=int(self.solution.nfev),
                                family="generalized Michie/lowered-isothermal (LIMEPY family)",
                                density_closure=force_error)
        if max(force_error, self.diagnostics["potential_closure"]) > n.closure_tolerance or mass_error > n.mass_tolerance:
            raise DFConvergenceError(f"lowered-isothermal closure failed: {self.diagnostics}")

    def _W(self, r):
        r = np.asarray(r, float)
        x = np.clip(r/self.radius_scale, self.solution.t[0], self.solution.t[-1])
        W = self.solution.sol(x.ravel())[0].reshape(r.shape)
        return np.where(r < self.r_t, np.maximum(W, 0.), 0.)

    def enclosed_stellar_mass(self, r):
        r = np.asarray(r, float)
        x = np.clip(r/self.radius_scale, self.solution.t[0], self.solution.t[-1])
        mass = self.solution.sol(x.ravel())[1].reshape(r.shape)*self.mass_scale
        return np.where(r < self.solution.t[0]*self.radius_scale,
                        4*np.pi*self.config.rho0*r**3/3, mass)

    def intrinsic_moments(self, r, order=None):
        r = np.asarray(r, float)
        if np.any(~np.isfinite(r)) or np.any(r < 0):
            raise ValueError("radii must be finite and nonnegative")
        c = self.config
        values = self.amplitude*lowered_moments(self._W(r), r, c.velocity_scale,
                                                c.anisotropy_radius, c.truncation_index,
                                                order or c.numerics.velocity_nodes)
        rho, pr, pt = np.moveaxis(values, -1, 0)
        beta = 1-np.divide(pt, pr, out=np.ones_like(pr), where=pr > 0)
        return dict(rho=rho, radial_pressure=pr, tangential_pressure=pt, beta=beta)

    def _check_inside_tracer(self, r):
        r = np.asarray(r)
        if np.any(~np.isfinite(r)) or np.any(r <= 0) or np.any(r >= self.r_t):
            raise ValueError("projected radii must lie inside the stellar boundary")

    def projected_moments(self, R, *, order=None, velocity_order=None):
        R = np.atleast_1d(np.asarray(R, float))
        self._check_inside_tracer(R)
        c = self.config
        x, w = _quadrature(order or c.numerics.projection_nodes)
        umax = np.arccosh(self.r_t/R)
        r = R[:, None]*np.cosh(umax[:, None]*x)
        weight = 2*r*umax[:, None]*w
        moments = self.intrinsic_moments(r, velocity_order)
        rho, pr, pt = (moments[k] for k in ("rho", "radial_pressure", "tangential_pressure"))
        q = (R[:, None]/r)**2
        return dict(Sigma=np.sum(weight*rho, axis=1),
                    los=np.sum(weight*((1-q)*pr+q*pt), axis=1),
                    pmr=np.sum(weight*(q*pr+(1-q)*pt), axis=1), pmt=np.sum(weight*pt, axis=1))

    def direct_projected_moments(self, R):
        """Refined direct energy/velocity integration (not AGAMA projection)."""
        n = self.config.numerics
        return self.projected_moments(R, order=2*n.projection_nodes, velocity_order=2*n.velocity_nodes)

    def phase_space_density(self, posvel):
        points = np.asarray(posvel, float)
        if points.ndim != 2 or points.shape[1] != 6 or np.any(~np.isfinite(points)):
            raise ValueError("expected finite N x 6 phase-space coordinates")
        c = self.config
        energy = self.potential.potential(points[:, :3])+.5*np.sum(points[:, 3:]**2, axis=1)
        l2 = np.sum(np.cross(points[:, :3], points[:, 3:])**2, axis=1)
        return self.amplitude*np.exp(-l2/(2*c.anisotropy_radius**2*c.velocity_scale**2))* (
            lowered_exponential(c.truncation_index, (self.energy_cutoff-energy)/c.velocity_scale**2))
