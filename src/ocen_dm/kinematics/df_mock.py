"""Independent, positive energy-based equilibrium mocks for DF capacity tests.

Each component has f(Q) = A expm1(Q/s^2) for Q>0, zero otherwise, with
Q = Phi(rt)-E-L^2/(2 ra^2).  All gravitating components enter a fresh spherical
Poisson solve. This is an Osipkov--Merritt lowered-isothermal construction,
not AGAMA's action DoublePowerLaw family and not a Jeans-generated mock.
It establishes collisionless equilibrium, not stability or tidal survival.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import hyp1f1

from .positive_df import agama_pc


def king_integrals(w):
    """Dimensionless density and one-component pressure, stable at the edge.

rho = A (2 pi s^2)^1.5 k0(W); pr = A (2 pi s^2)^1.5 s^2 k1(W).
    """
    w = np.maximum(np.asarray(w, float), 0.)
    return (8/(15*np.sqrt(np.pi))*w**2.5*hyp1f1(1., 3.5, w),
            16/(105*np.sqrt(np.pi))*w**3.5*hyp1f1(1., 4.5, w))


@dataclass(frozen=True)
class MockPopulation:
    name: str
    central_fraction: float
    velocity_ratio: float
    anisotropy_radius_pc: float | None = None
    light_per_mass: float = 1.

    def __post_init__(self):
        values = (self.central_fraction, self.velocity_ratio, self.light_per_mass)
        if not all(np.isfinite(values)) or min(values[:2]) <= 0 or values[2] < 0:
            raise ValueError("invalid mock population")
        if self.anisotropy_radius_pc is not None and (
                not np.isfinite(self.anisotropy_radius_pc) or self.anisotropy_radius_pc <= 0):
            raise ValueError("anisotropy radius must be positive or None (isotropic)")


class LoweredIsothermalMock:
    """Self-consistent mass components; luminosity can have different weights."""

    def __init__(self, populations, W0=7., radius_pc=3., velocity_kms=15., rtol=2e-10):
        self.populations = tuple(populations)
        if not self.populations or not np.isclose(sum(p.central_fraction for p in populations), 1):
            raise ValueError("central density fractions must sum to one")
        if min(W0, radius_pc, velocity_kms, rtol) <= 0 or not all(
                np.isfinite([W0, radius_pc, velocity_kms, rtol])):
            raise ValueError("invalid mock scale")
        if not any(p.light_per_mass > 0 for p in populations):
            raise ValueError("need at least one luminous population")
        self.W0, self.radius_pc, self.velocity_kms = W0, radius_pc, velocity_kms
        self.G = agama_pc().G
        self.rho0 = 9*velocity_kms**2/(4*np.pi*self.G*radius_pc**2)
        self.mass_scale = 4*np.pi*self.rho0*radius_pc**3
        self.fractions = np.array([p.central_fraction for p in populations])
        self.s2 = np.array([p.velocity_ratio**2 for p in populations])
        self.ra = np.array([np.inf if p.anisotropy_radius_pc is None else
                            p.anisotropy_radius_pc for p in populations])
        self.light = np.array([p.light_per_mass for p in populations])
        self.central_integral = king_integrals(W0/self.s2)[0]

        def rhs(x, y):
            a = 1+(x*radius_pc/self.ra)**2
            rho = self.fractions*king_integrals(y[0]/self.s2)[0]/self.central_integral/a
            return np.r_[-9*sum(y[1:])/x**2, x*x*rho]

        def edge(x, y):
            return y[0]
        edge.terminal = True
        edge.direction = -1
        x0 = 1e-6
        y0 = np.r_[W0-1.5*x0*x0, self.fractions*x0**3/3]
        self.solution = solve_ivp(rhs, (x0, 1e5), y0, events=edge, dense_output=True,
                                  rtol=rtol, atol=rtol*1e-3, max_step=2.)
        if not self.solution.success or not len(self.solution.t_events[0]):
            raise ValueError("mock Poisson solution did not reach a finite edge")
        self.r_t = float(self.solution.t_events[0][0]*radius_pc)
        self.masses = self.solution.y_events[0][0, 1:]*self.mass_scale
        self.total_mass = float(sum(self.masses))

    def _state(self, r):
        r = np.asarray(r, float)
        if np.any(~np.isfinite(r)) or np.any(r < 0):
            raise ValueError("mock radii must be finite and nonnegative")
        x = np.clip(r/self.radius_pc, 1e-6, self.r_t/self.radius_pc)
        state = self.solution.sol(x.ravel()).reshape((1+len(self.populations),)+r.shape)
        state[0] = np.where(r < self.r_t, state[0], 0.)
        tiny = r/self.radius_pc < 1e-6
        state[0] = np.where(tiny, self.W0-1.5*(r/self.radius_pc)**2, state[0])
        for i, fraction in enumerate(self.fractions):
            state[i+1] = np.where(tiny, fraction*(r/self.radius_pc)**3/3, state[i+1])
        return state

    def intrinsic_components(self, r):
        r = np.asarray(r, float)
        psi = self._state(r)[0]
        a = 1+(r[..., None]/self.ra)**2
        k0, k1 = king_integrals(psi[..., None]/self.s2)
        factor = self.rho0*self.fractions/self.central_integral/a
        rho = factor*k0
        pr = factor*k1*self.s2*self.velocity_kms**2
        return np.stack((rho, pr, pr/a), axis=-1)  # (..., population, rho/pr/pt)

    def enclosed_mass(self, r):
        return self._state(r)[1:].sum(axis=0)*self.mass_scale

    def potential_value(self, r):
        r = np.asarray(r, float)
        # Inside the truncation radius Phi = Phi(rt)-Psi; outside is Keplerian.
        return np.where(r <= self.r_t,
                        -self.G*self.total_mass/self.r_t-self._state(r)[0]*self.velocity_kms**2,
                        -self.G*self.total_mass/np.maximum(r, self.r_t))

    def phase_space_components(self, r, vr, vt):
        r, vr, vt = np.broadcast_arrays(r, vr, vt)
        psi = self._state(r)[0]*self.velocity_kms**2
        q = psi[..., None]-.5*vr[..., None]**2-.5*vt[..., None]**2*(1+(r[..., None]/self.ra)**2)
        s2 = self.s2*self.velocity_kms**2
        A = self.rho0*self.fractions/(self.central_integral*(2*np.pi*s2)**1.5)
        return A*np.expm1(np.maximum(q/s2, 0.))

    def projected_moments(self, R, order=160):
        R = np.atleast_1d(np.asarray(R, float))
        if np.any(~np.isfinite(R)) or np.any(R <= 0) or np.any(R >= self.r_t):
            raise ValueError("mock projection requires 0 < R < truncation radius")
        x, w = np.polynomial.legendre.leggauss(order)
        umax = np.arccosh(self.r_t/R)
        u = umax[:, None]*(x+1)/2
        r = R[:, None]*np.cosh(u)
        weight = r*umax[:, None]*w
        intrinsic = np.einsum("rqic,i->rqc", self.intrinsic_components(r), self.light)
        rho, pr, pt = np.moveaxis(intrinsic, -1, 0)
        q = (R[:, None]/r)**2
        return dict(Sigma=np.sum(weight*rho, axis=1),
                    los=np.sum(weight*((1-q)*pr+q*pt), axis=1),
                    pmr=np.sum(weight*(q*pr+(1-q)*pt), axis=1),
                    pmt=np.sum(weight*pt, axis=1))

    def _check_inside_tracer(self, r):
        if np.any(np.asarray(r) <= 0) or np.any(np.asarray(r) >= self.r_t):
            raise ValueError("observations extend beyond mock truncation radius")

    def agama_potential(self, nodes=400):
        agama = agama_pc()
        density = agama.Density(density=lambda xyz: self.intrinsic_components(
            np.linalg.norm(xyz, axis=1))[..., 0].sum(axis=-1), symmetry="s")
        return agama.Potential(type="Multipole", density=density, lmax=0,
                               gridSizeR=nodes, rmin=1e-5, rmax=2*self.r_t)

    def to_dict(self):
        return dict(family="Osipkov-Merritt lowered isothermal", W0=self.W0,
                    radius_pc=self.radius_pc, velocity_kms=self.velocity_kms,
                    populations=[asdict(p) for p in self.populations],
                    truncation_radius_pc=self.r_t, total_mass=self.total_mass,
                    component_masses=self.masses.tolist(), central_density=self.rho0,
                    rotating=False, stability_tested=False)


def challenge_mock(name):
    if name == "single_isotropic":
        return LoweredIsothermalMock([MockPopulation("stars", 1., 1.)])
    if name not in ("mixed_no_dm", "mixed_with_dm"):
        raise ValueError(f"unknown mock case {name}")
    dm = .015 if name == "mixed_with_dm" else 0.
    populations = [MockPopulation("compact_stars", .52-dm, .90, 35., 1.),
                   MockPopulation("extended_stars", .28, 1.12, 120., .7),
                   MockPopulation("dark_remnants", .20, .65, None, 0.)]
    if dm:
        populations.append(MockPopulation("extended_dark_component", dm, 1.5, None, 0.))
    return LoweredIsothermalMock(populations)
