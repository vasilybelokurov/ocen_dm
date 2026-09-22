"""Regularized spherical exponential action DFs (pc, Msun, km/s).

The contour prescription follows Binney (2026, MNRAS 549, stag854).
The *exact spherical* frequency ratio is the default.  Every equilibrium
iteration rebuilds the potential-dependent contour map and normalization.
This module does not alter the older DoublePowerLaw models or their schema.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from functools import cached_property

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.special import gammainccinv, gammaln

from .positive_df import (
    DFConvergenceError, DFNumerics, PositiveDFModel, _density_interpolator,
    _finite, _xyz, agama_pc, spherical_velocity_moments,
)


@dataclass(frozen=True)
class ExponentialParameters:
    J0: float = 300.
    alpha: float = 1.5
    b_out: float = .5
    J_a: float = 300.
    # Optional second transition: b_out is then the intermediate amplitude.
    b_outer: float | None = None
    J_outer: float | None = None

    def __post_init__(self):
        _finite(**{k: v for k, v in asdict(self).items() if v is not None})
        if min(self.J0, self.J_a) <= 0 or self.alpha <= .5:
            raise ValueError("J0 and J_a must be positive; alpha > 1/2 ensures a smooth joint origin")
        if (self.b_outer is None) != (self.J_outer is None):
            raise ValueError("second transition needs both b_outer and J_outer")
        if self.J_outer is not None and self.J_outer <= self.J_a:
            raise ValueError("second transition requires J_outer > J_a")
        # Computational guard, not an inference prior; prevents exponential overflow.
        if max(abs(self.b_out), abs(self.b_outer or 0.)) > 8:
            raise ValueError("anisotropy amplitude exceeds the supported numerical range")

    def bias(self, q):
        q = np.asarray(q, float)
        s = lambda a: (q/np.hypot(q, a))**2
        result = self.b_out*s(self.J_a)
        if self.b_outer is not None:
            result += (self.b_outer-self.b_out)*s(self.J_outer)
        return result


@dataclass(frozen=True)
class ContourNumerics:
    action_nodes: int = 192
    circularity_nodes: int = 65
    normalization_order: int = 160
    ode_rtol: float = 2e-9
    normalization_tolerance: float = 2e-5
    frequency_mode: str = "exact"

    def __post_init__(self):
        for k in ("action_nodes", "circularity_nodes", "normalization_order"):
            if not isinstance(getattr(self, k), int) or getattr(self, k) < 24:
                raise ValueError(f"invalid {k}")
        if not 0 < self.ode_rtol < .01 or not 0 < self.normalization_tolerance < .01:
            raise ValueError("invalid contour tolerances")
        if self.frequency_mode not in ("exact", "epicycle"):
            raise ValueError("frequency_mode must be exact or epicycle")


class SphericalFrequencyRatio:
    """Omega_r/Omega_t in a spherical, finite-central-potential model.

    ActionMapper returns spherical frequencies without numerical differencing
    of the Hamiltonian.  The optional epicycle prescription is an explicitly
    different, approximate reference, evaluated pointwise along every contour.
    """
    def __init__(self, potential, mode="exact"):
        agama = agama_pc()
        if not np.isfinite(potential.potential([0., 0., 0.])):
            raise ValueError("regularized finite-centre DF does not support a point mass")
        if mode not in ("exact", "epicycle"):
            raise ValueError("unknown frequency reference")
        self.potential, self.mode, self.agama = potential, mode, agama
        self.mapper = agama.ActionMapper(potential)

    def __call__(self, jr, angular):
        jr, angular = np.broadcast_arrays(jr, angular)
        shape = jr.shape
        jr, angular = jr.ravel(), angular.ravel()
        if self.mode == "exact":
            ratio = np.full(len(jr), 2.)
            radial = angular == 0
            # Direct period integrals lose precision at a coalescing pair of
            # turning points. The epicycle limit has O(Jr/L) error here.
            circular = (~radial) & (jr < 1e-5*angular)
            general = ~(radial | circular)
            if np.any(circular):
                r = np.atleast_1d(self.potential.Rcirc(L=angular[circular]))
                force, deriv = self.potential.eval(_xyz(r), acc=True, der=True)
                omega2 = -force[:, 0]/r
                ratio[circular] = np.sqrt((-deriv[:, 0]+3*omega2)/omega2)
            aa = np.zeros((int(general.sum()), 6))
            aa[:, 0], aa[:, 1], aa[:, 3:] = jr[general], angular[general], 1.
            if len(aa):
                xv = self.mapper(aa)
                # Interpolated Hamiltonian derivatives may turn negative on
                # very eccentric weakly bound orbits. Use direct integration.
                _, omega = self.agama.actions(self.potential, xv, frequencies=True)
                ratio[general] = omega[:, 0]/omega[:, 1]
        else:
            ell = angular+1.4*jr
            r = np.atleast_1d(self.potential.Rcirc(L=np.maximum(ell, 1e-100)))
            force, deriv = self.potential.eval(_xyz(r), acc=True, der=True)
            omega2 = -force[:, 0]/r
            kappa2 = -deriv[:, 0]+3*omega2
            c = np.divide(angular, jr+angular, out=np.zeros_like(jr), where=jr+angular > 0)
            ratio = 1/(.5+c*(np.sqrt(omega2/kappa2)-.5))
        # Analytic radial-orbit limit, including the joint origin by continuity.
        ratio[angular == 0] = 2.
        if np.any(~np.isfinite(ratio)) or np.any(ratio <= 0):
            bad = (~np.isfinite(ratio)) | (ratio <= 0)
            raise DFConvergenceError(f"invalid spherical frequencies at actions "
                                     f"{np.column_stack((jr[bad], angular[bad]))[:3].tolist()}")
        return ratio.reshape(shape)


class _ContourSpline:
    """Tensor cubic interpolation, with characteristic-derived c derivatives.

    The radial boundary obeys h_c=-h_x/2 exactly for the interpolated h=ln Lc,
    not merely at the tabulated nodes.  Thus interpolation preserves the
    low-L derivative ratio. No clipping or repair of negative DFs is used.
    """
    def __init__(self, x, c, h, slopes):
        ch = CubicHermiteSpline(c, h.T, slopes.T, axis=0)
        self.coeff = CubicSpline(x, np.moveaxis(ch.c, -1, 0), axis=0).c
        self.x, self.c = x, c
        # Rebuild the first Hermite segment with d0=-h0'/2 as a polynomial
        # identity in x. h0/h1 and d1 are the existing x-spline polynomials.
        hx = CubicSpline(x, h, axis=0)
        h0, h1 = hx.c[:, :, 0], hx.c[:, :, 1]
        d0 = np.zeros_like(h0)
        d0[1:] = -.5*np.arange(3, 0, -1)[:, None]*h0[:-1]
        d1 = CubicSpline(x, slopes[:, 1]).c
        dc = c[1]-c[0]
        self.coeff[:, :, 3, 0] = h0
        self.coeff[:, :, 2, 0] = d0
        self.coeff[:, :, 1, 0] = 3*(h1-h0)/dc**2-(2*d0+d1)/dc
        self.coeff[:, :, 0, 0] = 2*(h0-h1)/dc**3+(d0+d1)/dc**2

    def __call__(self, x, c):
        x, c = np.broadcast_arrays(x, c)
        ix = np.clip(np.searchsorted(self.x, x)-1, 0, len(self.x)-2)
        ic = np.clip(np.searchsorted(self.c, c)-1, 0, len(self.c)-2)
        dx, dc = x-self.x[ix], c-self.c[ic]
        result = np.zeros_like(x)
        for i in range(4):
            row = np.zeros_like(x)
            for j in range(4):
                row = row*dc+self.coeff[i, ix, j, ic]
            result = result*dx+row
        return result


class RegularizedExponentialDF:
    """Callable AGAMA-compatible stellar DF normalized in action space.

    q=Jr+L, c=L/q. A contour satisfies
    d ln(q)/dc = (g-1)/(c+g(1-c)). We integrate from c=1 (q=Lc),
    then invert the monotone chart onto a rectangular (ln q,c) grid.
    """
    def __init__(self, mass, parameters=None, *, potential=None, frequency_ratio=None,
                 numerics=None):
        if not np.isfinite(mass) or mass <= 0:
            raise ValueError("stellar DF mass must be positive")
        self.mass = mass
        self.parameters = parameters or ExponentialParameters()
        self.numerics = numerics or ContourNumerics()
        if (potential is None) == (frequency_ratio is None):
            raise ValueError("supply exactly one potential or test frequency reference")
        self.frequency_ratio = (frequency_ratio if frequency_ratio is not None else
                                SphericalFrequencyRatio(potential, self.numerics.frequency_mode))
        p, n = self.parameters, self.numerics
        # Generous finite map domain. Outside it, integrate the contour directly;
        # never extrapolate a spline into unsupported actions or truncate the DF.
        tmax = gammainccinv(3/p.alpha, 1e-14)
        padding = 1.5+max(abs(p.b_out), abs(p.b_outer or 0.))
        self.x = np.linspace(np.log(p.J0)-14., np.log(p.J0)+max(12., np.log(tmax)/p.alpha+padding),
                             n.action_nodes)
        self.c = np.linspace(0., 1., n.circularity_nodes)**2
        s = np.linspace(self.x[0]-padding, self.x[-1]+padding, n.action_nodes+64)

        def rhs(c, logq):
            return self._flow(np.exp(logq), c)

        solution = solve_ivp(rhs, (1., 0.), s, t_eval=self.c[::-1],
                             rtol=n.ode_rtol, atol=n.ode_rtol*.1)
        if not solution.success:
            raise DFConvergenceError("contour integration failed")
        chart = solution.y[:, ::-1]
        if np.any(np.diff(chart, axis=0) <= 0):
            raise DFConvergenceError("contours crossed or lost monotonicity")
        h = np.column_stack([CubicSpline(chart[:, i], s)(self.x)
                             for i in range(len(self.c))])
        hx = CubicSpline(self.x, h, axis=0)(self.x, 1)
        slopes = -self._flow(np.exp(self.x[:, None]), self.c[None, :])*hx
        self.spline = _ContourSpline(self.x, self.c, h, slopes)
        self.integral = self.normalization_integral(n.normalization_order)
        refined = self.normalization_integral(2*n.normalization_order)
        self.normalization_error = abs(self.integral/refined-1)
        if not np.isfinite(refined) or refined <= 0 or self.normalization_error > n.normalization_tolerance:
            raise DFConvergenceError(f"action normalization error {self.normalization_error:g}")
        self.integral = refined
        self.amplitude = mass/((2*np.pi)**3*self.integral)

    def g(self, jr, angular):
        jr, angular = np.broadcast_arrays(jr, angular)
        q = jr+angular
        c = np.divide(angular, q, out=np.zeros_like(q, dtype=float), where=q > 0)
        return self.frequency_ratio(jr, angular)*np.exp(-self.parameters.bias(q)*np.sin(np.pi*c/2))

    def _flow(self, q, c):
        g = self.g(q*(1-c), q*c)
        return (g-1)/(c+g*(1-c))

    def direct_contour_label(self, jr, angular):
        """Independent endpoint integration, also used outside the lookup grid."""
        jr, angular = np.broadcast_arrays(np.asarray(jr, float), np.asarray(angular, float))
        q = jr+angular
        if np.any(~np.isfinite(q)) or np.any(jr < 0) or np.any(angular < 0):
            raise ValueError("actions must be finite and nonnegative")
        shape = q.shape
        positive = q.ravel() > 0
        answer = np.zeros(q.size)
        if not np.any(positive):
            return answer.reshape(shape)
        q0 = q.ravel()[positive]
        c0 = angular.ravel()[positive]/q0

        def rhs(t, logq):
            c = c0+(1-c0)*t
            return (1-c0)*self._flow(np.exp(logq), c)

        solution = solve_ivp(rhs, (0., 1.), np.log(q0), rtol=self.numerics.ode_rtol,
                             atol=self.numerics.ode_rtol*.1)
        if not solution.success:
            raise DFConvergenceError("direct contour integration failed")
        answer[positive] = np.exp(solution.y[:, -1])
        return answer.reshape(shape)

    def contour_label(self, jr, angular):
        jr, angular = np.broadcast_arrays(np.asarray(jr, float), np.asarray(angular, float))
        q = jr+angular
        if np.any(~np.isfinite(q)) or np.any(jr < 0) or np.any(angular < 0):
            raise ValueError("actions must be finite and nonnegative")
        result = np.zeros(q.shape)
        positive = q > 0
        x = np.log(q[positive])
        inside = (x >= self.x[0]) & (x <= self.x[-1])
        values = np.empty(len(x))
        values[inside] = np.exp(self.spline(x[inside], angular[positive][inside]/q[positive][inside]))
        if np.any(~inside):
            values[~inside] = self.direct_contour_label(jr[positive][~inside], angular[positive][~inside])
        result[positive] = values
        return result

    def normalization_integral(self, order, extend=0.):
        """2 int c q^3 exp[-(Lc/J0)^alpha] dln(q) dc; tail checks can extend domain."""
        z, w = np.polynomial.legendre.leggauss(order)
        lo, hi = self.x[0]-extend, self.x[-1]+extend
        x, c = lo+(z+1)*(hi-lo)/2, (z+1)/2
        q = np.exp(x[:, None])
        lc = self.contour_label(q*(1-c), q*c)
        values = 2*c*q**3*np.exp(-(lc/self.parameters.J0)**self.parameters.alpha)
        return float(np.einsum("i,ij,j->", w, values, w)*(hi-lo)/4)

    def __call__(self, actions):
        actions = np.asarray(actions, float)
        if actions.shape[-1:] != (3,) or np.any(~np.isfinite(actions)):
            raise ValueError("expected finite (...,3) actions")
        if np.any(actions[..., :2] < 0):
            raise ValueError("Jr and Jz must be nonnegative")
        angular = actions[..., 1]+abs(actions[..., 2])
        # Lc >= L because g>0. Avoid computing contours whose DF is already
        # guaranteed below floating-point range, without imposing a tail cut.
        cutoff = (np.log(self.amplitude)-np.log(np.nextafter(0., 1.))+2)**(1/self.parameters.alpha)
        active = angular < self.parameters.J0*cutoff
        result = np.zeros(angular.shape)
        lc = self.contour_label(actions[..., 0][active], angular[active])
        with np.errstate(over="ignore"):
            result[active] = self.amplitude*np.exp(-(lc/self.parameters.J0)**self.parameters.alpha)
        return result

    def totalMass(self):
        """Fresh quadrature, not a return of the requested mass."""
        return (2*np.pi)**3*self.amplitude*self.normalization_integral(
            2*self.numerics.normalization_order)


@dataclass(frozen=True)
class PrescribedMatter:
    rho20: float = 0.
    r_s: float = 30.
    gamma: float = 0.
    r_t: float = 1000.
    M_rem: float = 0.
    a_rem: float = 3.

    def __post_init__(self):
        _finite(**asdict(self))
        if min(self.rho20, self.M_rem) < 0 or min(self.r_s, self.r_t, self.a_rem) <= 0:
            raise ValueError("invalid prescribed density mass/scale")
        if self.gamma not in (0., 1.):
            raise ValueError("choose gamma=0 (core) or gamma=1 (cusp) as separate models")

    def halo_density(self, r):
        r = np.asarray(r, float)
        if np.any(~np.isfinite(r)) or np.any(r < 0):
            raise ValueError("halo radii must be finite and nonnegative")
        if not self.rho20:
            return np.zeros_like(r)
        x, x20 = r/self.r_s, 20/self.r_s
        with np.errstate(divide="ignore", over="ignore"):
            logshape = (-self.gamma*np.log(np.maximum(x, np.finfo(float).tiny)/x20)
                        +(self.gamma-3)*(np.log1p(x)-np.log1p(x20))
                        -(r/self.r_t)**2+(20/self.r_t)**2)
            density = self.rho20*np.exp(logshape)
        return np.where((r == 0) & (self.gamma == 1), np.inf, density)

    def potentials(self, numerics):
        agama = agama_pc()
        result = []
        if self.M_rem:
            result.append(agama.Potential(type="Plummer", mass=self.M_rem, scaleRadius=self.a_rem))
        if self.rho20:
            density = agama.Density(density=lambda xyz: self.halo_density(np.linalg.norm(xyz, axis=1)),
                                    symmetry="s")
            result.append(agama.Potential(type="Multipole", density=density, lmax=0,
                                          gridSizeR=numerics.potential_nodes,
                                          rmin=numerics.r_min, rmax=numerics.r_max))
        return result


@dataclass(frozen=True)
class RegularizedDFConfig:
    M_star: float = 3e6
    stellar: ExponentialParameters = field(default_factory=ExponentialParameters)
    matter: PrescribedMatter = field(default_factory=PrescribedMatter)
    numerics: DFNumerics = field(default_factory=DFNumerics)
    contours: ContourNumerics = field(default_factory=ContourNumerics)
    distance_kpc: float = 5.43
    seed_radius: float = 5.

    def __post_init__(self):
        if not all(np.isfinite([self.M_star, self.distance_kpc, self.seed_radius])) or min(
                self.M_star, self.distance_kpc, self.seed_radius) <= 0:
            raise ValueError("stellar mass, distance and seed radius must be positive")

    def to_dict(self):
        return dict(schema_version=1, family="regularized_exponential", **asdict(self))

    @classmethod
    def from_dict(cls, payload):
        row = dict(payload)
        if row.pop("schema_version") != 1 or row.pop("family") != "regularized_exponential":
            raise ValueError("unsupported regularized DF schema")
        for k, constructor in (("stellar", ExponentialParameters), ("matter", PrescribedMatter),
                               ("numerics", DFNumerics), ("contours", ContourNumerics)):
            row[k] = constructor(**row[k])
        return cls(**row)


def _resolved_density(agama, radii, density):
    """Discard only underflow-dominated grid tails before a log interpolation.

    The continuation is a fitted declining power law in an already negligible
    tail, not a floor in f. Stellar mass closure independently checks its mass.
    """
    if np.any(~np.isfinite(density)) or np.any(density < 0) or not np.any(density > 0):
        raise DFConvergenceError("invalid stellar density")
    valid = density > density.max()*1e-30
    last = np.flatnonzero(valid)[-1]
    if last < 10 or not np.all(valid[:last+1]):
        raise DFConvergenceError("stellar grid does not resolve a contiguous central density")
    return _density_interpolator(agama, radii[:last+1], density[:last+1])


class RegularizedDFModel(PositiveDFModel):
    """Self-gravitating stars plus prescribed dark densities; no central BH.

    Convergence compares stellar density where rho > 1e-12 peak OR the
    stellar mass per log-radius exceeds 1e-8 Mstar, plus total force everywhere.
    This avoids a meaningless relative error on exponential numerical zeros.
    """
    def __init__(self, config: RegularizedDFConfig, *, initial_stellar_potential=None, progress=None):
        self.config = config
        agama, n = agama_pc(), config.numerics
        self.static_potentials = config.matter.potentials(n)
        star = (initial_stellar_potential if initial_stellar_potential is not None else
                agama.Potential(type="Plummer", mass=config.M_star, scaleRadius=config.seed_radius))
        self.potential = agama.Potential(star, *self.static_potentials)
        r = np.geomspace(n.r_min, n.r_max, n.potential_nodes)
        points = _xyz(r)
        previous = None
        self.iterations = []
        for i in range(n.max_iterations):
            df = RegularizedExponentialDF(config.M_star, config.stellar, potential=self.potential,
                                          numerics=config.contours)
            af = agama.ActionFinder(self.potential)
            density = spherical_velocity_moments(self.potential, df, af, r, n.velocity_nodes)[:, 0]
            self.stellar_density = _resolved_density(agama, r, density)
            star = agama.Potential(type="Multipole", density=self.stellar_density, lmax=0,
                                   gridSizeR=n.potential_nodes, rmin=n.r_min, rmax=n.r_max)
            potential = agama.Potential(star, *self.static_potentials)
            force = -potential.force(points)[:, 0]
            active = self._active(r, density)
            if previous is not None:
                old_force, old_density = previous
                mask = active | self._active(r, old_density)
                change = max(float(np.max(abs(force/old_force-1))),
                             float(np.max(abs(density[mask]/old_density[mask]-1))))
                self.iterations.append(change)
                if progress is not None:
                    progress(dict(iteration=i+1, max_relative_change=change))
                if change < n.iteration_tolerance:
                    self.potential = potential
                    break
            self.potential = potential
            previous = force, density
        else:
            raise DFConvergenceError(f"regularized DF equilibrium did not converge in {n.max_iterations} iterations")
        # Rebuild once more in the final potential: never retain the seed's map.
        self.df = RegularizedExponentialDF(config.M_star, config.stellar, potential=self.potential,
                                           numerics=config.contours)
        self.af = agama.ActionFinder(self.potential)
        self.stellar_potential = star
        self.galaxy = agama.GalaxyModel(self.potential, self.df, self.af)
        check_r = np.sqrt(r[:-1]*r[1:])
        direct = spherical_velocity_moments(self.potential, self.df, self.af, check_r,
                                            2*n.velocity_nodes)[:, 0]
        active = self._active(check_r, direct)
        reference = self.stellar_density.density(_xyz(check_r))
        closure = float(np.max(abs(direct[active]/reference[active]-1)))
        mass = float(star.totalMass())
        self.diagnostics = dict(iterations=i+1, last_change=self.iterations[-1],
                                density_closure=closure, stellar_mass=mass,
                                mass_error=abs(mass/config.M_star-1), df_mass=self.df.totalMass(),
                                normalization_error=self.df.normalization_error,
                                frequency_mode=config.contours.frequency_mode,
                                frequency_map_rebuilt=True,
                                closure_radius_pc=[float(check_r[active][0]), float(check_r[active][-1])])
        if closure > n.closure_tolerance or self.diagnostics["mass_error"] > n.mass_tolerance:
            raise DFConvergenceError(f"regularized DF closure failed: {self.diagnostics}")

    def _active(self, r, rho):
        return (rho > rho.max()*1e-12) | (4*np.pi*r**3*rho > self.config.M_star*1e-8)

    @cached_property
    def moment_table(self):
        n = self.config.numerics
        r = np.geomspace(n.r_min, n.r_max, n.moment_nodes)
        values = spherical_velocity_moments(self.potential, self.df, self.af, r, n.velocity_nodes)
        if np.any(~np.isfinite(values)) or np.any(values < 0):
            raise DFConvergenceError("invalid regularized DF moments")
        # Analytic logarithmic continuation only after all three moments have
        # fallen below 1e-30 of their peaks. Keeps legacy projection quadrature
        # well-defined over its full radius grid without altering the DF itself.
        for j in range(3):
            last = np.flatnonzero(values[:, j] > values[:, j].max()*1e-30)[-1]
            slope = np.log(values[last, j]/values[last-1, j])/np.log(r[last]/r[last-1])
            if slope >= -3:
                raise DFConvergenceError("moment grid does not reach a convergent stellar tail")
            # Stop continuation at the smallest representable positive value;
            # this is only an interpolation aid far below any resolved mass.
            logtail = np.log(values[last, j])+slope*np.log(r[last+1:]/r[last])
            values[last+1:, j] = np.exp(np.maximum(logtail, np.log(np.finfo(float).tiny)))
        rho, pr, pt = values.T
        return dict(r=r, rho=rho, radial_pressure=pr, tangential_pressure=pt, beta=1-pt/pr)
