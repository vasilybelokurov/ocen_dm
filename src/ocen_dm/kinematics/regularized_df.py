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


class _SphericalOrbitQuadrature:
    """Direct radial-orbit quadratures in a spherical potential (``self.potential``).

    Validated against analytic Kepler and isochrone actions (~1e-11) and against
    numerical orbit integration (~5e-11). Used where AGAMA's spherical action
    routines are not accurate enough: with an exact central point mass they are
    off by ~1e-4 for near-radial orbits (also in pure Kepler) and can fail for
    weakly bound ones (Jr = 0 at 1e4 pc).
    """
    def _phi_radial(self, r):
        r = np.asarray(r, float)
        return self.potential.potential(np.column_stack((r.ravel(), 0*r.ravel(), 0*r.ravel()))).reshape(r.shape)

    def _circular_radius(self, L, iterations=400):
        """R_circ(L) from r^3 dPhi/dr = L^2 by bisection in log r (monotone for any
        nonnegative spherical density)."""
        L = np.asarray(L, float)
        lo, hi = np.full_like(L, np.log(1e-30)), np.full_like(L, np.log(1e30))
        for _ in range(iterations):
            mid = .5*(lo+hi)
            r = np.exp(mid)
            force = self.potential.force(np.column_stack((r, 0*r, 0*r)))[:, 0]
            inside = -force*r**3 < L*L
            lo, hi = np.where(inside, mid, lo), np.where(inside, hi, mid)
            if np.all(hi-lo < 1e-14):
                break
        return np.exp(.5*(lo+hi))

    def _turning_points(self, energy, L, r0=None, iterations=200):
        """Peri- and apocentre by bisection in log r, bracketed outwards from an interior radius.

        Rows where r0 is missing or not strictly inside the orbit use R_circ(L),
        which lies inside every non-circular orbit; an orbit with F(R_circ) <= 0
        is circular to rounding and returns r_p = r_a = R_circ.
        """
        energy, L = np.asarray(energy, float), np.asarray(L, float)
        if np.any(~np.isfinite(energy)) or np.any(energy >= 0) or np.any(~(L > 0)):
            raise DFConvergenceError("turning points need bound orbits with L > 0")

        def F(r):          # > 0 between the turning points
            return 2*(energy-self._phi_radial(r))*r*r-L*L

        r0 = self._circular_radius(L) if r0 is None else np.array(r0, float)
        outside = ~(F(r0) > 0)
        if np.any(outside):
            r0[outside] = self._circular_radius(L[outside])

        def bisect(lo, hi, inside_is_hi):
            for _ in range(iterations):
                mid = np.sqrt(lo*hi)
                positive = F(mid) > 0
                if inside_is_hi:
                    hi, lo = np.where(positive, mid, hi), np.where(positive, lo, mid)
                else:
                    lo, hi = np.where(positive, mid, lo), np.where(positive, hi, mid)
                if np.all(hi/lo-1 < 1e-15):
                    break
            return np.sqrt(lo*hi)
        r_out = r0*2.
        for _ in range(1000):
            grow = F(r_out) > 0
            if not np.any(grow):
                break
            r_out = np.where(grow, r_out*2, r_out)
        if np.any(F(r_out) > 0):
            raise DFConvergenceError("apocentre not bracketed in radial-orbit quadrature")
        r_in = r0*.5
        for _ in range(1000):       # F -> -L^2 < 0 as r -> 0 (point-mass potential ~ -GM/r)
            shrink = F(r_in) > 0
            if not np.any(shrink):
                break
            r_in = np.where(shrink, r_in*.5, r_in)
        if np.any(F(r_in) > 0):
            raise DFConvergenceError("pericentre not bracketed in radial-orbit quadrature")
        return bisect(r_in, r0, inside_is_hi=True), bisect(r0, r_out, inside_is_hi=False)

    def _quadratures_at(self, energy, L, rp, ra, nodes):
        z, w = np.polynomial.legendre.leggauss(nodes)
        theta = .5*np.pi*(z+1)
        sin = np.sin(theta)[None, :]
        quad = lambda f: .5*np.pi*np.sum(w[None, :]*f, axis=1)

        def radicand(r, rad):
            peak = np.max(np.abs(rad), axis=1, keepdims=True)
            if np.any(rad < -1e-6*peak):
                raise DFConvergenceError("negative radicand inside the orbit (non-monotone potential?)")
            return np.maximum(rad, 0.)
        up, ua = 1/rp, 1/ra
        um, du = .5*(up+ua), .5*(up-ua)
        u = um[:, None]+du[:, None]*np.cos(theta)[None, :]
        rad = radicand(1/u, 2*(energy[:, None]-self._phi_radial(1/u))-(L[:, None]*u)**2)
        with np.errstate(divide="ignore", invalid="ignore"):
            dphi = 2*quad(np.where(rad > 0, L[:, None]*du[:, None]*sin/np.sqrt(rad), 0.))
        half = .5*np.log(ra/rp)
        r = rp[:, None]*np.exp(2*half[:, None]*np.sin(theta/2)[None, :]**2)
        pr2 = radicand(r, 2*(energy[:, None]-self._phi_radial(r))-(L[:, None]/r)**2)
        with np.errstate(divide="ignore", invalid="ignore"):
            inv = np.where(pr2 > 0, r*half[:, None]*sin/np.sqrt(pr2), 0.)
        jr = quad(r*half[:, None]*sin*np.sqrt(pr2))/np.pi
        djr = quad(inv)/np.pi
        return dphi, jr, djr

    def _radial_quadratures(self, energy, L, r0=None, nodes=96, rtol=(1e-6, 1e-9, 1e-6), max_nodes=1536):
        """Delta phi, Jr and dJr/dE = T_r/(2 pi), vectorised over orbits.

        Delta phi uses u = 1/r = u_mid + du cos(theta) (pericentre-weighted; exactly
        2 pi in pure Kepler); Jr and T_r use ln r = ln r_p + Delta (1-cos theta), which
        resolves the outer orbit that dominates them at extreme eccentricity
        (r_p/r_a ~ 1e-20 around a point mass). Both substitutions remove the
        inverse-square-root singularities at the turning points. The node count is
        doubled per orbit until successive orders agree to rtol (Delta phi, Jr, T_r; None = not required);
        Delta phi and T_r have a ~1e-7 rounding floor on near-circular orbits.
        """
        energy, L = np.asarray(energy, float), np.asarray(L, float)
        rp, ra = self._turning_points(energy, L, r0)
        out = np.array(self._quadratures_at(energy, L, rp, ra, nodes))
        todo = np.arange(len(energy))
        while len(todo):
            if 2*nodes > max_nodes:
                raise DFConvergenceError(f"radial-orbit quadrature unconverged at {nodes} nodes")
            nodes *= 2
            finer = np.array(self._quadratures_at(energy[todo], L[todo], rp[todo], ra[todo], nodes))
            # Jr is judged relative to Jr + L (the DF depends on the actions through
            # it): near-circular orbits have a ~1e-9 relative rounding floor in Jr alone.
            scale = np.abs(finer)
            scale[1] += L[todo]
            with np.errstate(divide="ignore", invalid="ignore"):
                change = np.abs(finer-out[:, todo])/scale
            change = np.where(finer == out[:, todo], 0., change)
            out[:, todo] = finer
            tol = np.array([np.inf if t is None else t for t in rtol])[:, None]   # None: not required
            todo = todo[np.any(~(change <= tol), axis=0)]
        return out[0], out[1], out[2]

    def _apsidal_ratio(self, energy, L, r0=None):
        """Omega_r/Omega_t = 2 pi/Delta phi, Delta phi = 2 int L du/sqrt(2(E-Phi(1/u))-L^2 u^2).

        E is taken from the phase-space point mapped from the actions (or from
        inverting Jr(E, L) where the mapper fails); r0 is a radius on the orbit.
        """
        dphi = self._radial_quadratures(energy, L, r0, rtol=(1e-5, None, None))[0]
        if np.any(~np.isfinite(dphi)) or np.any(dphi <= 0):
            raise DFConvergenceError("invalid apsidal angle")
        return 2*np.pi/dphi

    def _energy_from_actions(self, jr, L, r_circ=None, tolerance=1e-12, iterations=300):
        """E with Jr(E, L) = jr: safeguarded Newton in w = 1/sqrt(-2E).

        Jr(w) is exactly linear in Kepler (Jr = GM w - L) and monotone in general;
        the bracket starts at the circular orbit (Jr = 0) and is open above.
        """
        jr, L = np.asarray(jr, float), np.asarray(L, float)
        r_circ = self._circular_radius(L) if r_circ is None else np.asarray(r_circ, float)
        e_circ = self._phi_radial(r_circ)+.5*(L/r_circ)**2
        lo, hi = 1/np.sqrt(-2*e_circ), np.full_like(e_circ, np.inf)
        gm = -self._phi_radial(r_circ)*r_circ
        w = lo+jr/gm
        previous = np.full_like(w, np.inf)
        for _ in range(iterations):
            energy = -.5/w**2
            _, j, dj = self._radial_quadratures(energy, L, r_circ, rtol=(None, 1e-9, None))   # T_r only steers Newton
            resid = j-jr
            if np.all(np.abs(resid) <= tolerance*jr):
                break
            lo, hi = np.where(resid < 0, w, lo), np.where(resid < 0, hi, w)
            with np.errstate(divide="ignore", invalid="ignore"):
                step = w-resid/(dj/w**3)                     # dJr/dw = dJr/dE / w^3
            # Newton unless it leaves the bracket or the residual failed to halve
            # (Jr(E) steepens where the orbit leaves the point mass's sphere of influence).
            newton = (step > lo) & (step < hi) & (np.abs(resid) < .5*previous)
            fallback = np.where(np.isfinite(hi), np.sqrt(lo*hi), 2*lo)
            w = np.where(newton, step, fallback)
            previous = np.where(newton, np.abs(resid), np.inf)
        energy = -.5/w**2
        if np.any(np.abs(self._radial_quadratures(energy, L, r_circ, rtol=(None, 1e-9, None))[1]/jr-1) > 1e-9):
            raise DFConvergenceError("radial-action inversion failed")
        return energy


class ExactSphericalActions(_SphericalOrbitQuadrature):
    """AGAMA-compatible action finder by direct quadrature: N x 6 points -> (Jr, Jz, Jphi).

    For spherical potentials with an exact central point mass (see
    _SphericalOrbitQuadrature); several times slower than agama.ActionFinder.
    """
    def __init__(self, potential, chunk=8192):
        self.potential, self.agama, self.chunk = potential, agama_pc(), chunk

    def __call__(self, points):
        p = np.asarray(points, float)
        if p.ndim != 2 or p.shape[1] != 6 or np.any(~np.isfinite(p)):
            raise ValueError("expected a finite N x 6 array")
        x, v = p[:, :3], p[:, 3:]
        r = np.linalg.norm(x, axis=1)
        energy = self._phi_radial(r)+.5*np.sum(v*v, axis=1)
        if np.any(~(energy < 0)):
            raise DFConvergenceError("actions requested for unbound or non-finite points")
        lvec = np.cross(x, v)
        lz, L = lvec[:, 2], np.linalg.norm(lvec, axis=1)
        # Radial rows: Jr is continuous in L, so a relative 1e-14 floor changes it negligibly.
        L_eff = np.maximum(L, 1e-14*r*np.linalg.norm(v, axis=1))
        jr = np.empty(len(p))
        for s in range(0, len(p), self.chunk):
            k = slice(s, s+self.chunk)
            jr[k] = self._radial_quadratures(energy[k], L_eff[k], r[k].copy(), rtol=(None, 1e-9, None))[1]
        return np.column_stack((jr, L-np.abs(lz), lz))


class GuardedActionFinder:
    """agama.ActionFinder for point-mass models, with exact quadrature where its errors could add mass.

    Compared with exact actions, AGAMA's moves the kinematic bin predictions of point-mass models by
    <= 0.01 sigma (bin/compare_df_actions.py), but for some DFs it returns Jr ~ 0 for slow orbits far
    outside the cluster, where the true DF has underflowed: spurious mass in the outermost grid cells.
    Points beyond r_guard whose AGAMA actions give f > 0, and any non-finite AGAMA actions, are
    recomputed exactly. Not caught:
    an AGAMA overestimate of Jr out there, which could only remove mass where f is negligible anyway.
    """
    def __init__(self, potential, df, r_guard=300.):
        self.af, self.exact = agama_pc().ActionFinder(potential), ExactSphericalActions(potential)
        self.df, self.r_guard = df, r_guard

    def __call__(self, points):
        points = np.asarray(points, float)
        actions = np.array(self.af(points))
        broken = ~np.all(np.isfinite(actions), axis=1) | np.any(actions[:, :2] < 0, axis=1)
        if np.any(broken):          # AGAMA can also return NaN actions around a point mass
            actions[broken] = self.exact(points[broken])
        far = np.flatnonzero(np.linalg.norm(points[:, :3], axis=1) > self.r_guard)
        if len(far):
            suspect = far[self.df(actions[far]) > 0]
            if len(suspect):
                actions[suspect] = self.exact(points[suspect])
        return actions


class SphericalFrequencyRatio(_SphericalOrbitQuadrature):
    """Omega_r/Omega_t in a spherical model, with a finite centre or an exact central point mass.

    ActionMapper returns spherical frequencies without numerical differencing
    of the Hamiltonian.  The optional epicycle prescription is an explicitly
    different, approximate reference, evaluated pointwise along every contour.

    The radial-orbit limit (L -> 0) depends on the centre: an orbit through a
    finite (harmonic) core advances by pi per radial period, so the ratio is 2;
    with a central point mass the near-radial orbit swings by 2 pi around it at
    pericentre, so the ratio is 1, at every energy.  ``radial_limit`` holds it.
    """
    def __init__(self, potential, mode="exact"):
        agama = agama_pc()
        if mode not in ("exact", "epicycle"):
            raise ValueError("unknown frequency reference")
        self.point_mass = not np.isfinite(potential.potential([0., 0., 0.]))
        if self.point_mass and mode != "exact":
            raise ValueError("the epicycle frequency reference assumes a finite centre; use mode='exact'")
        self.radial_limit = 1. if self.point_mass else 2.
        self.potential, self.mode, self.agama = potential, mode, agama
        self.mapper = agama.ActionMapper(potential)

    def __call__(self, jr, angular):
        jr, angular = np.broadcast_arrays(jr, angular)
        shape = jr.shape
        jr, angular = jr.ravel(), angular.ravel()
        if self.mode == "exact":
            ratio = np.full(len(jr), self.radial_limit)
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
            # Spherical frequencies depend on Jr and L, not the orbital plane.
            # Use Jz=0, Jphi=L: the polar representative Jz=L, Jphi=0 can
            # hit R=0 in AGAMA's map and divide by R in its velocity transform.
            aa[:, 0], aa[:, 2], aa[:, 3:] = jr[general], angular[general], 1.
            if len(aa):
                xv = self.mapper(aa)
                if self.point_mass:
                    # AGAMA's frequencies are inaccurate for highly eccentric orbits
                    # around a point mass (e.g. ratio 0.69 instead of 1 in pure Kepler
                    # at L/Jr = 1e-3); integrate the apsidal angle directly instead.
                    # The mapper also returns NaN, or unbound points (E >= 0), at large
                    # Jr (>~1e6, apocentres >~1e8 pc) around a point mass: invert Jr(E, L)
                    # there (exact, but ~50x slower than the mapper).
                    L = angular[general]
                    finite = np.all(np.isfinite(xv), axis=1)
                    r0 = np.where(finite, np.linalg.norm(np.where(finite[:, None], xv, 1.)[:, :3], axis=1), 1.)
                    energy = np.where(finite, self._phi_radial(r0)+.5*np.sum(np.where(finite[:, None], xv, 0.)[:, 3:]**2, axis=1), 0.)
                    lost = ~finite | (energy >= 0)
                    target = jr[general]
                    dphi = np.empty(len(L))
                    kept = np.flatnonzero(~lost)
                    if len(kept):
                        # Also reject finite mapped orbits whose own Jr misses the request by
                        # > 1% of Jr + L (seen: 20x too small at ~1e3 pc); its usual ~1e-3
                        # errors move the frequency ratio by < 3e-5.
                        d, j, _ = self._radial_quadratures(energy[kept], L[kept], r0[kept], rtol=(1e-5, 1e-9, None))
                        dphi[kept] = d
                        lost[kept[np.abs(j-target[kept]) > 1e-2*(target[kept]+L[kept])]] = True
                    if np.any(lost):
                        r0[lost] = self._circular_radius(L[lost])
                        energy[lost] = self._energy_from_actions(target[lost], L[lost], r0[lost])
                        dphi[lost] = self._radial_quadratures(energy[lost], L[lost], r0[lost], rtol=(1e-5, None, None))[0]
                    if np.any(~np.isfinite(dphi)) or np.any(dphi <= 0):
                        raise DFConvergenceError("invalid apsidal angle")
                    ratio[general] = 2*np.pi/dphi
                elif np.any(~np.isfinite(xv)):
                    raise DFConvergenceError("nonfinite mapped orbit in spherical frequency calculation")
                else:
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
        # Analytic radial-orbit limit (2 for a finite centre, 1 with a point mass),
        # including the joint origin by continuity.
        ratio[angular == 0] = self.radial_limit
        if np.any(~np.isfinite(ratio)) or np.any(ratio <= 0):
            bad = (~np.isfinite(ratio)) | (ratio <= 0)
            raise DFConvergenceError(f"invalid spherical frequencies at actions "
                                     f"{np.column_stack((jr[bad], angular[bad]))[:3].tolist()}")
        return ratio.reshape(shape)

class _ContourSpline:
    """Tensor cubic interpolation, with characteristic-derived c derivatives.

    The radial boundary obeys h_c=-k h_x exactly for the interpolated h=ln Lc,
    not merely at the tabulated nodes, with k=(g0-1)/g0 from the radial-orbit
    frequency ratio g0 (k=1/2 for a finite centre, 0 with a central point mass).
    Thus interpolation preserves the low-L derivative ratio. No clipping or
    repair of negative DFs is used.
    """
    def __init__(self, x, c, h, slopes, radial_factor=.5):
        ch = CubicHermiteSpline(c, h.T, slopes.T, axis=0)
        self.coeff = CubicSpline(x, np.moveaxis(ch.c, -1, 0), axis=0).c
        self.x, self.c = x, c
        # Rebuild the first Hermite segment with d0=-h0'/2 as a polynomial
        # identity in x. h0/h1 and d1 are the existing x-spline polynomials.
        hx = CubicSpline(x, h, axis=0)
        h0, h1 = hx.c[:, :, 0], hx.c[:, :, 1]
        d0 = np.zeros_like(h0)
        d0[1:] = -radial_factor*np.arange(3, 0, -1)[:, None]*h0[:-1]
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
        g0 = float(np.ravel(self.frequency_ratio(np.array([1.]), np.array([0.])))[0])
        self.spline = _ContourSpline(self.x, self.c, h, slopes, radial_factor=(g0-1)/g0)
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
        # Also Lc >= q exp(1-e^B): the frequency ratio is >= 1 for any nonnegative
        # spherical density, |bias| <= B, so the contour flow is >= 1-e^B. This keeps
        # near-escape orbits around a point mass (Jr ~ 1e12+) out of the action mapper.
        p = self.parameters
        bound = abs(p.b_out)+(abs(p.b_outer-p.b_out) if p.b_outer is not None else 0.)
        q = actions[..., 0]+angular
        cutoff = (np.log(self.amplitude)-np.log(np.nextafter(0., 1.))+2)**(1/p.alpha)
        active = (angular < p.J0*cutoff) & (q*np.exp(1-np.exp(bound)) < p.J0*cutoff)
        result = np.zeros(angular.shape)
        lc = self.contour_label(actions[..., 0][active], angular[active])
        with np.errstate(over="ignore"):
            result[active] = self.amplitude*np.exp(-(lc/p.J0)**p.alpha)
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
    M_bh: float = 0.          # exact central point mass (Msun); 0 = none

    def __post_init__(self):
        _finite(**asdict(self))
        if min(self.rho20, self.M_rem, self.M_bh) < 0 or min(self.r_s, self.r_t, self.a_rem) <= 0:
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
        if self.M_bh:
            result.append(agama.Potential(type="Plummer", mass=self.M_bh, scaleRadius=0.))  # exact point mass
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
    """Self-gravitating stars plus prescribed dark densities and an optional central point mass.

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
            af = self._action_finder(self.potential, df)
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
        self.af = self._action_finder(self.potential, self.df)
        self.stellar_potential = star
        # GalaxyModel (direct_projected_moments, validation only) needs AGAMA's own finder.
        self.galaxy = agama.GalaxyModel(self.potential, self.df, agama.ActionFinder(self.potential))
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

    def _action_finder(self, potential, df):
        if self.config.numerics.exact_actions:
            return ExactSphericalActions(potential)
        if self.config.matter.M_bh > 0:
            return GuardedActionFinder(potential, df)
        return agama_pc().ActionFinder(potential)

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
