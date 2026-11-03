"""Model-side centre marginalisation: data frozen about AvdM10, model centred at c.

Replaces the data-side construction of :mod:`ocen_dm.centre_marg` (re-binning the stars about
every trial centre), which averages likelihoods of *different* data summaries and whose
centre weights collapse onto one or two lucky re-binnings (effective sample size 1.3-3.3 of
36-100 nodes; docs/reviews/{codex,fable}_recentre_plan_20261103.md).

Here every binned product (counts, HST PM and MUSE dispersions, their errors) stays exactly as
measured about AvdM10. A model whose centre sits at c = (dRA*, dDec) [arcsec] predicts, for a
fixed bin k,

* counts: the model surface density averaged over the bin's *covered sky*,
  G_k(c) = <Sigma(|p - c|)>_p / <Sigma(|p|)>_p  (p: sky points of the bin, area weighted);
* HST PM dispersions: the per-star mean of the model variance **projected on the data's
  radial/tangential axes** (defined about AvdM10), with psi_i the angle between the star's
  radial direction about AvdM10 and about c,
  radial:     q_i = cos^2 psi_i s_R^2(r_ic) + sin^2 psi_i s_T^2(r_ic)
  tangential: q_i = sin^2 psi_i s_R^2(r_ic) + cos^2 psi_i s_T^2(r_ic),
  F_k(c) = <q(c)>_i / <q(0)>_i;
* MUSE LOS dispersion: F_k(c) = <s_los^2(r_ic)>_i / <s_los^2(r_i0)>_i over the bin's stars.

The factors multiply the frozen problem's own bin averages (DFJointProblem.evaluate(factors=)),
so c = 0 reproduces the frozen likelihood exactly and only the *change* with c uses the star /
sky-point averages. The streaming terms (MUSE v_rot^2/2, HST and Gaia external rotation) are
kept fixed: they belong to the fixed data. Gaia kinematic bins (300-2400") are not moved; the
effect of a few-arcsec shift there is second order (tested).

ln L_m(theta) = logsumexp_k [ln w_k + ln L(theta, c_k | d)]  -- a marginal likelihood of fixed data.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator

from .centre_marg import GAIA_NAME, HST_NAME, _hst_catalogue
from .kinematics.likelihood import ARCSEC_PER_RAD

__all__ = ["ModelMoments", "PointSet", "CentreGeometry", "projected_variance", "effective_sample_size"]

HST_PMR, HST_PMT, MUSE = "hst_pm_radial_ours", "hst_pm_tangential_ours", "muse_los_dispersion"


class ModelMoments:
    """Projected moments of a model tabulated on a log grid in R (arcsec), PCHIP in log-log.

    ``__call__(r)`` returns Sigma (arbitrary units) and the projected variances s2r, s2t, s2los
    (the model's own units: per-unit-Sigma second moments; only ratios are used). Radii below the
    grid are clamped to its first point.
    """

    def __init__(self, model, distance_kpc: float, r_min_arcsec: float = 1e-3,
                 r_max_arcsec: float = 4000., n: int = 900):
        pc = distance_kpc*1e3/ARCSEC_PER_RAD
        R = np.geomspace(r_min_arcsec, r_max_arcsec, n)
        m = model.projected_moments(R*pc)
        self._set(R, m)

    def _set(self, R, m):
        self.table = dict(R=np.asarray(R, float), **{k: np.asarray(m[k], float) for k in ("Sigma", "pmr", "pmt", "los")})
        self._lo = np.log(R[0])
        x = np.log(R)
        self._f = {"Sigma": PchipInterpolator(x, np.log(m["Sigma"])),
                   "s2r": PchipInterpolator(x, np.log(m["pmr"]/m["Sigma"])),
                   "s2t": PchipInterpolator(x, np.log(m["pmt"]/m["Sigma"])),
                   "s2los": PchipInterpolator(x, np.log(m["los"]/m["Sigma"]))}

    @classmethod
    def from_table(cls, table):
        """Rebuild from ``self.table`` (cached projected moments: R [arcsec], Sigma, pmr, pmt, los)."""
        self = cls.__new__(cls)
        self._set(table["R"], table)
        return self

    @classmethod
    def from_functions(cls, **funcs):
        """Analytic moments (tests, mocks): callables of r [arcsec] for Sigma, s2r, s2t, s2los."""
        self = cls.__new__(cls)
        self._lo = None
        self._funcs = funcs
        return self

    def __call__(self, r, keys=("Sigma", "s2r", "s2t", "s2los")):
        r = np.asarray(r, float)
        if self._lo is None:
            return {k: np.asarray(self._funcs[k](r), float) for k in keys}
        x = np.maximum(np.log(np.maximum(r, 1e-300)), self._lo)
        return {k: np.exp(self._f[k](x)) for k in keys}


@dataclass
class PointSet:
    """Points (stars or sky samples) of several bins: positions about AvdM10 [arcsec], bin index, weight."""

    x: np.ndarray
    y: np.ndarray
    idx: np.ndarray
    w: np.ndarray
    n_bins: int

    def mean(self, values):
        num = np.bincount(self.idx, self.w*values, minlength=self.n_bins)
        den = np.bincount(self.idx, self.w, minlength=self.n_bins)
        return num/den

    @classmethod
    def from_positions(cls, x, y, r_lower, r_upper, w=None):
        """Assign points to bins [r_lower, r_upper) by their radius about AvdM10 (outside: dropped)."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        w = np.ones_like(x) if w is None else np.asarray(w, float)
        R = np.hypot(x, y)
        lo, hi = np.asarray(r_lower, float), np.asarray(r_upper, float)
        if np.any(hi[:-1] > lo[1:]):
            raise ValueError("bins overlap")
        k = np.searchsorted(lo, R, side="right")-1
        ok = (k >= 0) & (R < hi[np.clip(k, 0, len(hi)-1)])
        return cls(x[ok], y[ok], k[ok], w[ok], len(lo))

    @classmethod
    def annuli(cls, r_lower, r_upper, n_r=24, n_phi=256):
        """Full annuli: Gauss-Legendre in R (weight R dR) x uniform phi (area weights)."""
        u, wu = np.polynomial.legendre.leggauss(n_r)
        phi = (np.arange(n_phi)+.5)/n_phi*2*np.pi
        xs, ys, ks, ws = [], [], [], []
        for k, (lo, hi) in enumerate(zip(r_lower, r_upper)):
            lo = max(lo, 1e-3*hi)
            r = lo+(u+1)/2*(hi-lo)
            wr = wu/2*(hi-lo)*r
            xs.append((r[:, None]*np.sin(phi)[None, :]).ravel())
            ys.append((r[:, None]*np.cos(phi)[None, :]).ravel())
            ws.append(np.repeat(wr, n_phi)*(2*np.pi/n_phi))
            ks.append(np.full(n_r*n_phi, k))
        return cls(np.concatenate(xs), np.concatenate(ys), np.concatenate(ks), np.concatenate(ws), len(r_lower))


class ContinuousSet:
    """Bins with full azimuthal coverage, the selected stars replaced by a smooth density.

    Averaging the model over the *discrete* stars fails near a central point mass: the model
    variance grows as 1/r_c, and a trial centre landing within ~0.03" of one star lets that star
    dominate its bin (ln L jumps by 20-40 between neighbouring centres for the 4e4 Msun model).
    The expectation over a smooth star density is finite (the 1/r cusp is integrable in 2-D), and
    our estimator does not respond to such stars anyway (mock test). Here the star surface density
    is a power law n(R0) ~ R0^-gamma (one ML slope for all these bins' stars), and the bin mean
    <f>_k(c) = int n f dA / int n dA is computed in polar coordinates (rho, theta) about c, where
    dA = rho drho dtheta removes the singularity. Composite Gauss-Legendre in rho, uniform theta.
    """

    def __init__(self, r_lower, r_upper, radii, n_panel=48, n_gl=4, n_theta=192):
        self.lo, self.hi = np.asarray(r_lower, float), np.asarray(r_upper, float)
        # one slope over all the bins (per-bin fits to 50-100 stars are noisy and hit the bounds)
        g = self._fit_gamma(np.concatenate(radii), self.lo.min(), self.hi.max())
        self.gamma = np.full(len(self.lo), g)
        self.n_panel, self.n_theta = n_panel, n_theta
        self._u, self._wu = np.polynomial.legendre.leggauss(n_gl)
        self.n_bins = len(self.lo)

    @staticmethod
    def _fit_gamma(r, lo, hi):
        """ML power-law slope of the surface density in [lo, hi): p(R) ~ R^(1-gamma)."""
        from scipy.optimize import minimize_scalar
        r = np.asarray(r, float)
        lo_ = max(lo, 1e-3*hi)
        lr = np.log(r).sum()

        def nll(g):
            a = 2.-g
            norm = np.log(hi/lo_) if abs(a) < 1e-8 else np.log(abs((hi**a-lo_**a)/a))
            return -((1-g)*lr-len(r)*norm)
        return float(minimize_scalar(nll, bounds=(-3., 1.9), method="bounded").x)

    def _points(self, k, c):
        """Quadrature points (x, y) about AvdM10 and weights n(R0) rho drho dtheta for bin k, centre c."""
        lo, hi = max(self.lo[k], 1e-3*self.hi[k]), self.hi[k]
        cr = float(np.hypot(*c))
        r1, r2 = max(0., lo-cr), hi+cr
        # panels: uniform in rho, plus refinement near rho = 0 (where the model cusp sits)
        e = np.unique(np.r_[np.linspace(r1, r2, self.n_panel+1), np.geomspace(max(r1, 1e-4), r2, self.n_panel+1)]) \
            if r1 == 0. else np.linspace(r1, r2, self.n_panel+1)
        if r1 == 0. and e[0] > 0:
            e = np.r_[0., e]
        a, b = e[:-1], e[1:]
        rho = (a[:, None]+(self._u[None, :]+1)/2*(b-a)[:, None]).ravel()
        wr = ((b-a)[:, None]/2*self._wu[None, :]).ravel()*rho
        th = (np.arange(self.n_theta)+.5)/self.n_theta*2*np.pi
        x = c[0]+(rho[:, None]*np.sin(th)[None, :]).ravel()
        y = c[1]+(rho[:, None]*np.cos(th)[None, :]).ravel()
        R0 = np.hypot(x, y)
        w = np.repeat(wr, self.n_theta)*(2*np.pi/self.n_theta)
        inside = (R0 >= self.lo[k]) & (R0 < hi)
        w = np.where(inside, w*np.maximum(R0, 1e-12)**(-self.gamma[k]), 0.)
        return x[inside], y[inside], w[inside]

    def mean(self, func, c):
        """Per-bin mean of func(x, y, c) (vectorised over points) under the smooth star density."""
        out = np.empty(self.n_bins)
        for k in range(self.n_bins):
            x, y, w = self._points(k, c)
            out[k] = np.sum(w*func(x, y, c))/np.sum(w)
        return out


def projected_variance(mom, x, y, c, component):
    """Model variance of each star along the data's axes (about AvdM10), model centred at c.

    component: 'r' or 't' (PM, uses s2r/s2t) or 'los'.
    """
    dx, dy = x-c[0], y-c[1]
    rc = np.hypot(dx, dy)
    if component == "los":
        return mom(rc, ("s2los",))["s2los"]
    m = mom(rc, ("s2r", "s2t"))
    r0 = np.hypot(x, y)
    cos = (x*dx+y*dy)/np.maximum(r0*rc, 1e-300)
    cos2 = np.clip(cos*cos, 0., 1.)
    if component == "r":
        return cos2*m["s2r"]+(1-cos2)*m["s2t"]
    return (1-cos2)*m["s2r"]+cos2*m["s2t"]


class CentreGeometry:
    """Fixed point sets of the frozen products and the per-bin model factors for a centre c."""

    def __init__(self, data, counts, hst_stars=None, muse_stars=None, r_cont=15.):
        prof = {p.name: p for p in data.profiles}
        self.sets, self.kind, self.cont = {}, {}, {}
        if hst_stars is None:
            from .kinematics.hst_profile import load_hst_sample
            s = load_hst_sample(require_flag=True)
            hst_stars = (s.r_arcsec*np.sin(s.phi), s.r_arcsec*np.cos(s.phi))
        if muse_stars is None:
            c = _hst_catalogue()
            muse_stars = (c["lx"], c["ly"])
        for name, comp in ((HST_PMR, "r"), (HST_PMT, "t")):
            if name in prof:
                p = prof[name]
                self.sets[name] = PointSet.from_positions(*hst_stars, p.r_lower, p.r_upper)
                self.kind[name] = comp
        if MUSE in prof:
            p = prof[MUSE]
            self.sets[MUSE] = PointSet.from_positions(*muse_stars, p.r_lower, p.r_upper)
            self.kind[MUSE] = "los"
        for cp in counts:
            if cp.name.startswith(HST_NAME):
                self.sets[cp.name] = self._hst_count_points(cp)
            elif cp.name.startswith(GAIA_NAME):
                self.sets[cp.name] = PointSet.annuli(cp.r_lower, cp.r_upper, n_r=32, n_phi=512)
            else:
                continue
            self.kind[cp.name] = "Sigma"
        for name, ps in self.sets.items():
            if np.any(np.bincount(ps.idx, minlength=ps.n_bins) == 0):
                raise ValueError(f"{name}: a bin has no points")
        # kinematic bins starting inside r_cont: smooth density (see ContinuousSet)
        for name in [n for n, k in self.kind.items() if k != "Sigma"]:
            p, ps = prof[name], self.sets[name]
            m = int(np.count_nonzero(p.r_lower < r_cont))
            if m:
                radii = [np.hypot(ps.x, ps.y)[ps.idx == k] for k in range(m)]
                self.cont[name] = ContinuousSet(p.r_lower[:m], p.r_upper[:m], radii)

    @staticmethod
    def _hst_count_points(cp, full_inside=30.):
        """Covered sky of each HST count bin: full annuli inside ``full_inside`` (as the builder
        assumes), the occupied-cell point set (0.3" sub-cells of 3" cells) outside."""
        cat = _hst_catalogue()
        inner = cp.r_upper <= full_inside
        ann = PointSet.annuli(cp.r_lower[inner], cp.r_upper[inner], n_r=24, n_phi=256)
        k0 = int(inner.sum())
        cells = PointSet.from_positions(cat["px"], cat["py"], cp.r_lower[~inner], cp.r_upper[~inner],
                                        w=np.full(len(cat["px"]), cat["dA"]))
        return PointSet(np.r_[ann.x, cells.x], np.r_[ann.y, cells.y], np.r_[ann.idx, cells.idx+k0],
                        np.r_[ann.w, cells.w], cp.n)

    def _means(self, mom, name, c):
        ps, kind = self.sets[name], self.kind[name]
        if kind == "Sigma":
            v = mom(np.hypot(ps.x-c[0], ps.y-c[1]), ("Sigma",))["Sigma"]
            return ps.mean(v)
        out = ps.mean(projected_variance(mom, ps.x, ps.y, c, kind))
        if name in self.cont:
            cs = self.cont[name]
            out[:cs.n_bins] = cs.mean(lambda x, y, cc: projected_variance(mom, x, y, cc, kind), c)
        return out

    def reference(self, mom):
        """Per-bin means at c = 0 (cache once per model)."""
        return {n: self._means(mom, n, (0., 0.)) for n in self.sets}

    def factors(self, mom, c, ref=None):
        ref = self.reference(mom) if ref is None else ref
        return {n: self._means(mom, n, c)/ref[n] for n in self.sets}

    def prediction(self, mom, name, c):
        """Absolute per-bin mean (variance for kinematics, Sigma for counts) at centre c (mocks)."""
        return self._means(mom, name, c)


def effective_sample_size(lnl, weights):
    """Kish ESS of the normalised centre contributions a_k = w_k L_k / sum_j w_j L_j."""
    lnl = np.asarray(lnl, float)
    a = np.asarray(weights, float)*np.exp(lnl-lnl.max())
    a = a/a.sum()
    return float(1/np.sum(a*a)), float(a.max())
