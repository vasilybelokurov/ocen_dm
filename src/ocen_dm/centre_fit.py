"""Joint measurement of omega Cen's centre from star counts and proper-motion dispersions (2026-10-20).

Design: docs/reviews/codex_centre_20261020.md and docs/reviews/fable_centre_20261020.md. Positions are offsets (arcsec)
from the adopted Anderson & van der Marel (2010) centre: x = dRA* (east +), y = dDec (north +).

Counts (positions of a magnitude-limited sample in a fixed circular window W of radius R_max about the origin):
    lambda(x) = C(x) Sigma(x),   Sigma = (1 + m^2/a^2)^(-b),   m^2 = u^2 + v^2/q^2  (u, v rotated by theta about c),
    C(x) = exp(beta_c D(x))  (completeness tied to the depth proxy D; beta_c free),
    optional lopsided term  [1 + A1 (R/100") cos(phi - phi1)]  (a radial law different from a translation's),
    ln L = sum_i ln lambda(x_i) - N ln int_W lambda dA   (conditional on N).
Dispersions (locally corrected PMs: zero mean everywhere, so only second moments carry information):
    v_i ~ N(0, Sigma_i),  Sigma_i = R(phi_i) diag(s_R^2, s_T^2) R(phi_i)^T + f_i^2 diag(e_ra^2, e_dec^2),
    s^2 = s0^2 (1 + m^2/rc^2)^(-g) (1 + eta (mag - 21))^2,  s_R^2 = s^2 (1 + kappa w),  s_T^2 = s^2 (1 - kappa w),
    w = R^2 / (R^2 + 200"^2)  (radial anisotropy growing outwards; its axis points at the centre),
    f_i = f0 exp(gamma D_i)  (PM-error scale tied to depth).
The two are combined with a shared centre: ln L = ln L_counts + ln L_disp (velocities conditional on positions).
Depth proxy D: per 20" cell, ln(median PM error of F625W 20-21 stars / its median over cells); 0 in cells with < 30 such
stars. It is built from errors, not counts, so it does not reuse the counted stars (corr. -0.58 with the faint/bright
ratio at 100-200").
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

CELL = 20.0


@dataclass
class Stars:
    x: np.ndarray
    y: np.ndarray
    mag: np.ndarray
    pmx: np.ndarray | None = None
    pmy: np.ndarray | None = None
    ex: np.ndarray | None = None
    ey: np.ndarray | None = None
    hq: np.ndarray | None = None

    def subset(self, k):
        return Stars(*[None if v is None else v[k] for v in (self.x, self.y, self.mag, self.pmx, self.pmy, self.ex, self.ey, self.hq)])


def load_omegacat(rmax=210.):
    """oMEGACat VI stars with F625W within rmax of the adopted centre."""
    from astropy.table import Table
    from .paths import raw_dir
    t = Table.read(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")
    xp, yp, m = (np.asarray(t[c], float) for c in ("x", "y", "f625w"))
    ok = np.isfinite(xp) & np.isfinite(yp) & np.isfinite(m)
    x, y = -0.04*(xp-15000.), 0.04*(yp-15000.)
    ok &= np.hypot(x, y) < rmax
    g = lambda c: np.asarray(t[c], float)[ok]
    return Stars(x[ok], y[ok], m[ok], g("pmra_corrected"), g("pmdec_corrected"), g("pmra_corrected_err"), g("pmdec_corrected_err"),
                 np.asarray(t["selection_hq_astrometry_and_membership"]).astype(bool)[ok])


class DepthMap:
    """Cell lookup of the depth proxy D (see module docstring)."""

    def __init__(self, s: Stars, cell=CELL, half=220., mag=(20., 21.), nmin=30):
        self.cell, self.half = cell, half
        n = int(round(2*half/cell)); self.D = np.zeros((n, n))
        e = np.hypot(s.ex, s.ey)/np.sqrt(2); k = (s.mag >= mag[0]) & (s.mag < mag[1]) & np.isfinite(e)
        ix, iy = self._idx(s.x[k], s.y[k]); med = np.full((n, n), np.nan)
        order = np.lexsort((iy, ix)); ix, iy, ee = ix[order], iy[order], e[k][order]
        key = ix*n+iy; u, st, cnt = np.unique(key, return_index=True, return_counts=True)
        for kk, a, c in zip(u, st, cnt):
            if c >= nmin:
                med[kk//n, kk % n] = np.median(ee[a:a+c])
        ref = np.nanmedian(med); self.D = np.where(np.isfinite(med), np.log(med/ref), 0.)

    def _idx(self, x, y):
        n = self.D.shape[0]
        return (np.clip(((x+self.half)//self.cell).astype(int), 0, n-1), np.clip(((y+self.half)//self.cell).astype(int), 0, n-1))

    def __call__(self, x, y):
        return self.D[self._idx(np.asarray(x), np.asarray(y))]


def _ell_m2(x, y, cx, cy, q, th):
    dx, dy = x-cx, y-cy; c, s = np.cos(th), np.sin(th)
    u, v = dx*c+dy*s, -dx*s+dy*c
    return u*u+v*v/(q*q)


class CountModel:
    """Unbinned Poisson likelihood of positions in the window R < rmax (see module docstring)."""
    names = ("cx", "cy", "log_a", "log_b", "q", "theta", "beta_c", "A1", "phi1")

    def __init__(self, s: Stars, depth: DepthMap, rmax, lopsided=False, grid=1.0, mask_depth=None):
        k = np.hypot(s.x, s.y) < rmax
        g = np.arange(-rmax+grid/2, rmax, grid); X, Y = np.meshgrid(g, g); inside = np.hypot(X, Y) < rmax
        self.gx, self.gy, self.dA = X[inside], Y[inside], grid*grid
        self.Dg = depth(self.gx, self.gy); x, y = s.x[k], s.y[k]; Ds = depth(x, y)
        if mask_depth is not None:                                     # drop cells whose depth deviates by more than mask_depth
            keep, keepg = np.abs(Ds) <= mask_depth, np.abs(self.Dg) <= mask_depth
            x, y, Ds = x[keep], y[keep], Ds[keep]; self.gx, self.gy, self.Dg = self.gx[keepg], self.gy[keepg], self.Dg[keepg]
        self.x, self.y, self.Ds, self.lopsided = x, y, Ds, lopsided

    def _lam(self, x, y, D, p):
        cx, cy, la, lb, q, th, bc, A1, ph = p
        m2 = _ell_m2(x, y, cx, cy, q, th); out = np.exp(-np.exp(lb)*np.log1p(m2/np.exp(2*la))+bc*D)
        if self.lopsided:
            R = np.hypot(x-cx, y-cy); out = out*(1+A1*(R/100.)*np.cos(np.arctan2(y-cy, x-cx)-ph))
        return out

    def loglike(self, p):
        p = np.asarray(p, float)
        if not (0.7 < p[4] < 1.3 and -3 < p[3] < 3 and 1. < p[2] < 9. and abs(p[7]) < 0.3):
            return -np.inf
        ls, lg = self._lam(self.x, self.y, self.Ds, p), self._lam(self.gx, self.gy, self.Dg, p)
        if np.any(ls <= 0) or np.any(lg <= 0):
            return -np.inf
        return float(np.sum(np.log(ls))-len(self.x)*np.log(np.sum(lg)*self.dA))


class DispersionModel:
    """Per-star Gaussian likelihood of the locally corrected PMs (see module docstring)."""
    names = ("cx", "cy", "log_s0", "log_rc", "g", "q_s", "theta_s", "eta", "kappa", "log_f0", "gamma")

    def __init__(self, s: Stars, depth: DepthMap, rmax, mag=(-99, 99)):
        k = (np.hypot(s.x, s.y) < rmax) & s.hq & (s.mag >= mag[0]) & (s.mag < mag[1]) & np.isfinite(s.pmx) & np.isfinite(s.pmy) \
            & np.isfinite(s.ex) & np.isfinite(s.ey)
        self.s = s.subset(k); self.D = depth(self.s.x, self.s.y)

    def loglike(self, p):
        cx, cy, ls0, lrc, g, qs, ths, eta, kap, lf0, gam = p
        if not (0.7 < qs < 1.3 and abs(kap) < 0.9 and 0. < g < 3. and 2. < lrc < 8. and abs(eta) < 0.3):
            return -np.inf
        s = self.s; dx, dy = s.x-cx, s.y-cy; R2 = dx*dx+dy*dy; R = np.sqrt(R2)+1e-9
        m2 = _ell_m2(s.x, s.y, cx, cy, qs, ths)
        sig2 = np.exp(2*ls0)*(1+m2/np.exp(2*lrc))**(-g)*(1+eta*(s.mag-21.))**2
        w = R2/(R2+200.**2); sR2, sT2 = sig2*(1+kap*w), sig2*(1-kap*w)
        c, sn = dx/R, dy/R; f2 = np.exp(2*lf0+2*gam*self.D)
        a = sR2*c*c+sT2*sn*sn+f2*s.ex**2; d = sR2*sn*sn+sT2*c*c+f2*s.ey**2; b = (sR2-sT2)*c*sn
        det = a*d-b*b
        if np.any(det <= 0):
            return -np.inf
        q = (d*s.pmx**2-2*b*s.pmx*s.pmy+a*s.pmy**2)/det
        return float(-0.5*np.sum(q+np.log(det)))


def _start_counts(cx=0., cy=0.):
    return np.array([cx, cy, np.log(150.), 0.3, 1.0, 0.5, 0., 0., 0.])


def _start_disp(cx=0., cy=0.):
    return np.array([cx, cy, np.log(0.8), np.log(150.), 0.4, 1.0, 0.5, 0.0, 0.1, 0.0, 0.0])


@dataclass
class FitResult:
    params: dict
    centre: tuple
    centre_err: tuple
    loglike: float
    converged: bool
    message: str


def _hessian_centre_err(f, p, idx=(0, 1), h=0.05):
    """Errors on the centre from the numerical Hessian of -lnL over ALL parameters (central differences)."""
    n = len(p); H = np.zeros((n, n)); f0 = f(p)
    steps = np.full(n, 1e-3); steps[list(idx)] = h
    for i in range(n):
        for j in range(i, n):
            ei, ej = np.zeros(n), np.zeros(n); ei[i], ej[j] = steps[i], steps[j]
            H[i, j] = H[j, i] = (f(p+ei+ej)-f(p+ei-ej)-f(p-ei+ej)+f(p-ei-ej))/(4*steps[i]*steps[j])
    C = np.linalg.pinv(H, rcond=1e-10)                              # pseudo-inverse: unconstrained directions drop out
    return tuple(float(np.sqrt(C[i, i])) if C[i, i] > 0 else np.nan for i in idx)


def fit(counts: CountModel | None = None, disp: DispersionModel | None = None, start=(0., 0.), errors=True, maxiter=3000):
    """Maximise lnL_counts + lnL_disp with a shared centre. Returns a FitResult (centre error from the Hessian)."""
    pc = _start_counts(*start) if counts is not None else np.zeros(0)
    pd = _start_disp(*start)[2:] if disp is not None else np.zeros(0)
    nc = len(pc)
    fixed = set()
    if counts is not None:
        if not counts.lopsided:
            fixed |= {7, 8}
        if not np.any(counts.Ds) and not np.any(counts.Dg):        # flat depth map: beta_c is unconstrained
            fixed.add(6)
    free_c = [i for i in range(nc) if i not in fixed]
    disp_fixed = disp is not None and not np.any(disp.D)           # flat depth map: gamma is unconstrained
    if disp_fixed:
        pd = pd[:-1]
    z0 = np.r_[start, pc[[i for i in free_c if i > 1]] if counts is not None else [], pd]

    def unpack(z):
        cx, cy = z[0], z[1]; k = 2; out_c = out_d = None
        if counts is not None:
            q = pc.copy(); q[0], q[1] = cx, cy
            for i in [i for i in free_c if i > 1]:
                q[i] = z[k]; k += 1
            out_c = q
        if disp is not None:
            out_d = np.r_[cx, cy, z[k:], [0.] if disp_fixed else []]
        return out_c, out_d

    def nll(z):
        a, b = unpack(z); tot = 0.
        if a is not None:
            tot += counts.loglike(a)
        if b is not None:
            tot += disp.loglike(b)
        return -tot if np.isfinite(tot) else 1e30

    steps = np.full(len(z0), 0.1); steps[:2] = 1.0
    simplex = np.vstack([z0]+[z0+d for d in np.diag(steps)])
    r = minimize(nll, z0, method="Nelder-Mead", options=dict(maxiter=maxiter, xatol=1e-4, fatol=1e-4, initial_simplex=simplex, adaptive=True))
    r2 = minimize(nll, r.x, method="Powell", options=dict(maxiter=maxiter, xtol=1e-4, ftol=1e-7))
    best = r2 if r2.fun <= r.fun else r
    a, b = unpack(best.x); params = {}
    if a is not None:
        params.update({f"counts.{n}": float(v) for n, v in zip(CountModel.names, a)})
    if b is not None:
        params.update({f"disp.{n}": float(v) for n, v in zip(DispersionModel.names, b)})
    err = _hessian_centre_err(nll, best.x) if errors else (np.nan, np.nan)
    return FitResult(params, (float(best.x[0]), float(best.x[1])), err, float(-best.fun), bool(best.success),
                     f"NM: {r.message}; Powell: {r2.message}")


def symmetry_centre(x, y, r_in, r_out, half=6., nw=16):
    """Trial centre minimising the opposite-wedge chi^2 (8 pairs) in the annulus r_in-r_out; stars kept to r_out + half + 1
    so that no displaced annulus is truncated (the 2026-10-20 first version cut at r_out + 5). Coarse 0.5" grid over +-half,
    then 0.1" over +-0.6", parabolic refinement."""
    k = np.hypot(x, y) < r_out+half+1.; x, y = x[k], y[k]

    def grid(c0, h, step):
        g0 = np.arange(-h, h+step/2, step); gx, gy = c0[0]+g0, c0[1]+g0; chi = np.empty((len(g0), len(g0)))
        for i, cx in enumerate(gx):
            dx = x-cx
            for j, cy in enumerate(gy):
                dy = y-cy; R = np.hypot(dx, dy); kk = (R >= r_in) & (R < r_out)
                w = np.floor((np.arctan2(dy[kk], dx[kk])+np.pi)/(2*np.pi)*nw).astype(int) % nw
                n = np.bincount(w, minlength=nw); a, b = n[:nw//2], n[nw//2:]
                chi[i, j] = np.sum((a-b)**2/np.maximum(a+b, 1))
        i, j = np.unravel_index(np.argmin(chi), chi.shape)

        def para(v, kk, g):
            if 0 < kk < len(v)-1:
                d = v[kk-1]-2*v[kk]+v[kk+1]
                return g[kk]+(0.5*(v[kk-1]-v[kk+1])/d*step if d > 0 else 0.)
            return g[kk]
        edge = i in (0, len(g0)-1) or j in (0, len(g0)-1)
        return np.array([para(chi[:, j], i, gx), para(chi[i, :], j, gy)]), edge
    c, edge1 = grid(np.zeros(2), half, 0.5)
    c, edge2 = grid(c, 0.6, 0.1)
    return c, bool(edge1)
