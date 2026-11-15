#!/usr/bin/env python3
"""Unbinned fit of omega Cen's projected shape to star positions: flattened MGE with q' a smooth function of the
Gaussian width (docs/JAM_FIT_PLAN.md step 4c; 2026-11-15).

Model (per sample s, sky arcsec, x east, y north; position angle east of north):
  Sigma(x, y) = sum_k w_k / (2 pi s_k^2) exp(-(x'^2 + y'^2/q_k^2) / (2 sigma_maj,k^2)),  sigma_maj,k = s_k/sqrt(q_k),
  i.e. jam_mock_rotation_test.flatten (total of each Gaussian = w_k, circularised scale s_k), x' along the major axis
  at PA, y' along the minor axis, about a per-sample centre offset (dx_s, dy_s) (absorbs m = 1 / centring).
  q_k = q(log s_k): linear interpolation between NODES (log s), constant beyond the end nodes; q in [0.3, 1.5]
  (q > 1 = elongated along PA + 90; NO q > cos i limit -- consistency with i = 44 is checked afterwards).
  lambda_s = A_s Sigma (HST: A_HST free) or Sigma c_G(R) + B(x, y) (Gaia; c_G the radial Gaia completeness of the count
  product). Weights w_k >= 0 (linear, so a Gaussian can switch off; log-weights crawled towards -inf), free.
  Field B = exp(lnF + g_x x/1000 + g_y y/1000), FIXED from measurements (measure_field; 2026-11-15, after Codex + Fable
  reviews): tilt (g_x, g_y) and, for the all-star sample, level from the PM-REJECTED stars at 1000-rmax" (real field
  stars at the same sky position); level for the PM-member sample from field stars in a PM control shell
  (chi2 9-36, parallax cut) at 1000-rmax", scaled by PM area (/3). --field free fits lnF, g_x, g_y (diagnostic only).
Likelihood: Poisson point process, ln L = sum_i ln lambda(x_i) - int lambda dA over each sample's annulus
  (polar quadrature, 320 log-R x 256 angles; full coverage: HST 60-250", Gaia 300-2400").
Fit: L-BFGS-B (ftol = 0: the default relative ftol stops when dF < ~5e-4 at |F| ~ 2e5, i.e. short) from several
  starts to get close, then <= 15 eigenvalue-modified Newton steps with the exact JAX Hessian (built from Hessian-vector
  products with lax.map: 1 GB instead of 5 GB) over the parameters not on a bound; convergence = predicted Newton gain
  and last dlnL < 1e-3. Errors: inverse Hessian (pseudo-inverse, flagged, if not positive definite); isophote errors
  from exact derivatives. JAX is held to one CPU thread with NPROC=1 (XLA flags do not work on macOS; verified).
Output: isophotal q'(R) (axis ratio of the model's isophotes, b/a, at major-axis radius R), q_k, PA, profile check.

Usage: PYTHONPATH=src:bin python bin/fit_shape_counts.py [--pa free|fixed] [--hst-mag 19] [--gaia-mag 17]
           [--gaia all|members] [--rmax-gaia 2400] [--field fixed|free] [--field-scale 1] [--tag T]
Output: results/shape/shape_fit_<tag>.json, plots/shape_fit_<tag>.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
os.environ.setdefault("NPROC", "1")              # one CPU thread for XLA (the only setting that works on macOS)
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
from scipy.optimize import minimize

import shape_counts as SC
from ocen_dm.paths import processed_dir

SIG = np.array([75., 105., 150., 210., 300., 420., 590., 830., 1160., 1620.])      # circularised widths (arcsec)
NODES = np.log(np.array([100., 250., 600., 1400.]))                                 # q' nodes in log width
NR, NPHI = 320, 256
WS = 1000.                                                                          # weight scale: w_k = WS v_k


def q_of(qn, s):
    return jnp.interp(jnp.log(s), jnp.asarray(NODES), qn)


def sigma_xy(x, y, w, qn, pa):
    """Model surface density (stars/arcsec^2) at (x, y); w = Gaussian totals (linear)."""
    a = jnp.radians(pa)
    xp = x*jnp.sin(a)+y*jnp.cos(a)                 # along the major axis (PA east of north)
    yp = -x*jnp.cos(a)+y*jnp.sin(a)
    q = q_of(qn, jnp.asarray(SIG))
    smaj2 = SIG**2/q
    return jnp.sum(w/(2*jnp.pi*SIG**2)*jnp.exp(-(xp[..., None]**2+(yp[..., None]/q)**2)/(2*smaj2)), axis=-1)


def polar_grid(rmin, rmax):
    e = np.geomspace(rmin, rmax, NR+1); r = np.sqrt(e[1:]*e[:-1])
    ph = (np.arange(NPHI)+.5)/NPHI*2*np.pi
    dA = (np.pi*(e[1:]**2-e[:-1]**2)/NPHI)[:, None]*np.ones(NPHI)[None, :]
    R, P = np.meshgrid(r, ph, indexing="ij")
    return (R*np.sin(P)).ravel(), (R*np.cos(P)).ravel(), dA.ravel(), R.ravel()


def measure_field(gaia_mag=17., rmax=2400., rmin=1000.):
    """Field model from the Gaia catalogue itself (no free field parameters in the shape fit).
    Tilt and all-star level: unbinned log-linear plane fitted to the PM-REJECTED stars at rmin-rmax" (field-dominated;
    cluster stars failing the cut are a few per cent there). Member level: stars passing the parallax cut with PM
    chi2 in [9, 36) (a PM control shell with 3x the PM area of the member window chi2 < 9; cluster leakage
    exp(-4.5) ~ 1%), divided by 3. Returns dict(all=(lnF, gx, gy), members=(lnF, gx, gy), n_rej, n_shell)."""
    from astropy.table import Table
    from scipy.optimize import minimize as _min
    t = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    x, y = SC.gnomonic(t["ra"], t["dec"]); R = np.hypot(x, y); g = np.asarray(t["g_mag"], float)
    k = np.isfinite(g) & (g < gaia_mag) & (R >= rmin) & (R < rmax)
    mem = SC.member_cut(t)
    dx = np.asarray(t["pmra"])-SC.MU_SYS[0]; dy = np.asarray(t["pmdec"])-SC.MU_SYS[1]
    chi2 = dx**2/(SC.S_INT**2+np.asarray(t["pmra_error"])**2)+dy**2/(SC.S_INT**2+np.asarray(t["pmdec_error"])**2)
    plx = np.abs(np.asarray(t["parallax"])-1/SC.OCEN_DISTANCE_KPC) < 3*np.asarray(t["parallax_error"])
    xr, yr = x[k & ~mem]/1000, y[k & ~mem]/1000
    e = np.geomspace(rmin, rmax, 201)/1000; ph = (np.arange(256)+.5)/256*2*np.pi; rc = np.sqrt(e[1:]*e[:-1])
    Rg, Pg = np.meshgrid(rc, ph, indexing="ij"); dA = (np.pi*(e[1:]**2-e[:-1]**2)/256)[:, None]*np.ones(256)
    X, Y = Rg*np.sin(Pg), Rg*np.cos(Pg)
    f = lambda p: -(np.sum(p[0]+p[1]*xr+p[2]*yr)-np.sum(np.exp(p[0]+p[1]*X+p[2]*Y)*dA))
    r = _min(f, [np.log(len(xr)/dA.sum()), 0., 0.], method="BFGS")
    lnF_all = r.x[0]-np.log(1e6)                                  # per arcsec^2 (fit was per (1000")^2)
    gx, gy = r.x[1], r.x[2]
    n_shell = int((k & plx & (chi2 >= 9) & (chi2 < 36)).sum())
    # member level: shell count / 3, spread with the same tilt over the same annulus
    norm = np.sum(np.exp(gx*X+gy*Y)*dA)*1e6
    lnF_mem = np.log(n_shell/3/norm)
    err = np.sqrt(np.diag(r.hess_inv))
    return dict(all=(float(lnF_all), float(gx), float(gy)), members=(float(lnF_mem), float(gx), float(gy)),
                tilt_err=[float(err[1]), float(err[2])], n_rej=int(len(xr)), n_shell=n_shell, rmin=rmin, rmax=rmax)


class Fit:
    def __init__(self, hst, gaia, pa_fixed=None, field=None):
        """field: (lnF, gx, gy) fixed field (per arcsec^2, per 1000"), or None = free (diagnostic)."""
        self.hst, self.gaia, self.pa_fixed = hst, gaia, pa_fixed
        self.field_fixed = None if field is None else tuple(float(v) for v in field)
        self.gh = polar_grid(hst["rmin"], hst["rmax"]); self.gg = polar_grid(gaia["rmin"], gaia["rmax"])
        c = json.loads((processed_dir()/"kinematics"/"ocen_counts_gaia_g17_c.json").read_text())
        rc = np.sqrt(np.array(c["r_lower"])*np.array(c["r_upper"]))
        self.cg = lambda R: jnp.interp(R, jnp.asarray(rc), jnp.asarray(c["completeness"]))
        self.names = [f"w{int(s)}" for s in SIG]+[f"q{int(round(np.exp(n)))}" for n in NODES]+\
            (["pa"] if pa_fixed is None else [])+["lnA_hst", "dxh", "dyh", "dxg", "dyg"]+\
            (["lnF", "gx", "gy"] if field is None else [])
        self.n = len(self.names)
        self._f = jax.jit(self.neglnl); self._g = jax.jit(jax.grad(self.neglnl))

        def hess(p):                                    # Hessian from HVPs, one column at a time (memory ~1 GB)
            hvp = lambda v: jax.jvp(jax.grad(self.neglnl), (p,), (v,))[1]
            return jax.lax.map(hvp, jnp.eye(self.n))
        self._h = jax.jit(hess)

    def unpack(self, p):
        d = dict(zip(self.names, [p[i] for i in range(self.n)]))
        lw = WS*jnp.stack([d[f"w{int(s)}"] for s in SIG]); qn = jnp.stack([d[k] for k in self.names if k.startswith("q")])
        pa = d.get("pa", self.pa_fixed)
        return d, lw, qn, pa

    def field(self, d, x, y):
        lnF, gx, gy = self.field_fixed if self.field_fixed is not None else (d["lnF"], d["gx"], d["gy"])
        return jnp.exp(lnF+gx*x/1000+gy*y/1000)

    def neglnl(self, p):
        d, lw, qn, pa = self.unpack(p)
        h, g = self.hst, self.gaia
        A = jnp.exp(d["lnA_hst"])
        lh = A*sigma_xy(h["x"]-d["dxh"], h["y"]-d["dyh"], lw, qn, pa)
        xg, yg, dAg, Rg = self.gg; xh, yh, dAh, _ = self.gh
        ih = jnp.sum(A*sigma_xy(xh-d["dxh"], yh-d["dyh"], lw, qn, pa)*dAh)
        Rs = jnp.hypot(g["x"], g["y"])
        lg = sigma_xy(g["x"]-d["dxg"], g["y"]-d["dyg"], lw, qn, pa)*self.cg(Rs)+self.field(d, g["x"], g["y"])
        ig = jnp.sum((sigma_xy(xg-d["dxg"], yg-d["dyg"], lw, qn, pa)*self.cg(Rg)+self.field(d, xg, yg))*dAg)
        return -(jnp.sum(jnp.log(lh))-ih+jnp.sum(jnp.log(lg))-ig)

    def bounds(self):
        b = []
        for k in self.names:
            b.append((0., 100.) if k.startswith("w") else (0.3, 1.5) if k.startswith("q") else (40., 170.) if k == "pa"
                     else (-5., 5.) if k == "lnA_hst" else (-15., 15.) if k[0] == "d" else (-25., 2.) if k == "lnF"
                     else (-2., 2.))                    # gx, gy (free field only)
        return b

    def start(self, rng=None):
        p = {k: 0. for k in self.names}
        R = np.hypot(self.gaia["x"], self.gaia["y"])
        for s in SIG:
            p[f"w{int(s)}"] = len(self.gaia["x"])/len(SIG)*4/WS
        for k in self.names:
            if k.startswith("q"):
                p[k] = 0.88
        p["pa"] = 104.2; p["lnA_hst"] = np.log(10.); p["lnF"] = np.log(max(np.sum(R > 2000), 1)/(np.pi*(2400**2-2000**2))*0.5)
        x = np.array([p[k] for k in self.names])
        if rng is not None:
            x = x*np.exp(rng.normal(0, 0.3, self.n)*np.array([k.startswith("w") for k in self.names]))+\
                rng.normal(0, 0.3, self.n)*np.array([0.1 if k.startswith("q") else 10 if k == "pa" else 0 for k in self.names])
        lo, hi = np.array(self.bounds()).T
        return np.clip(x, lo+1e-6, hi-1e-6)

    def fit(self, nstart=3, seed=1, newton_max=15):
        """L-BFGS-B (ftol = 0) from nstart starts, then eigenvalue-modified Newton on the best (see module doc)."""
        lo, hi = np.array(self.bounds()).T
        f = lambda p: float(self._f(jnp.asarray(p))); g = lambda p: np.asarray(self._g(jnp.asarray(p)))
        rng = np.random.default_rng(seed); best = None
        iw = np.array([k.startswith("w") for k in self.names])
        # stage A: weights in log space (well scaled from a poor start; linear weights stalled after 1 iteration from
        # 2 of 3 starts); stage B: linear weights from there (a weight can reach 0), ftol = 0
        def z2x(z):
            x = z.copy(); x[iw] = np.exp(z[iw]); return x
        loA, hiA = lo.copy(), hi.copy(); loA[iw], hiA[iw] = -12., np.log(hi[iw])
        fA = lambda z: f(z2x(z))
        gA = lambda z: g(z2x(z))*np.where(iw, z2x(z), 1.)
        for i in range(nstart):
            x0 = self.start(None if i == 0 else rng); z0 = x0.copy(); z0[iw] = np.log(x0[iw])
            rA = minimize(fA, z0, jac=gA, method="L-BFGS-B", bounds=list(zip(loA, hiA)), options=dict(maxiter=5000))
            r = minimize(f, z2x(rA.x), jac=g, method="L-BFGS-B", bounds=list(zip(lo, hi)),
                         options=dict(maxiter=20000, maxfun=40000, ftol=0., gtol=1e-8, maxcor=30))
            print(f"   start {i}: -lnL {rA.fun:.3f} (log-w, {rA.nit} it) -> {r.fun:.3f} (linear w, {r.nit} it)", flush=True)
            if best is None or r.fun < best.fun:
                best = r
        x, fx = best.x.copy(), float(best.fun); f_lbfgs = fx
        for it in range(newton_max):
            free = (x > lo+1e-6*(hi-lo)) & (x < hi-1e-6*(hi-lo)) | (g(x)*np.sign(np.where(x <= lo+1e-6*(hi-lo), -1, 1)) > 0)
            H = np.asarray(self._h(jnp.asarray(x)))[np.ix_(free, free)]; gr = g(x)[free]
            ev, V = np.linalg.eigh(H); evm = np.maximum(np.abs(ev), 1e-8*np.abs(ev).max())
            step = -(V/evm)@(V.T@gr); gain = 0.5*float(gr@(V/evm)@(V.T@gr))
            t, ok = 1., False
            for _ in range(30):
                xn = x.copy(); xn[free] = np.clip(x[free]+t*step, lo[free], hi[free]); fn = f(xn)
                if fn < fx:
                    ok = True; break
                t *= 0.5
            df = fx-fn if ok else 0.
            if ok:
                x, fx = xn, fn
            if gain < 1e-3 and df < 1e-3:
                break
        self.conv = dict(f_lbfgs=f_lbfgs, dlnL_newton=f_lbfgs-fx, newton_it=it+1, pred_gain=gain, last_dlnL=df,
                         converged=bool(gain < 1e-3 and df < 1e-3))
        print(f"   Newton: d(-lnL) {fx-f_lbfgs:+.4f} in {it+1} it; predicted gain {gain:.2g}; converged {self.conv['converged']}",
              flush=True)
        best.x, best.fun = x, fx
        return best, self.covariance(x)

    def covariance(self, x):
        """Inverse Hessian over the parameters not on a bound (bound parameters: zero variance, listed). If the Hessian
        is not positive definite, pseudo-inverse over the positive eigen-directions (anti-conservative: flagged)."""
        lo, hi = np.array(self.bounds()).T
        free = (x > lo+1e-6*(hi-lo)) & (x < hi-1e-6*(hi-lo))
        H = np.asarray(self._h(jnp.asarray(x)))[np.ix_(free, free)]
        ev, V = np.linalg.eigh(H)
        cov = np.zeros((self.n, self.n))
        if ev.min() > 0:
            cf = np.linalg.inv(H)
        else:
            keep = ev > 1e-9*ev.max(); cf = (V[:, keep]/ev[keep])@V[:, keep].T
        cov[np.ix_(free, free)] = cf
        self.hess_min_eig = float(ev.min()); self.hess_pd = bool(ev.min() > 0)
        self.at_bound = [k for k, fr in zip(self.names, free) if not fr]
        if self.at_bound or not self.hess_pd:
            print(f"   Hessian: positive definite {self.hess_pd} (min eig {ev.min():.3g}); at bound: {self.at_bound}", flush=True)
        return cov


def isophote_q(w, qn, pa, R):
    """Axis ratio b/a of the model isophote through major-axis radius R (pa irrelevant: evaluated in the model frame)."""
    a = np.radians(pa); out = []
    for r in R:
        s0 = float(sigma_xy(jnp.array([r*np.sin(a)]), jnp.array([r*np.cos(a)]), w, qn, pa)[0])
        lo, hi = 0.2*r, 1.6*r
        for _ in range(60):
            m = .5*(lo+hi)
            v = float(sigma_xy(jnp.array([-m*np.cos(a)]), jnp.array([m*np.sin(a)]), w, qn, pa)[0])
            lo, hi = (m, hi) if v > s0 else (lo, m)
        out.append(.5*(lo+hi)/r)
    return np.array(out)


def isophote_with_errors(F, best, cov, Rq):
    """Isophotal q'(Rq) with delta-method errors from EXACT derivatives (implicit function theorem; JAX): the minor-axis
    radius b of the isophote through major-axis radius R solves S(b e_min; p) = S(R e_maj; p), so
    db/dp = (dS_maj/dp - dS_min/dp) / (dS_min/db) and dq/dp = (db/dp)/R. (Finite differences with steps 0.1 sigma_p were
    unreliable along the near-flat weight directions.) Returns q, err, pa, pa_err."""
    x = jnp.asarray(best.x)
    d, lw, qn, pa = F.unpack(x)
    qiso = isophote_q(lw, qn, float(pa), Rq)

    def S(p, r, minor):
        dd, l2, q2, p2 = F.unpack(p); a = jnp.radians(p2)
        xs, ys = (-r*jnp.cos(a), r*jnp.sin(a)) if minor else (r*jnp.sin(a), r*jnp.cos(a))
        return sigma_xy(jnp.atleast_1d(xs), jnp.atleast_1d(ys), l2, q2, p2)[0]
    gS = jax.jit(jax.grad(S), static_argnums=2); dSdr = jax.jit(jax.grad(S, argnums=1), static_argnums=2)
    J = np.zeros((len(Rq), F.n))
    for i, (r, q) in enumerate(zip(Rq, qiso)):
        b = q*r
        J[i] = (np.asarray(gS(x, r, False))-np.asarray(gS(x, b, True)))/float(dSdr(x, b, True))/r
    qerr = np.sqrt(np.clip(np.einsum("ri,ij,rj->r", J, cov, J), 0, None))
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    pae = float(err[F.names.index("pa")]) if "pa" in F.names else 0.
    return qiso, qerr, float(pa), pae


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pa", choices=("free", "fixed"), default="free")
    ap.add_argument("--hst-mag", type=float, default=19.)
    ap.add_argument("--gaia-mag", type=float, default=17.)
    ap.add_argument("--gaia", choices=("members", "all"), default="members")
    ap.add_argument("--rmax-gaia", type=float, default=2400.)
    ap.add_argument("--field", choices=("fixed", "free"), default="fixed")
    ap.add_argument("--field-scale", type=float, default=1., help="multiply the fixed field level (sensitivity test)")
    ap.add_argument("--nstart", type=int, default=3)
    ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    tag = a.tag or (f"pa{a.pa}_h{a.hst_mag:g}_g{a.gaia_mag:g}{a.gaia}_r{a.rmax_gaia:g}" + ("_fieldfree" if a.field == "free" else "")
                    + (f"_fx{a.field_scale:g}" if a.field_scale != 1 else ""))
    hst = SC.load_hst(a.hst_mag); gaia = SC.load_gaia(a.gaia_mag, rmax=a.rmax_gaia, members=a.gaia == "members")
    print(f"{hst['name']}: {len(hst['x'])} stars; {gaia['name']}: {len(gaia['x'])} stars", flush=True)
    t0 = time.perf_counter()
    fm = measure_field(a.gaia_mag, a.rmax_gaia)
    fld = fm[a.gaia]; fld = (fld[0]+np.log(a.field_scale), fld[1], fld[2])
    print(f"field ({a.gaia}): lnF {fld[0]:.3f} /arcsec^2, gx {fld[1]:+.3f} gy {fld[2]:+.3f} per 1000\" "
          f"(rejected stars {fm['n_rej']}, PM-shell stars {fm['n_shell']})" + ("; FREE in the fit" if a.field == "free" else ""),
          flush=True)
    F = Fit(hst, gaia, pa_fixed=None if a.pa == "free" else SC.PA_KIN, field=None if a.field == "free" else fld)
    best, cov = F.fit(a.nstart)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    d, lw, qn, pa = F.unpack(jnp.asarray(best.x))
    Rq = np.geomspace(70, 2300, 40)
    qiso, qerr, pa, _ = isophote_with_errors(F, best, cov, Rq)
    qk = np.asarray(q_of(qn, jnp.asarray(SIG)))
    # JAM export (van de Ven et al. 2006 prescription): q'_k = isophotal axis ratio at R = s_k; smooth, unlike the
    # per-Gaussian node values, which trade off between overlapping Gaussians
    qk_iso = np.interp(np.log(SIG), np.log(Rq), qiso)
    out = dict(tag=tag, jax=jax.__version__, samples=[hst["name"], gaia["name"]], n_stars=[len(hst["x"]), len(gaia["x"])], neglnl=float(best.fun),
               names=F.names, x=best.x.tolist(), err=err.tolist(), pa=pa, pa_err=float(err[F.names.index("pa")]) if "pa" in F.names else 0.,
               sigma_circ_arcsec=SIG.tolist(), q_k=qk.tolist(), q_k_iso=qk_iso.tolist(), q_nodes_width=np.exp(NODES).tolist(),
               q_nodes=np.asarray(qn).tolist(), q_nodes_err=[float(err[F.names.index(k)]) for k in F.names if k.startswith("q")],
               weights=np.asarray(lw).tolist(), field=dict(fixed=a.field == "fixed", lnF_gx_gy=list(fld), measurement=fm),
               conv=F.conv, hess_pd=F.hess_pd, at_bound=F.at_bound, R_iso=Rq.tolist(), q_iso=qiso.tolist(), q_iso_err=qerr.tolist(),
               elapsed_s=round(time.perf_counter()-t0, 1))
    (ROOT/"results/shape").mkdir(parents=True, exist_ok=True)
    (ROOT/f"results/shape/shape_fit_{tag}.json").write_text(json.dumps(out, indent=1))
    print(f"{tag}: -lnL {best.fun:.2f}; PA {pa:.1f}+-{out['pa_err']:.1f}; q nodes "
          + ", ".join(f"{w:.0f}\": {q:.3f}+-{e:.3f}" for w, q, e in zip(out['q_nodes_width'], out['q_nodes'], out['q_nodes_err']))
          + f"; ({out['elapsed_s']:.0f} s)")
    for r, q, e in zip(Rq[::4], qiso[::4], qerr[::4]):
        print(f"   R {r:6.0f}\"  q_iso {q:.3f} +- {e:.3f}")


if __name__ == "__main__":
    main()
