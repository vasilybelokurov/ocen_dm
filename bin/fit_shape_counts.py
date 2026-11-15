#!/usr/bin/env python3
"""Unbinned fit of omega Cen's projected shape to star positions: flattened MGE with q' a smooth function of the
Gaussian width (docs/JAM_FIT_PLAN.md step 4c; 2026-11-15).

Model (per sample s, sky arcsec, x east, y north; position angle east of north):
  Sigma(x, y) = sum_k w_k / (2 pi s_k^2) exp(-(x'^2 + y'^2/q_k^2) / (2 sigma_maj,k^2)),  sigma_maj,k = s_k/sqrt(q_k),
  i.e. jam_mock_rotation_test.flatten (total of each Gaussian = w_k, circularised scale s_k), x' along the major axis
  at PA, y' along the minor axis, about a per-sample centre offset (dx_s, dy_s) (absorbs m = 1 / centring).
  q_k = q(log s_k): linear interpolation between NODES (log s), constant beyond the end nodes; q in [0.3, 1.5]
  (q > 1 = elongated along PA + 90; NO q > cos i limit -- consistency with i = 44 is checked afterwards).
  lambda_s = A_s Sigma (HST: A_HST free) or Sigma c_G(R) + F (1 + g_x x/1000 + g_y y/1000) (Gaia; uniform + gradient
  field, F >= 0; c_G the radial Gaia completeness of the count product). Weights w_k > 0 (log), free.
Likelihood: Poisson point process, ln L = sum_i ln lambda(x_i) - int lambda dA over each sample's annulus
  (polar quadrature, 320 log-R x 256 angles; full coverage: HST 60-250", Gaia 300-2400").
Fit: L-BFGS-B with JAX gradients from several starts; errors from the inverse Hessian.
Output: isophotal q'(R) (axis ratio of the model's isophotes, b/a, at major-axis radius R), q_k, PA, profile check.

Usage: PYTHONPATH=src:bin python bin/fit_shape_counts.py [--pa free|fixed] [--hst-mag 19] [--gaia-mag 17]
           [--gaia all|members] [--rmax-gaia 2400] [--quad-field] [--tag T]
Output: results/shape/shape_fit_<tag>.json, plots/shape_fit_<tag>.png
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
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


def q_of(qn, s):
    return jnp.interp(jnp.log(s), jnp.asarray(NODES), qn)


def sigma_xy(x, y, lw, qn, pa):
    """Model surface density (stars/arcsec^2 per unit total) at (x, y)."""
    a = jnp.radians(pa)
    xp = x*jnp.sin(a)+y*jnp.cos(a)                 # along the major axis (PA east of north)
    yp = -x*jnp.cos(a)+y*jnp.sin(a)
    q = q_of(qn, jnp.asarray(SIG))
    smaj2 = SIG**2/q
    w = jnp.exp(lw)
    return jnp.sum(w/(2*jnp.pi*SIG**2)*jnp.exp(-(xp[..., None]**2+(yp[..., None]/q)**2)/(2*smaj2)), axis=-1)


def polar_grid(rmin, rmax):
    e = np.geomspace(rmin, rmax, NR+1); r = np.sqrt(e[1:]*e[:-1])
    ph = (np.arange(NPHI)+.5)/NPHI*2*np.pi
    dA = (np.pi*(e[1:]**2-e[:-1]**2)/NPHI)[:, None]*np.ones(NPHI)[None, :]
    R, P = np.meshgrid(r, ph, indexing="ij")
    return (R*np.sin(P)).ravel(), (R*np.cos(P)).ravel(), dA.ravel(), R.ravel()


class Fit:
    def __init__(self, hst, gaia, pa_fixed=None, quad_field=False):
        self.hst, self.gaia, self.pa_fixed, self.quad = hst, gaia, pa_fixed, quad_field
        self.gh = polar_grid(hst["rmin"], hst["rmax"]); self.gg = polar_grid(gaia["rmin"], gaia["rmax"])
        c = json.loads((processed_dir()/"kinematics"/"ocen_counts_gaia_g17_c.json").read_text())
        rc = np.sqrt(np.array(c["r_lower"])*np.array(c["r_upper"]))
        self.cg = lambda R: jnp.interp(R, jnp.asarray(rc), jnp.asarray(c["completeness"]))
        self.names = [f"lw{int(s)}" for s in SIG]+[f"q{int(round(np.exp(n)))}" for n in NODES]+\
            (["pa"] if pa_fixed is None else [])+["lnA_hst", "dxh", "dyh", "dxg", "dyg", "lnF", "gx", "gy"]+\
            (["gxx", "gxy", "gyy"] if quad_field else [])
        self.n = len(self.names)
        self._f = jax.jit(self.neglnl); self._g = jax.jit(jax.grad(self.neglnl)); self._h = jax.jit(jax.hessian(self.neglnl))

    def unpack(self, p):
        d = dict(zip(self.names, [p[i] for i in range(self.n)]))
        lw = jnp.stack([d[f"lw{int(s)}"] for s in SIG]); qn = jnp.stack([d[k] for k in self.names if k.startswith("q")])
        pa = d.get("pa", self.pa_fixed)
        return d, lw, qn, pa

    def field(self, d, x, y):
        f = 1+d["gx"]*x/1000+d["gy"]*y/1000
        if self.quad:
            f = f+(d["gxx"]*x**2+d["gxy"]*x*y+d["gyy"]*y**2)/1e6
        return jnp.exp(d["lnF"])*jnp.clip(f, 1e-6, None)

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
            b.append((-5., 15.) if k.startswith("lw") else (0.3, 1.5) if k.startswith("q") else (40., 170.) if k == "pa"
                     else (-5., 5.) if k == "lnA_hst" else (-15., 15.) if k[0] == "d" else (-12., 2.) if k == "lnF"
                     else (-2., 2.))
        return b

    def start(self, rng=None):
        p = {k: 0. for k in self.names}
        R = np.hypot(self.gaia["x"], self.gaia["y"])
        for s in SIG:
            p[f"lw{int(s)}"] = np.log(len(self.gaia["x"])/len(SIG)*4)
        for k in self.names:
            if k.startswith("q"):
                p[k] = 0.88
        p["pa"] = 104.2; p["lnA_hst"] = np.log(10.); p["lnF"] = np.log(max(np.sum(R > 2000), 1)/(np.pi*(2400**2-2000**2))*0.5)
        x = np.array([p[k] for k in self.names])
        if rng is not None:
            x = x+rng.normal(0, 0.3, self.n)*np.array([1 if k.startswith("lw") else 0.1 if k.startswith("q")
                                                         else 10 if k == "pa" else 0 for k in self.names])
        lo, hi = np.array(self.bounds()).T
        return np.clip(x, lo+1e-6, hi-1e-6)

    def fit(self, nstart=4, seed=1):
        rng = np.random.default_rng(seed); best = None
        for i in range(nstart):
            x0 = self.start(None if i == 0 else rng)
            r = minimize(lambda p: float(self._f(jnp.asarray(p))), x0, jac=lambda p: np.asarray(self._g(jnp.asarray(p))),
                         method="L-BFGS-B", bounds=self.bounds(), options=dict(maxiter=5000, maxfun=20000))
            print(f"   start {i}: -lnL {r.fun:.3f} ({r.nit} it; {r.message})", flush=True)
            if best is None or r.fun < best.fun:
                best = r
        H = np.asarray(self._h(jnp.asarray(best.x)))
        try:
            cov = np.linalg.inv(H)
        except np.linalg.LinAlgError:
            cov = np.full((self.n, self.n), np.nan)
        return best, cov


def isophote_q(lw, qn, pa, R):
    """Axis ratio b/a of the model isophote through major-axis radius R (pa irrelevant: evaluated in the model frame)."""
    a = np.radians(pa); out = []
    for r in R:
        s0 = float(sigma_xy(jnp.array([r*np.sin(a)]), jnp.array([r*np.cos(a)]), lw, qn, pa)[0])
        lo, hi = 0.2*r, 1.6*r
        for _ in range(60):
            m = .5*(lo+hi)
            v = float(sigma_xy(jnp.array([-m*np.cos(a)]), jnp.array([m*np.sin(a)]), lw, qn, pa)[0])
            lo, hi = (m, hi) if v > s0 else (lo, m)
        out.append(.5*(lo+hi)/r)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pa", choices=("free", "fixed"), default="free")
    ap.add_argument("--hst-mag", type=float, default=19.)
    ap.add_argument("--gaia-mag", type=float, default=17.)
    ap.add_argument("--gaia", choices=("members", "all"), default="members")
    ap.add_argument("--rmax-gaia", type=float, default=2400.)
    ap.add_argument("--quad-field", action="store_true")
    ap.add_argument("--nstart", type=int, default=4)
    ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    tag = a.tag or f"pa{a.pa}_h{a.hst_mag:g}_g{a.gaia_mag:g}{a.gaia}_r{a.rmax_gaia:g}" + ("_quad" if a.quad_field else "")
    hst = SC.load_hst(a.hst_mag); gaia = SC.load_gaia(a.gaia_mag, rmax=a.rmax_gaia, members=a.gaia == "members")
    print(f"{hst['name']}: {len(hst['x'])} stars; {gaia['name']}: {len(gaia['x'])} stars", flush=True)
    t0 = time.perf_counter()
    F = Fit(hst, gaia, pa_fixed=None if a.pa == "free" else SC.PA_KIN, quad_field=a.quad_field)
    best, cov = F.fit(a.nstart)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    d, lw, qn, pa = F.unpack(jnp.asarray(best.x))
    pa = float(pa)
    Rq = np.geomspace(70, 2300, 40)
    qiso = isophote_q(lw, qn, pa, Rq)
    # isophote errors: linear propagation (delta method) of the inverse-Hessian covariance, J = d q_iso / d p by
    # central differences with steps 0.1 sigma_p (random draws from an ill-conditioned covariance are not reliable)
    J = np.zeros((len(Rq), F.n))
    for i in range(F.n):
        h = 0.1*err[i] if err[i] > 0 else 1e-4
        for sgn in (1, -1):
            x = best.x.copy(); x[i] += sgn*h
            dd, l2, q2, p2 = F.unpack(jnp.asarray(x))
            J[:, i] += sgn*isophote_q(l2, q2, float(p2), Rq)/(2*h)
    qerr = np.sqrt(np.clip(np.einsum("ri,ij,rj->r", J, cov, J), 0, None))
    qk = np.asarray(q_of(qn, jnp.asarray(SIG)))
    # JAM export (van de Ven et al. 2006 prescription): q'_k = isophotal axis ratio at R = s_k; smooth, unlike the
    # per-Gaussian node values, which trade off between overlapping Gaussians
    qk_iso = np.interp(np.log(SIG), np.log(Rq), qiso)
    out = dict(tag=tag, samples=[hst["name"], gaia["name"]], n_stars=[len(hst["x"]), len(gaia["x"])], neglnl=float(best.fun),
               names=F.names, x=best.x.tolist(), err=err.tolist(), pa=pa, pa_err=float(err[F.names.index("pa")]) if "pa" in F.names else 0.,
               sigma_circ_arcsec=SIG.tolist(), q_k=qk.tolist(), q_k_iso=qk_iso.tolist(), q_nodes_width=np.exp(NODES).tolist(),
               q_nodes=np.asarray(qn).tolist(), q_nodes_err=[float(err[F.names.index(k)]) for k in F.names if k.startswith("q")],
               weights=np.exp(np.asarray(lw)).tolist(), R_iso=Rq.tolist(), q_iso=qiso.tolist(), q_iso_err=qerr.tolist(),
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
