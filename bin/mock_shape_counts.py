#!/usr/bin/env python3
"""Fake-cluster (mock) tests of the unbinned shape fit (bin/fit_shape_counts.py; docs/JAM_FIT_PLAN.md step 4d; 2026-11-15).

Each mock draws star positions from a known model through the SAME annuli, Gaia radial completeness c_G(R) and field
model as the real fit, with the real star numbers (Poisson), then refits with fit_shape_counts.Fit and compares the
isophotal q'(R) with the truth's isophotal q'(R) (same isophote routine).

Truth = a real fit result (default the main fit) modified per test:
  fit     : the real best fit itself (flattened, q' 0.81-1.0)
  round   : all q' = 1
  step    : q'_k = 0.82 for s_k <= 590", 1.0 for s_k >= 830" (sharper than the fitted node model can follow)
  steep   : real fit, outer two Gaussians (1160", 1620") x 0.5 (steeper outer profile)
  offset  : real fit, true centre at (+5", -3") in both samples
  m1      : real fit, HST density x (1 + 0.015 cos(phi - 290 deg)) (lopsided, not a pure centre shift)
  curved  : truth from the Gaia-all-stars fit plus a quadratic field term (gxx 0.15, gyy -0.1, gxy 0.1 per 1e6 arcsec^2);
            fitted with the linear field (as the default fit)
Not done here: PM-selection injection (needs a PM model of the selection efficiency; open).

Pass (plan): at each R, |mean(q_fit - q_true)| < 0.01 and the pull (q_fit - q_true)/err has s.d. ~1 (within 1 sigma).

Usage: PYTHONPATH=src:bin python bin/mock_shape_counts.py --test round --n 20 [--seed0 0] [--nstart 2]
Output: results/shape/mocks/mock_<test>.json (one record per mock; appended, so runs can be split by seed)
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import fit_shape_counts as FS
import jax.numpy as jnp
import numpy as np

import shape_counts as SC

ROOT = Path(__file__).resolve().parents[1]
RQ = np.geomspace(70, 2300, 40)
IDX = [int(np.argmin(abs(RQ-r))) for r in (100, 420, 860, 1760)]
TESTS = ("fit", "round", "step", "steep", "offset", "m1", "curved")
REAL = {"members": "pafree_h19_g17members_r2400", "all": "pafree_h19_g17all_r2400"}


def sigma_k(x, y, lw, qk, pa, c=(0., 0.)):
    """Flattened MGE with per-Gaussian q'_k (numpy; same convention as fit_shape_counts.sigma_xy)."""
    a = np.radians(pa); x = x-c[0]; y = y-c[1]
    xp = x*np.sin(a)+y*np.cos(a); yp = -x*np.cos(a)+y*np.sin(a)
    s = FS.SIG; smaj2 = s**2/qk
    return np.sum(np.exp(lw)/(2*np.pi*s**2)*np.exp(-(xp[..., None]**2+(yp[..., None]/qk)**2)/(2*smaj2)), axis=-1)


def iso_q(lw, qk, pa, R):
    """Isophotal b/a at major-axis radius R for a per-Gaussian model (bisection, as fit_shape_counts.isophote_q)."""
    a = np.radians(pa); out = []
    for r in R:
        s0 = sigma_k(np.array([r*np.sin(a)]), np.array([r*np.cos(a)]), lw, qk, pa)[0]
        lo, hi = 0.2*r, 1.6*r
        for _ in range(60):
            m = .5*(lo+hi)
            v = sigma_k(np.array([-m*np.cos(a)]), np.array([m*np.sin(a)]), lw, qk, pa)[0]
            lo, hi = (m, hi) if v > s0 else (lo, m)
        out.append(.5*(lo+hi)/r)
    return np.array(out)


def draw(lam, rmin, rmax, rng, nr=1200, nphi=720):
    """Poisson point process with intensity lam(x, y) over an annulus: cell sampling, uniform in area within cells."""
    e = np.geomspace(rmin, rmax, nr+1); ph = np.arange(nphi+1)/nphi*2*np.pi
    rc = np.sqrt(e[1:]*e[:-1]); pc = .5*(ph[1:]+ph[:-1])
    R, P = np.meshgrid(rc, pc, indexing="ij")
    dA = (np.pi*(e[1:]**2-e[:-1]**2)/nphi)[:, None]*np.ones(nphi)[None, :]
    mu = lam(R*np.sin(P), R*np.cos(P))*dA
    n = rng.poisson(mu.sum())
    cell = rng.choice(mu.size, n, p=(mu/mu.sum()).ravel())
    i, j = np.unravel_index(cell, mu.shape)
    r = np.sqrt(rng.uniform(e[i]**2, e[i+1]**2)); p = rng.uniform(ph[j], ph[j+1])
    return r*np.sin(p), r*np.cos(p)


def truth(test, F):
    """Truth parameters (dict) for a test, from a real fit result."""
    d = json.loads((ROOT/f"results/shape/shape_fit_{REAL['all' if test == 'curved' else 'members']}.json").read_text())
    p = dict(zip(d["names"], d["x"]))
    lw = np.array([p[f"lw{int(s)}"] for s in FS.SIG]); qk = np.array(d["q_k"]); pa = d["pa"]
    t = dict(lw=lw, qk=qk, pa=pa, A=np.exp(p["lnA_hst"]), ch=(p["dxh"], p["dyh"]), cg=(p["dxg"], p["dyg"]),
             F=np.exp(p["lnF"]), gx=p["gx"], gy=p["gy"], quad=(0., 0., 0.), m1=0.)
    if test == "round":
        t["qk"] = np.ones_like(qk)
    elif test == "step":
        t["qk"] = np.where(FS.SIG <= 590, 0.82, 1.0)
    elif test == "steep":
        t["lw"] = lw+np.log(np.where(FS.SIG >= 1160, 0.5, 1.0))
    elif test == "offset":
        t["ch"] = t["cg"] = (5., -3.)
    elif test == "m1":
        t["m1"] = 0.015
    elif test == "curved":
        t["quad"] = (0.15, 0.1, -0.1)
    return t


def make_mock(t, hst, gaia, F, rng):
    def lam_h(x, y):
        s = t["A"]*sigma_k(x, y, t["lw"], t["qk"], t["pa"], t["ch"])
        return s*(1+t["m1"]*np.cos(np.arctan2(x, y)-np.radians(290.)))

    def lam_g(x, y):
        gxx, gxy, gyy = t["quad"]
        f = 1+t["gx"]*x/1000+t["gy"]*y/1000+(gxx*x**2+gxy*x*y+gyy*y**2)/1e6
        return sigma_k(x, y, t["lw"], t["qk"], t["pa"], t["cg"])*np.asarray(F.cg(np.hypot(x, y)))+t["F"]*np.clip(f, 1e-6, None)
    xh, yh = draw(lam_h, hst["rmin"], hst["rmax"], rng)
    xg, yg = draw(lam_g, gaia["rmin"], gaia["rmax"], rng)
    return dict(hst, x=jnp.asarray(xh), y=jnp.asarray(yh)), dict(gaia, x=jnp.asarray(xg), y=jnp.asarray(yg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", required=True, choices=TESTS)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--nstart", type=int, default=2)
    a = ap.parse_args()
    hst = SC.load_hst(19.); gaia = SC.load_gaia(17., rmax=2400., members=a.test != "curved")
    F0 = FS.Fit(hst, gaia)                                      # for c_G only
    t = truth(a.test, F0)
    q_true = iso_q(t["lw"], t["qk"], t["pa"], RQ)
    out = ROOT/f"results/shape/mocks/mock_{a.test}.json"; out.parent.mkdir(parents=True, exist_ok=True)
    recs = json.loads(out.read_text()) if out.exists() else []
    for k in range(a.seed0, a.seed0+a.n):
        if any(r["seed"] == k for r in recs):
            continue
        t0 = time.perf_counter()
        rng = np.random.default_rng([TESTS.index(a.test), k])
        h, g = make_mock(t, hst, gaia, F0, rng)
        F = FS.Fit(h, g)
        best, cov = F.fit(a.nstart, seed=k)
        q, e, pa, pae = FS.isophote_with_errors(F, best, cov, RQ)
        recs.append(dict(seed=k, n=[int(len(h["x"])), int(len(g["x"]))], neglnl=float(best.fun), pa=pa, pa_err=pae,
                         q_iso=q.tolist(), q_iso_err=e.tolist(), elapsed_s=round(time.perf_counter()-t0, 1)))
        out.write_text(json.dumps(recs, indent=0))
        meta = out.with_suffix(".truth.json")
        meta.write_text(json.dumps(dict(test=a.test, R=RQ.tolist(), q_true=q_true.tolist(), pa_true=float(t["pa"]))))
        print(f"{a.test} seed {k}: PA {pa:.1f}+-{pae:.1f}; dq at 100/420/860/1760\": "
              + " ".join(f"{q[i]-q_true[i]:+.3f}({e[i]:.3f})" for i in IDX) + f" ({recs[-1]['elapsed_s']:.0f} s)",
              flush=True)


if __name__ == "__main__":
    main()
