#!/usr/bin/env python3
"""Mock calibration of the joint centre fit (2026-10-20; bin/fit_centre.py, src/ocen_dm/centre_fit.py).

Truth = the real-data joint best fit (results/plot_data/centre_fit.json, 'JOINT ...'), with its centre moved by an injected
offset (0, 0), (+2, 0), (0, +2)" (x = dRA* east +).
  counts: N(F625W 19-22, R < 200") positions drawn from lambda = C(x) Sigma(x) of the truth (rejection on the 1" window
          grid, uniform jitter inside the cell), so the depth-tied completeness pattern is reproduced;
  PMs:    the real hq stars' positions, magnitudes and quoted errors; velocities drawn from the truth's anisotropic
          dispersion model with the truth's error scale f(D);
  stress set ('depth-misscaled'): the drawn PM noise uses f(D) x (1 + 0.15 sign(D)) -- a depth-correlated error mis-scale
          that the fitted form f0 exp(gamma D) can only partly absorb.
Each mock is refitted jointly from the adopted origin. Reported per set: mean bias, scatter, mean Hessian error, and the
pull distribution (fit - truth)/err (should have mean 0, sd 1).
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src python bin/centre_mocks.py --n 15 --workers 4 [--stress]
Output: results/plot_data/centre_mocks.json
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
OFFSETS = ((0., 0.), (2., 0.), (0., 2.))


def _truth():
    j = json.loads((ROOT/"results/plot_data/centre_fit.json").read_text())["fits"]
    return [v for k, v in j.items() if k.startswith("JOINT")][0]["params"]


def one(args):
    seed, off, stress = args
    from ocen_dm.centre_fit import CountModel, DepthMap, DispersionModel, Stars, _ell_m2, fit, load_omegacat
    rng = np.random.default_rng(seed); P = _truth()
    s = load_omegacat(); depth = DepthMap(s)
    pc = np.array([P[f"counts.{n}"] for n in CountModel.names]); pd = np.array([P[f"disp.{n}"] for n in DispersionModel.names])
    cx, cy = pc[0]+off[0], pc[1]+off[1]; pc[:2] = pd[:2] = (cx, cy)
    # counts mock
    mid = (s.mag >= 19) & (s.mag < 22) & (np.hypot(s.x, s.y) < 200.); n = int(mid.sum())
    cm = CountModel(s.subset(mid), depth, 200.); lam = cm._lam(cm.gx, cm.gy, cm.Dg, pc); pr = lam/lam.sum()
    idx = rng.choice(len(pr), n, p=pr); mx = cm.gx[idx]+rng.uniform(-0.5, 0.5, n); my = cm.gy[idx]+rng.uniform(-0.5, 0.5, n)
    k = np.hypot(mx, my) < 200.; sc = Stars(mx[k], my[k], np.full(k.sum(), 20.5))
    # PM mock on the real hq stars
    dmr = DispersionModel(s, depth, 200.); h = dmr.s
    ls0, lrc, g, qs, ths, eta, kap, lf0, gam = pd[2:]
    dx, dy = h.x-cx, h.y-cy; R2 = dx*dx+dy*dy; R = np.sqrt(R2)+1e-9; m2 = _ell_m2(h.x, h.y, cx, cy, qs, ths)
    sig2 = np.exp(2*ls0)*(1+m2/np.exp(2*lrc))**(-g)*(1+eta*(h.mag-21.))**2; w = R2/(R2+200.**2)
    vR, vT = rng.normal(0, np.sqrt(sig2*(1+kap*w))), rng.normal(0, np.sqrt(sig2*(1-kap*w)))
    c, sn = dx/R, dy/R; D = depth(h.x, h.y); f = np.exp(lf0+gam*D)*((1+0.15*np.sign(D)) if stress else 1.)
    vx = vR*c-vT*sn+rng.normal(0, f*h.ex); vy = vR*sn+vT*c+rng.normal(0, f*h.ey)
    sd = Stars(h.x, h.y, h.mag, vx, vy, h.ex, h.ey, np.ones(len(h.x), bool))
    r = fit(counts=CountModel(sc, depth, 200.), disp=DispersionModel(sd, depth, 200.))
    return dict(seed=seed, offset=off, stress=stress, truth=(cx, cy), centre=r.centre, err=r.centre_err, converged=r.converged)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=15); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--stress", action="store_true"); a = ap.parse_args()
    jobs = [(1000*i+j+(50000 if a.stress else 0), off, a.stress) for i, off in enumerate(OFFSETS) for j in range(a.n)]
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(one, jobs))
    out = dict(runs=res, summary={})
    for off in OFFSETS:
        rr = [r for r in res if tuple(r["offset"]) == off]
        d = np.array([np.subtract(r["centre"], r["truth"]) for r in rr]); e = np.array([r["err"] for r in rr]); pull = d/e
        out["summary"][str(off)] = dict(n=len(rr), bias=d.mean(0).tolist(), scatter=d.std(0, ddof=1).tolist(), mean_err=np.nanmean(e, 0).tolist(),
                                       pull_mean=np.nanmean(pull, 0).tolist(), pull_sd=np.nanstd(pull, 0, ddof=1).tolist(),
                                       converged=int(sum(r["converged"] for r in rr)))
        sm = out["summary"][str(off)]
        print(f"offset {off}: bias ({sm['bias'][0]:+.2f}, {sm['bias'][1]:+.2f}) scatter ({sm['scatter'][0]:.2f}, {sm['scatter'][1]:.2f}) "
              f"Hessian err ({sm['mean_err'][0]:.2f}, {sm['mean_err'][1]:.2f}) pull sd ({sm['pull_sd'][0]:.2f}, {sm['pull_sd'][1]:.2f})", flush=True)
    tag = "_stress" if a.stress else ""
    (ROOT/f"results/plot_data/centre_mocks{tag}.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
