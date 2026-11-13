#!/usr/bin/env python3
"""Step 2 of the flattened-mock plan (2026-11-13): tune a flattened, rotating JAM model to omega Cen's 89 binned data.

Purpose: a realistic fake omega Cen (flattened, rotating, whose binned kinematics look like the real data) to test the
spherical Jeans pipeline (step 5). NOT a measurement of omega Cen's shape: q' and i are assumed per run.

Model: the K1 (no-DM) mass model -- stars (count MGE x M/L), Plummer remnants (NNLS MGE), point BH -- flattened by
jam_mock_rotation_test.jam_models (sigma_major = sigma/sqrt(q'), mass preserved; remnant q' separately), turnover
beta(r) as JAM's beta(r) (spherical-aligned ellipsoid, constant in theta), inclination i, distance D.
Free: M_star, M_rem, a_rem, M_bh, beta_0, beta_mid, beta_inf, r_beta, delta_r_beta, D (the K1 family's parameters and
prior ranges; D has its Gaussian prior).

Prediction per bin: the estimator rules recombine to the bin's weighted mean TOTAL second moment of the component
(HST: local mean removed + <vbar^2> added back; MUSE: ring sinusoid; Gaia: ring mean), independent of how rotation
splits <v_phi^2>. So predicted sigma = sqrt(<m2>_bin - streaming2), with streaming2 the per-estimator rotation term of
the REAL data, corrected (2026-11-13): HST R and T = ring <vbar^2> of the measured rotation field (MeanField, mas/yr at
the measurement distance), MUSE = published v_rot^2/2, Gaia = <mu_T>^2 (VB21, as in the snapshot). JAM second moments
are computed on a polar grid in (x', y') (36 log radii 2-2600" x 24 angles in [0, pi), point symmetry) and
interpolated to the sample points (real HST/MUSE star positions; Sigma-weighted Gaia rings). Likelihood: the
pipeline's split normal; objective = lnL + distance prior. Negative/NaN second moments -> rejected.

Usage: PYTHONPATH=src:bin python bin/tune_jam_to_data.py --qobs 0.90 [--qrem Q] [--inc 44] [--nang 31] [--tag T]
       --check: compare grid interpolation with direct evaluation at the start point, then exit.
Output: results/rotation/jam_tuned_<tag>.json (git-ignored).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

import jampy as jam
import jam_mock_rotation_test as J
from ocen_dm.kinematics import FitProblem, KinematicData
from ocen_dm.kinematics.likelihood import ProfileLikelihood
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from plot_jeans_centre2 import map_point

GR = np.geomspace(2., 2600., 36)
GPHI = (np.arange(24)+.5)/24*np.pi
D_MEAS = 5.43            # distance at which the PM rotation field was converted to km/s


class Setup:
    def __init__(self, nang=31, nlos=150, per_bin=250):
        run = load_run_metadata(ROOT/"results/fits/jeans_K1_turnover_counts_c")
        self.pb0 = problem_for(run)
        self.fam = self.pb0.family
        self.x_map, _, _ = map_point("jeans_K1_turnover_counts_c", self.fam)
        self.nang, self.nlos = nang, nlos
        rng = np.random.default_rng(1)
        pos = J.star_positions(self.pb0.data, rng, per_bin=per_bin)
        a = np.radians(J.PA)
        self.pts = {}
        for p in self.pb0.data.profiles:
            x, y, b, nd = pos[p.name]
            xp, yp = x*np.sin(a)+y*np.cos(a), -x*np.cos(a)+y*np.sin(a)
            self.pts[p.name] = (xp, yp, b, nd)
        self.data = self._corrected_data()

    def _corrected_data(self):
        """Real data with the corrected rotation terms (HST R and T from the measured rotation field)."""
        mf = J.MeanField(1.0)
        out = []
        for p in self.pb0.data.profiles:
            if p.instrument == "HST":
                xp, yp, b, _ = self.pts[p.name]
                v = mf(xp, yp)
                R = np.hypot(xp, yp)
                e = (np.array([xp/R, yp/R]) if p.kind == "pmr" else np.array([-yp/R, xp/R]))
                m1 = (e[0]*v[0]+e[1]*v[1])/(4.74047*D_MEAS)
                st2 = np.array([np.mean(m1[b == k]**2) for k in range(p.n)])
                out.append(replace(p, streaming2=st2))
            else:
                out.append(p)
        return KinematicData(tuple(out))

    def grid_moments(self, lum, pot, mbh, D, anis, inc):
        xg = (GR[:, None]*np.cos(GPHI)[None, :]).ravel(); yg = (GR[:, None]*np.sin(GPHI)[None, :]).ravel()
        beta = lambda r, th: (anis.beta(r), 0.)
        o = jam.axi.proj(lum, pot, inc, mbh, D/1e3, xg, yg, beta=beta, ml=1, moment="zz", interp=False,
                         plot=False, quiet=True, nang=self.nang, nlos=self.nlos)
        v2 = o.vel2
        comps = {k: v2[i, j].reshape(len(GR), len(GPHI)) for k, (i, j) in
                 dict(xx=(0, 0), yy=(1, 1), xy=(0, 1), zz=(2, 2)).items()}
        return comps

    def interp(self, comps, xp, yp):
        """Second-moment components at points; point symmetry (x,y)->(-x,-y) maps phi into [0, pi)."""
        r = np.clip(np.hypot(xp, yp), GR[0], GR[-1]); ph = np.mod(np.arctan2(yp, xp), np.pi)
        phig = np.r_[GPHI[-1]-np.pi, GPHI, GPHI[0]+np.pi]
        q = np.column_stack([np.log(r), ph])
        out = {}
        for k, g in comps.items():
            gg = np.c_[g[:, -1:], g, g[:, :1]]
            out[k] = RegularGridInterpolator((np.log(GR), phig), gg)(q)
        return out

    def predict(self, theta, qobs, qrem, inc, direct=False):
        jeans, D, _ = self.fam.build(theta)
        lum, pot, mbh = J.jam_models(jeans, D, qobs, qrem)
        comps = None if direct else self.grid_moments(lum, pot, mbh, D, jeans.anisotropy, inc)
        pred = {}
        for p in self.data.profiles:
            xp, yp, b, nd = self.pts[p.name]
            if direct:
                beta = lambda r, th: (jeans.anisotropy.beta(r), 0.)
                o = jam.axi.proj(lum, pot, inc, mbh, D/1e3, xp, yp, beta=beta, ml=1, interp=False, plot=False,
                                 quiet=True, nang=41, nlos=400)
                m = dict(xx=o.vel2[0, 0], yy=o.vel2[1, 1], xy=o.vel2[0, 1], zz=o.vel2[2, 2])
            else:
                m = self.interp(comps, xp, yp)
            R = np.hypot(xp, yp)
            if p.kind == "los":
                m2 = m["zz"]
            else:
                e = np.array([xp/R, yp/R]) if p.kind == "pmr" else np.array([-yp/R, xp/R])
                m2 = e[0]**2*m["xx"]+2*e[0]*e[1]*m["xy"]+e[1]**2*m["yy"]
            w = J.weights(lum, xp, yp, nd)
            tot = np.array([np.sum(w[b == k]*m2[b == k])/np.sum(w[b == k]) for k in range(p.n)])
            if p.kind != "los":
                tot = tot/(4.74047*D)**2
            if np.any(~np.isfinite(tot)) or np.any(tot <= 0):
                raise ValueError("non-positive second moment")
            st = 0 if p.streaming2 is None else p.streaming2
            if np.any(tot-st <= 0):
                raise ValueError("streaming exceeds second moment")
            pred[p.name] = np.sqrt(tot-st)
        return pred, D

    def lnlike(self, pred):
        return {p.name: float(np.sum(ProfileLikelihood._split_normal_lnlike(pred[p.name], p)))
                for p in self.data.profiles}

    def stats(self, pred):
        out = {}
        for p in self.data.profiles:
            err = np.where(pred[p.name] > p.value, p.err_hi, p.err_lo)
            d = (pred[p.name]-p.value)/err
            out[p.name] = dict(mean=float(d.mean()), rms=float(np.sqrt(np.mean(d**2))), chi2=float(np.sum(d**2)))
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qobs", type=float, required=True)
    ap.add_argument("--qrem", type=float, default=None)
    ap.add_argument("--inc", type=float, default=J.INC)
    ap.add_argument("--nang", type=int, default=31)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--maxfev", type=int, default=4000)
    args = ap.parse_args()
    qrem = args.qobs if args.qrem is None else args.qrem
    tag = args.tag or f"q{args.qobs:.2f}_r{qrem:.2f}_i{args.inc:g}"
    S = Setup(nang=args.nang)
    fam = S.fam
    names = list(fam.names)
    prior = {p.name: p.prior for p in fam.parameters}
    lo = np.array([prior[n].transform(1e-6) for n in names]); hi = np.array([prior[n].transform(1-1e-6) for n in names])
    logp = {"M_star", "M_rem", "M_bh", "a_rem", "r_beta", "delta_r_beta"}
    to_u = lambda x: np.array([np.log10(v) if n in logp else v for n, v in zip(names, x)])
    from_u = lambda u: np.array([10**v if n in logp else v for n, v in zip(names, u)])
    dprior = prior["distance"]

    def objective(u, record=None):
        x = np.clip(from_u(u), lo, hi)
        try:
            pred, D = S.predict(dict(zip(names, x)), args.qobs, qrem, args.inc)
        except (ValueError, FloatingPointError, ZeroDivisionError):
            return 1e30
        ln = sum(S.lnlike(pred).values())-0.5*((D-dprior.mu)/dprior.sigma)**2
        pen = np.sum(np.clip(from_u(u)-hi, 0, None)**2+np.clip(lo-from_u(u), 0, None)**2)*1e6
        return -ln+pen

    x0 = S.x_map.copy()
    t0 = time.perf_counter()
    pred0, _ = S.predict(dict(zip(names, x0)), args.qobs, qrem, args.inc)
    print(f"{tag}: one evaluation {time.perf_counter()-t0:.2f} s; start (spherical MAP) per-dataset (pred-data)/err:")
    for k, v in S.stats(pred0).items():
        print(f"   {k:28s} mean {v['mean']:+6.2f} rms {v['rms']:5.2f}")
    if args.check:
        pd, _ = S.predict(dict(zip(names, x0)), args.qobs, qrem, args.inc, direct=True)
        for p in S.data.profiles:
            print(f"   grid vs direct {p.name:28s} max |rel| {np.max(np.abs(pred0[p.name]/pd[p.name]-1)):.4f}, "
                  f"max |diff|/err {np.max(np.abs(pred0[p.name]-pd[p.name])/p.err_lo):.3f}")
        return
    starts = [x0]
    x1 = x0.copy(); x1[names.index("beta_inf")] = -0.4; starts.append(x1)
    best = None
    for s in starts:
        r = minimize(objective, to_u(s), method="Powell", options=dict(maxfev=args.maxfev, xtol=1e-4, ftol=1e-7))
        print(f"   start beta_inf {s[names.index('beta_inf')]:+.2f}: -lnpost {r.fun:.3f} ({r.nfev} calls, "
              f"{time.perf_counter()-t0:.0f} s)", flush=True)
        if best is None or r.fun < best.fun:
            best = r
    xb = np.clip(from_u(best.x), lo, hi)
    pred, D = S.predict(dict(zip(names, xb)), args.qobs, qrem, args.inc)
    st = S.stats(pred)
    out = dict(tag=tag, created=time.strftime("%Y-%m-%dT%H:%M:%S"), qobs=args.qobs, qrem=qrem, inc=args.inc,
               nang=args.nang, names=names, x=xb.tolist(), neg_lnpost=float(best.fun),
               lnL=S.lnlike(pred), stats=st, pred={k: v.tolist() for k, v in pred.items()},
               spherical_map=dict(zip(names, map(float, x0))))
    (ROOT/"results/rotation").mkdir(parents=True, exist_ok=True)
    (ROOT/f"results/rotation/jam_tuned_{tag}.json").write_text(json.dumps(out, indent=1))
    print(f"{tag}: best -lnpost {best.fun:.3f}; parameters (tuned | spherical MAP):")
    for n, v, v0 in zip(names, xb, x0):
        print(f"   {n:14s} {v:12.5g} | {v0:12.5g}")
    for k, v in st.items():
        print(f"   {k:28s} mean {v['mean']:+6.2f} rms {v['rms']:5.2f} chi2 {v['chi2']:.1f}")


if __name__ == "__main__":
    main()
