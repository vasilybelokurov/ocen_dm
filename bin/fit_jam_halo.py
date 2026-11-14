#!/usr/bin/env python3
"""Flattened JAM fit with an optional NFW halo and the flattening q' as a free parameter (2026-11-14).

Model: tune_jam_to_data.Setup (all stars, precomputed bin operators, corrected rotation terms and errors) with
jam_halo_models.jam_models (fast remnant MGE; round truncated-NFW halo). Families: K1 (no DM) or NFW (the production
K2_nfw parameters M_dm_100, r_s with their priors). q' (stars = remnants) is free in [0.75, 1.0] unless --qobs is
given; inclination fixed (--inc, default 44). Objective: ln L + distance prior (as the tuner); Powell from several
starts (2, run in parallel: 2 cores per fit), each restarted until the gain is < 0.05; the best is kept.
Default --method de (2026-11-14, Powell got stuck far from the known truth on the mocks): differential evolution over
the prior box, population evaluated on --workers processes, then a Powell polish.
--fix-mdm M fixes M_dm_100 = M (profile likelihood point; r_s stays free).

Data: --data real (the corrected real data) or --data <jam_tuned tag> (noiseless mock: the saved tuned prediction),
or --data <tag> with --inject-mdm M --inject-rs R: the tuned truth plus an NFW halo, predicted with this model
(mock not tuned to the data; the truth then has M_dm_100 = M, r_s = R).

Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/fit_jam_halo.py --data q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1 \
           --family NFW [--fix-mdm 3e5] [--qobs 0.9] [--maxfev 3000]
Output: results/rotation/jamfit_<data>[_inj..]_<family>[_q..][_mdm..].json
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
from scipy.optimize import minimize

import jam_halo_models as H
from ocen_dm.kinematics import KinematicData
from tune_jam_to_data import Setup

RUNS = {"K1": "jeans_K1_turnover_counts_c", "NFW": "jeans_K2_nfw_turnover_counts_c"}
LOGP = {"M_star", "M_rem", "M_bh", "a_rem", "r_beta", "delta_r_beta", "M_dm_100", "r_s"}
Q_RANGE = (0.75, 1.0)


class Problem:
    """Data, free parameters and objective for one configuration (built once per process)."""

    def __init__(self, a):
        self.a = a
        self.S = S = Setup(run=RUNS[a.family], models=H.jam_models)
        names = list(S.fam.names)
        prior = {p.name: p.prior for p in S.fam.parameters}
        self.dprior = prior["distance"]
        self.truth, self.tag = None, a.data
        if a.data != "real":
            T = json.loads((ROOT/f"results/rotation/jam_tuned_{a.data}.json").read_text())
            self.truth = dict(zip(T["names"], T["x"]), qobs=T["qobs"], qrem=T["qrem"], inc=T["inc"])
            if a.inject_mdm is None:
                pred = {k: np.array(v) for k, v in T["pred"].items()}
            else:
                SN = S if a.family == "NFW" else Setup(run=RUNS["NFW"], models=H.jam_models)
                self.truth.update(M_dm_100=a.inject_mdm, r_s=a.inject_rs)
                pred, _ = SN.predict({n: self.truth[n] for n in SN.fam.names}, T["qobs"], T["qrem"], T["inc"])
                self.tag += f"_injM{a.inject_mdm:.3g}_rs{a.inject_rs:g}"
            S.data = KinematicData(tuple(replace(p, value=pred[p.name]) for p in S.data.profiles))
        self.qfree = a.qobs is None
        free = [n for n in names if not (n == "M_dm_100" and a.fix_mdm is not None)]
        lo = {n: prior[n].transform(1e-6) for n in names}; hi = {n: prior[n].transform(1-1e-6) for n in names}
        lo["qobs"], hi["qobs"] = Q_RANGE
        self.pars = free+(["qobs"] if self.qfree else [])
        self.lo = np.array([lo[n] for n in self.pars]); self.hi = np.array([hi[n] for n in self.pars])
        self.base = dict(zip(names, S.x_map))

    def to_u(self, x):
        return np.array([np.log10(v) if n in LOGP else v for n, v in zip(self.pars, x)])

    def from_u(self, u):
        return np.array([10**v if n in LOGP else v for n, v in zip(self.pars, u)])

    def unpack(self, x):
        d = dict(zip(self.pars, np.clip(x, self.lo, self.hi)))
        q = d.pop("qobs", self.a.qobs)
        if self.a.fix_mdm is not None:
            d["M_dm_100"] = self.a.fix_mdm
        return d, q

    def objective(self, u):
        x = self.from_u(u)
        th, q = self.unpack(x)
        try:
            pred, D = self.S.predict(th, q, q, self.a.inc)
        except (ValueError, FloatingPointError, ZeroDivisionError, AssertionError):
            return 1e30
        ln = sum(self.S.lnlike(pred).values())-0.5*((D-self.dprior.mu)/self.dprior.sigma)**2
        return -ln+np.sum(np.clip(x-self.hi, 0, None)**2+np.clip(self.lo-x, 0, None)**2)*1e6

    def starts(self):
        """Spherical MAP of the family; q' 0.90 and 0.97 if free; halo 1e5 (with q' 0.90) and 1e6 (with 0.97)."""
        out = []
        qs = (0.90, 0.97) if self.qfree else (self.a.qobs, self.a.qobs)
        ms = (1e5, 1e6) if (self.a.family == "NFW" and self.a.fix_mdm is None) else (None, None)
        for q0, m0 in zip(qs, ms):
            s = dict(self.base, qobs=q0)
            if m0 is not None:
                s["M_dm_100"], s["r_s"] = m0, 300.
            if self.a.family == "NFW" and self.a.fix_mdm is not None:
                s["r_s"] = 300. if q0 == qs[0] else 100.
            out.append(np.clip(np.array([s[n] for n in self.pars]), self.lo, self.hi))
        return out


_P = None


def _init(a):
    global _P
    _P = Problem(a)


def _run_start(k):
    P, t0 = _P, time.perf_counter()
    s = P.starts()[k]
    u, fprev, nf, hist = P.to_u(s), P.objective(P.to_u(s)), 0, []
    for rd in range(P.a.max_rounds):
        r = minimize(P.objective, u, method="Powell", options=dict(maxfev=P.a.maxfev, xtol=1e-4, ftol=1e-7))
        u, nf = r.x, nf+r.nfev; hist.append(float(r.fun))
        print(f"   start {k} round {rd+1}: -lnpost {r.fun:.3f} ({r.nfev} calls, {time.perf_counter()-t0:.0f} s)",
              flush=True)
        if fprev-r.fun < 0.05:
            break
        fprev = r.fun
    return dict(start=s.tolist(), x=P.from_u(u).tolist(), f=float(r.fun), rounds=hist, nfev=nf)


def _obj(u):
    return _P.objective(u)


def _run_de(a, P):
    """Differential evolution in the optimiser coordinates (log10 for scale parameters) over the prior box, the
    population evaluated in parallel on a.workers processes; then a serial Powell polish of the best member."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    from scipy.optimize import differential_evolution
    t0 = time.perf_counter()
    bounds = list(zip(P.to_u(P.lo), P.to_u(P.hi)))
    hist = []

    def cb(intermediate_result):
        hist.append(float(intermediate_result.fun))
        if len(hist) % 10 == 0:
            print(f"   DE gen {len(hist)}: best -lnpost {hist[-1]:.3f} ({time.perf_counter()-t0:.0f} s)", flush=True)

    with ProcessPoolExecutor(a.workers, mp_context=get_context("spawn"), initializer=_init, initargs=(a,)) as ex:
        r = differential_evolution(_obj, bounds, popsize=a.popsize, maxiter=a.maxiter, tol=1e-8, atol=a.de_atol,
                                   mutation=(0.5, 1.0), recombination=0.7, seed=a.seed, init="sobol",
                                   updating="deferred", workers=ex.map, polish=False, callback=cb)
    print(f"   DE done: -lnpost {r.fun:.3f} after {r.nit} generations, {r.nfev} calls "
          f"({time.perf_counter()-t0:.0f} s); {r.message}", flush=True)
    p = minimize(P.objective, r.x, method="Powell", options=dict(maxfev=a.maxfev, xtol=1e-4, ftol=1e-8))
    print(f"   Powell polish: -lnpost {p.fun:.3f} ({p.nfev} calls, {time.perf_counter()-t0:.0f} s)", flush=True)
    best = p if p.fun < r.fun else r
    return [dict(method="differential_evolution", x=P.from_u(best.x).tolist(), f=float(best.fun), de_f=float(r.fun),
                 polish_f=float(p.fun), generations=int(r.nit), nfev=int(r.nfev+p.nfev), history=hist,
                 popsize=a.popsize, seed=a.seed)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--family", choices=sorted(RUNS), required=True)
    ap.add_argument("--qobs", type=float, default=None, help="fix q' (default: free)")
    ap.add_argument("--inc", type=float, default=44.0)
    ap.add_argument("--fix-mdm", type=float, default=None)
    ap.add_argument("--inject-mdm", type=float, default=None)
    ap.add_argument("--inject-rs", type=float, default=300.)
    ap.add_argument("--maxfev", type=int, default=3000)
    ap.add_argument("--max-rounds", type=int, default=5)
    ap.add_argument("--method", choices=("de", "powell"), default="de")
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--popsize", type=int, default=15, help="DE population = popsize x n_parameters")
    ap.add_argument("--maxiter", type=int, default=400, help="DE generations")
    ap.add_argument("--de-atol", type=float, default=0.02, help="DE stops when the population's -lnpost s.d. < this")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    t0 = time.perf_counter()
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    P = Problem(a)
    if a.method == "de":
        runs = _run_de(a, P)
    else:
        with ProcessPoolExecutor(2, mp_context=get_context("spawn"), initializer=_init, initargs=(a,)) as ex:
            runs = list(ex.map(_run_start, range(2)))           # the two starts run in parallel (2 cores)
    best = min(runs, key=lambda r: r["f"])
    th, q = P.unpack(np.array(best["x"]))
    pred, D = P.S.predict(th, q, q, a.inc)
    st, lnL = P.S.stats(pred), P.S.lnlike(pred)
    suffix = (f"_q{a.qobs:.2f}" if not P.qfree else "")+(f"_mdm{a.fix_mdm:.3g}" if a.fix_mdm is not None else "")
    out = dict(data=P.tag, family=a.family, inc=a.inc, q_free=P.qfree, fix_mdm=a.fix_mdm, theta=th, qobs=q, D=D,
               neg_lnpost=best["f"], lnL_total=float(sum(lnL.values())), lnL=lnL, stats=st,
               chi2_total=float(sum(v["chi2"] for v in st.values())), pred={k: v.tolist() for k, v in pred.items()},
               data_values={p.name: p.value.tolist() for p in P.S.data.profiles}, starts=runs, truth=P.truth,
               created=time.strftime("%Y-%m-%dT%H:%M:%S"), elapsed_s=round(time.perf_counter()-t0, 1))
    path = ROOT/f"results/rotation/jamfit_{P.tag}_{a.family}{suffix}.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"{path.name}: -lnpost {best['f']:.3f} (starts: {[round(r['f'], 2) for r in runs]}); chi2 "
          f"{out['chi2_total']:.1f}; q' {q:.3f}; M_dm_100 {th.get('M_dm_100', 0):.3g}; M_bh {th['M_bh']:.0f}; "
          f"D {D:.3f}; ({out['elapsed_s']/60:.0f} min)", flush=True)


if __name__ == "__main__":
    main()
