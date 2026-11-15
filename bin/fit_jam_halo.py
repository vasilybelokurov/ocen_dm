#!/usr/bin/env python3
"""Flattened JAM fit with an optional NFW halo and the flattening q' as a free parameter (2026-11-14).

Model: tune_jam_to_data.Setup (all stars, precomputed bin operators, corrected rotation terms and errors) with
jam_halo_models.jam_models (fast remnant MGE; round truncated-NFW halo). Families: K1 (no DM) or NFW (the production
K2_nfw parameters M_dm_100, r_s with their priors). q' (stars = remnants) is free in [0.75, 1.0] unless --qobs is
given; inclination fixed (--inc, default 44). Objective: ln L + distance prior (as the tuner).
Default --method lsq (2026-11-14): bounded least squares (capfit) on the whitened residuals from several starts in
parallel on --workers processes. Starts: the family's standard starts, any --start-from fits, and the best points of a
Sobol screen of the prior box, picked best-first with a minimum separation (--min-sep, normalised box); every start
must give a valid model. Each run is restarted from its end point (fresh Jacobian) until the gain is < --lsq-gain
(at most --lsq-rounds). Other methods: --method de (differential evolution + Powell polish), --method powell.
--fix-mdm M fixes M_dm_100 = M (profile likelihood point; r_s stays free).
--rmax R fits only the bins with radius <= R arcsec (the model is still evaluated on all bins; 2026-11-15).

Data: --data real (the corrected real data) or --data <jam_tuned tag> (noiseless mock: the saved tuned prediction),
or --data <tag> with --inject-mdm M --inject-rs R: the tuned truth plus an NFW halo, predicted with this model
(mock not tuned to the data; the truth then has M_dm_100 = M, r_s = R).

Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/fit_jam_halo.py --data q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1 \
           --family NFW [--fix-mdm 3e5] [--qobs 0.9] [--maxfev 3000]
Output: results/rotation/jamfit_<data>[_inj..]_<family>[_q..][_mdm..][_lsq][_<label>].json; refuses to overwrite an
existing file unless --overwrite. Saves the arguments, the data (values, errors, rotation terms) and the git commit.
"""
from __future__ import annotations

import argparse
import json
import subprocess
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
from ocen_dm.kinematics.likelihood import ProfileLikelihood
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
                self.truth.update(M_dm_100=a.inject_mdm, r_s=a.inject_rs)
                pred = getattr(a, "inject_pred", None)            # computed once in the parent (workers reuse it)
                if pred is None:
                    SN = S if a.family == "NFW" else Setup(run=RUNS["NFW"], models=H.jam_models)
                    pred, _ = SN.predict({n: self.truth[n] for n in SN.fam.names}, T["qobs"], T["qrem"], T["inc"])
                    a.inject_pred = pred
                self.tag += f"_injM{a.inject_mdm:.3g}_rs{a.inject_rs:g}"
            S.data = KinematicData(tuple(replace(p, value=pred[p.name]) for p in S.data.profiles))
        self.qfree = a.qobs is None
        free = [n for n in names if not (n == "M_dm_100" and a.fix_mdm is not None)]
        lo = {n: prior[n].transform(1e-6) for n in names}; hi = {n: prior[n].transform(1-1e-6) for n in names}
        lo["qobs"], hi["qobs"] = Q_RANGE
        self.pars = free+(["qobs"] if self.qfree else [])
        self.lo = np.array([lo[n] for n in self.pars]); self.hi = np.array([hi[n] for n in self.pars])
        self.base = dict(zip(names, S.x_map))
        self.nfail = 0                                            # failed model evaluations (this process)
        rmax = np.inf if getattr(a, "rmax", None) is None else a.rmax
        self.keep = {p.name: np.asarray(p.r) <= rmax for p in S.data.profiles}   # bins entering the fit

    def lnlike(self, pred):
        """Split-normal lnL per dataset over the fitted bins."""
        return {p.name: float(np.sum(ProfileLikelihood._split_normal_lnlike(pred[p.name], p)[self.keep[p.name]]))
                for p in self.S.data.profiles}

    def stats(self, pred):
        out = {}
        for p in self.S.data.profiles:
            k = self.keep[p.name]
            d = ((pred[p.name]-p.value)/np.where(pred[p.name] > p.value, p.err_hi, p.err_lo))[k]
            out[p.name] = dict(mean=float(d.mean()), rms=float(np.sqrt(np.mean(d**2))), chi2=float(np.sum(d**2)),
                               n=int(k.sum()))
        return out

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
            self.nfail += 1
            return 1e30
        ln = sum(self.lnlike(pred).values())-0.5*((D-self.dprior.mu)/self.dprior.sigma)**2
        return -ln+np.sum(np.clip(x-self.hi, 0, None)**2+np.clip(self.lo-x, 0, None)**2)*1e6

    def residuals(self, u):
        """Whitened residuals: (model - data)/err (split-normal side) for the 89 bins, plus the distance-prior term.
        0.5 sum(r^2) = -lnpost + const (the split-normal normalisation is shared by both sides). u inside the box
        (bounded solvers); failed models return large finite residuals."""
        th, q = self.unpack(self.from_u(u))
        try:
            pred, D = self.S.predict(th, q, q, self.a.inc)
        except (ValueError, FloatingPointError, ZeroDivisionError, AssertionError):
            self.nfail += 1
            return np.full(self.nres, 1e3)
        r = [((pred[p.name]-p.value)/np.where(pred[p.name] > p.value, p.err_hi, p.err_lo))[self.keep[p.name]]
             for p in self.S.data.profiles]
        return np.r_[np.concatenate(r), (D-self.dprior.mu)/self.dprior.sigma]

    @property
    def nres(self):
        return sum(int(k.sum()) for k in self.keep.values())+1

    @property
    def lnnorm(self):
        """Constant: -lnpost = 0.5 sum(r^2) - lnnorm."""
        return float(sum(np.sum(np.log(np.sqrt(2/np.pi)/(p.err_lo+p.err_hi))[self.keep[p.name]])
                         for p in self.S.data.profiles))

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


def _screen(u):
    return _P.objective(u)


def _lsq_once(P, u0, lo, hi):
    a = P.a
    step = np.full(len(u0), a.lsq_step)
    if a.lsq_solver == "capfit":
        from capfit.capfit import capfit
        r = capfit(P.residuals, u0, bounds=(lo, hi), abs_step=step, x_scale="jac", ftol=a.lsq_ftol, xtol=1e-8,
                   max_nfev=a.lsq_maxfev)
    else:
        from scipy.optimize import least_squares
        r = least_squares(P.residuals, u0, bounds=(lo, hi), method="trf", x_scale="jac", diff_step=None,
                          ftol=a.lsq_ftol, xtol=1e-8, gtol=1e-10, max_nfev=a.lsq_maxfev,
                          jac=lambda u: _fd_jac(P.residuals, u, step, lo, hi))
    return r.x, int(r.nfev), int(r.status)


def _lsq_start(job):
    """Bounded least squares (capfit or scipy TRF) from u0 in the optimiser coordinates, restarted from its end point
    (fresh Jacobian) until the gain in -lnpost is < a.lsq_gain or a.lsq_rounds runs are done."""
    k, u0 = job
    P, a, t0 = _P, _P.a, time.perf_counter()
    lo, hi = P.to_u(P.lo), P.to_u(P.hi)
    u0 = np.clip(u0, lo+1e-6*(hi-lo), hi-1e-6*(hi-lo))
    nfail0, f0 = P.nfail, P.objective(u0)
    u, fprev, nf, rounds = u0, f0, 0, []
    for _ in range(a.lsq_rounds):
        u, n, status = _lsq_once(P, u, lo, hi)
        f = P.objective(u); nf += n
        rounds.append(dict(f=float(f), status=status, nfev=n))
        if fprev-f < a.lsq_gain:
            break
        fprev = f
    nfail = P.nfail-nfail0
    print(f"   start {k}: -lnpost {f0:.2f} -> {' -> '.join(f'{r['f']:.3f}' for r in rounds)} ({nf} calls + jacobians, "
          f"status {status}, {nfail} failed models, {time.perf_counter()-t0:.0f} s)", flush=True)
    return dict(method=f"lsq_{a.lsq_solver}", start=P.from_u(u0).tolist(), x=P.from_u(u).tolist(), f=float(f),
                f_start=float(f0), status=status, nfev=nf, rounds=rounds, n_failed_models=nfail,
                valid=bool(f < 1e29), seconds=round(time.perf_counter()-t0, 1))


def _fd_jac(fun, u, step, lo, hi):
    """Forward differences, stepping inwards at the upper bound."""
    f0, cols = fun(u), []
    for i in range(len(u)):
        h = step[i] if u[i]+step[i] <= hi[i] else -step[i]
        v = u.copy(); v[i] += h
        cols.append((fun(v)-f0)/h)
    return np.array(cols).T


def _run_lsq(a, P):
    """Screen a Sobol design over the box (plus the family's standard starts), then run bounded least squares from
    the --lsq-starts best valid points, in parallel on a.workers processes."""
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    from scipy.stats import qmc
    t0 = time.perf_counter()
    lo, hi = P.to_u(P.lo), P.to_u(P.hi)
    design = qmc.scale(qmc.Sobol(len(lo), seed=a.seed).random(a.screen), lo, hi) if a.screen else np.empty((0, len(lo)))
    std = np.array([P.to_u(s) for s in P.starts()])
    extra = [P.to_u(np.array([dict(zip(*_xs(j)))[n] for n in P.pars])) for j in (a.start_from or [])]
    cand = np.vstack([std, *([np.array(extra)] if extra else []), design])
    with ProcessPoolExecutor(a.workers, mp_context=get_context("spawn"), initializer=_init, initargs=(a,)) as ex:
        f = np.array(list(ex.map(_screen, cand, chunksize=4)))
        ngiven = len(std)+len(extra)
        bad = [i for i in range(ngiven) if f[i] >= 1e29]
        print(f"   screened {len(cand)} points: {np.sum(f < 1e29)} valid, best -lnpost {f.min():.2f}; "
              f"invalid standard/given starts dropped: {bad} ({time.perf_counter()-t0:.0f} s)", flush=True)
        keep, mind = _pick_starts(f, (cand-lo)/(hi-lo), ngiven, a.lsq_starts, a.min_sep)
        print(f"   {len(keep)} starts; min separation (normalised box) {mind:.2f}", flush=True)
        runs = list(ex.map(_lsq_start, [(k, cand[i]) for k, i in enumerate(keep)]))
    fs = sorted(r["f"] for r in runs)
    print(f"   least squares done: best -lnpost {fs[0]:.3f}; all {[round(v, 2) for v in fs]} "
          f"({time.perf_counter()-t0:.0f} s)", flush=True)
    return runs


def _pick_starts(f, z, ngiven, n, min_sep):
    """Indices of up to n valid candidates: the given ones (first ngiven, if valid), then the rest best-first,
    accepted only at normalised distance >= min_sep from every accepted start; if too few, filled best-first.
    Returns (indices, smallest pairwise separation among them)."""
    keep = [i for i in range(ngiven) if f[i] < 1e29]
    order = [i for i in np.argsort(f) if i >= ngiven and f[i] < 1e29]
    for i in order:
        if len(keep) >= n:
            break
        if not keep or np.min(np.linalg.norm(z[keep]-z[i], axis=1)) >= min_sep:
            keep.append(i)
    for i in order:
        if len(keep) >= n:
            break
        if i not in keep:
            keep.append(i)
    d = [np.linalg.norm(z[i]-z[j]) for k, i in enumerate(keep) for j in keep[k+1:]]
    return keep, (min(d) if d else np.nan)


def _git_commit():
    try:
        h = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "bin", "src"], cwd=ROOT, capture_output=True,
                               text=True).stdout.strip()
        return h.strip()+("+dirty" if dirty else "")
    except OSError:
        return None


def _xs(path):
    """(names, values) of a previous jamfit JSON's best point (theta plus qobs)."""
    d = json.loads(Path(path).read_text())
    th = dict(d["theta"], qobs=d["qobs"])
    return list(th), list(th.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--family", choices=sorted(RUNS), required=True)
    ap.add_argument("--qobs", type=float, default=None, help="fix q' (default: free)")
    ap.add_argument("--inc", type=float, default=44.0)
    ap.add_argument("--fix-mdm", type=float, default=None)
    ap.add_argument("--rmax", type=float, default=None, help="fit only bins with R <= rmax (arcsec)")
    ap.add_argument("--inject-mdm", type=float, default=None)
    ap.add_argument("--inject-rs", type=float, default=300.)
    ap.add_argument("--maxfev", type=int, default=3000)
    ap.add_argument("--max-rounds", type=int, default=5)
    ap.add_argument("--method", choices=("lsq", "de", "powell"), default="lsq")
    ap.add_argument("--lsq-solver", choices=("capfit", "trf"), default="capfit")
    ap.add_argument("--lsq-starts", type=int, default=6, help="number of least-squares starts (incl. standard ones)")
    ap.add_argument("--screen", type=int, default=128, help="Sobol points screened to pick starts")
    ap.add_argument("--start-from", nargs="*", default=None, help="previous jamfit JSONs used as extra starts")
    ap.add_argument("--lsq-step", type=float, default=1e-3, help="finite-difference step in optimiser coordinates")
    ap.add_argument("--lsq-ftol", type=float, default=1e-6)
    ap.add_argument("--lsq-maxfev", type=int, default=None)
    ap.add_argument("--lsq-rounds", type=int, default=3, help="max restarts of each least-squares run")
    ap.add_argument("--lsq-gain", type=float, default=0.02, help="restart while -lnpost improves by more than this")
    ap.add_argument("--min-sep", type=float, default=0.25, help="min start separation in the normalised box")
    ap.add_argument("--label", default="", help="suffix for the output file name")
    ap.add_argument("--overwrite", action="store_true", help="allow overwriting an existing output file")
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--popsize", type=int, default=15, help="DE population = popsize x n_parameters")
    ap.add_argument("--maxiter", type=int, default=400, help="DE generations")
    ap.add_argument("--de-atol", type=float, default=0.02, help="DE stops when the population's -lnpost s.d. < this")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    args = dict(vars(a))
    t0 = time.perf_counter()
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    P = Problem(a)
    suffix = (f"_q{a.qobs:.2f}" if not P.qfree else "")+(f"_mdm{a.fix_mdm:.3g}" if a.fix_mdm is not None else "") \
        + (f"_rmax{a.rmax:g}" if a.rmax is not None else "")+("_lsq" if a.method == "lsq" else "") \
        + (f"_{a.label}" if a.label else "")
    path = ROOT/f"results/rotation/jamfit_{P.tag}_{a.family}{suffix}.json"
    if path.exists() and not a.overwrite:
        sys.exit(f"{path.name} exists; use --label or --overwrite")
    if a.method == "lsq":
        runs = _run_lsq(a, P)
    elif a.method == "de":
        runs = _run_de(a, P)
    else:
        with ProcessPoolExecutor(2, mp_context=get_context("spawn"), initializer=_init, initargs=(a,)) as ex:
            runs = list(ex.map(_run_start, range(2)))           # the two starts run in parallel (2 cores)
    best = min(runs, key=lambda r: r["f"])
    th, q = P.unpack(np.array(best["x"]))
    pred, D = P.S.predict(th, q, q, a.inc)
    st, lnL = P.stats(pred), P.lnlike(pred)
    out = dict(data=P.tag, family=a.family, inc=a.inc, q_free=P.qfree, fix_mdm=a.fix_mdm, theta=th, qobs=q, D=D,
               neg_lnpost=best["f"], lnL_total=float(sum(lnL.values())), lnL=lnL, stats=st,
               chi2_total=float(sum(v["chi2"] for v in st.values())), pred={k: v.tolist() for k, v in pred.items()},
               data_values={p.name: p.value.tolist() for p in P.S.data.profiles}, starts=runs, truth=P.truth,
               data_used={p.name: dict(r=p.r.tolist(), value=p.value.tolist(), err_lo=np.asarray(p.err_lo).tolist(),
                                  err_hi=np.asarray(p.err_hi).tolist(),
                                  streaming2=None if p.streaming2 is None else np.asarray(p.streaming2).tolist())
                     for p in P.S.data.profiles},
               rotation_terms="one field (fit_rotation_field cylindrical MAP) for HST, MUSE, Gaia (2026-11-15)",
               args=args, git=_git_commit(), lnnorm=P.lnnorm,
               created=time.strftime("%Y-%m-%dT%H:%M:%S"), elapsed_s=round(time.perf_counter()-t0, 1))
    path.write_text(json.dumps(out, indent=1))
    print(f"{path.name}: -lnpost {best['f']:.3f} (starts: {[round(r['f'], 2) for r in runs]}); chi2 "
          f"{out['chi2_total']:.1f}; q' {q:.3f}; M_dm_100 {th.get('M_dm_100', 0):.3g}; M_bh {th['M_bh']:.0f}; "
          f"D {D:.3f}; ({out['elapsed_s']/60:.0f} min)", flush=True)


if __name__ == "__main__":
    main()
