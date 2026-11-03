#!/usr/bin/env python3
"""Model-side centre marginalisation of the saved fixed-BH DF fits (src/ocen_dm/centre_model.py).

Data frozen about AvdM10 (the runs' own snapshot); the model centre c moves. For each centre of
the Gauss-Hermite grids the saved fits (no BH, 8,200, 4e4 Msun; not refitted) are scored with
DFJointProblem.evaluate(factors=CentreGeometry.factors(...)) and

    ln L_m = logsumexp_k [ln w_k + ln L(theta, c_k | d)].

Tests: T2 (c = 0 with unit factors = frozen loglike_profiled, <= 1e-9), T4 (n = 10 and 14 vs 6:
|change of -2 Delta ln L_m| < 0.1), ESS of the centre contributions per model (> 10), Gaia PM
bins (not moved): change of their prediction at |c| = 6" < 1e-2 of the bin error.
Central count deficit: prior-predictive P(N_inner2 <= n_obs), n_obs fixed (AvdM10 bins).

Usage: python bin/centre_model_marginal.py [--label 20261103]   (1 core; mostly the 3 model builds)
Output: results/plot_data/centre_model_marginal_<label>.json, plots/centre_model_marginal_<label>.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"bin"))

import numpy as np
from scipy.stats import poisson

from centre_marginal import RUNS, load_run, scores
from ocen_dm import centre_marg as cm
from ocen_dm.centre_model import CentreGeometry, ModelMoments, PointSet, effective_sample_size, projected_variance
from ocen_dm.kinematics.df_fit import DFJointProblem, build_df_model
from ocen_dm.kinematics.likelihood import ARCSEC_PER_RAD

GRIDS = {"s1_n6": (1., 6), "s2_n6": (2., 6), "s3_n6": (3., 6), "s2_n10": (2., 10), "s2_n14": (2., 14)}


def gaia_pm_check(data, mom, c=(6., 0.)):
    """Change of the Gaia PM bin predictions (not moved) for a shifted model centre, in bin errors."""
    out = {}
    for p in data.profiles:
        if not p.name.startswith("gaia_edr3_ours"):
            continue
        ps = PointSet.annuli(p.r_lower, p.r_upper, n_r=32, n_phi=512)
        comp = "r" if p.kind == "pmr" else "t"
        S0 = mom(np.hypot(ps.x, ps.y), ("Sigma",))["Sigma"]
        S1 = mom(np.hypot(ps.x-c[0], ps.y-c[1]), ("Sigma",))["Sigma"]
        f = (ps.mean(S1*projected_variance(mom, ps.x, ps.y, c, comp))/ps.mean(S1)) / \
            (ps.mean(S0*projected_variance(mom, ps.x, ps.y, (0., 0.), comp))/ps.mean(S0))
        out[p.name] = float(np.max(np.abs((np.sqrt(f)-1)*p.value/p.err_lo)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="20261103")
    args = ap.parse_args()
    t_start = time.time()
    runs = {k: load_run(v) for k, v in RUNS.items()}
    shas = {json.dumps(r[3]["data_snapshot"]["sha256"])+json.dumps(r[3]["counts_sha256"]) for r in runs.values()}
    if len(shas) != 1:
        raise SystemExit("the three runs do not share the same frozen data")
    _, data, counts, _ = runs["no BH"]
    geom = CentreGeometry(data, counts)
    print("point sets:", {k: int(len(v.x)) for k, v in geom.sets.items()}, flush=True)
    hst = next(c for c in counts if c.name.startswith(cm.HST_NAME))
    n_obs2 = int(np.sum(hst.counts[:2]))

    grids = {g: cm.centre_grid(*v) for g, v in GRIDS.items()}
    centres = sorted({(0., 0.)} | {tuple(np.round(x, 12)) for C, _ in grids.values() for x in C})
    rows = {c: {} for c in centres}
    frozen, tests = {}, {"T2": {}, "gaia_pm_shift6_over_err": {}}
    for k, (cfg, d, cps, _) in runs.items():
        t0 = time.time()
        model = build_df_model(cfg)
        prob = DFJointProblem(d, None, cps)
        ev0 = prob.evaluate(cfg, model)
        frozen[k] = scores(ev0)
        pc = cfg.distance_kpc*1e3/ARCSEC_PER_RAD
        mom = ModelMoments(model, cfg.distance_kpc, r_min_arcsec=max(cfg.numerics.r_min*1.0001/pc, 1e-3))
        ref = geom.reference(mom)
        f0 = geom.factors(mom, (0., 0.), ref)
        tests["T2"][k] = abs(prob.evaluate(cfg, model, factors=f0)["loglike_profiled"]-ev0["loglike_profiled"])
        tests["gaia_pm_shift6_over_err"][k] = gaia_pm_check(d, mom)
        print(f"{k}: built in {time.time()-t0:.0f}s; T2 {tests['T2'][k]:.2e}", flush=True)
        t1 = time.time()
        for c in centres:
            ev = prob.evaluate(cfg, model, factors=geom.factors(mom, c, ref))
            s = scores(ev)
            s["P_inner2"] = float(poisson.cdf(n_obs2, s["mu_inner2"]))
            rows[c][k] = s
        print(f"{k}: {len(centres)} centres in {time.time()-t1:.0f}s", flush=True)

    out = dict(runs=RUNS, frozen=frozen, tests=tests, n_obs_inner2=n_obs2, grids={})
    names = list(runs)
    for g, (C, W) in grids.items():
        L = np.array([[rows[tuple(np.round(x, 12))][k]["lnl"] for x in C] for k in names])
        Lm = cm.marginal_loglike(L, W)
        res = {"lnl_m": dict(zip(names, map(float, Lm)))}
        res["dm2"] = {k: float(-2*(Lm[j]-Lm[0])) for j, k in enumerate(names)}
        dd = L-L[0]
        res["per_centre_dm2"] = {k: dict(mean=float(-2*W@dd[j]), sd=float(2*np.sqrt(W@(dd[j]-W@dd[j])**2)),
                                         min=float(-2*dd[j].max()), max=float(-2*dd[j].min()))
                                 for j, k in enumerate(names)}
        res["ess"], res["post_centre"] = {}, {}
        for j, k in enumerate(names):
            ess, amax = effective_sample_size(L[j], W)
            a = W*np.exp(L[j]-L[j].max()); a /= a.sum()
            mc = a @ C
            res["ess"][k] = dict(ess=ess, max_weight=amax, n=len(W))
            res["post_centre"][k] = dict(mean=mc.tolist(), sd=np.sqrt(a @ (C-mc)**2).tolist())
        res["lnl_spread_sd"] = {k: float(np.sqrt(W@(L[j]-W@L[j])**2)) for j, k in enumerate(names)}
        res["P_inner2_prior_predictive"] = {
            k: float(W@np.array([rows[tuple(np.round(x, 12))][k]["P_inner2"] for x in C])) for k in names}
        out["grids"][g] = res
    tests["T4"] = {f"{g}_minus_s2_n6": {k: out["grids"][g]["dm2"][k]-out["grids"]["s2_n6"]["dm2"][k] for k in names}
                   for g in ("s2_n10", "s2_n14")}
    tests["T4_pass"] = bool(all(abs(v) < .1 for d in tests["T4"].values() for v in d.values()))
    tests["ESS_min_s2_n10"] = min(out["grids"]["s2_n10"]["ess"][k]["ess"] for k in names)
    out["rows"] = [dict(c=list(c), models=rows[c]) for c in centres]
    out["elapsed_s"] = time.time()-t_start
    path = ROOT/f"results/plot_data/centre_model_marginal_{args.label}.json"
    path.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k not in ("rows", "frozen")}, indent=1))
    print("written", path.relative_to(ROOT))


if __name__ == "__main__":
    main()
