#!/usr/bin/env python3
"""Jeans count-tracer fits with a sampled model centre (normal(0, 2") per axis) vs the fixed-centre twins (2026-11-13).

Fixed centre: jeans_K1_turnover_counts_c, jeans_K2_{nfw,cored}_turnover_counts_c (2026-10-16).
Free centre:  the same labels with suffix _centre2 (bin/run_jeans_counts_tracer.py --centre-sigma 2).
Left: kinematic data with the maximum-a-posteriori (MAP) predictions of K1 (fixed centre), K1 and NFW (free centre),
with residuals. MAP = the sampled point (ultranest/chains/weighted_post.txt) maximising ln L + ln prior, the prior
density taken in the sampled coordinates (flat in log10 for loguniform, flat for uniform, Gaussian for normal), so
only the normal priors (distance, centre) enter. Not the maximum-likelihood point: with a free centre that lies
7.6-9.1" south, outside the validated range of the centre model. The free-centre curves include the model-centre
factors (FitProblem.predict); the replayed ln L is checked against the stored value. Right top: posterior of the model centre offset (dRA*, dDec) about AvdM10 for the three
free-centre fits, with the 1- and 2-sigma prior circles. Right bottom: M_bh posteriors.
Output: plots/jeans_centre2_20261113.png, results/plot_data/jeans_centre2_20261113.json.
Usage: PYTHONPATH=src python bin/plot_jeans_centre2.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.report import load_run_metadata, problem_for

CURVES = [("jeans_K1_turnover_counts_c", "no DM, fixed centre", "0.55", "--"),
          ("jeans_K1_turnover_counts_c_centre2", "no DM, free centre", "C3", "-"),
          ("jeans_K2_nfw_turnover_counts_c_centre2", "NFW, free centre", "C0", "-")]
FREE = [("jeans_K1_turnover_counts_c_centre2", "no DM", "C3"),
        ("jeans_K2_nfw_turnover_counts_c_centre2", "NFW", "C0"),
        ("jeans_K2_cored_turnover_counts_c_centre2", "cored", "C2")]
SIGMA = 2.0


def map_point(label, family):
    """Sampled point maximising ln L + ln prior (sampled coordinates); returns (x, lnL, lnprior)."""
    d = ROOT/"results/fits"/label/"ultranest/chains/weighted_post.txt"
    names = open(d).readline().split()
    rows = np.loadtxt(d, skiprows=1)
    cols = {n: rows[:, i] for i, n in enumerate(names)}
    lnprior = np.zeros(len(rows))
    for p in family.parameters:
        if p.prior.kind in ("normal", "truncnormal"):
            lnprior += -0.5*((cols[p.name]-p.prior.mu)/p.prior.sigma)**2
    k = int(np.argmax(cols["logl"]+lnprior))
    return np.array([cols[n][k] for n in family.names]), float(cols["logl"][k]), float(lnprior[k])


def main():
    runs = {}
    for label, *_ in CURVES:
        run = load_run_metadata(ROOT/"results/fits"/label)
        problem = problem_for(run)                          # replays ln L at the ML point
        s = run["summary"]
        x, lnl, lnp = map_point(label, problem.family)
        assert np.isclose(problem.loglike_vector(x), lnl, atol=1e-6), label
        chi2 = problem.chi2(x)
        runs[label] = dict(problem=problem, pred=problem.predict(x), chi2=chi2, summary=s, x=x, lnl=lnl,
                           chi2_total=sum(v[0] for v in chi2.values()),
                           centre={n: float(x[problem.family.names.index(n)]) for n in ("centre_x", "centre_y")
                                   if n in problem.family.names})
    data = runs[CURVES[0][0]]["problem"].data
    profiles = {p.name: p for p in data.profiles}
    fig = plt.figure(figsize=(13, 7.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[2.2, 1], hspace=0.06, wspace=0.34)
    gs2 = fig.add_gridspec(2, 4, hspace=0.35, wspace=0.34)
    order = [("muse_los_dispersion", 0), ("hst_pm_radial_ours", 1), ("gaia_edr3_ours_radial", 1),
             ("hst_pm_tangential_ours", 2), ("gaia_edr3_ours_tangential", 2)]
    axes = {c: (fig.add_subplot(gs[0, c]), fig.add_subplot(gs[1, c])) for c in range(3)}
    for name, c in order:
        p = profiles[name]
        a, r = axes[c]
        mk = "s" if "gaia" in name else "o"
        a.errorbar(p.r, p.value, [p.err_lo, p.err_hi], fmt=mk, color="k", ms=3, lw=0.7,
                   mfc="none" if "gaia" in name else "k")
        for label, _, col, ls in CURVES:
            pred = runs[label]["pred"][name]
            a.plot(p.r, pred, ls, color=col, lw=1.2)
            e = np.where(pred > p.value, p.err_hi, p.err_lo)
            r.plot(p.r, (p.value-pred)/e, mk+ls, color=col, ms=2.5, lw=0.8)
        for q in (a, r):
            q.set_xscale("log")
        a.set_xticklabels([])
        r.axhline(0, color="k", lw=0.5); r.set_ylim(-4, 4); r.set_xlabel("projected radius [arcsec]")
    axes[0][0].set_ylabel(r"$\sigma_{\rm los}$ [km/s]")
    axes[1][0].set_ylabel(r"$\sigma_{\rm R}$ [mas/yr] (HST ●, Gaia □)")
    axes[2][0].set_ylabel(r"$\sigma_{\rm T}$ [mas/yr] (HST ●, Gaia □)")
    axes[0][1].set_ylabel("(data − model)/err")
    for label, name, col, ls in CURVES:
        m = runs[label]
        axes[0][0].plot([], [], ls, color=col, label=f"{name}: MAP χ² {m['chi2_total']:.1f}/89, ln Z {m['summary']['logz']:.1f}")
    fig.legend(loc="upper center", bbox_to_anchor=(0.42, 0.96), ncol=3, fontsize=8, frameon=False)

    out = dict(prior_sigma_arcsec=SIGMA, models={})
    ac, ab = fig.add_subplot(gs2[0, 3]), fig.add_subplot(gs2[1, 3])
    t = np.linspace(0, 2*np.pi, 200)
    for k in (1, 2):
        ac.plot(k*SIGMA*np.cos(t), k*SIGMA*np.sin(t), color="k", lw=0.6, ls=":" if k == 2 else "-")
    for label, name, col in FREE:
        post = Table.read(ROOT/"results/fits"/label/"posterior.ecsv")
        cx, cy = np.asarray(post["centre_x"]), np.asarray(post["centre_y"])
        w = np.asarray(post["weight"]) if "weight" in post.colnames else np.ones(len(cx))
        idx = np.random.default_rng(1).choice(len(cx), 1500, p=w/w.sum())
        ac.scatter(cx[idx], cy[idx], s=1.5, color=col, alpha=0.35, lw=0)
        ac.plot(np.median(cx[idx]), np.median(cy[idx]), "x", color=col, ms=7, mew=1.5, label=name)
        ab.hist(np.asarray(post["M_bh"])[idx]/1e4, bins=30, histtype="step", color=col, density=True, label=f"{name}, free")
        s = json.loads((ROOT/"results/fits"/label/"summary.json").read_text())
        out["models"][label] = dict(logz=s["logz"], logzerr=s["logzerr"], lnL_max=s["lnL_max"],
                                    chi2={k: v["chi2"] for k, v in s["chi2_ml"].items()},
                                    parameters={k: {q: v[q] for q in ("ml", "p16", "p50", "p84")}
                                                for k, v in s["parameters"].items()})
    old = Table.read(ROOT/"results/fits/jeans_K1_turnover_counts_c/posterior.ecsv")
    ab.hist(np.asarray(old["M_bh"])/1e4, bins=30, histtype="step", color="0.55", ls="--", density=True,
            label="no DM, fixed")
    ac.set_aspect("equal"); ac.set_xlim(-7, 7); ac.set_ylim(-7, 7)
    ac.axhline(0, color="k", lw=0.3); ac.axvline(0, color="k", lw=0.3)
    ac.set_xlabel(r"model centre $\Delta$RA* [arcsec]"); ac.set_ylabel(r"$\Delta$Dec [arcsec]")
    ac.text(0.03, 0.97, "about AvdM10; circles: prior 1σ, 2σ", transform=ac.transAxes, fontsize=6.5, va="top")
    ac.legend(fontsize=7, loc="lower right", frameon=False)
    ab.set_xlabel(r"$M_{\rm BH}$ [$10^4$ M$_\odot$]"); ab.set_ylabel("posterior density"); ab.legend(fontsize=6.5, frameon=False)
    fig.savefig(ROOT/"plots/jeans_centre2_20261113.png", dpi=140, bbox_inches="tight")
    (ROOT/"results/plot_data").mkdir(exist_ok=True)
    (ROOT/"results/plot_data/jeans_centre2_20261113.json").write_text(json.dumps(out, indent=1))
    for label, m in runs.items():
        out["models"].setdefault(label, {})["map"] = dict(lnL=m["lnl"], chi2_total=m["chi2_total"], centre=m["centre"],
                                                         chi2={k: v[0] for k, v in m["chi2"].items()})
        print(label, "MAP lnL %.2f chi2 %.2f" % (m["lnl"], m["chi2_total"]), m["centre"],
              {k: round(v[0], 2) for k, v in m["chi2"].items()})
    (ROOT/"results/plot_data/jeans_centre2_20261113.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
