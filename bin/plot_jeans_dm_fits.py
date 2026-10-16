#!/usr/bin/env python3
"""Dark-matter fits of the Jeans reference models, both tracers (2026-10-16).

Runs (results/fits/): J-Trager rung2_K1_turnover, rung2_K2_{nfw,cored}_turnover; J-counts jeans_K1_turnover_counts_c,
jeans_K2_{nfw,cored}_turnover_counts_c. Same 89-bin data, priors (halo: M_dm(<100 pc) log-uniform 1e3-3e7 Msun, r_s
log-uniform 3-1000 pc; NFW gamma = 1 or cored gamma = 0, taper 1 kpc) and sampler; the two families differ only in the
tracer MGE.
Figures:
  plots/jeans_dm_fits_20261016_mass.png   (a) M_DM(<r), posterior median and 16-84% (300 saved draws per run);
                                          (b) DM share M_DM(<r)/M_tot(<r); (c) total M(<r), all six runs.
  plots/jeans_dm_fits_20261016_post.png   (a) posterior of log M_DM(<100 pc) against the (flat-in-log) prior;
                                          (b) profile-likelihood envelope: max ln L of the sampler's points in bins of
                                              log M_DM(<100 pc), relative to the no-DM ln L_max of the same tracer;
                                          (c, d) J-counts posteriors of M_DM(<100 pc) against M_rem and M_star.
  plots/jeans_dm_fits_20261016_kin.png    data - model residuals of the J-counts no-DM, NFW and cored best fits
                                          (exact likelihood-bin predictions; chi2 replayed).
Numbers: results/plot_data/jeans_dm_fits_20261016.json.
Usage: PYTHONPATH=src python bin/plot_jeans_dm_fits.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.report import load_run_metadata, problem_for

RUNS = [("rung2_K1_turnover", "J-Trager, no DM", "C0", ":"), ("rung2_K2_nfw_turnover", "J-Trager + NFW", "C0", "-"),
        ("rung2_K2_cored_turnover", "J-Trager + cored", "C0", "--"),
        ("jeans_K1_turnover_counts_c", "J-counts, no DM", "C3", ":"), ("jeans_K2_nfw_turnover_counts_c", "J-counts + NFW", "C3", "-"),
        ("jeans_K2_cored_turnover_counts_c", "J-counts + cored", "C3", "--")]
NODM = {"C0": "rung2_K1_turnover", "C3": "jeans_K1_turnover_counts_c"}


def load(label):
    d = ROOT/"results/fits"/label
    s = json.loads((d/"summary.json").read_text())
    prof = np.load(d/"profiles.npz")
    wp = np.loadtxt(d/"ultranest/chains/weighted_post.txt", skiprows=1)
    names = open(d/"ultranest/chains/weighted_post.txt").readline().split()
    return dict(summary=s, prof=prof, post=Table.read(d/"posterior.ecsv"), dead=dict(zip(names, wp.T)))


def main():
    R = {lab: load(lab) for lab, *_ in RUNS}
    out = dict(runs={})
    for lab, name, col, ls in RUNS:
        s = R[lab]["summary"]
        out["runs"][lab] = dict(name=name, logz=s["logz"], logzerr=s["logzerr"], lnL_max=s["lnL_max"], chi2_ml=s["chi2_ml_total"],
                                parameters={k: {q: v[q] for q in ("ml", "p16", "p50", "p84")} for k, v in s["parameters"].items()})
    # ---- mass figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for lab, name, col, ls in RUNS:
        p = R[lab]["prof"]; r = p["r"]
        if "K2" in lab:
            lo, med, hi = np.percentile(p["M_dm"], [16, 50, 84], axis=0)
            ax[0].plot(r, med, ls, color=col, lw=1.3, label=name); ax[0].fill_between(r, lo, hi, color=col, alpha=.12 if ls == "-" else .07)
            f = p["M_dm"]/p["M_total"]
            lo, med, hi = np.percentile(f, [16, 50, 84], axis=0)
            ax[1].plot(r, med, ls, color=col, lw=1.3, label=name); ax[1].fill_between(r, lo, hi, color=col, alpha=.10)
        ax[2].plot(r, np.median(p["M_total"], axis=0), ls, color=col, lw=1.3, label=name)
    for a in ax:
        a.set_xscale("log"); a.set_xlabel("r [pc]")
    ax[0].set_yscale("log"); ax[0].set_ylim(1e2, 3e7); ax[0].set_ylabel(r"$M_{\rm DM}(<r)$ [M$_\odot$]")
    ax[1].set_ylim(0, 1); ax[1].set_ylabel(r"$M_{\rm DM}(<r)/M_{\rm tot}(<r)$")
    ax[2].set_yscale("log"); ax[2].set_ylim(1e3, 2e7); ax[2].set_ylabel(r"$M_{\rm tot}(<r)$ [M$_\odot$] (median)")
    for a, t in zip(ax, "abc"):
        a.axvspan(0.1, 63, color="0.93", zorder=0); a.text(.03, .94, t, transform=a.transAxes, fontweight="bold", va="top")
    ax[0].legend(fontsize=7, loc="lower right"); ax[2].legend(fontsize=7, loc="lower right")
    ax[1].text(.97, .94, "shaded: radii with kinematic data\nbands: 16–84%", transform=ax[1].transAxes, ha="right", va="top", fontsize=7)
    fig.savefig(ROOT/"plots/jeans_dm_fits_20261016_mass.png", dpi=140, bbox_inches="tight"); plt.close(fig)
    # ---- posterior figure
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.0))
    bins = np.linspace(3, np.log10(3e7), 36)
    for lab, name, col, ls in RUNS:
        if "K2" not in lab:
            continue
        x = np.log10(np.asarray(R[lab]["post"]["M_dm_100"], float))
        ax[0].hist(x, bins, density=True, histtype="step", color=col, ls=ls, lw=1.4, label=name)
        dd = R[lab]["dead"]; lx = np.log10(dd["M_dm_100"]); ll = dd["logl"]
        ref = R[NODM[col]]["summary"]["lnL_max"]
        c, env = [], []
        for a, b in zip(bins[:-1], bins[1:]):
            k = (lx >= a) & (lx < b)
            if k.sum():
                c.append(.5*(a+b)); env.append(ll[k].max()-ref)
        ax[1].plot(c, env, ls, color=col, lw=1.4, marker=".", label=name)
        out["runs"][lab]["profile_envelope"] = dict(log10_Mdm100=c, dlnL_vs_noDM=env)
    ax[0].axhline(1/(bins[-1]-bins[0]), color="k", lw=.8, ls=":", label="prior")
    ax[0].set_xlabel(r"$\log_{10} M_{\rm DM}(<100\,{\rm pc})$ [M$_\odot$]"); ax[0].set_ylabel("posterior density"); ax[0].legend(fontsize=7)
    ax[1].axhline(0, color="k", lw=.6); ax[1].set_ylim(-8, 5)
    ax[1].set_xlabel(r"$\log_{10} M_{\rm DM}(<100\,{\rm pc})$"); ax[1].set_ylabel(r"max $\ln L$ − no-DM $\ln L_{\max}$ (same tracer)")
    ax[1].legend(fontsize=7)
    for j, lab in enumerate(("jeans_K2_nfw_turnover_counts_c", "jeans_K2_cored_turnover_counts_c")):
        t = R[lab]["post"]
        ax[2].scatter(np.log10(t["M_dm_100"]), np.log10(t["M_rem"]), s=1, alpha=.15, color=("C3", "C1")[j], rasterized=True,
                      label=("J-counts + NFW", "J-counts + cored")[j])
        ax[3].scatter(np.log10(t["M_dm_100"]), np.log10(t["M_star"]), s=1, alpha=.15, color=("C3", "C1")[j], rasterized=True)
    ax[2].set_xlabel(r"$\log_{10} M_{\rm DM}(<100\,{\rm pc})$"); ax[2].set_ylabel(r"$\log_{10} M_{\rm rem}$")
    ax[3].set_xlabel(r"$\log_{10} M_{\rm DM}(<100\,{\rm pc})$"); ax[3].set_ylabel(r"$\log_{10} M_\star$")
    leg = ax[2].legend(fontsize=7, markerscale=8)
    for a, t in zip(ax, "abcd"):
        a.text(.03, .94, t, transform=a.transAxes, fontweight="bold", va="top")
    fig.savefig(ROOT/"plots/jeans_dm_fits_20261016_post.png", dpi=140, bbox_inches="tight"); plt.close(fig)
    # ---- kinematic residuals (J-counts trio)
    trio = [r for r in RUNS if r[0].endswith("counts_c")]
    preds = {}
    for lab, *_ in trio:
        run = load_run_metadata(ROOT/"results/fits"/lab)
        prob = problem_for(run)
        x = np.array([run["summary"]["parameters"][n]["ml"] for n in prob.family.names])
        preds[lab] = (prob, prob.predict(x), prob.chi2(x))
        assert np.isclose(sum(v[0] for v in preds[lab][2].values()), run["summary"]["chi2_ml_total"], atol=.011)
    data = preds[trio[0][0]][0].data
    fig, ax = plt.subplots(1, 5, figsize=(18, 3.6), sharey=True)
    for i, p in enumerate(data.profiles):
        for (lab, name, col, ls), mk in zip(trio, ("o", "s", "^")):
            pr = preds[lab][1][p.name]
            e = np.where(pr > p.value, p.err_hi, p.err_lo)
            c2 = preds[lab][2][p.name][0]
            ax[i].plot(p.r, (p.value-pr)/e, marker=mk, ls=ls, ms=3, lw=.9, color={":": "0.4", "-": "C3", "--": "C1"}[ls],
                       label=f"{name}: χ² {c2:.1f}")
        ax[i].axhline(0, color="k", lw=.5); ax[i].set_xscale("log"); ax[i].set_title(f"{p.name} ({p.n} bins)", fontsize=8)
        ax[i].set_xlabel("R [arcsec]"); ax[i].legend(fontsize=6)
    ax[0].set_ylabel("(data − model)/err"); ax[0].set_ylim(-4, 4)
    fig.savefig(ROOT/"plots/jeans_dm_fits_20261016_kin.png", dpi=140, bbox_inches="tight"); plt.close(fig)
    (ROOT/"results/plot_data/jeans_dm_fits_20261016.json").write_text(json.dumps(out, indent=1))
    for lab, o in out["runs"].items():
        env = o.get("profile_envelope")
        print(lab, round(o["logz"], 2), round(o["lnL_max"], 2), round(o["chi2_ml"], 1),
              "" if env is None else f"envelope max {max(env['dlnL_vs_noDM']):+.2f} at 10^{env['log10_Mdm100'][int(np.argmax(env['dlnL_vs_noDM']))]:.2f}")


if __name__ == "__main__":
    main()
