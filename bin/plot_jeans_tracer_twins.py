#!/usr/bin/env python3
"""J-Trager vs J-counts: the two Jeans reference models that differ only in their tracer (2026-10-16).

J-Trager = rung2_K1_turnover (tracer and stellar-mass shape: MGE of the composite HST-count/Trager light profile);
J-counts = jeans_K1_turnover_counts_c (MGE fitted to the completeness-corrected HST + Gaia counts with the tied
amplitude ratio). Same 89-bin data snapshot, parameters, priors and sampler. Curves are the exact likelihood-bin
predictions at each maximum-likelihood point (chi2 replayed against the saved value).
Panels: the four kinematic datasets (LOS, HST PM radial/tangential, Gaia PM radial/tangential) with residuals; the
two tracer surface-density shapes (normalised to the same HST count amplitude); enclosed mass (posterior median).
Output: plots/jeans_tracer_twins_20261016.png, results/plot_data/jeans_tracer_twins_20261016.json.
Usage: PYTHONPATH=src python bin/plot_jeans_tracer_twins.py
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

from ocen_dm.kinematics.report import load_run_metadata, problem_for

FITS = [("rung2_K1_turnover", "J-Trager (composite light profile)", "C0"),
        ("jeans_K1_turnover_counts_c", "J-counts (HST + Gaia counts)", "C3")]
PANELS = [("muse_los_dispersion", r"$\sigma_{\rm los}$ [km/s]"), ("hst_pm_radial_ours", r"HST $\sigma_{\rm R}$ [mas/yr]"),
          ("hst_pm_tangential_ours", r"HST $\sigma_{\rm T}$ [mas/yr]"), ("gaia_edr3_ours_radial", r"Gaia $\sigma_{\rm R}$ [mas/yr]"),
          ("gaia_edr3_ours_tangential", r"Gaia $\sigma_{\rm T}$ [mas/yr]")]


def main():
    runs = []
    for label, name, col in FITS:
        run = load_run_metadata(ROOT/"results/fits"/label)
        problem = problem_for(run)
        s = run["summary"]
        x = np.array([s["parameters"][n]["ml"] for n in problem.family.names])
        pred, chi2 = problem.predict(x), problem.chi2(x)
        assert np.isclose(sum(v[0] for v in chi2.values()), s["chi2_ml_total"], atol=.011)
        prof = np.load(ROOT/"results/fits"/label/"profiles.npz")
        runs.append(dict(label=label, name=name, col=col, problem=problem, pred=pred, chi2=chi2, summary=s, prof=prof,
                         mge=problem.family.mge_fit))
    data = runs[0]["problem"].data
    fig = plt.figure(figsize=(13, 7.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[2.2, 1], hspace=0.06, wspace=0.32)
    gs2 = fig.add_gridspec(2, 4, hspace=0.35, wspace=0.32)
    out = dict(models={})
    profiles = {p.name: p for p in data.profiles}
    order = [("muse_los_dispersion", 0), ("hst_pm_radial_ours", 1), ("gaia_edr3_ours_radial", 1),
             ("hst_pm_tangential_ours", 2), ("gaia_edr3_ours_tangential", 2)]
    axes = {}
    for col in range(3):
        axes[col] = (fig.add_subplot(gs[0, col]), fig.add_subplot(gs[1, col]))
    for name, col in order:
        p = profiles[name]
        a, r = axes[col]
        mk = "s" if "gaia" in name else "o"
        a.errorbar(p.r, p.value, [p.err_lo, p.err_hi], fmt=mk, color="k", ms=3, lw=0.7, mfc="none" if "gaia" in name else "k")
        for m in runs:
            a.plot(p.r, m["pred"][name], "-", color=m["col"], lw=1.2)
            e = np.where(m["pred"][name] > p.value, p.err_hi, p.err_lo)
            r.plot(p.r, (p.value-m["pred"][name])/e, mk+"-", color=m["col"], ms=2.5, lw=0.8)
        for q in (a, r):
            q.set_xscale("log")
        r.axhline(0, color="k", lw=0.5); r.set_ylim(-4, 4); r.set_xlabel("projected radius [arcsec]")
    axes[0][0].set_ylabel(r"$\sigma_{\rm los}$ [km/s]"); axes[1][0].set_ylabel(r"$\sigma_{\rm R}$ [mas/yr] (HST ●, Gaia □)")
    axes[2][0].set_ylabel(r"$\sigma_{\rm T}$ [mas/yr] (HST ●, Gaia □)"); axes[0][1].set_ylabel("(data − model)/err")
    for m in runs:
        tot = sum(v[0] for v in m["chi2"].values())
        axes[0][0].plot([], [], color=m["col"], label=f"{m['name']}: χ² {tot:.1f}/89, ln L$_{{\\max}}$ {m['summary']['lnL_max']:.1f}")
    fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=2, fontsize=8, frameon=False)
    # tracer shapes and enclosed mass
    at = fig.add_subplot(gs2[0, 3]); am = fig.add_subplot(gs2[1, 3])
    Rg = np.geomspace(1, 3000, 300)
    ref = runs[0]["mge"].surface_intensity(Rg)
    for m in runs:
        at.plot(Rg, m["mge"].surface_intensity(Rg)/ref, color=m["col"], lw=1.2)
        Ms = np.percentile(m["prof"]["M_total"], 50, axis=0)
        am.plot(m["prof"]["r"], Ms, color=m["col"], lw=1.2)
        am.plot(m["prof"]["r"], np.percentile(m["prof"]["M_stars"], 50, axis=0), "--", color=m["col"], lw=0.8)
        am.plot(m["prof"]["r"], np.percentile(m["prof"]["M_remnants"], 50, axis=0), ":", color=m["col"], lw=0.8)
        out["models"][m["label"]] = dict(chi2={k: v[0] for k, v in m["chi2"].items()}, lnL_max=m["summary"]["lnL_max"],
                                         logz=m["summary"]["logz"], parameters={k: {q: v[q] for q in ("ml", "p16", "p50", "p84")}
                                                                                for k, v in m["summary"]["parameters"].items()})
    at.set_xscale("log"); at.axhline(1, color="k", lw=0.5); at.set_ylim(0.5, 1.6)
    at.set_xlabel("R [arcsec]"); at.set_ylabel("tracer Σ / J-Trager tracer Σ\n(each MGE normalised to unit total)")
    am.set_xscale("log"); am.set_yscale("log"); am.set_xlabel("r [pc]"); am.set_ylabel(r"$M(<r)$ [M$_\odot$]")
    am.text(0.03, 0.9, "solid total, dashed stars, dotted remnants\n(posterior medians)", transform=am.transAxes, fontsize=6.5, va="top")
    am.set_ylim(1e3, 5e6)
    fig.savefig(ROOT/"plots/jeans_tracer_twins_20261016.png", dpi=140, bbox_inches="tight")
    (ROOT/"results/plot_data/jeans_tracer_twins_20261016.json").write_text(json.dumps(out, indent=1))
    for m in runs:
        print(m["label"], {k: round(v[0], 2) for k, v in m["chi2"].items()})


if __name__ == "__main__":
    main()
