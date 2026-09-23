#!/usr/bin/env python3
"""Central kinematics: is there evidence for a point mass beyond the current compact DF?

Compares, bin by bin inside R_max (default 2 pc), the observed HST PM (radial,
tangential) and MUSE LOS dispersions with (a) the best two-transition compact
DF fit (no black hole) and (b) the rung-2 Jeans turnover fits, which include a
free point mass (ML M_bh ~ 4.4e4 Msun). Predictions are the saved likelihood-bin
values: the DF ones from results/plot_data/<batch>_best_fits.json and the Jeans
ones from results/plot_data/rung0_rung1_rung2_no_dm_vs_dm.ecsv. Also reports the
chi2 of the central bins alone for each model.

Usage: python bin/plot_central_bh_evidence.py
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
PC_PER_ARCSEC = 5.43e3*np.pi/(180*3600)
R_MAX = 2.0
DATASETS = [("hst_pm_radial_ours", r"HST $\sigma_{\rm pm,R}$ [mas/yr]"),
            ("hst_pm_tangential_ours", r"HST $\sigma_{\rm pm,T}$ [mas/yr]"),
            ("muse_los_dispersion", r"MUSE $\sigma_{\rm los}$ [km/s]")]
DF_FITS = {"DF two transitions, no BH, no halo": ("dftwo_20260923_best_fits", "no_halo_start0")}
JEANS = {"Jeans turnover + BH, no DM": "rung2_K1_turnover",
         "Jeans turnover + BH, cored DM": "rung2_K2_cored_turnover"}


def chi2(value, model, lo, hi):
    err = np.where(model > value, hi, lo)
    return (value-model)/err


def main():
    table = Table.read(ROOT/"results/plot_data/rung0_rung1_rung2_no_dm_vs_dm.ecsv")
    bh = {run: json.loads((ROOT/"results/fits"/run/"summary.json").read_text())["parameters"]["M_bh"]
          for run in JEANS.values()}
    fig, axes = plt.subplots(2, 3, figsize=(14, 6.5), sharex=True, constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[3, 1.4]))
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), R_max_pc=R_MAX,
                  jeans_M_bh={k: {q: v[q] for q in ("ml", "p16", "p50", "p84")} for k, v in bh.items()},
                  central_chi2={}, bins={})
    for col, (name, ylabel) in enumerate(DATASETS):
        rows = table[table["dataset"] == name]
        R = np.asarray(rows["radius_arcsec"])*PC_PER_ARCSEC
        keep = R <= R_MAX
        v, lo, hi = (np.asarray(rows[k])[keep] for k in ("value", "err_lo", "err_hi"))
        ax, rx = axes[0, col], axes[1, col]
        ax.errorbar(R[keep], v, yerr=[lo, hi], fmt="o", color="black", ms=4, label="observed")
        record["bins"][name] = dict(R_pc=R[keep].tolist(), value=v.tolist())
        models = {}
        for label, (plotdata, fit) in DF_FITS.items():
            pred = json.loads((ROOT/"results/plot_data"/(plotdata+".json")).read_text())["fits"][fit][name]
            models[label] = np.asarray(pred)[keep]
        for label, run in JEANS.items():
            models[label] = np.asarray(rows[run])[keep]
        for (label, m), style in zip(models.items(), ("-", "--", ":")):
            z = chi2(v, m, lo, hi)
            ax.plot(R[keep], m, ls=style, lw=2, label=label)
            rx.plot(R[keep], z, ls=style, marker=".", lw=1)
            record["central_chi2"].setdefault(label, {})[name] = dict(chi2=float(z @ z), n=int(keep.sum()))
            record["bins"][name][label] = m.tolist()
        rx.axhline(0, color="black", lw=.7)
        ax.set(xscale="log", ylabel=ylabel)
        rx.set(xscale="log", xlabel="projected radius [pc]", ylabel="(data-model)/err")
        if col == 0:
            ax.legend(fontsize=7)
    text = "; ".join(f"{k}: M_bh ML {v['ml']:.2g} (p16-p84 {v['p16']:.2g}-{v['p84']:.2g})" for k, v in bh.items())
    fig.suptitle(f"Central kinematics (R <= {R_MAX} pc). Jeans point masses: {text}", fontsize=9)
    plot = ROOT/"plots"/"central_bh_evidence_20260923.png"
    fig.savefig(plot, dpi=160)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/"central_bh_evidence_20260923.json").write_text(json.dumps(record))
    for label, per in record["central_chi2"].items():
        total = sum(t["chi2"] for t in per.values()); n = sum(t["n"] for t in per.values())
        print(f"{label:38s} central chi2 {total:6.1f} / {n}  ", {k.split('_')[0]+'_'+k.split('_')[-2]: round(t['chi2'], 1) for k, t in per.items()})
    print(plot)


if __name__ == "__main__":
    main()
