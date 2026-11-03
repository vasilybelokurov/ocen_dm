#!/usr/bin/env python3
"""Figures for bin/centre_marginal.py (docs/CENTRE_MARGINALISATION.md).

1. plots/centre_marginal_<label>.png -- (a, b) per-centre -2 Delta lnL of the BH models relative to
   no BH on the n = 10, sigma_c = 2" grid, with the 1 and 2 sigma_c prior circles; (c) marginal
   -2 Delta lnL against sigma_c (points: marginal; bars: weighted s.d. over centres) with the
   fixed-centre value; (d) per-centre split of -2 Delta lnL (BH 4e4) into datasets.
2. plots/centre_marginal_<label>_data.png -- inner data against centre (grey: every grid centre;
   black: AvdM10) with the three models (fixed at their saved fits) and residuals.
Usage: python bin/centre_marginal_plot.py [--label 20261103]
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COL = {"no BH": "C0", "BH 8200": "C1", "BH 40000": "C3"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="20261103")
    a = ap.parse_args()
    J = json.loads((ROOT/f"results/plot_data/centre_marginal_{a.label}.json").read_text())
    rows = J["rows"]
    idx = {tuple(r["c"]): r for r in rows}
    import sys
    sys.path.insert(0, str(ROOT/"src"))
    from ocen_dm.centre_marg import centre_grid
    C, W = centre_grid(2., 10)
    get = lambda c: idx[tuple(np.round(c, 12))]

    fig, ax = plt.subplots(1, 4, figsize=(20, 4.8), constrained_layout=True)
    for j, k in enumerate(("BH 8200", "BH 40000")):
        d = np.array([-2*(get(c)["models"][k]["lnl"]-get(c)["models"]["no BH"]["lnl"]) for c in C])
        sc = ax[j].scatter(C[:, 0], C[:, 1], c=d, s=20+800*W, cmap="coolwarm", vmin=-np.abs(d).max(), vmax=np.abs(d).max())
        for s in (2., 4.):
            ax[j].add_patch(plt.Circle((0, 0), s, fill=False, ls=":", color="k"))
        ax[j].plot(0, 0, "k+", ms=12)
        ax[j].plot(-1.95, 1.12, "g*", ms=12, label="joint centre fit")
        ax[j].plot(-1.31, 1.44, "mx", ms=9, label="oMEGACat VI kinematic")
        ax[j].set(xlabel="dRA* [arcsec]", ylabel="dDec [arcsec]", aspect="equal",
                  title=f"{k} - no BH: -2 dlnL per centre (size ~ weight)")
        ax[j].legend(fontsize=7, loc="lower left")
        fig.colorbar(sc, ax=ax[j], shrink=.8)
    sig = [1., 2., 3.]
    for k in ("BH 8200", "BH 40000"):
        m = [J["grids"][f"s{int(s)}_n6"]["lnl_dm2"][k]["marginal"] for s in sig]
        sd = [J["grids"][f"s{int(s)}_n6"]["lnl_dm2"][k]["sd"] for s in sig]
        ax[2].errorbar(sig, m, sd, marker="o", capsize=3, color=COL[k], label=f"{k} (n = 6)")
        ax[2].plot([2.], [J["grids"]["s2_n10"]["lnl_dm2"][k]["marginal"]], "s", mfc="none", color=COL[k], ms=9,
                   label=f"{k} (n = 10)")
        f0 = -2*(J["frozen"][k]["lnl"]-J["frozen"]["no BH"]["lnl"])
        ax[2].plot([0], [f0], "D", color=COL[k], label=f"{k} at AvdM10 only")
    ax[2].axhline(0, color="k", lw=.8)
    ax[2].axhline(5, color="0.5", ls="--", lw=.8)
    ax[2].set(xlabel="prior width sigma_c [arcsec]", ylabel="-2 dlnL (BH - no BH), marginal",
              title="marginal over the centre (bars: s.d. over centres)")
    ax[2].legend(fontsize=7)
    keys = [k for k in rows[0]["models"]["no BH"] if k.startswith("lnl:")]
    D = {k: np.array([-2*(get(c)["models"]["BH 40000"][k]-get(c)["models"]["no BH"][k]) for c in C]) for k in keys}
    lab = [k[4:].replace("ocen_counts_", "counts ").replace("_ours", "") for k in keys]
    ax[3].boxplot([D[k] for k in keys], vert=False, labels=lab, whis=(0, 100))
    ax[3].axvline(0, color="k", lw=.8)
    ax[3].set(xlabel="-2 dlnL (BH 4e4 - no BH) per dataset", title="split over datasets, all n = 10 centres")
    fig.savefig(ROOT/f"plots/centre_marginal_{a.label}.png", dpi=130)

    # ------------------------------------------------------------------ data
    from ocen_dm.kinematics.run_io import read_data_snapshot
    od = ROOT/"results/df/dftwo_counts_area_20261013/observed"
    rec = json.loads((od/"problem.json").read_text())
    prof = {p.name: p for p in read_data_snapshot(od, rec["data_snapshot"]).profiles}
    cnt = json.loads((od/"ocen_counts_hst_f625w19.json").read_text())
    lo, hi = np.array(cnt["r_lower"])[:14], np.array(cnt["r_upper"])[:14]
    rc, area = np.sqrt(lo*hi), np.array(cnt["area_arcsec2"])[:14]
    fig, ax = plt.subplots(2, 4, figsize=(20, 7.5), sharex="col", gridspec_kw=dict(height_ratios=(2.2, 1)),
                           constrained_layout=True)
    r0 = idx[(0., 0.)]
    # counts
    for c in C:
        ax[0, 0].plot(rc, np.array(get(c)["counts_hst"][:14])/area, color="0.8", lw=.6)
    ax[0, 0].errorbar(rc, np.array(r0["counts_hst"][:14])/area, np.sqrt(np.array(r0["counts_hst"][:14]))/area,
                      fmt="ko", ms=3, capsize=2, label="AvdM10")
    for k in COL:
        mu = np.array(r0["models"][k]["mu_hst"])
        ax[0, 0].plot(rc, mu/area, color=COL[k], label=k)
        ax[1, 0].plot(rc, (np.array(r0["counts_hst"][:14])-mu)/np.sqrt(mu), "o-", ms=3, color=COL[k])
    spread = np.array([[(get(c)["counts_hst"][i]-r0["models"]["no BH"]["mu_hst"][i])/np.sqrt(r0["models"]["no BH"]["mu_hst"][i])
                        for i in range(14)] for c in C])
    ax[1, 0].fill_between(rc, spread.min(0), spread.max(0), color="0.85", label="no BH, all centres")
    ax[0, 0].set(xscale="log", yscale="log", ylabel="F625W<19 stars per arcsec$^2$", title="HST counts 2-59\"")
    ax[1, 0].set(xscale="log", xlabel="R [arcsec]", ylabel="(N - mu)/sqrt(mu)")
    ax[0, 0].legend(fontsize=7)
    # kinematics
    for j, (name, key, n) in enumerate((("hst_pm_radial_ours", "hst_pmr", 6), ("hst_pm_tangential_ours", "hst_pmt", 6),
                                        ("muse_los_dispersion", "muse", 3))):
        p = prof[name]
        r = p.r[:n]
        V = np.array([get(c)[key] for c in C])
        for v in V:
            ax[0, j+1].plot(r, v, color="0.8", lw=.6)
        ax[0, j+1].errorbar(r, p.value[:n], [p.err_lo[:n], p.err_hi[:n]], fmt="ko", ms=3, capsize=2, label="AvdM10")
        ax[0, j+1].errorbar(r*1.02, W@V, np.sqrt(W@(V-W@V)**2), fmt="s", color="g", mfc="none", ms=4, capsize=2,
                            label="prior-weighted mean +- s.d. over centres")
        for k in COL:
            pr = np.array(r0["models"][k][f"pred:{name}"])[:n]
            ax[0, j+1].plot(r, pr, color=COL[k], label=k)
            ax[1, j+1].plot(r, (p.value[:n]-pr)/p.err_lo[:n], "o-", ms=3, color=COL[k])
        prn = np.array(r0["models"]["no BH"][f"pred:{name}"])[:n]
        res = (V-prn)/p.err_lo[:n]
        ax[1, j+1].fill_between(r, res.min(0), res.max(0), color="0.85")
        ax[0, j+1].set(title=name+(" (published + our shift)" if j == 2 else ""),
                       ylabel="km/s" if j == 2 else "mas/yr")
        ax[1, j+1].set(xlabel="R [arcsec]", ylabel="(data - model)/err")
        ax[0, j+1].legend(fontsize=6)
    for a_ in ax[1]:
        a_.axhline(0, color="k", lw=.8)
    fig.suptitle("Inner data against the centre: grey = each of the 100 grid centres (sigma_c = 2\"), black = AvdM10; "
                 "models fixed at their saved fits; residual bands = no-BH model over all centres")
    fig.savefig(ROOT/f"plots/centre_marginal_{a.label}_data.png", dpi=130)
    print("written plots")


if __name__ == "__main__":
    main()
