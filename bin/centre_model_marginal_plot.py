#!/usr/bin/env python3
"""Figure for bin/centre_model_marginal.py: -2 Delta ln L maps over the model centre and the marginal vs sigma_c.

Usage: python bin/centre_model_marginal_plot.py [--label 20261103] -> plots/centre_model_marginal_<label>.png
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="20261103")
    a = ap.parse_args()
    d = json.loads((ROOT/f"results/plot_data/centre_model_marginal_{a.label}.json").read_text())
    names = list(d["runs"])
    C = np.array([r["c"] for r in d["rows"]])
    L = {k: np.array([r["models"][k]["lnl"] for r in d["rows"]]) for k in names}
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
    for j, k in enumerate(names):
        z = -2*(L[k]-L[k].max())
        sc = ax[j].scatter(C[:, 0], C[:, 1], c=np.minimum(z, 20), s=14, cmap="viridis_r")
        for s in (1, 2, 3):
            ax[j].add_patch(plt.Circle((0, 0), 2*s, fill=False, lw=.6, ls=":", color="k"))
        ax[j].plot(0, 0, "r+", ms=10)
        ax[j].set_title(f"{k}: -2(lnL - max) over model centre")
        ax[j].set_xlabel("dRA* [arcsec]"); ax[j].set_ylabel("dDec [arcsec]"); ax[j].set_aspect("equal")
        plt.colorbar(sc, ax=ax[j], fraction=.046)
    s_c = {"s1_n6": 1, "s2_n14": 2, "s3_n6": 3}
    for k, col in zip(names[1:], ("C0", "C3")):
        x = [s_c[g] for g in s_c]
        ax[3].plot(x, [d["grids"][g]["dm2"][k] for g in s_c], "o-", color=col, label=f"{k}: all bins")
        ax[3].plot(x, [d["grids"][g]["dm2_drop2"][k] for g in s_c], "s--", color=col, label=f"{k}: without 2 innermost PM bins")
        ax[3].plot([0], [d["grids"]["s2_n6"]["dm2_at_avdm10"][k]], "*", color=col, ms=12)
        ax[3].plot([0], [d["grids"]["s2_n6"]["dm2_drop2_at_avdm10"][k]], "*", mfc="none", color=col, ms=12)
    ax[3].axhline(0, color="k", lw=.5)
    ax[3].set_xlabel("centre prior sigma_c per axis [arcsec] (0 = AvdM10 fixed)")
    ax[3].set_ylabel("marginal -2 Delta lnL vs no BH")
    ax[3].legend(fontsize=7)
    fig.suptitle("Model-side centre marginalisation (data fixed about AvdM10; saved fits, not refitted)")
    fig.tight_layout()
    out = ROOT/f"plots/centre_model_marginal_{a.label}.png"
    fig.savefig(out, dpi=110)
    print("wrote", out.relative_to(ROOT))


if __name__ == "__main__":
    main()
