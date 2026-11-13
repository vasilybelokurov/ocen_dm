#!/usr/bin/env python3
"""Posterior anisotropy profiles beta(r) of Jeans runs with the turnover anisotropy.

beta = 1 - sigma_t^2/sigma_r^2 from TurnoverAnisotropy for each posterior sample (equal-weight posterior.ecsv);
median and 16-84 / 2.5-97.5 per cent bands vs 3D radius r [pc]. Grey ticks on the top axis: projected radii of the
data bins (arcsec converted at D = 5.43 kpc; projected R, so only indicative of the 3D radii probed). Dotted lines:
the prior edges beta_0 <= -0.5 (An & Evans 2006 bound with a central point mass) and beta_inf >= -1.
Usage: PYTHONPATH=src python bin/plot_jeans_beta_profile.py LABEL [LABEL ...]
Output: plots/jeans_beta_profile_20261113.png
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.anisotropy import TurnoverAnisotropy
from ocen_dm.kinematics.report import data_for, load_run_metadata

NAMES = {"jeans_K1_turnover_counts_c": ("no DM, fixed centre", "0.5"),
         "jeans_K1_turnover_counts_c_centre2": ("no DM, free centre", "C3"),
         "jeans_K2_nfw_turnover_counts_c_centre2": ("NFW, free centre", "C0"),
         "jeans_K2_cored_turnover_counts_c_centre2": ("cored, free centre", "C2")}
PC_PER_ARCSEC = 5430/206264.806


def main(labels):
    r = np.geomspace(0.01, 200, 300)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for label in labels:
        post = Table.read(ROOT/"results/fits"/label/"posterior.ecsv")
        b = np.array([TurnoverAnisotropy(*(float(row[k]) for k in ("beta_0", "beta_mid", "beta_inf", "r_beta",
                                                                    "delta_r_beta"))).beta(r) for row in post])
        name, col = NAMES.get(label, (label, None))
        q = np.percentile(b, [2.5, 16, 50, 84, 97.5], axis=0)
        line, = ax.plot(r, q[2], color=col, lw=1.6, label=name)
        ax.fill_between(r, q[1], q[3], color=line.get_color(), alpha=0.25, lw=0)
        ax.fill_between(r, q[0], q[4], color=line.get_color(), alpha=0.10, lw=0)
        print(label, "beta at 0.1/1/10/30/100 pc:",
              [f"{np.interp(x, r, q[2]):.2f} [{np.interp(x, r, q[1]):.2f},{np.interp(x, r, q[3]):.2f}]"
               for x in (0.1, 1, 10, 30, 100)])
    run = load_run_metadata(ROOT/"results/fits"/labels[0])
    data, _ = data_for(run["summary"], run_dir=run["dir"])
    R = np.concatenate([p.r for p in data.profiles])*PC_PER_ARCSEC
    ax.plot(R, np.full(len(R), 0.97), "|", color="0.4", ms=8, transform=ax.get_xaxis_transform())
    ax.axhline(0, color="k", lw=0.5)
    ax.axhline(-0.5, color="k", lw=0.6, ls=":")
    ax.axhline(-1, color="k", lw=0.6, ls=":")
    ax.text(0.012, -0.47, r"prior edge $\beta_0\leq-0.5$ (centre)", fontsize=7)
    ax.text(25, -0.97, r"prior edge $\beta_\infty\geq-1$ (outskirts)", fontsize=7)
    ax.set_xscale("log"); ax.set_xlim(r[0], r[-1]); ax.set_ylim(-1.1, 0.6)
    ax.set_xlabel("3D radius r [pc]"); ax.set_ylabel(r"$\beta = 1-\sigma_t^2/\sigma_r^2$")
    ax.set_title("Jeans anisotropy profiles: posterior median, 68% and 95% bands; ticks = data bin radii (projected)",
                 fontsize=8)
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(0, 0.93))
    out = ROOT/"plots/jeans_beta_profile_20261113.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main(sys.argv[1:])
