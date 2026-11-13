#!/usr/bin/env python3
"""Corner plot of a Jeans posterior (equal-weight samples in results/fits/<label>/posterior.ecsv).

Masses and scale radii are shown as log10; the maximum-likelihood point (summary.json "ml") is marked in red.
Usage: PYTHONPATH=src python bin/plot_jeans_corner.py LABEL [LABEL ...]
Output: plots/corner_<label>.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import corner
import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
LOG = {"M_star", "M_rem", "M_bh", "M_dm_100", "r_s", "a_rem", "r_beta"}
TEX = {"M_star": r"$\log M_\star$", "M_rem": r"$\log M_{\rm rem}$", "a_rem": r"$\log a_{\rm rem}$",
       "M_bh": r"$\log M_{\rm BH}$", "beta_0": r"$\beta_0$", "beta_mid": r"$\beta_{\rm mid}$",
       "beta_inf": r"$\beta_\infty$", "r_beta": r"$\log r_\beta$", "delta_r_beta": r"$\Delta r_\beta$",
       "M_dm_100": r"$\log M_{\rm DM}(<100\,{\rm pc})$", "r_s": r"$\log r_s$", "distance": r"$D$ [kpc]",
       "centre_x": r"$\Delta$RA* ['']", "centre_y": r"$\Delta$Dec ['']"}


def main(labels):
    for label in labels:
        d = ROOT/"results/fits"/label
        post = Table.read(d/"posterior.ecsv")
        s = json.loads((d/"summary.json").read_text())
        names = [n for n in post.colnames if n in s["parameters"]]
        f = lambda n, v: np.log10(v) if n in LOG else v
        x = np.column_stack([f(n, np.asarray(post[n], float)) for n in names])
        ml = [f(n, s["parameters"][n]["ml"]) for n in names]
        fig = corner.corner(x, labels=[TEX.get(n, n) for n in names], truths=ml, truth_color="C3",
                            quantiles=[.16, .5, .84], show_titles=True, title_fmt=".2f",
                            title_kwargs=dict(fontsize=8), label_kwargs=dict(fontsize=9), bins=30, smooth=1.0)
        fig.suptitle(f"{label}: ln Z = {s['logz']:.2f} ± {s['logzerr']:.2f}, ln L$_{{\\max}}$ = {s['lnL_max']:.2f}, "
                     f"{len(post)} samples; red = maximum likelihood", fontsize=12, y=1.01)
        out = ROOT/"plots"/f"corner_{label}.png"
        fig.savefig(out, dpi=90, bbox_inches="tight")
        print(out)


if __name__ == "__main__":
    main(sys.argv[1:])
