#!/usr/bin/env python3
"""Action coverage of the data: q = J_r + L of DF tracers versus projected radius.

Samples the selected fit's DF (tracer = stellar mass DF; mass follows light),
computes actions in its converged potential, and shows p(q | R_proj) with the
5/50/95% quantiles and the anisotropy transition scales of chosen fits. The
bias of the regularized DF is a function of q, so transition scales outside the
sampled q range are not identified. Projection assumes spherical symmetry.

Usage: python bin/plot_action_coverage.py results/df/observed_df_wide_20260923 observed_no_halo_start0
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict

# Transition scales of the subset fits (docs/DF_ANISOTROPY_DIAGNOSTICS.md).
MARKS = {"J_a, HST+MUSE fit (32)": 32., "J_a, Gaia fit (1665)": 1665.}
BINS = [(.1, 1.), (1., 3.), (3., 9.5), (9.4, 20.), (20., 63.)]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batch", type=Path)
    parser.add_argument("job")
    parser.add_argument("--n", type=int, default=200000)
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()
    s = json.loads((args.batch/"fits"/args.job/"summary.json").read_text())
    model = build_df_model(model_config_from_dict(s["best"]["config"]))
    np.random.seed(args.seed)  # AGAMA sampling draws from numpy's global state
    posvel, _ = model.galaxy.sample(args.n)
    act = model.af(posvel)
    q = act[:, 0]+act[:, 1]+np.abs(act[:, 2])
    R = np.hypot(posvel[:, 0], posvel[:, 1])
    edges = np.geomspace(.05, 80., 41)
    centers = np.sqrt(edges[1:]*edges[:-1])
    quant = np.full((len(centers), 3), np.nan)
    for i in range(len(centers)):
        k = (R >= edges[i]) & (R < edges[i+1])
        if k.sum() >= 50:
            quant[i] = np.percentile(q[k], [5, 50, 95])
    bins = {f"{lo}-{hi} pc": np.percentile(q[(R > lo) & (R < hi)], [5, 50, 95]).tolist() for lo, hi in BINS}

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    keep = (R > .05) & (R < 80) & (q > 1)
    ax.hist2d(np.log10(R[keep]), np.log10(q[keep]), bins=(60, 60), cmap="Greys", cmin=1)
    ax.plot(np.log10(centers), np.log10(quant[:, 1]), color="tab:blue", lw=2, label="median q")
    ax.fill_between(np.log10(centers), np.log10(quant[:, 0]), np.log10(quant[:, 2]), color="tab:blue",
                    alpha=.2, label="5-95% of q")
    for (label, value), ls in zip(MARKS.items(), ("--", ":")):
        ax.axhline(np.log10(value), color="tab:red", ls=ls, label=label)
    for x, text in ((np.log10(9.5), "HST/MUSE | Gaia"), (np.log10(63.), "Gaia edge")):
        ax.axvline(x, color="tab:green", lw=1)
        ax.text(x, ax.get_ylim()[0]+.05, text, rotation=90, fontsize=7, va="bottom", ha="right")
    ax.set_ylim(top=max(ax.get_ylim()[1], np.log10(max(MARKS.values()))+.15))
    ax.set(xlabel=r"$\log_{10}$ projected radius [pc]", ylabel=r"$\log_{10}\,q = J_r+L$ [pc km s$^{-1}$]",
           title=f"Action coverage of the data ({args.job.removeprefix('observed_')}, {args.n} tracers)")
    ax.legend(fontsize=7, loc="upper left")
    name = "action_coverage_"+args.batch.resolve().name
    plot = ROOT/"plots"/(name+".png")
    fig.savefig(plot, dpi=160)
    record = dict(batch=str(args.batch.resolve()), job=args.job, n=args.n, seed=args.seed,
                  created_utc=datetime.now(timezone.utc).isoformat(), R_centers_pc=centers.tolist(),
                  q_quantiles_5_50_95=quant.tolist(), data_bins=bins,
                  plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/(name+".json")).write_text(json.dumps(record))
    print(json.dumps(bins, indent=1))
    print(plot)


if __name__ == "__main__":
    main()
