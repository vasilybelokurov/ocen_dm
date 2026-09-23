#!/usr/bin/env python3
"""Structure of a fitted compact DF: intrinsic moments, bias b(q) vs data coverage, f(J_r, L).

Panels: (a) intrinsic sigma_r and sigma_t (one component) and beta(r);
(b) the anisotropy bias b(q), q = J_r + L, over the q distributions of tracers
in the HST/MUSE (R <= 9.5 pc) and Gaia (9.4-63 pc) projected ranges, sampled
from the DF; (c) log10 f in the (J_r, L) plane, to check for artificial
ridges from the transitions. Reads saved fits; builds one refined model.

Usage: python bin/plot_df_structure.py results/df/dftwo_20260923 observed_no_halo_start0
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

from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batch", type=Path)
    parser.add_argument("job")
    parser.add_argument("--n", type=int, default=200000)
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()
    s = json.loads((args.batch/"fits"/args.job/"summary.json").read_text())
    config = refined_config(model_config_from_dict(s["best"]["config"]))
    model = build_df_model(config)
    p = config.stellar

    r = np.geomspace(.05, 80., 200)
    mom = model.intrinsic_moments(r)
    sig_r = np.sqrt(mom["radial_pressure"]/mom["rho"])
    sig_t = np.sqrt(mom["tangential_pressure"]/mom["rho"])

    np.random.seed(args.seed)
    posvel, _ = model.galaxy.sample(args.n)
    act = model.af(posvel)
    q_s = act[:, 0]+act[:, 1]+np.abs(act[:, 2])
    R = np.hypot(posvel[:, 0], posvel[:, 1])
    inner, outer = q_s[(R > .1) & (R <= 9.5)], q_s[(R > 9.4) & (R < 63)]

    q = np.geomspace(1., 5000., 400)
    bias = model.df.parameters.bias(q)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    ax = axes[0]
    ax.plot(r, sig_r, label=r"$\sigma_r$")
    ax.plot(r, sig_t, label=r"$\sigma_t$ (one component)")
    ax.set(xscale="log", xlabel="r [pc]", ylabel="intrinsic dispersion [km/s]", title="(a) intrinsic moments")
    tw = ax.twinx()
    tw.plot(r, mom["beta"], color="black", ls="--", label=r"$\beta$")
    tw.axhline(0, color="grey", lw=.5)
    tw.set_ylabel(r"$\beta=1-\sigma_t^2/\sigma_r^2$")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = tw.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=7, loc="lower left")

    ax = axes[1]
    edges = np.geomspace(5., 3000., 60)
    ax.hist(inner, bins=edges, density=True, alpha=.4, color="tab:blue", label="tracers, R <= 9.5 pc (HST/MUSE)")
    ax.hist(outer, bins=edges, density=True, alpha=.4, color="tab:green", label="tracers, 9.4-63 pc (Gaia)")
    ax.set(xscale="log", xlabel=r"$q=J_r+L$ [pc km/s]", ylabel="p(q)", title="(b) bias b(q) vs sampled actions")
    tw = ax.twinx()
    tw.plot(q, bias, color="tab:red", lw=2, label="b(q)")
    for v, lab in ((p.J_a, r"$J_a$"), (p.J_outer, r"$J_{\rm outer}$")):
        if v:
            tw.axvline(v, color="tab:red", ls=":", lw=1)
            tw.text(v, 0, lab, color="tab:red", fontsize=8, rotation=90, va="bottom")
    tw.axhline(0, color="grey", lw=.5)
    tw.set_ylabel("b(q)  (b > 0: radial bias)")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = tw.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=7, loc="upper left")

    ax = axes[2]
    jr, L = np.meshgrid(np.geomspace(1., 3000., 160), np.geomspace(1., 3000., 160))
    grid = np.column_stack((jr.ravel(), np.zeros(jr.size), L.ravel()))  # J_z = 0, J_phi = L
    f = np.asarray(model.df(grid)).reshape(jr.shape)
    with np.errstate(divide="ignore"):
        logf = np.log10(np.where(f > 0, f, np.nan))
    top = np.nanmax(logf)
    cs = ax.contourf(np.log10(jr), np.log10(L), logf, levels=np.linspace(top-12, top, 25), cmap="viridis")
    ax.contour(np.log10(jr), np.log10(L), logf, levels=np.linspace(top-12, top, 13), colors="white", linewidths=.4)
    for v in (p.J_a, p.J_outer):
        if v:  # the line q = J_r + L = v, drawn where both actions exceed 1 pc km/s
            jline = np.geomspace(1., v-1., 200)
            ax.plot(np.log10(jline), np.log10(v-jline), color="tab:red", ls=":", lw=1.5)
    ax.set(xlim=(0, np.log10(3000.)), ylim=(0, np.log10(3000.)))
    fig.colorbar(cs, ax=ax, label=r"$\log_{10} f$")
    ax.set(xlabel=r"$\log_{10} J_r$ [pc km/s]", ylabel=r"$\log_{10} L$ [pc km/s]",
           title=r"(c) $f(J_r, L)$; dotted: $q = J_a, J_{\rm outer}$")
    fig.suptitle(f"Compact DF structure: {args.batch.resolve().name} / {args.job}")
    name = f"df_structure_{args.batch.resolve().name}_{args.job}"
    plot = ROOT/"plots"/(name+".png")
    fig.savefig(plot, dpi=160)
    record = dict(batch=str(args.batch.resolve()), job=args.job, created_utc=datetime.now(timezone.utc).isoformat(),
                  r_pc=r.tolist(), sigma_r=sig_r.tolist(), sigma_t=sig_t.tolist(), beta=mom["beta"].tolist(),
                  q=q.tolist(), bias=bias.tolist(),
                  q_inner_5_50_95=np.percentile(inner, [5, 50, 95]).tolist(),
                  q_outer_5_50_95=np.percentile(outer, [5, 50, 95]).tolist(),
                  plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/(name+".json")).write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
