#!/usr/bin/env python3
"""Star-count tracer profiles (Poisson) against compact-DF model predictions.

Left: HST F625W<19 and Gaia G<17 surface number densities with Poisson errors,
and each model's prediction mu_i/A_i (amplitude and field profiled per dataset).
Right: signed Poisson deviance residuals per bin (squared sum = deviance).
Models are given as BATCH::JOB::LABEL and rebuilt at refined resolution.

Usage: python bin/plot_count_profiles.py --name counts_vs_models_20260923 \
    "results/df/dftwo_20260923::observed_no_halo_start0::two transitions, no halo"
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
from ocen_dm.kinematics.counts import CountProfile
from ocen_dm.kinematics.df_fit import DFJointProblem, build_df_model, model_config_from_dict
from ocen_dm.kinematics.run_io import read_data_snapshot

PRODUCTS = ("ocen_counts_hst_f625w19", "ocen_counts_gaia_g17")
PC = 5.43e3*np.pi/(180*3600)


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--name", required=True)
    parser.add_argument("--data", default="results/df/dftwo_20260923/observed", help="directory with the kinematic snapshot")
    parser.add_argument("specs", nargs="*", help="BATCH::JOB::LABEL")
    args = parser.parse_args()
    counts = [CountProfile.from_dict(read(ROOT/"data/processed/kinematics"/f"{n}.json")) for n in PRODUCTS]
    d = ROOT/args.data
    data = read_data_snapshot(d, read(d/"problem.json")["data_snapshot"])
    problem = DFJointProblem(data, None, counts)
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[3, 1.5]))
    ax, rx = axes
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), products=PRODUCTS, models={})
    for c, marker in zip(counts, ("o", "s")):
        dens = c.counts/c.area_arcsec2*3600.
        err = np.sqrt(np.maximum(c.counts, 1))/c.area_arcsec2*3600.
        ax.errorbar(c.r_median*PC, dens, err, fmt=marker, color="black", ms=4, mfc="white", label=f"{c.name} (N={int(c.counts.sum())})")
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, spec in enumerate(args.specs):
        batch, job, label = spec.split("::")
        s = read(ROOT/batch/"fits"/job/"summary.json")
        config = refined_config(model_config_from_dict(s["best"]["config"]))
        ev = problem.evaluate(config, build_df_model(config))
        text = []
        for c, ls in zip(counts, ("-", "--")):
            out = ev["counts"][c.name]
            ax.plot(c.r_median*PC, out["mu"]/c.area_arcsec2*3600., ls=ls, color=colors[i], lw=1.5,
                    label=f"{label}: {c.name.split('_')[2]} deviance {out['deviance']:.0f}/{c.n}"
                          + (f", field {out['field_density_per_arcmin2']:.2f}/arcmin$^2$" if c.fit_field else ""))
            rx.plot(c.r_median*PC, out["residual"], ls=ls, marker=".", color=colors[i], lw=.8)
            text.append(f"{c.name}: deviance {out['deviance']:.1f}/{c.n}")
            record["models"].setdefault(label, {})[c.name] = dict(mu=out["mu"].tolist(), residual=np.asarray(out["residual"]).tolist(),
                                                                  deviance=out["deviance"], amplitude=out["amplitude"], field=out["field"])
        print(label, "; ".join(text), flush=True)
    for y in (-2, 0, 2):
        rx.axhline(y, color="black" if y == 0 else "0.7", lw=.7)
    ax.set(xscale="log", yscale="log", ylabel=r"stars per arcmin$^2$", title="Poisson star-count tracer profiles vs models")
    rx.set(xscale="log", xlabel="projected radius [pc]", ylabel="deviance residual")
    ax.legend(fontsize=7)
    plot = ROOT/"plots"/(args.name+".png")
    fig.savefig(plot, dpi=150)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/(args.name+".json")).write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
