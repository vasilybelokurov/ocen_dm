#!/usr/bin/env python3
"""Plot observed Omega Cen data against compact-DF fits from an observed-data batch.

Rebuilds each finished fit's selected model at refined resolution (about a
minute per fit) and plots surface brightness, LOS and PM dispersions with the
adopted errors, plus normalized residuals. This is a fit-quality figure; it
does not display or certify a stellar/dark mass decomposition.

Usage: python bin/plot_observed_df_fits.py results/df/observed_df_20260923
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
from ocen_dm.kinematics.df_fit import DFJointProblem, PhotometricData, model_config_from_dict
from ocen_dm.kinematics.run_io import read_data_snapshot

PANELS = [("phot", r"$\mu$ [mag]", "surface brightness"), ("los", r"$\sigma_{\rm los}$ [km s$^{-1}$]", "line of sight"),
          ("pmr", r"$\sigma_{\rm pm,R}$ [mas yr$^{-1}$]", "PM radial"),
          ("pmt", r"$\sigma_{\rm pm,T}$ [mas yr$^{-1}$]", "PM tangential")]


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batch", type=Path)
    parser.add_argument("--jobs", nargs="*", help="plot only these job ids (default: all finished)")
    parser.add_argument("--suffix", default="_fits", help="output name suffix")
    args = parser.parse_args()
    out = args.batch.resolve()
    manifest = read(out/"batch.json")
    d = out/"observed"
    problem = DFJointProblem(read_data_snapshot(d, read(d/"problem.json")["data_snapshot"]),
                             PhotometricData.from_dict(read(d/"photometry.json")))
    pc = manifest["fixed"]["distance_kpc"]*1e3*np.pi/(180*3600)
    fits = []
    for job in manifest["jobs"]:
        if args.jobs and job["id"] not in args.jobs:
            continue
        f = out/"fits"/job["id"]/"summary.json"
        if f.exists() and "best" in read(f):  # infeasible starts have no model
            s = read(f)
            ev = problem.evaluate(refined_config(model_config_from_dict(s["best"]["config"])))
            fits.append((job["id"].removeprefix("observed_"), s, ev))
            print(job["id"], "chi2_kin %.1f" % ev["chi2_kinematic"], flush=True)
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5), sharex="col", constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[3, 1.4]))
    record = dict(batch=str(out), created_utc=datetime.now(timezone.utc).isoformat(),
                  source_commit=manifest["source_commit"], fits={})
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for col, (kind, ylabel, title) in enumerate(PANELS):
        ax, rx = axes[0, col], axes[1, col]
        if kind == "phot":
            ph = problem.photometry
            ax.errorbar(ph.r_arcsec*pc, ph.mu, yerr=ph.sigma_mag, fmt=".", color="0.4", ms=3, lw=.6, label="observed")
            ax.invert_yaxis()
        for p in (p for p in problem.data.profiles if p.kind == kind):
            ax.errorbar(p.r*pc, p.value, yerr=[p.err_lo, p.err_hi], fmt=".", color="0.4", ms=3, lw=.6,
                        label="observed" if p is next(q for q in problem.data.profiles if q.kind == kind) else None)
        for i, (name, s, ev) in enumerate(fits):
            gates = s["gates"]
            label = f"{name}: $\\chi^2_{{kin}}/N$={gates['chi2_kin_per_point']:.2f}"
            if kind == "phot":
                ph = problem.photometry
                order = np.argsort(ph.r_arcsec)
                ax.plot(ph.r_arcsec[order]*pc, ev["photometry"]["prediction"][order], color=colors[i], lw=1, label=label)
                rx.plot(ph.r_arcsec*pc, ph.mu-ev["photometry"]["prediction"], ".", color=colors[i], ms=3)
                rx.set_ylabel("data - model [mag]")
                record["fits"].setdefault(name, {})["photometry"] = ev["photometry"]["prediction"].tolist()
                continue
            for j, p in enumerate(q for q in problem.data.profiles if q.kind == kind):
                m = ev["predictions"][p.name]
                err = np.where(m > p.value, p.err_hi, p.err_lo)
                ax.plot(p.r*pc, m, color=colors[i], lw=1, label=label if j == 0 else None)
                rx.plot(p.r*pc, (p.value-m)/err, ".-", color=colors[i], ms=3, lw=.5)
                record["fits"].setdefault(name, {})[p.name] = m.tolist()
            rx.set_ylabel("(data - model) / error")
        rx.axhline(0, color="black", lw=.8)
        ax.set(xscale="log", ylabel=ylabel, title=title)
        rx.set(xscale="log", xlabel="Projected radius [pc]")
        ax.legend(fontsize=6)
    fig.suptitle(f"Compact regularized DF fits to observed Omega Cen data ({out.name}); adopted errors")
    plot = ROOT/"plots"/(out.name+args.suffix+".png")
    fig.savefig(plot, dpi=160)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/(out.name+args.suffix+".json")).write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
