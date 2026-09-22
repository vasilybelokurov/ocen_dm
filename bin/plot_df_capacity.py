#!/usr/bin/env python3
"""Plot saved independent DF capacity tests; PNGs and separate JSON provenance."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
os.environ.setdefault("OMP_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ocen_dm.kinematics.df_mock import LoweredIsothermalMock, MockPopulation
from ocen_dm.kinematics.run_io import sha256

LABELS = {"single_isotropic":"Single isotropic population",
          "mixed_no_dm":"Stellar mixture + remnants",
          "mixed_with_dm":"Stellar mixture + remnants + dark halo"}
COLORS = ["#2166ac", "#e08214", "#1b7837"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, nargs="+", required=True)
    parser.add_argument("--controls", type=Path, nargs="*", default=[], help="wider-bound runs for comparison")
    parser.add_argument("--prefix", default="df_capacity_20260922")
    args = parser.parse_args()
    prefix = args.prefix
    if Path(prefix).name != prefix:
        parser.error("prefix must be a filename stem")
    plotdir, datadir = ROOT/"plots", ROOT/"results/plot_data"
    plotdir.mkdir(exist_ok=True)
    datadir.mkdir(parents=True, exist_ok=True)
    records = []
    for run in args.runs:
        summary = json.loads((run/"summary.json").read_text())
        if summary["status"] != "completed":
            raise ValueError(f"run incomplete: {run}")
        truth = json.loads((run/"truth.json").read_text())
        obs = json.loads((run/"mock_observations.json").read_text())
        profiles = json.loads((run/"data_snapshot.json").read_text())["profiles"]
        photo = json.loads((run/"photometry_template.json").read_text())
        case = summary["case"]
        record = dict(run=str(run.resolve()), case=case, summary=summary, truth=truth, observations=obs,
                      profiles=profiles, photometry=photo,
                      source_sha256={f:sha256(run/f) for f in ("summary.json", "truth.json",
                                     "mock_observations.json", "data_snapshot.json", "photometry_template.json", "run.json")})
        records.append(record)
        fig = plt.figure(figsize=(15, 10), constrained_layout=True)
        grid = fig.add_gridspec(4, 3, height_ratios=[2.3, 1, 2.3, 1])
        start = 0
        for index in range(6):
            row, col = divmod(index, 3)
            ax, residual_ax = fig.add_subplot(grid[2*row, col]), fig.add_subplot(grid[2*row+1, col])
            if index < 5:
                p = profiles[index]
                r = np.array(p["r"])
                n = len(r)
                title = p["name"].replace("_ours", "").replace("_", " ")
                ylabel = "Dispersion [km/s]" if p["kind"] == "los" else "Dispersion [mas/yr]"
            else:
                r = np.array(photo["r_arcsec"])
                n = len(r)
                title, ylabel = "Photometric shape", "Relative surface brightness [mag]"
            sl = slice(start, start+n)
            y, err = np.array(obs["truth_values"])[sl], np.array(obs["errors"])[sl]
            offset = 15-y[0] if index == 5 else 0.
            ax.errorbar(r, y+offset, yerr=err, fmt=".", color="black", ms=3,
                        alpha=.75, label="Independent mock; adopted error bars", zorder=1)
            for entry, color in zip(summary["ladder"], COLORS):
                resid = np.array(entry["residual_sigma"])[sl]
                pred = y+resid*err  # includes the profiled photometric zero point
                ndf = entry["components"]
                ax.plot(r, pred+offset, color=color, lw=1.5, label=f"{ndf} positive DF"+("s" if ndf > 1 else ""))
                residual_ax.plot(r, resid, color=color, lw=1.3)
            residual_ax.axhspan(-.3, .3, color="grey", alpha=.13)
            residual_ax.axhline(0, color="grey", lw=.6)
            residual_ax.axhline(1, color="grey", lw=.5, ls=":")
            residual_ax.axhline(-1, color="grey", lw=.5, ls=":")
            ax.set(title=title, ylabel=ylabel, xscale="log")
            residual_ax.set(xscale="log", xlabel="Projected radius [arcsec]", ylabel="Residual / error")
            ax.tick_params(labelbottom=False)
            if index == 5:
                ax.invert_yaxis()
            if index == 0:
                ax.legend(fontsize=7.5)
            start += n
        fig.suptitle(f"{LABELS[case]} — known potential, noise-free mock", fontsize=14)
        fig.savefig(plotdir/f"{prefix}_{case}.png", dpi=170)
        plt.close(fig)

    controls = {}
    for run in args.controls:
        summary = json.loads((run/"summary.json").read_text())
        if summary["status"] != "completed":
            raise ValueError(f"control incomplete: {run}")
        controls[summary["case"]] = dict(run=str(run.resolve()), summary=summary,
                                         source_sha256={f:sha256(run/f) for f in ("summary.json", "run.json")})
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), constrained_layout=True)
    markers = ["o", "s", "^"]
    for case_index, (record, marker) in enumerate(zip(records, markers)):
        entries = record["summary"]["ladder"]
        for ax, metric in zip(axes, ("score", "rms_sigma", "max_sigma")):
            y = [e["validation"]["refined_metrics"][metric] for e in entries]
            color = f"C{case_index}"
            ax.plot([e["components"] for e in entries], y, color=color, marker=marker, label=LABELS[record["case"]])
            if record["case"] in controls:
                wide = controls[record["case"]]["summary"]["ladder"]
                ax.plot([e["components"] for e in wide],
                        [e["validation"]["refined_metrics"][metric] for e in wide],
                        color=color, ls="--", marker=marker, markerfacecolor="white", alpha=.8)
            for e, v in zip(entries, y):
                if not e["validation"]["passed"]:
                    ax.scatter(e["components"], v, marker="x", color="red", s=80)
            ax.set(xlabel="Number of fitted positive DFs", xticks=[1,2,3], yscale="log")
    axes[0].set(ylabel="Q = sum of squared normalized residuals", title="Joint approximation error")
    axes[1].set(ylabel="RMS residual / error", title="Average discrepancy")
    axes[2].set(ylabel="Maximum |residual| / error", title="Worst bin")
    axes[1].axhline(.1, color="k", ls=":", label="Predefined target")
    axes[2].axhline(.3, color="k", ls=":")
    axes[0].legend(fontsize=8)
    axes[0].text(.03, .04, "Dashed / open symbols: wider search bounds", transform=axes[0].transAxes, fontsize=8)
    axes[1].legend(fontsize=8)
    fig.suptitle("Independent DF capacity challenge — doubled-resolution evaluation; no mass inference")
    fig.savefig(plotdir/f"{prefix}_summary.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(len(records), 3, figsize=(14, 3.7*len(records)), squeeze=False, constrained_layout=True)
    for row, record in enumerate(records):
        t = record["truth"]
        mock = LoweredIsothermalMock([MockPopulation(**p) for p in t["populations"]],
                                    W0=t["W0"], radius_pc=t["radius_pc"], velocity_kms=t["velocity_kms"])
        r = np.geomspace(.03, .99*mock.r_t, 220)
        components = mock.intrinsic_components(r)
        states = mock._state(r)[1:]*mock.mass_scale
        for i, pop in enumerate(mock.populations):
            label = pop.name.replace("_", " ")
            axes[row, 0].plot(r, components[:, i, 0], label=label)
            axes[row, 1].plot(r, states[i], label=label)
            if pop.light_per_mass:
                axes[row, 2].plot(r, 1-components[:, i, 2]/components[:, i, 1], label=label)
        axes[row, 0].set(xscale="log", yscale="log", ylabel="Mass density [M☉ pc⁻³]", ylim=(1e-3, 1e4), title=LABELS[record["case"]])
        axes[row, 1].set(xscale="log", yscale="log", ylabel="Enclosed mass [M☉]", ylim=(1e2, 1e7))
        axes[row, 2].set(xscale="log", ylabel="Stellar anisotropy β", ylim=(-.05, 1.05))
        for ax in axes[row]:
            ax.axvline(20, color="grey", ls=":", lw=.7)
            ax.set_xlabel("Radius [pc]")
        axes[row, 0].legend(fontsize=8)
        axes[row, 2].legend(fontsize=8)
        record["intrinsic_truth"] = dict(r_pc=r.tolist(), components_rho_pr_pt=components.tolist(), enclosed_component_masses=states.tolist())
    fig.suptitle("Known self-consistent mocks: positive energy-based DFs, distinct gravitating populations")
    fig.savefig(plotdir/f"{prefix}_truth.png", dpi=160)
    plt.close(fig)

    payload = dict(description="Fixed-potential noise-free capacity challenge; diagnostic errors, no evidence or mass constraints",
                   generator_sha256=sha256(Path(__file__)), cases=records, wider_bound_controls=controls)
    (datadir/f"{prefix}.json").write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")
    print(f"Saved {len(records)+2} PNG figures and {datadir/f'{prefix}.json'}")


if __name__ == "__main__":
    main()
