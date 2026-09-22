#!/usr/bin/env python3
"""Reproduce the numerical-validation and optimizer plots in agama_df_models.tex.

Reads saved results only. PNGs go to plots/, plotted arrays/provenance to
results/plot_data/. No equilibrium builds or fits are launched.
"""
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-nbody-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PLOTS = ROOT/"plots"
TABLES = ROOT/"results/plot_data"


def read(path):
    return json.loads(path.read_text())


def main():
    plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "legend.fontsize": 9,
                         "savefig.dpi": 200, "axes.spines.top": False,
                         "axes.spines.right": False})
    names = ["validation_no_dm_20260922", "validation_cored_dm_20260922",
             "validation_cusped_dm_20260922", "validation_two_component_20260922",
             "pilot_no_dm_20260922"]
    labels = ["Single DF\nNo DM", "Single DF\nCore", "Single DF\nCusp",
              "Two DFs\nNo DM", "Single DF\nPilot"]
    sources = [ROOT/"results/df"/name/"summary.json" for name in names]
    summaries = [read(p) for p in sources]
    metrics = {
        "density_closure": [d["diagnostics"]["density_closure"] for d in summaries],
        "mass_closure": [d["diagnostics"]["mass_error"] for d in summaries],
        "direct_projection": [max(d["validation"]["direct_projection_relative_error"].values())
                              for d in summaries],
        "dispersion_shift_in_data_sigma": [max(d["validation"]["refinement_max_shift_in_data_sigma"].values())
                                           for d in summaries],
    }
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2), constrained_layout=True)
    for j, (key, label) in enumerate([
            ("density_closure", "DF density closure"),
            ("mass_closure", "Stellar mass closure"),
            ("direct_projection", "Direct AGAMA projection")]):
        axes[0].semilogy(x+(j-1)*.16, np.array(metrics[key])*100, "os^"[j],
                         ms=7, label=label)
    axes[0].axhline(.5, color="0.4", ls="--", lw=1.2, label="Acceptance: 0.5%")
    axes[0].set(ylabel="Maximum fractional discrepancy [%]", ylim=(8e-5, 2),
                title="Density, gravity and projection checks")
    axes[0].legend(loc="center", bbox_to_anchor=(.52, .40), ncol=2,
                   fontsize=8, framealpha=.95)
    axes[1].semilogy(x, metrics["dispersion_shift_in_data_sigma"], "o", ms=8, color="C3",
                     label="Largest shift among 89 bins")
    axes[1].axhline(.1, color="0.4", ls="--", lw=1.2, label="Acceptance: 0.1 data errors")
    axes[1].set(ylabel="Dispersion shift / observational error", ylim=(2e-4, .25),
                title="Complete model at doubled resolution")
    axes[1].legend(loc="lower right", framealpha=.95)
    for ax in axes:
        ax.set_xticks(x, labels)
        ax.grid(axis="y", alpha=.18)
    fig.savefig(PLOTS/"df_numerical_validation.png")
    plt.close(fig)

    pilot = summaries[-1]
    history_path = ROOT/"results/df/pilot_no_dm_20260922/evaluations.jsonl"
    history = [json.loads(line) for line in history_path.read_text().splitlines()]
    calls = [0]
    q = [pilot["initial_objective"]]
    kin = [pilot["initial_chi2_kinematic"]]
    photo = [q[0]-kin[0]]
    best = q[0]
    rejected = []
    for row in history:
        calls.append(row["evaluation"])
        if "rejected" in row:
            rejected.append(dict(evaluation=row["evaluation"], reason=row["rejected"]))
        if row.get("objective", np.inf) < best:
            best = row["objective"]
            q.append(best)
            kin.append(row["chi2_kinematic"])
            photo.append(row["chi2_photometric"])
        else:
            q.append(q[-1]); kin.append(kin[-1]); photo.append(photo[-1])
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.1), constrained_layout=True)
    finite = [r for r in history if "objective" in r]
    axes[0].semilogy([r["evaluation"] for r in finite], [r["objective"] for r in finite],
                     ".", color="0.7", ms=5, label="Evaluated models")
    axes[0].step(calls, q, where="post", color="C0", lw=2, label="Best joint objective so far")
    for row in rejected:
        axes[0].axvline(row["evaluation"], ymin=.0, ymax=.035, color="C3", lw=2)
    axes[0].set(ylabel=r"$Q=\chi^2_{\rm kin}+\chi^2_{\rm phot}$", title="Bounded local optimization")
    axes[0].legend(loc="upper right")
    axes[1].step(calls, kin, where="post", label=r"$\chi^2_{\rm kin}$", color="C1", lw=2)
    axes[1].step(calls, photo, where="post", label=r"$\chi^2_{\rm phot}$", color="C2", lw=2)
    axes[1].set(yscale="log", ylabel=r"Contribution to $Q$",
                title="Contributions of the best joint model")
    axes[1].legend(loc="upper right")
    for ax in axes:
        ax.set(xlabel="Objective evaluation (0 = starting model)", xlim=(0, 99))
        ax.axvline(96, color="0.4", ls="--", lw=1)
        ax.grid(alpha=.18)
    fig.suptitle("No-DM pilot: 96-evaluation limit reached; optimizer unconverged", fontsize=12)
    fig.savefig(PLOTS/"df_pilot_optimization.png")
    plt.close(fig)
    sources.append(history_path)
    record = dict(
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sources},
        models=names, labels=labels, numerical_metrics=metrics,
        optimization=dict(evaluation=calls, best_joint_objective=q,
                          best_joint_kinematic_contribution=kin,
                          best_joint_photometric_contribution=photo, rejected=rejected))
    (TABLES/"df_writeup_diagnostics.json").write_text(json.dumps(record, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
