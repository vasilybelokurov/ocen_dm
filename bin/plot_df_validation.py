#!/usr/bin/env python3
"""Plot the density and anisotropy of the four unfitted DF validation examples."""
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    records = []
    families = [("no_dm", "Single stellar DF; no DM"),
                ("cored_dm", "Single stellar DF; cored DM"),
                ("cusped_dm", "Single stellar DF; cusped DM"),
                ("two_component", "Two positive stellar DFs; no DM")]
    for index, (key, label) in enumerate(families):
        name = f"validation_{key}_20260922"
        source = ROOT/"results/plot_data"/f"df_{name}.json"
        row = json.loads(source.read_text())
        summary = json.loads((ROOT/"results/df"/name/"summary.json").read_text())
        r = np.array(row["r_pc"])
        keep = r <= 100
        density = np.array(row["intrinsic"]["rho"])/summary["model"]["M_star"]
        beta = np.array(row["intrinsic"]["beta"])
        style = dict(color=f"C{index}", lw=2.5 if index == 3 else 1.7, label=label)
        axes[0].loglog(r[keep], density[keep], **style)
        axes[1].semilogx(r[keep], beta[keep], **style)
        records.append(dict(label=label, source=str(source.relative_to(ROOT)),
                            r_pc=r[keep].tolist(), density_per_stellar_mass=density[keep].tolist(),
                            beta=beta[keep].tolist(), numerical_validation=summary["validation"]))
    axes[0].set(xlabel="Radius [pc]", ylabel="Stellar density / total stellar mass [pc⁻³]")
    axes[1].set(xlabel="Radius [pc]", ylabel="Anisotropy β(r)")
    axes[1].axhline(0, color="0.6", ls=":", lw=1)
    axes[1].text(.04, .94, "β > 0: radial preference\nβ < 0: tangential preference",
                 transform=axes[1].transAxes, va="top", fontsize=9)
    axes[0].legend(fontsize=8, loc="lower left")
    for ax in axes:
        ax.grid(alpha=.15)
    fig.suptitle("Positive AGAMA DFs: example equilibria, not fits to ω Cen")
    fig.savefig(ROOT/"plots/df_validation_profiles.png", dpi=180)
    plt.close(fig)
    (ROOT/"results/plot_data/df_validation_profiles.json").write_text(
        json.dumps(dict(description="Unfitted code-validation examples", models=records),
                   indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
