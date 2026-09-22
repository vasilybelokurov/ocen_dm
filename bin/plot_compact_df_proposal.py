#!/usr/bin/env python3
"""Analytic illustrations for the proposed compact DF model; no fitted results."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    plots = root / "plots"
    data_dir = root / "results" / "plot_data"
    plots.mkdir(exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 200})

    fig, ax = plt.subplots(figsize=(10.3, 4.4))
    ax.set(xlim=(0, 10), ylim=(0, 4.35))
    ax.axis("off")
    colors = ["#e0ecf5", "#eee6f4", "#faecd7"]

    def box(x, y, w, h, label, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                     boxstyle="round,pad=0.06,rounding_size=0.09",
                     facecolor=color, edgecolor="#425568", linewidth=1.1))
        ax.text(x+w/2, y+h/2, label, ha="center", va="center", fontsize=11)

    def arrow(start, end, **kw):
        ax.annotate("", xy=end, xytext=start,
                    arrowprops=dict(arrowstyle="->", color="#425568", lw=1.4, **kw))

    labels = ["Regularized stellar DF\n"+r"$f_\star(J_r,L)$",
              "Prescribed DM density\n"+r"$\rho_{\rm DM}(r)$",
              "Optional remnant density\n"+r"$\rho_{\rm rem}(r)$"]
    for x, label, color in zip([.3, 3.65, 7.0], labels, colors):
        box(x, 2.87, 2.65, .85, label, color)
        arrow((x+1.325, 2.85), (5, 2.35))
    box(2.68, 1.50, 4.64, .83,
        "Common potential and component densities\n"+r"$\Phi_{\rm tot}\;\leftrightarrow\;\rho_\star+\rho_{\rm DM}+\rho_{\rm rem}$",
        "#edf1f3")
    ax.plot([7.37, 9.85, 9.85, 1.625], [1.93, 1.93, 4.14, 4.14],
            color="#425568", lw=1.4)
    arrow((1.625, 4.14), (1.625, 3.76))
    ax.text(5.7, 3.88, "Update actions, contours, normalization and density", ha="center", fontsize=10)
    ax.text(8.72, 1.40, "Dark profiles fixed\nwithin each trial", ha="center", fontsize=9)
    arrow((5, 1.46), (5, 1.01))
    box(1.86, .14, 6.28, .85,
        "Project the observed stellar tracer\nPhotometry + LOS velocities + both PM components",
        "#e5f0e7")
    fig.subplots_adjust(left=.01, right=.99, top=.99, bottom=.02)
    fig.savefig(plots / "compact_df_proposal_architecture.png", bbox_inches="tight")
    plt.close(fig)

    # Show input functions only: these are not equilibrium or beta(r) results.
    c = np.linspace(0., 1., 401)
    q = np.geomspace(.01, 100., 601)
    palette = ["#0077aa", "#bb5500", "#7657a6"]
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 3.7), layout="constrained")
    contour_slopes = {}
    for b, color in zip([-.8, 0., .8], palette):
        slope = 2 * np.exp(-b * np.sin(np.pi*c/2))
        contour_slopes[str(b)] = slope.tolist()
        axes[0].plot(c, slope, color=color, label=rf"$b={b:g}$")
    axes[0].axhline(1., color=".4", ls="--", label=r"Linear precursor: $\eta=1$")
    axes[0].set(xlabel=r"$c=L/(L+J_r)$", ylabel=r"$g=f_{J_r}/f_L$",
                title=r"Regularity with illustrative $g_H=2$", xlim=(0, 1), ylim=(.6, 4.7))
    transition = lambda x: x**2/(1+x**2)
    b_inputs = {
        "one_radial": .8*transition(q),
        "one_tangential": -.6*transition(q),
        "two_transition": .8*transition(q/.5) + (-.6-.8)*transition(q/5.)
    }
    for name, color, label in [
        ("one_radial", palette[0], r"One: $b_{\rm out}=0.8,\ J_a=J_0$"),
        ("one_tangential", palette[1], r"One: $b_{\rm out}=-0.6,\ J_a=J_0$"),
        ("two_transition", palette[2], r"Two: $J_1=0.5J_0,\ J_2=5J_0$")]:
        axes[1].plot(q, b_inputs[name], color=color, label=label)
    axes[1].axhline(0., color=".5", lw=.8)
    axes[1].set(xscale="log", xlim=(.01, 100), ylim=(-.7, .9),
                xlabel=r"$q/J_{0,\star}$", ylabel=r"Input $b(q)$",
                title=r"Action transitions; not $\beta(r)$")
    for ax in axes:
        ax.legend(frameon=False, fontsize=8.5)
        ax.grid(alpha=.18)
    fig.savefig(plots / "compact_df_proposal_regularization.png")
    plt.close(fig)

    s = np.linspace(0, 10, 601)
    alpha_values = [.5, 1., 2.]
    ell = np.linspace(0, 4, 401)
    eta_values = [.5, 1., 2.]
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 3.7), layout="constrained")
    palette = ["#0077aa", "#bb5500", "#7657a6"]
    profiles, angular = {}, {}
    for alpha, color in zip(alpha_values, palette):
        value = np.exp(-s**alpha)
        profiles[str(alpha)] = value.tolist()
        axes[0].plot(s, value, color=color, label=rf"$\alpha={alpha:g}$")
    axes[0].set(xlabel=r"$s=(J_r+\eta L)/J_{0,\star}$", ylabel=r"$f_\star/f_\star(0,0)$",
                yscale="log", ylim=(1.e-4, 1.2), xlim=(0, 10), title="Action-space fall-off")
    for eta, color in zip(eta_values, palette):
        value = np.exp(-(1+eta*ell))
        angular[str(eta)] = value.tolist()
        axes[1].plot(ell, value, color=color, label=rf"$\eta={eta:g}$")
    axes[1].set(xlabel=r"$L/J_{0,\star}$", ylabel=r"$f_\star/f_\star(0,0)$",
                yscale="log", ylim=(1.e-4, 1.2), xlim=(0, 4),
                title=r"Orbital weighting ($\alpha=1$, $J_r/J_{0,\star}=1$)")
    for ax in axes:
        ax.legend(frameon=False)
        ax.grid(alpha=.18)
    fig.savefig(plots / "compact_df_proposal_action_shapes.png")
    plt.close(fig)
    data = dict(scope="Analytic input-function and diagnostic-precursor illustrations; no contour DF implementation, equilibrium, data fit, beta(r), or mass recovery is claimed.",
                architecture="Recommended model: regularized stellar DF plus prescribed DM and optional remnant densities. Update the stellar actions, frequency map, contour labels, normalization and density during potential iteration.",
                regularization=dict(c=c.tolist(), illustrative_gH=2., g=contour_slopes,
                                    q_over_J0=q.tolist(), n=2,
                                    b_inputs={key:value.tolist() for key,value in b_inputs.items()},
                                    two_transition=dict(b1=.8,b2=-.6,J1_over_J0=.5,J2_over_J0=5.)),
                formula="f/f(0,0) = exp(-((Jr+eta*L)/J0)**alpha)",
                scaled_action=s.tolist(), alpha_profiles=profiles,
                scaled_angular_momentum=ell.tolist(), eta_profiles=angular,
                eta_panel_fixed_alpha=1., eta_panel_fixed_Jr_over_J0=1.,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (data_dir / "compact_df_proposal.json").write_text(json.dumps(data, indent=2)+"\n")


if __name__ == "__main__":
    main()
