"""Data-versus-model figures for the Jeans fits."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ..kinematics.fit import FitProblem
from ..kinematics.likelihood import ARCSEC_PER_RAD
from .style import apply as apply_style, SERIES, COLOR_FIELD as FIELD_GREY

_LABEL = {"los": r"$\sigma_{\rm LOS}$ [km/s]", "pmr": r"$\sigma_{\mu,R}$ [mas/yr]",
          "pmt": r"$\sigma_{\mu,T}$ [mas/yr]", "pmc": r"$\sigma_{\mu}$ [mas/yr]"}


def plot_profile_fit(problem: FitProblem, x: np.ndarray, path: Path, title: str = "",
                     samples: np.ndarray | None = None, n_band: int = 100) -> Path:
    """One panel per dataset: data with asymmetric errors, model curve, residuals.

    ``samples`` (optional, shape (n, n_params)) adds a 16-84 per cent band from a
    random subset of posterior samples.
    """
    apply_style()
    fam, like = problem.family, problem.likelihood
    profiles = problem.data.profiles
    n = len(profiles)
    fig, axes = plt.subplots(2, n, figsize=(3.6 * n, 6.2), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1.2], "hspace": 0.05, "wspace": 0.3})
    axes = np.atleast_2d(axes).reshape(2, n)
    theta = fam.to_dict(x)
    jeans, D, scales = fam.build(theta)
    pc_per_arcsec = D * 1e3 / ARCSEC_PER_RAD
    pred = like.predict(jeans, D, scales)

    rng = np.random.default_rng(0)
    sub = None
    if samples is not None and len(samples):
        sub = samples[rng.choice(len(samples), size=min(n_band, len(samples)), replace=False)]

    for i, p in enumerate(profiles):
        ax, rax = axes[0, i], axes[1, i]
        # smooth curve: point-evaluated model over the data range
        R = np.geomspace(0.8 * p.r.min(), 1.2 * p.r.max(), 80)
        m = jeans.dispersions_kms(R * pc_per_arcsec)
        s2 = 0.5 * (m["pmr"] ** 2 + m["pmt"] ** 2) if p.kind == "pmc" else m[p.kind] ** 2
        curve = np.sqrt(s2) * (1.0 if p.kind == "los" else 1.0 / (4.740470463533348 * D)) * scales.get(p.instrument, 1.0)
        if sub is not None:
            band = []
            for xs in sub:
                js, Ds, scs = fam.build(fam.to_dict(xs))
                ms = js.dispersions_kms(R * Ds * 1e3 / ARCSEC_PER_RAD)
                s2s = 0.5 * (ms["pmr"] ** 2 + ms["pmt"] ** 2) if p.kind == "pmc" else ms[p.kind] ** 2
                band.append(np.sqrt(s2s) * (1.0 if p.kind == "los" else 1.0 / (4.740470463533348 * Ds)) * scs.get(p.instrument, 1.0))
            lo, hi = np.percentile(band, [16, 84], axis=0)
            ax.fill_between(R, lo, hi, color=SERIES[0], alpha=0.25, lw=0)
        ax.plot(R, curve, color=SERIES[0], lw=2, label="model (curve)")
        xerr = None if not p.has_edges else [p.r - p.r_lower, p.r_upper - p.r]
        ax.errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], xerr=xerr, fmt="o", ms=4, color="k",
                    ecolor=FIELD_GREY, elinewidth=1, capsize=0, label="data")
        ax.plot(p.r, pred[p.name], "_", color=SERIES[1], ms=10, mew=2, label="model (bin mean)")
        ax.set_xscale("log"); ax.set_ylabel(_LABEL[p.kind])
        ax.set_title(f"{p.name}  ({p.instrument}, n={p.n})", fontsize=10)
        err = np.where(pred[p.name] > p.value, p.err_hi, p.err_lo)
        res = (p.value - pred[p.name]) / err
        rax.axhline(0, color=FIELD_GREY, lw=1); rax.axhspan(-1, 1, color=FIELD_GREY, alpha=0.2, lw=0)
        rax.plot(p.r, res, "o", ms=4, color="k")
        rax.set_ylim(-max(4, 1.1 * np.abs(res).max()), max(4, 1.1 * np.abs(res).max()))
        rax.set_xlabel("R [arcsec]"); rax.set_ylabel(r"(data$-$model)/$\sigma$")
        rax.text(0.03, 0.85, r"$\chi^2$=%.0f" % np.sum(res ** 2), transform=rax.transAxes, fontsize=9)
        if i == 0:
            ax.legend(fontsize=8, loc="lower left")
    if title:
        fig.suptitle(title, y=0.995)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_posterior_profiles(run_dirs: dict[str, Path], path: Path, title: str = "") -> Path:
    """Enclosed-mass, dark fraction and circular-speed bands from ``profiles.npz`` files.

    Parameters
    ----------
    run_dirs : dict
        ``label -> results/fits/<label>`` for each run to overlay (K1, K2 ...).
    """
    apply_style()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for i, (label, d) in enumerate(run_dirs.items()):
        prof = np.load(Path(d) / "profiles.npz")
        r = prof["r"]; color = SERIES[i % len(SERIES)]
        for ax, key, ylabel, log in ((axes[0], "M_total", r"$M(<r)$ [M$_\odot$]", True),
                                     (axes[1], "f_dm", r"$f_{\rm DM}(<r)$", False),
                                     (axes[2], "v_circ", r"$v_{\rm circ}$ [km/s]", False)):
            lo, mid, hi = np.percentile(prof[key], [16, 50, 84], axis=0)
            ax.fill_between(r, lo, hi, color=color, alpha=0.25, lw=0)
            ax.plot(r, mid, color=color, lw=2, label=label)
            ax.set_xscale("log"); ax.set_xlabel("r [pc]"); ax.set_ylabel(ylabel)
            if log:
                ax.set_yscale("log")
        if "M_dm" in prof.files:
            lo, mid, hi = np.percentile(prof["M_dm"], [16, 50, 84], axis=0)
            axes[0].fill_between(r, np.maximum(lo, 1), hi, color=color, alpha=0.12, lw=0, hatch="//")
            axes[0].plot(r, np.maximum(mid, 1), color=color, lw=1.2, ls="--", label=f"{label}: DM only")
    axes[0].set_ylim(1e3, None); axes[0].legend(fontsize=8)
    for R_hst, R_gaia in ((9.0, 53.0),):                      # data extents in pc at 5.43 kpc
        for ax in axes:
            ax.axvline(R_hst, color=FIELD_GREY, lw=1, ls=":"); ax.axvline(R_gaia, color=FIELD_GREY, lw=1, ls=":")
    axes[1].text(9.2, 0.95, "HST/MUSE edge", fontsize=8, color=FIELD_GREY, transform=axes[1].get_xaxis_transform())
    axes[1].text(55, 0.95, "Gaia edge", fontsize=8, color=FIELD_GREY, transform=axes[1].get_xaxis_transform())
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path
