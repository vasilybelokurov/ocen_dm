"""Where the kinematic constraints come from, and whether the outer tracers are trustworthy.

Two figures:

``plot_constraint_map`` -- every dispersion dataset on one axis in km/s, the K1 (no dark
matter) and K2 (dark halo) model curves over them, the model ratio and the data residuals
underneath, and the number of tracer stars per bin. It answers, directly from the data,
over what radii the two families are distinguishable at all.

``plot_outer_tracer_audit`` -- the audit of the outer Gaia sample: how many secure members
there are at each radius, how their per-star errors compare with the dispersion they are
supposed to measure, and how the measured dispersion moves under membership, quality,
magnitude and error-calibration cuts.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from ..kinematics.likelihood import BinnedProfile, KinematicData, ProfileLikelihood
from ..kinematics.outer_profile import KMS_PER_MASYR_KPC, binned_dispersion, load_members
from ..kinematics.report import _family_for
from ..paths import processed_dir, results_dir
from . import style
from .style import add_pc_axis

__all__ = ["plot_constraint_map", "plot_outer_tracer_audit", "our_outer_profile", "OUTER_EDGES"]

#: log-spaced annuli for our own outer measurement (arcsec)
OUTER_EDGES = np.geomspace(300.0, 2400.0, 11)
#: Vasiliev & Baumgardt quality bit: the stars their own profiles are built from
QUALITY_BIT = 2


def our_outer_profile(edges: np.ndarray = OUTER_EDGES, prob_min: float = 0.9,
                      quality_mask: int | None = QUALITY_BIT, **kwargs) -> Table:
    """Our error-deconvolved PM dispersion profile from the member catalogue."""
    return binned_dispersion(load_members(), edges, prob_min=prob_min, quality_mask=quality_mask, **kwargs)


def _model(label: str):
    """``(jeans, distance_kpc, scales, family)`` of a finished run's best sample."""
    s = json.loads((results_dir() / "fits" / label / "summary.json").read_text())
    fam = _family_for(s)
    x = np.array([s["parameters"][n]["ml"] for n in fam.names])
    jeans, D, scales = fam.build(fam.to_dict(x))
    return jeans, D, scales, fam


def _sigma_1d_kms(jeans, D: float, R_arcsec: np.ndarray, kind: str = "pm") -> np.ndarray:
    m = jeans.dispersions_kms(np.asarray(R_arcsec, float) * D * 1e3 / 206264.806)
    return np.sqrt(0.5 * (m["pmr"] ** 2 + m["pmt"] ** 2)) if kind == "pm" else m["los"]


def _datasets_in_kms(D: float) -> list[dict]:
    """Every dispersion dataset converted to a 1-D dispersion in km/s."""
    k = KMS_PER_MASYR_KPC * D
    out = []
    hst_r = Table.read(processed_dir() / "kinematics" / "omegacat_vi_pm_radial.ecsv")
    hst_t = Table.read(processed_dir() / "kinematics" / "omegacat_vi_pm_tangential.ecsv")
    sig = np.sqrt(0.5 * (np.asarray(hst_r["sigma_pmr"]) ** 2 + np.asarray(hst_t["sigma_pmt"]) ** 2)) * k
    err = 0.5 * (np.asarray(hst_r["sigma_pmr_err_lo"]) + np.asarray(hst_r["sigma_pmr_err_hi"])) * k / np.sqrt(2)
    out.append(dict(name="HST PM (oMEGACat)", r=np.asarray(hst_r["r_median"]), sigma=sig, err=err,
                    n=np.asarray(hst_r["n_stars"]), kind="pm", color=style.SERIES[0], marker="o"))
    los = Table.read(processed_dir() / "kinematics" / "omegacat_vi_los_dispersion.ecsv")
    out.append(dict(name="MUSE line of sight", r=np.asarray(los["r_median"]), sigma=np.asarray(los["sigma_los"]),
                    err=0.5 * (np.asarray(los["sigma_los_err_lo"]) + np.asarray(los["sigma_los_err_hi"])),
                    n=np.asarray(los["n_stars"]), kind="los", color=style.SERIES[1], marker="s"))
    dr2 = Table.read(processed_dir() / "kinematics" / "baumgardt2019_ocen_pm_dispersion.ecsv")
    out.append(dict(name="Gaia DR2 (Baumgardt+ 2019)", r=np.asarray(dr2["r"]), sigma=np.asarray(dr2["sigma_pm"]) * k,
                    err=0.5 * (np.asarray(dr2["sigma_pm_err_lo"]) + np.asarray(dr2["sigma_pm_err_hi"])) * k,
                    n=None, kind="pm", color=style.SERIES[2], marker="^"))
    ours = our_outer_profile()
    out.append(dict(name="Gaia EDR3, our measurement", r=np.asarray(ours["r_median"]),
                    sigma=np.asarray(ours["sigma_pm"]) * k, err=np.asarray(ours["sigma_pmr_err"]) * k,
                    n=np.asarray(ours["n_stars"]), kind="pm", color=style.INK, marker="D",
                    edges=(np.asarray(ours["r_lower"]), np.asarray(ours["r_upper"])), table=ours))
    return out


def plot_constraint_map(path: Path | str = "plots/constraint_map.png",
                        k1: str = "K1_noDM_composite", k2: str = "K2_cored_composite") -> Path:
    """All dispersion data in km/s with the K1 and K2 curves, their ratio, and the tracer counts."""
    style.apply()
    j1, D, s1, _ = _model(k1)
    j2, D2, s2, _ = _model(k2)
    data = _datasets_in_kms(D)
    R = np.geomspace(2.0, 2400.0, 220)

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 11), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.5, 1.4], "hspace": 0.07})
    ax, axr, axn = axes

    for kind, ls, lab in (("pm", "-", "proper motion"), ("los", "--", "line of sight")):
        ax.plot(R, _sigma_1d_kms(j1, D, R, kind), color=style.SERIES[0], lw=2, ls=ls, alpha=0.9,
                label=f"K1 no dark matter, {lab}")
        ax.plot(R, _sigma_1d_kms(j2, D2, R, kind), color=style.SERIES[1], lw=2, ls=ls, alpha=0.9,
                label=f"K2 cored halo, {lab}")
    for d in data:
        ax.errorbar(d["r"], d["sigma"], yerr=d["err"], fmt=d["marker"], ms=4, color=d["color"],
                    ecolor=d["color"], elinewidth=1, capsize=0, lw=0, label=d["name"], zorder=5)
    ax.set_yscale("log"); ax.set_ylabel("1-D velocity dispersion  [km/s]")
    ax.legend(fontsize=8, ncol=2, loc="lower left")
    ax.set_title("What constrains the mass, and where", fontsize=11)
    add_pc_axis(ax, D)

    # residuals of every dataset against K1, and the K2/K1 model ratio
    ratio = _sigma_1d_kms(j2, D2, R) / _sigma_1d_kms(j1, D, R)
    axr.plot(R, 100 * (ratio - 1), color=style.SERIES[1], lw=2, label="K2 / K1 model (proper motion)")
    axr.axhline(0, color=style.INK_SECONDARY, lw=0.8)
    for d in data:
        model = _sigma_1d_kms(j1, D, d["r"], d["kind"])
        scale = s1.get("MUSE" if d["kind"] == "los" else "GaiaDR2" if "DR2" in d["name"] else "GaiaEDR3", 1.0) \
            if d["name"] != "HST PM (oMEGACat)" else 1.0
        axr.errorbar(d["r"], 100 * (d["sigma"] / (model * scale) - 1), yerr=100 * d["err"] / (model * scale),
                     fmt=d["marker"], ms=4, color=d["color"], ecolor=d["color"], elinewidth=1, capsize=0, lw=0)
    axr.set_ylabel("deviation from K1  [%]\n(instrument scales applied)")
    axr.set_ylim(-25, 40); axr.legend(fontsize=8, loc="upper left")

    for d in data:
        if d["n"] is None:
            continue
        axn.step(d["r"], d["n"], where="mid", color=d["color"], lw=1.6, label=d["name"])
    axn.set_yscale("log"); axn.set_ylabel("tracer stars per bin")
    axn.set_xlabel("R  [arcsec]"); axn.set_xscale("log")
    axn.legend(fontsize=8, loc="upper left")

    for a in axes:
        a.axvspan(1250, 2400, color=style.SERIES[1], alpha=0.06, lw=0)
    axr.text(1300, 32, "K1 and K2 differ by more than 3 %", fontsize=8, color=style.SERIES[1])
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_outer_tracer_audit(path: Path | str = "plots/outer_tracer_audit.png") -> Path:
    """Are the distant tracers secure, and is their dispersion resolved rather than deconvolved?"""
    style.apply()
    sample = load_members()
    edges = OUTER_EDGES
    variants = [
        ("P > 0.9, quality flag (reference)", dict(prob_min=0.9, quality_mask=QUALITY_BIT), style.INK, "o"),
        ("P > 0.5", dict(prob_min=0.5, quality_mask=QUALITY_BIT), style.SERIES[0], "v"),
        ("P > 0.99", dict(prob_min=0.99, quality_mask=QUALITY_BIT), style.SERIES[2], "^"),
        ("G < 18.5 (errors 0.13 mas/yr)", dict(prob_min=0.9, quality_mask=QUALITY_BIT, g_range=(0, 18.5)),
         style.SERIES[1], "s"),
        ("no quality cut", dict(prob_min=0.9), style.COLOR_FIELD, "x"),
        ("errors inflated 20 %", dict(prob_min=0.9, quality_mask=QUALITY_BIT, err_scale=1.2),
         style.SERIES_EXTRA, "d"),
    ]
    tables = [(lab, binned_dispersion(sample, edges, **kw), c, m) for lab, kw, c, m in variants]
    ref = tables[0][1]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    (a_n, a_err), (a_sig, a_rat) = axes

    for i, (lab, t, c, m) in enumerate(tables):
        a_n.step(t["r_median"], t["n_stars"], where="mid", color=c, lw=2.4 if i == 0 else 1.4,
                 ls="-" if i == 0 else (0, (4, 2)), label=lab)
    a_n.set_yscale("log"); a_n.set_xscale("log"); a_n.set_ylabel("members per annulus")
    a_n.set_xlabel("R  [arcsec]"); a_n.legend(fontsize=7.5); a_n.set_title("How many distant tracers", fontsize=10)
    add_pc_axis(a_n, 5.43)

    a_err.plot(ref["r_median"], ref["median_err"], "o-", color=style.COLOR_FIELD, lw=1.6,
               label="median per-star error, all members")
    bright = [t for lab, t, _, _ in tables if lab.startswith("G < 18.5")][0]
    a_err.plot(bright["r_median"], bright["median_err"], "s-", color=style.SERIES[1], lw=1.6,
               label="median per-star error, G < 18.5")
    a_err.errorbar(ref["r_median"], ref["sigma_pm"], yerr=ref["sigma_pmr_err"], fmt="o", color=style.INK,
                   ms=4, lw=0, ecolor=style.INK, elinewidth=1, label="measured dispersion")
    a_err.set_xscale("log"); a_err.set_yscale("log"); a_err.set_xlabel("R  [arcsec]")
    a_err.set_ylabel("mas / yr"); a_err.legend(fontsize=8)
    a_err.set_title("Is the dispersion resolved or deconvolved?", fontsize=10)
    add_pc_axis(a_err, 5.43)

    for lab, t, c, m in tables:
        a_sig.errorbar(t["r_median"], t["sigma_pm"], yerr=t["sigma_pmr_err"], fmt=m, ms=4, color=c,
                       ecolor=c, elinewidth=1, lw=0, label=lab)
    pub = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    a_sig.plot(pub["r"], pub["sigma_pm"], color=style.SERIES[0], lw=1.5, ls="--",
               label="Vasiliev & Baumgardt 2021 published")
    a_sig.set_xscale("log"); a_sig.set_yscale("log"); a_sig.set_xlabel("R  [arcsec]")
    a_sig.set_ylabel("1-D PM dispersion  [mas/yr]"); a_sig.legend(fontsize=7.5)
    a_sig.set_title("The measurement under different cuts", fontsize=10)
    add_pc_axis(a_sig, 5.43)

    for lab, t, c, m in tables[1:]:
        a_rat.plot(t["r_median"], 100 * (np.asarray(t["sigma_pm"]) / np.asarray(ref["sigma_pm"]) - 1), m + "-",
                   ms=4, color=c, lw=1.2, label=lab)
    a_rat.axhline(0, color=style.INK, lw=1.2)
    a_rat.fill_between(ref["r_median"], -100 * ref["sigma_pmr_err"] / ref["sigma_pm"],
                       100 * ref["sigma_pmr_err"] / ref["sigma_pm"], color=style.COLOR_FIELD, alpha=0.3, lw=0,
                       label="statistical error of the reference")
    a_rat.set_xscale("log"); a_rat.set_xlabel("R  [arcsec]")
    a_rat.set_ylabel("change vs reference  [%]"); a_rat.legend(fontsize=7.5)
    a_rat.set_title("Robustness of the outer dispersion", fontsize=10)
    add_pc_axis(a_rat, 5.43)

    fig.suptitle("Audit of the distant tracers (Vasiliev & Baumgardt 2021 EDR3 members)", y=1.0)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path
