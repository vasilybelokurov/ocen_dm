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

__all__ = ["plot_constraint_map", "plot_outer_tracer_audit", "plot_contamination_model", "our_outer_profile",
           "our_mixture_profile", "OUTER_EDGES"]

#: log-spaced annuli for our own outer measurement (arcsec)
OUTER_EDGES = np.geomspace(300.0, 2400.0, 11)
#: Vasiliev & Baumgardt quality bit: the stars their own profiles are built from
QUALITY_BIT = 2


_TRACER_CACHE: dict = {}


def _composite_tracer(distance_kpc: float = 5.43):
    """The composite (star-count) tracer MGE, for the line-of-sight depth term."""
    if distance_kpc not in _TRACER_CACHE:
        from ..light_model import build_stellar_mge, fit_mge_projected, load_tracer_profile
        fit = fit_mge_projected(load_tracer_profile("composite"), sigma_range_arcsec=(7.0, 3000.0))
        _TRACER_CACHE[distance_kpc] = build_stellar_mge(fit, distance_kpc, 1.0)
    return _TRACER_CACHE[distance_kpc]


def our_mixture_profile(edges: np.ndarray = OUTER_EDGES, distance_kpc: float = 5.43,
                        g_max: float = np.inf, field_sigma_min: float = 1.5) -> Table:
    """Contamination-modelled PM dispersion: free cluster + two-Gaussian field mixture per
    annulus over ALL quality stars (no membership probability used), exact systemic field,
    depth term included. Columns as :func:`our_outer_profile` plus ``f_field``."""
    from ..kinematics.outer_profile import mixture_dispersion_free
    from ..kinematics.perspective import depth_dispersion
    s = load_members(exact=True, distance_kpc=distance_kpc)
    q = s.select(((s.quality_flag & QUALITY_BIT) > 0) & (s.g_mag <= g_max))
    tracer = _composite_tracer(distance_kpc)
    phi_sys = np.arctan2(*q.mu_sys)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (q.r_arcsec >= lo) & (q.r_arcsec < hi)
        if m.sum() < 20:
            continue
        sd = depth_dispersion(tracer, q.r_arcsec[m] * distance_kpc * 1e3 / 206264.806, float(np.hypot(*q.mu_sys)), distance_kpc)
        dphi = q.phi[m] - phi_sys
        rr = mixture_dispersion_free(q.mu_r[m], q.err_r[m], sd**2 * np.cos(dphi) ** 2, field_sigma_min=field_sigma_min)
        tt = mixture_dispersion_free(q.mu_t[m], q.err_t[m], sd**2 * np.sin(dphi) ** 2, field_sigma_min=field_sigma_min)
        rows.append((lo, float(np.median(q.r_arcsec[m])), hi, int(m.sum()), rr["sigma"], rr["sigma_err"], rr["mean"],
                     tt["sigma"], tt["sigma_err"], tt["mean"], np.sqrt(0.5 * (rr["sigma"] ** 2 + tt["sigma"] ** 2)),
                     0.5 * np.hypot(rr["sigma_err"], tt["sigma_err"]), 0.5 * (rr["f"] + tt["f"]),
                     float(np.median(q.g_mag[m]))))
    return Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_stars", "sigma_pmr", "sigma_pmr_err", "mean_pmr",
                                   "sigma_pmt", "sigma_pmt_err", "mean_pmt", "sigma_pm", "sigma_pm_err", "f_field", "median_g"))


def our_outer_profile(edges: np.ndarray = OUTER_EDGES, prob_min: float = 0.9,
                      quality_mask: int | None = QUALITY_BIT, exact: bool = True, depth: bool = True,
                      distance_kpc: float = 5.43, **kwargs) -> Table:
    """Our error-deconvolved PM dispersion profile from the member catalogue.

    ``exact`` projects the systemic 3-D velocity onto each star's own tangent basis
    (perspective contraction and basis rotation removed exactly); ``depth`` removes the
    apparent dispersion from the unknown line-of-sight depth. Both default on; switching
    them off reproduces the naive constant-subtraction measurement for comparison.
    """
    sample = load_members(exact=exact, distance_kpc=distance_kpc)
    return binned_dispersion(sample, edges, prob_min=prob_min, quality_mask=quality_mask,
                             depth_tracer=_composite_tracer(distance_kpc) if depth else None,
                             distance_kpc=distance_kpc, **kwargs)


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


def _datasets_in_kms(D: float, contamination_modelled: bool = True) -> list[dict]:
    """Every dispersion dataset converted to a 1-D dispersion in km/s.

    ``contamination_modelled`` selects our Gaia measurement: the cluster+field mixture over
    all quality stars (True) or the P > 0.9 members with the quality flag (False).
    """
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
    if contamination_modelled:
        ours = our_mixture_profile()
        name, errcol = "Gaia EDR3, our measurement (field contamination modelled)", "sigma_pm_err"
        n = np.round(np.asarray(ours["n_stars"]) * (1 - np.asarray(ours["f_field"]))).astype(int)
    else:
        ours = our_outer_profile()
        name, errcol = "Gaia EDR3, our measurement (P > 0.9 members, no contamination model)", "sigma_pmr_err"
        n = np.asarray(ours["n_stars"])
    out.append(dict(name=name, r=np.asarray(ours["r_median"]), sigma=np.asarray(ours["sigma_pm"]) * k,
                    err=np.asarray(ours[errcol]) * k, n=n, kind="pm", color=style.INK, marker="D",
                    edges=(np.asarray(ours["r_lower"]), np.asarray(ours["r_upper"])), table=ours))
    return out


def plot_constraint_map(path: Path | str = "plots/constraint_map.png",
                        k1: str = "K1_noDM_composite", k2: str = "K2_cored_composite",
                        contamination_modelled: bool = True) -> Path:
    """All dispersion data in km/s with the K1 and K2 curves, their ratio, and the tracer counts."""
    style.apply()
    j1, D, s1, _ = _model(k1)
    j2, D2, s2, _ = _model(k2)
    data = _datasets_in_kms(D, contamination_modelled)
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
    ax.set_title("What constrains the mass, and where"
                 + ("  (Gaia: field contamination modelled)" if contamination_modelled else "  (Gaia: P > 0.9 members, no contamination model)"),
                 fontsize=11)
    add_pc_axis(ax, D)

    # residuals of every dataset against K1, and the K2/K1 model ratio. The Gaia points are
    # dispersions about the rotating mean, so the model second moment has mu_rot^2/2 removed
    # before the comparison (the same convention as the likelihood).
    from ..kinematics.likelihood import pm_rotation_curve
    ratio = _sigma_1d_kms(j2, D2, R) / _sigma_1d_kms(j1, D, R)
    axr.plot(R, 100 * (ratio - 1), color=style.SERIES[1], lw=2, label="K2 / K1 model (proper motion)")
    axr.axhline(0, color=style.INK_SECONDARY, lw=0.8)
    for d in data:
        model = _sigma_1d_kms(j1, D, d["r"], d["kind"])
        if "Gaia" in d["name"]:
            rot = pm_rotation_curve(d["r"]) * KMS_PER_MASYR_KPC * D
            model = np.sqrt(np.maximum(model**2 - 0.5 * rot**2, 1e-6))
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


def plot_outer_tracer_audit(path: Path | str = "plots/outer_tracer_audit.png", reference: str = "pcut") -> Path:
    """Are the distant tracers secure, and is their dispersion resolved rather than deconvolved?

    ``reference='pcut'``: the P > 0.9 + quality-flag measurement is the reference and the
    contamination-modelled one is a variant; ``'mixture'``: the other way round.
    """
    style.apply()
    edges = OUTER_EDGES
    variants = [
        ("P > 0.9, quality flag, exact projection (reference)", dict(prob_min=0.9), style.INK, "o"),
        ("P > 0.5", dict(prob_min=0.5), style.SERIES[0], "v"),
        ("P > 0.99", dict(prob_min=0.99), style.SERIES[2], "^"),
        ("G < 18.5 (errors 0.13 mas/yr)", dict(prob_min=0.9, g_range=(0, 18.5)), style.SERIES[1], "s"),
        ("no quality cut", dict(prob_min=0.9, quality_mask=None), style.COLOR_FIELD, "x"),
        ("errors inflated 20 %", dict(prob_min=0.9, err_scale=1.2), style.SERIES_EXTRA, "d"),
        ("naive: constant systemic PM, no depth term", dict(prob_min=0.9, exact=False, depth=False), "#7a5cc7", "P"),
    ]
    tables = [(lab, our_outer_profile(edges, **kw), c, m) for lab, kw, c, m in variants]
    mix = our_mixture_profile(edges)
    mix["sigma_pmr_err"] = mix["sigma_pm_err"]
    mix["n_stars"] = np.round(np.asarray(mix["n_stars"]) * (1 - np.asarray(mix["f_field"]))).astype(int)   # cluster stars, not all stars
    if reference == "mixture":
        pcut = tables[0]
        mix["median_err"] = np.interp(np.asarray(mix["r_median"]), np.asarray(pcut[1]["r_median"]), np.asarray(pcut[1]["median_err"]))
        tables[0] = ("field contamination modelled, no P cut (reference)", mix, style.INK, "o")
        tables.append(("P > 0.9, quality flag, no contamination model", pcut[1], "#b5175f", "*"))
    else:
        tables.append(("field contamination modelled (no P cut)", mix, "#b5175f", "*"))
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
        ref_at = np.interp(np.asarray(t["r_median"]), np.asarray(ref["r_median"]), np.asarray(ref["sigma_pm"]))
        a_rat.plot(t["r_median"], 100 * (np.asarray(t["sigma_pm"]) / ref_at - 1), m + "-",
                   ms=4, color=c, lw=1.2, label=lab)
    a_rat.axhline(0, color=style.INK, lw=1.2)
    a_rat.fill_between(ref["r_median"], -100 * ref["sigma_pmr_err"] / ref["sigma_pm"],
                       100 * ref["sigma_pmr_err"] / ref["sigma_pm"], color=style.COLOR_FIELD, alpha=0.3, lw=0,
                       label="statistical error of the reference")
    a_rat.set_xscale("log"); a_rat.set_xlabel("R  [arcsec]")
    a_rat.set_ylabel("change vs reference  [%]"); a_rat.legend(fontsize=7, ncol=2)
    a_rat.set_title("Robustness of the outer dispersion", fontsize=10)
    add_pc_axis(a_rat, 5.43)

    fig.suptitle("Audit of the distant tracers (Vasiliev & Baumgardt 2021 EDR3 members)"
                 + ("  -- reference: contamination modelled" if reference == "mixture" else "  -- reference: P > 0.9 members"), y=1.0)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_contamination_model(path: Path | str = "plots/contamination_model.png",
                             annuli: tuple[tuple[float, float], ...] = ((700.0, 1000.0), (1400.0, 1800.0), (1800.0, 2400.0)),
                             distance_kpc: float = 5.43) -> Path:
    """How the field is modelled: per-annulus PM histograms of all quality stars with the fitted
    cluster + field mixture, a zoom on the cluster peak, and the summary vs radius."""
    from ..kinematics.outer_profile import dispersion_ml, mixture_dispersion_free
    from ..kinematics.perspective import depth_dispersion
    style.apply()
    s = load_members(exact=True, distance_kpc=distance_kpc)
    q = s.select((s.quality_flag & QUALITY_BIT) > 0)
    tracer = _composite_tracer(distance_kpc)
    phi_sys = np.arctan2(*q.mu_sys)
    n_ann = len(annuli)
    fig, axes = plt.subplots(2, n_ann + 1, figsize=(4.2 * (n_ann + 1), 7.4))
    grid_wide = np.linspace(-15, 15, 601); grid_zoom = np.linspace(-1.5, 1.5, 601)

    def gauss(x, m, sd):
        return np.exp(-0.5 * ((x - m) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))

    for j, (lo, hi) in enumerate(annuli):
        m = (q.r_arcsec >= lo) & (q.r_arcsec < hi)
        v, e = q.mu_r[m], q.err_r[m]
        sd = depth_dispersion(tracer, q.r_arcsec[m] * distance_kpc * 1e3 / 206264.806, float(np.hypot(*q.mu_sys)), distance_kpc)
        fit = mixture_dispersion_free(v, e, sd**2 * np.cos(q.phi[m] - phi_sys) ** 2)
        w, fm, fs = fit["field"]; f = fit["f"]; n = m.sum()
        e2 = float(np.mean(e**2))
        for ax, grid, title in ((axes[0, j], grid_wide, f"{lo:.0f}-{hi:.0f} arcsec: all {n:,} quality stars"),
                                (axes[1, j], grid_zoom, "zoom on the cluster peak")):
            bins = np.linspace(grid[0], grid[-1], 121 if grid is grid_wide else 61)
            ax.hist(v, bins=bins, color=style.COLOR_FIELD, alpha=0.7, label="data (radial PM)")
            dx = bins[1] - bins[0]
            clu = n * (1 - f) * gauss(grid, fit["mean"], np.sqrt(fit["sigma"] ** 2 + e2)) * dx
            fld = n * f * (w * gauss(grid, fm[0], fs[0]) + (1 - w) * gauss(grid, fm[1], fs[1])) * dx
            ax.plot(grid, clu, color=style.SERIES[0], lw=2, label=r"cluster: $\sigma$ = %.3f, N = %.0f" % (fit["sigma"], n * (1 - f)))
            ax.plot(grid, fld, color=style.SERIES[1], lw=2, label="field: f = %.2f, $\sigma$ = %.1f, %.1f" % (f, fs[0], fs[1]))
            ax.plot(grid, clu + fld, color=style.INK, lw=1.2, ls="--", label="total")
            ax.set_xlabel(r"$\mu_R - \mu_{R,\rm sys}$  [mas/yr]"); ax.set_title(title, fontsize=9.5)
            if grid is grid_wide:
                ax.set_yscale("log"); ax.set_ylim(0.5, 3 * max(np.histogram(v, bins=bins)[0].max(), 1))
            else:
                # what a P-cut sample would contain: field stars under the cluster peak.
                # The curves are counts per histogram bin, so integrate with the grid step / bin width.
                inwin = np.abs(grid) < 3 * np.sqrt(fit["sigma"] ** 2 + e2)
                gstep = grid[1] - grid[0]
                n_under = fld[inwin].sum() * gstep / dx
                ax.fill_between(grid[inwin], 0, fld[inwin], color=style.SERIES[1], alpha=0.35, lw=0,
                                label="field under the peak (±3σ): %.0f stars = %.1f %% of the peak" % (n_under, 100 * n_under / max(clu[inwin].sum() * gstep / dx, 1)))
            ax.legend(fontsize=7)
        axes[0, j].set_ylabel("stars per bin"); axes[1, j].set_ylabel("stars per bin")

    # summary vs radius
    mix = our_mixture_profile(OUTER_EDGES); pcut = our_outer_profile(OUTER_EDGES)
    ax = axes[0, n_ann]
    ax.errorbar(pcut["r_median"], pcut["sigma_pm"], yerr=pcut["sigma_pmr_err"], fmt="o", color="#b5175f", ms=4, lw=0,
                ecolor="#b5175f", elinewidth=1, label="P > 0.9 members (no contamination model)")
    ax.errorbar(mix["r_median"], mix["sigma_pm"], yerr=mix["sigma_pm_err"], fmt="D", color=style.INK, ms=4, lw=0,
                ecolor=style.INK, elinewidth=1, label="cluster + field mixture, all stars")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("R  [arcsec]"); ax.set_ylabel("1-D PM dispersion  [mas/yr]")
    ax.legend(fontsize=7.5); ax.set_title("Effect on the dispersion profile", fontsize=9.5); add_pc_axis(ax, distance_kpc)
    ax = axes[1, n_ann]
    ax.plot(mix["r_median"], mix["f_field"], "o-", color=style.SERIES[1], lw=1.6, label="field fraction of all quality stars (mixture)")
    ax.plot(mix["r_median"], 100 * (np.asarray(pcut["sigma_pm"]) / np.asarray(mix["sigma_pm"]) - 1) / 100, "s--",
            color="#b5175f", lw=1.2, label=r"($\sigma_{P>0.9}/\sigma_{\rm mixture} - 1$)")
    ax.axhline(0, color=style.INK_SECONDARY, lw=0.8)
    ax.set_xscale("log"); ax.set_xlabel("R  [arcsec]"); ax.set_ylabel("fraction"); ax.legend(fontsize=7.5)
    ax.set_title("Field fraction and the bias of the P cut", fontsize=9.5); add_pc_axis(ax, distance_kpc)
    fig.suptitle("Field contamination of the Gaia EDR3 members: the mixture model", y=1.0)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path
