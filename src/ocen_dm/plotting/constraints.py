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
from ..selection.hst_gaia_match import OCEN_DEC, OCEN_RA
from ..paths import processed_dir, results_dir
from . import style
from .style import add_arcsec_axis, add_pc_axis

__all__ = ["plot_constraint_map", "plot_outer_tracer_audit", "plot_contamination_model", "plot_annulus_fits",
           "plot_method_comparison", "plot_residual_significance", "plot_dataset_step", "plot_offset_explained", "plot_datasets_unscaled", "plot_hst_gaia_star_by_star", "hst_gaia_excess_table", "hst_gaia_overlap_table", "plot_pm_datasets", "plot_periphery", "plot_periphery_density", "plot_extended_profile", "plot_master_datasets", "plot_hst_gaia_overlap", "hst_gaia_overlap_profile", "hst_gaia_overlap_tests",
           "fit_quality_table", "annulus_fits", "our_outer_profile", "our_mixture_profile", "OUTER_EDGES"]

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
                        g_max: float = np.inf, field_sigma_min: float = 1.5,
                        field: str = "2d", bw_method: float = 0.05) -> Table:
    """Contamination-modelled PM dispersion over ALL quality stars (no membership cut).

    ``field='2d'`` (default) scores every star against the **two-dimensional** empirical
    proper-motion density of an independent Gaia DR3 annulus outside the cluster, with the
    cluster a two-dimensional Gaussian whose axes follow each star's radial direction;
    ``field='template'`` uses the same field sample projected onto the radial/tangential
    directions (one dimension at a time), and ``field='free'`` fits a two-Gaussian field to
    the same stars. The last two are kept for comparison: projecting the field mixes its
    two unequal widths around an annulus and the two-Gaussian form cannot reproduce its
    heavy tails, so both misfit the observed distribution (JOURNAL 2026-09-18).
    """
    from ..kinematics.outer_profile import dispersion_2d, dispersion_with_field_template, mixture_dispersion_free
    from ..kinematics.perspective import depth_dispersion
    kdes = dens2d = None
    if field == "2d":
        from ..selection.field_template import field_density_2d
        dens2d = field_density_2d()
    elif field == "template":
        from ..selection.field_template import field_kde, load_field_template
        tpl = load_field_template()
        kdes = {c: field_kde(c, tpl, bw_method=bw_method) for c in ("r", "t")}
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
        if dens2d is not None:
            o = dispersion_2d(q, m, dens2d, depth_var=sd**2, field_at=q.absolute_pm)
            rows.append((lo, float(np.median(q.r_arcsec[m])), hi, int(m.sum()), o["sigma_r"], o["sigma_r_err"], o["mean_r"],
                         o["sigma_t"], o["sigma_t_err"], o["mean_t"], o["sigma"], o["sigma_err"], o["f"],
                         float(np.median(q.g_mag[m]))))
            continue
        if kdes is None:
            rr = mixture_dispersion_free(q.mu_r[m], q.err_r[m], sd**2 * np.cos(dphi) ** 2, field_sigma_min=field_sigma_min)
            tt = mixture_dispersion_free(q.mu_t[m], q.err_t[m], sd**2 * np.sin(dphi) ** 2, field_sigma_min=field_sigma_min)
        else:
            rr = dispersion_with_field_template(q.mu_r[m], q.err_r[m], kdes["r"], sd**2 * np.cos(dphi) ** 2)
            tt = dispersion_with_field_template(q.mu_t[m], q.err_t[m], kdes["t"], sd**2 * np.sin(dphi) ** 2)
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


def _datasets_in_kms(D: float, contamination_modelled: bool = True, streaming: str = "self") -> list[dict]:
    """Every dispersion dataset converted to a 1-D dispersion in km/s.

    ``contamination_modelled`` selects our Gaia measurement: the cluster+field mixture over
    all quality stars (True) or the P > 0.9 members with the quality flag (False).
    ``streaming`` selects how the rotation term is removed from the model for our own
    points: ``'self'`` uses the mean motions fitted in the same annuli as the dispersions,
    ``'published'`` uses the Vasiliev & Baumgardt rotation curve (the original treatment,
    kept for comparison -- it over-subtracts by 10-40 per cent near the rotation peak and
    manufactures a wiggle at 10-20 pc, JOURNAL 2026-09-18).
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
    # our own fit measures the mean motions as well, so the streaming term for these points
    # is measured rather than taken from the published rotation curve
    stream2 = None
    if streaming == "self" and "mean_pmt" in ours.colnames:
        stream2 = 0.5 * (np.asarray(ours["mean_pmr"]) ** 2 + np.asarray(ours["mean_pmt"]) ** 2) * k**2
    out.append(dict(name=name, r=np.asarray(ours["r_median"]), sigma=np.asarray(ours["sigma_pm"]) * k,
                    err=np.asarray(ours[errcol]) * k, n=n, kind="pm", color=style.INK, marker="D",
                    edges=(np.asarray(ours["r_lower"]), np.asarray(ours["r_upper"])), table=ours,
                    streaming2=stream2))
    return out


def plot_constraint_map(path: Path | str = "plots/constraint_map.png",
                        k1: str = "K1_noDM_composite", k2: str = "K2_cored_composite",
                        contamination_modelled: bool = True, streaming: str = "self") -> Path:
    """All dispersion data in km/s with the K1 and K2 curves, their ratio, and the tracer counts."""
    style.apply()
    j1, D, s1, _ = _model(k1)
    j2, D2, s2, _ = _model(k2)
    data = _datasets_in_kms(D, contamination_modelled, streaming)
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
                 + ("  (Gaia: field contamination modelled" if contamination_modelled else "  (Gaia: P > 0.9 members, no contamination model")
                 + ("; rotation from our own fit)" if streaming == "self" else "; rotation from the published curve)"),
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
        # published dispersions are measured about a rotating mean, so the model's second
        # moment must have the streaming term removed before comparison. Our own points carry
        # their measured mean motions; the others use the Vasiliev & Baumgardt rotation curve.
        if d.get("streaming2") is not None:
            model = np.sqrt(np.maximum(model**2 - d["streaming2"], 1e-6))
        elif "Gaia" in d["name"]:
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


def _model_histogram(sub, fit, dens2d, bins, depth_var, component: str = "r",
                     n_over: int = 12, seed: int = 0):
    """Expected counts per bin of the projected proper motion under the fitted 2-D model.

    The cluster part is integrated analytically for each star (its own error, its own depth
    term, the dispersion along the requested component); the field part is drawn from the
    empirical two-dimensional template **in absolute proper motion** and projected onto each
    star's own radial or tangential direction, which is what makes the projected field
    non-Gaussian.
    """
    from scipy.stats import norm

    from ..selection.field_template import load_field_template

    rng = np.random.default_rng(seed)
    n = len(sub.mu_r); f = fit["f"]
    radial = component == "r"
    sig = fit["sigma_r"] if radial else fit["sigma_t"]
    mean = fit["mean_r"] if radial else fit["mean_t"]
    err = sub.err_r if radial else sub.err_t
    # the depth term points along the systemic proper motion; project it onto this component
    mu_hat = np.asarray(sub.mu_sys, float); mu_hat = mu_hat / max(float(np.hypot(*mu_hat)), 1e-12)
    phi_sys = np.arctan2(mu_hat[0], mu_hat[1])
    proj_depth = np.cos(sub.phi - phi_sys) ** 2 if radial else np.sin(sub.phi - phi_sys) ** 2
    var = sig**2 + np.asarray(depth_var, float) * proj_depth + err**2
    cdf = norm.cdf((bins[None, :] - mean) / np.sqrt(var)[:, None])
    clu = (1 - f) * np.diff(cdf, axis=1).sum(axis=0)
    t = load_field_template(); keep = np.asarray(t["r_arcsec"], float) >= 3600.0
    fa = (np.asarray(t["mu_a"], float) + np.asarray(t["sys_a"], float))[keep]
    fd = (np.asarray(t["mu_d"], float) + np.asarray(t["sys_d"], float))[keep]
    n_f = int(round(f * n)) * n_over
    j = rng.choice(len(fa), size=n_f, replace=True); k = rng.choice(n, size=n_f, replace=True)
    # back to the cluster frame at the target star's position, then project
    ra_, rd_ = fa[j] - sub.sys_a[k], fd[j] - sub.sys_d[k]
    if radial:
        proj = ra_ * np.sin(sub.phi[k]) + rd_ * np.cos(sub.phi[k])
    else:
        proj = -ra_ * np.cos(sub.phi[k]) + rd_ * np.sin(sub.phi[k])
    fld = np.histogram(proj + rng.normal(0, err[k]), bins)[0] / n_over
    return clu, fld


def annulus_fits(edges: np.ndarray = OUTER_EDGES, distance_kpc: float = 5.43):
    """Fit every annulus and return ``(sample, [(lo, hi, subsample, fit, depth_var), ...])``."""
    from ..kinematics.outer_profile import dispersion_2d
    from ..kinematics.perspective import depth_dispersion
    from ..selection.field_template import field_density_2d

    s = load_members(exact=True, distance_kpc=distance_kpc)
    q = s.select((s.quality_flag & QUALITY_BIT) > 0)
    tracer = _composite_tracer(distance_kpc)
    dens2d = field_density_2d()
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (q.r_arcsec >= lo) & (q.r_arcsec < hi)
        if m.sum() < 20:
            continue
        sub = q.select(m)
        sd = depth_dispersion(tracer, sub.r_arcsec * distance_kpc * 1e3 / 206264.806,
                              float(np.hypot(*q.mu_sys)), distance_kpc)
        out.append((lo, hi, sub, dispersion_2d(q, m, dens2d, depth_var=sd**2, field_at=q.absolute_pm), sd**2))
    return q, dens2d, out


def fit_quality_table(edges: np.ndarray = OUTER_EDGES, distance_kpc: float = 5.43,
                      save: bool = True) -> Table:
    """Per-annulus parameters and goodness of fit of the 2-D cluster + field model.

    chi2 per bin of the projected histograms, over the full proper-motion range (where the
    field dominates) and over the cluster peak, for both components.
    """
    q, dens2d, fits = annulus_fits(edges, distance_kpc)
    rows = []
    for lo, hi, sub, fit, dvar in fits:
        chis = {}
        for comp in ("r", "t"):
            v = sub.mu_r if comp == "r" else sub.mu_t
            for tag, lim, nb in (("wide", (-15, 15), 121), ("peak", (-1.5, 1.5), 61)):
                bins = np.linspace(lim[0], lim[1], nb)
                data, _ = np.histogram(v, bins)
                clu, fld = _model_histogram(sub, fit, dens2d, bins, dvar, component=comp)
                model = clu + fld
                chis[f"chi2_{comp}_{tag}"] = float(np.sum((data - model) ** 2 / np.maximum(model, 1)) / (nb - 1))
        rows.append((lo, float(np.median(sub.r_arcsec)), hi, int(len(sub)), fit["f"], fit["n_cluster"],
                     fit["sigma_r"], fit["sigma_r_err"], fit["sigma_t"], fit["sigma_t_err"],
                     fit["mean_r"], fit["mean_t"], chis["chi2_r_wide"], chis["chi2_r_peak"],
                     chis["chi2_t_wide"], chis["chi2_t_peak"]))
    t = Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_stars", "f_field", "n_cluster",
                                "sigma_pmr", "sigma_pmr_err", "sigma_pmt", "sigma_pmt_err",
                                "mean_pmr", "mean_pmt", "chi2_r_wide", "chi2_r_peak",
                                "chi2_t_wide", "chi2_t_peak"))
    t.meta["description"] = ("2-D cluster + empirical-field fit per annulus; chi2 values are per bin of the "
                             "projected histogram, 'wide' = |mu| < 15 mas/yr (120 bins), 'peak' = |mu| < 1.5 (60 bins)")
    if save:
        path = processed_dir() / "kinematics" / "ocen_outer_fit_quality.ecsv"
        t.write(path, format="ascii.ecsv", overwrite=True)
    return t


def plot_annulus_fits(path: Path | str = "plots/outer_fit_annuli.png", component: str = "r",
                      edges: np.ndarray = OUTER_EDGES, distance_kpc: float = 5.43, n_col: int = 2) -> Path:
    """Data and fitted model in every annulus, with residuals and the chi2 of each panel."""
    style.apply()
    q, dens2d, fits = annulus_fits(edges, distance_kpc)
    label = r"\mu_R" if component == "r" else r"\mu_T"
    n = len(fits); n_row = int(np.ceil(n / n_col))
    fig = plt.figure(figsize=(6.6 * n_col, 3.4 * n_row))
    outer = fig.add_gridspec(n_row, n_col, hspace=0.42, wspace=0.22)
    for k, (lo, hi, sub, fit, dvar) in enumerate(fits):
        inner = outer[k // n_col, k % n_col].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax = fig.add_subplot(inner[0]); axr = fig.add_subplot(inner[1], sharex=ax)
        v = sub.mu_r if component == "r" else sub.mu_t
        bins = np.linspace(-15, 15, 121); mid = 0.5 * (bins[1:] + bins[:-1])
        data, _ = np.histogram(v, bins)
        clu, fld = _model_histogram(sub, fit, dens2d, bins, dvar, component=component)
        model = clu + fld
        ax.bar(mid, data, width=bins[1] - bins[0], color=style.COLOR_FIELD, alpha=0.75, label="data")
        ax.step(bins, np.append(clu, clu[-1]), where="post", color=style.SERIES[0], lw=1.6,
                label=r"cluster: $\sigma$ = %.3f $\pm$ %.3f, N = %.0f" % (
                    fit["sigma_r"] if component == "r" else fit["sigma_t"],
                    fit["sigma_r_err"] if component == "r" else fit["sigma_t_err"], fit["n_cluster"]))
        ax.step(bins, np.append(fld, fld[-1]), where="post", color=style.SERIES[1], lw=1.6,
                label="field: f = %.3f" % fit["f"])
        ax.step(bins, np.append(model, model[-1]), where="post", color=style.INK, lw=1.1, ls="--", label="total")
        ax.set_yscale("log"); ax.set_ylim(0.5, 3 * max(data.max(), 1))
        ax.set_ylabel("stars per bin"); ax.legend(fontsize=6.5, loc="upper right")
        ax.set_title(r"%.0f$-$%.0f arcsec  (%.1f$-$%.1f pc):  %s stars" % (
            lo, hi, lo * distance_kpc * 1e3 / 206264.806, hi * distance_kpc * 1e3 / 206264.806, f"{len(sub):,}"),
            fontsize=9.5)
        ax.tick_params(labelbottom=False)
        res = (data - model) / np.sqrt(np.maximum(model, 1))
        axr.bar(mid, res, width=bins[1] - bins[0], color=style.INK_SECONDARY)
        axr.axhline(0, color=style.INK, lw=0.8); axr.set_ylim(-4.5, 4.5)
        axr.set_xlabel((r"$\mu_R$" if component == "r" else r"$\mu_T$") + "  [mas/yr]  (systemic motion removed)")
        axr.set_ylabel(r"$\chi$", fontsize=8)
        peak = np.abs(mid) < 1.5
        axr.text(0.02, 0.78, r"$\chi^2$/bin = %.2f (all), %.2f (peak)" % (
            np.sum(res**2) / len(res), np.sum(res[peak] ** 2) / max(peak.sum(), 1)),
            transform=axr.transAxes, fontsize=7.5)
        axr.axvspan(-1.5, 1.5, color=style.SERIES[0], alpha=0.07, lw=0)
    fig.suptitle("Per-annulus fit, %s component: data, cluster + field model, residuals"
                 % ("radial" if component == "r" else "tangential"), y=0.995)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_contamination_model(path: Path | str = "plots/contamination_model.png",
                             annuli: tuple[tuple[float, float], ...] = ((700.0, 1000.0), (1400.0, 1800.0), (1800.0, 2400.0)),
                             distance_kpc: float = 5.43) -> Path:
    """How the field is modelled and how well it fits.

    Top row: the observed radial-PM histogram of all quality stars in each annulus with the
    fitted cluster and field components of the **two-dimensional** model; middle row: the
    zoom on the cluster peak; bottom row: residuals. Right column: the two-dimensional field
    density itself, the effect on the dispersion profile, and the resulting anisotropy.
    """
    from ..kinematics.outer_profile import dispersion_2d
    from ..kinematics.perspective import depth_dispersion
    from ..selection.field_template import field_density_2d, load_field_template
    style.apply()
    s = load_members(exact=True, distance_kpc=distance_kpc)
    q = s.select((s.quality_flag & QUALITY_BIT) > 0)
    tracer = _composite_tracer(distance_kpc)
    dens2d = field_density_2d()
    n_ann = len(annuli)
    fig = plt.figure(figsize=(4.6 * (n_ann + 1), 10.4))
    gs = fig.add_gridspec(3, n_ann + 1, height_ratios=[2.2, 2.2, 1.0], hspace=0.45, wspace=0.34)
    axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(n_ann + 1)] for i in range(3)])

    for j, (lo, hi) in enumerate(annuli):
        m = (q.r_arcsec >= lo) & (q.r_arcsec < hi)
        sub = q.select(m)
        sd = depth_dispersion(tracer, sub.r_arcsec * distance_kpc * 1e3 / 206264.806, float(np.hypot(*q.mu_sys)), distance_kpc)
        fit = dispersion_2d(q, m, dens2d, depth_var=sd**2)
        n = int(m.sum())
        for row, (lim, nb, title) in enumerate((((-15, 15), 121, f"{lo:.0f}-{hi:.0f} arcsec: all {n:,} quality stars"),
                                                ((-1.5, 1.5), 61, "zoom on the cluster peak"))):
            ax = axes[row, j]
            bins = np.linspace(lim[0], lim[1], nb); mid = 0.5 * (bins[1:] + bins[:-1])
            data, _ = np.histogram(sub.mu_r, bins)
            clu, fld = _model_histogram(sub, fit, dens2d, bins, sd**2, component="r")
            ax.bar(mid, data, width=bins[1] - bins[0], color=style.COLOR_FIELD, alpha=0.7, label="data (radial PM)")
            ax.step(bins, np.append(clu, clu[-1]), where="post", color=style.SERIES[0], lw=1.8,
                    label=r"cluster: $\sigma_R$ = %.3f, N = %.0f" % (fit["sigma_r"], fit["n_cluster"]))
            ax.step(bins, np.append(fld, fld[-1]), where="post", color=style.SERIES[1], lw=1.8,
                    label="field: f = %.2f (2-D template, projected)" % fit["f"])
            ax.step(bins, np.append(clu + fld, (clu + fld)[-1]), where="post", color=style.INK, lw=1.2, ls="--", label="total")
            ax.set_title(title, fontsize=9.5); ax.set_ylabel("stars per bin"); ax.legend(fontsize=6.8)
            if row == 0:
                ax.set_yscale("log"); ax.set_ylim(0.5, 3 * max(data.max(), 1))
            else:
                res = (data - clu - fld) / np.sqrt(np.maximum(clu + fld, 1))
                axr = axes[2, j]
                axr.bar(mid, res, width=bins[1] - bins[0], color=style.INK_SECONDARY)
                axr.axhline(0, color=style.INK, lw=0.8); axr.set_ylim(-4.5, 4.5)
                axr.set_xlabel(r"$\mu_R - \mu_{R,\rm sys}$  [mas/yr]"); axr.set_ylabel(r"(data$-$model)/$\sqrt{\rm model}$")
                axr.text(0.02, 0.82, r"$\chi^2$/bin = %.2f" % (np.sum(res**2) / len(res)), transform=axr.transAxes, fontsize=8.5)

    # the field template itself
    t = load_field_template(); keep = np.asarray(t["r_arcsec"], float) >= 3600.0
    ax = axes[0, n_ann]
    h = ax.hexbin(np.asarray(t["mu_a"])[keep], np.asarray(t["mu_d"])[keep], gridsize=60, extent=(-15, 15, -15, 15),
                  bins="log", cmap=style.SEQUENTIAL, mincnt=1)
    ax.plot(0, 0, "x", color=style.SERIES[1], ms=9, mew=2)
    ax.text(0.5, -2.2, "cluster sits here", color=style.SERIES[1], fontsize=7.5)
    ax.set_xlabel(r"$\mu_{\alpha*}$ residual [mas/yr]"); ax.set_ylabel(r"$\mu_\delta$ residual [mas/yr]")
    ax.set_title("the 2-D field template (%s stars, > 1 deg)" % f"{int(keep.sum()):,}", fontsize=9.5)
    fig.colorbar(h, ax=ax, label="field stars per cell", fraction=0.046, pad=0.03)

    mix = our_mixture_profile(OUTER_EDGES); pcut = our_outer_profile(OUTER_EDGES)
    ax = axes[1, n_ann]
    ax.errorbar(pcut["r_median"], pcut["sigma_pm"], yerr=pcut["sigma_pmr_err"], fmt="o", color="#b5175f", ms=4, lw=0,
                ecolor="#b5175f", elinewidth=1, label="P > 0.9 members (no field model)")
    ax.errorbar(mix["r_median"], mix["sigma_pm"], yerr=mix["sigma_pm_err"], fmt="D", color=style.INK, ms=4, lw=0,
                ecolor=style.INK, elinewidth=1, label="2-D field model, all stars")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("R  [arcsec]"); ax.set_ylabel("1-D PM dispersion  [mas/yr]")
    ax.legend(fontsize=7.5); ax.set_title("Effect on the dispersion profile", fontsize=9.5); add_pc_axis(ax, distance_kpc)
    ax = axes[2, n_ann]
    ratio = np.asarray(mix["sigma_pmt"]) / np.asarray(mix["sigma_pmr"])
    rerr = ratio * np.hypot(np.asarray(mix["sigma_pmt_err"]) / np.asarray(mix["sigma_pmt"]),
                            np.asarray(mix["sigma_pmr_err"]) / np.asarray(mix["sigma_pmr"]))
    ax.errorbar(mix["r_median"], ratio, yerr=rerr, fmt="o-", color=style.SERIES[2], ms=4, lw=1.4, ecolor=style.SERIES[2])
    ax.axhline(1.0, color=style.INK_SECONDARY, lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_xlabel("R  [arcsec]"); ax.set_ylabel(r"$\sigma_T / \sigma_R$")
    ax.text(0.03, 0.08, "radial", transform=ax.transAxes, fontsize=8, color=style.INK_SECONDARY)
    ax.text(0.03, 0.85, "tangential", transform=ax.transAxes, fontsize=8, color=style.INK_SECONDARY)
    ax.set_title("Anisotropy from the same fit", fontsize=9.5)

    fig.suptitle("Field contamination of the Gaia EDR3 outskirts: the two-dimensional model", y=0.995)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_method_comparison(path: Path | str = "plots/outer_method_comparison.png",
                           edges: np.ndarray = OUTER_EDGES, distance_kpc: float = 5.43) -> Path:
    """Membership cut versus cluster+field decomposition, component by component.

    Panels: the radial and tangential dispersion profiles measured both ways; their
    fractional difference against the statistical errors; and the anisotropy each method
    implies. The field fraction of each annulus runs along the top of the difference panel,
    since that is what sets how much the two can differ.
    """
    style.apply()
    mix = our_mixture_profile(edges, distance_kpc=distance_kpc)     # 2-D cluster + field, no P cut
    cut = our_outer_profile(edges)                                   # P > 0.9 members, error deconvolution
    r = np.asarray(mix["r_median"])
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4))
    (a_r, a_t), (a_d, a_b) = axes
    cols = {"cut": "#b5175f", "mix": style.INK}

    for ax, comp, name in ((a_r, "sigma_pmr", "radial"), (a_t, "sigma_pmt", "tangential")):
        ax.errorbar(np.asarray(cut["r_median"]), cut[comp], yerr=cut[f"{comp}_err"], fmt="o", ms=5,
                    color=cols["cut"], ecolor=cols["cut"], elinewidth=1.2, lw=0,
                    label="P > 0.9 members, errors deconvolved")
        ax.errorbar(r, mix[comp], yerr=mix[f"{comp}_err"], fmt="D", ms=5, color=cols["mix"],
                    ecolor=cols["mix"], elinewidth=1.2, lw=0, label="cluster + field decomposition (all stars)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("R  [arcsec]"); ax.set_ylabel(r"$\sigma_%s$  [mas/yr]" % ("R" if comp.endswith("pmr") else "T"))
        ax.set_title("%s component" % name, fontsize=10.5); ax.legend(fontsize=8)
        add_pc_axis(ax, distance_kpc)

    for comp, colour, marker, lab in (("sigma_pmr", style.SERIES[0], "o", r"$\sigma_R$"),
                                      ("sigma_pmt", style.SERIES[2], "^", r"$\sigma_T$")):
        d = np.asarray(mix[comp]) / np.asarray(cut[comp]) - 1.0
        e = np.hypot(np.asarray(mix[f"{comp}_err"]) / np.asarray(mix[comp]),
                     np.asarray(cut[f"{comp}_err"]) / np.asarray(cut[comp]))
        a_d.errorbar(r, 100 * d, yerr=100 * e, fmt=marker + "-", ms=5, color=colour, ecolor=colour,
                     elinewidth=1.2, lw=1.2, label=lab)
    a_d.axhline(0, color=style.INK, lw=1)
    a_d.set_xscale("log"); a_d.set_xlabel("R  [arcsec]")
    a_d.set_ylabel("decomposition $-$ membership cut  [%]")
    a_d.set_title("How much the two methods differ", fontsize=10.5); a_d.legend(fontsize=8, loc="upper left")
    add_pc_axis(a_d, distance_kpc)
    ax2 = a_d.twinx()
    ax2.plot(r, 100 * np.asarray(mix["f_field"]), ":", color=style.SERIES[1], lw=1.6)
    ax2.set_ylabel("field fraction of the annulus  [%]", color=style.SERIES[1])
    ax2.tick_params(axis="y", colors=style.SERIES[1]); ax2.set_ylim(0, 100); ax2.grid(False)

    for tab, colour, marker, lab in ((cut, cols["cut"], "o", "P > 0.9 members"),
                                     (mix, cols["mix"], "D", "decomposition")):
        ratio = np.asarray(tab["sigma_pmt"]) / np.asarray(tab["sigma_pmr"])
        err = ratio * np.hypot(np.asarray(tab["sigma_pmt_err"]) / np.asarray(tab["sigma_pmt"]),
                               np.asarray(tab["sigma_pmr_err"]) / np.asarray(tab["sigma_pmr"]))
        a_b.errorbar(np.asarray(tab["r_median"]), ratio, yerr=err, fmt=marker, ms=5, color=colour,
                     ecolor=colour, elinewidth=1.2, lw=0, label=lab)
    a_b.axhline(1.0, color=style.INK_SECONDARY, lw=0.9, ls="--")
    a_b.set_xscale("log"); a_b.set_xlabel("R  [arcsec]"); a_b.set_ylabel(r"$\sigma_T / \sigma_R$")
    a_b.set_title("Anisotropy implied by each method", fontsize=10.5); a_b.legend(fontsize=8)
    a_b.text(0.03, 0.08, "radial orbits", transform=a_b.transAxes, fontsize=8, color=style.INK_SECONDARY)
    a_b.text(0.03, 0.88, "tangential orbits", transform=a_b.transAxes, fontsize=8, color=style.INK_SECONDARY)
    add_pc_axis(a_b, distance_kpc)

    fig.suptitle("Membership cut versus cluster + background decomposition, Gaia EDR3 outskirts", y=0.995)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_dataset_step(path: Path | str = "plots/hst_gaia_step.png", k1: str = "K1_noDM_composite",
                      distance_kpc: float = 5.43) -> Path:
    """Where HST ends and Gaia begins: do the two datasets join smoothly?

    Both are compared with the **same unscaled** model, so the fitted Gaia scale cannot hide
    a mismatch. HST's outermost bin reaches 346 arcsec; Gaia's quality-flagged sample only
    becomes usable near 500 arcsec, so the two barely overlap and the model bridges the gap
    unconstrained.
    """
    import json

    from ..kinematics.report import _family_for
    style.apply()
    edges = np.array([300., 380, 460, 540, 630, 730, 850, 1000, 1200, 1500, 1900, 2400.])
    mix = our_mixture_profile(edges, distance_kpc=distance_kpc)
    s = load_members(exact=True, distance_kpc=distance_kpc)
    q = s.select((s.quality_flag & QUALITY_BIT) > 0)
    summary = json.loads((results_dir() / "fits" / k1 / "summary.json").read_text())
    fam = _family_for(summary)
    x = np.array([summary["parameters"][n]["ml"] for n in fam.names])
    jeans, D, _ = fam.build(fam.to_dict(x))

    hr = Table.read(processed_dir() / "kinematics" / "omegacat_vi_pm_radial.ecsv")
    ht = Table.read(processed_dir() / "kinematics" / "omegacat_vi_pm_tangential.ecsv")
    rh = np.asarray(hr["r_median"])
    hst = np.sqrt(0.5 * (np.asarray(hr["sigma_pmr"]) ** 2 + np.asarray(ht["sigma_pmt"]) ** 2))
    hst_err = 0.5 * (np.asarray(hr["sigma_pmr_err_lo"]) + np.asarray(hr["sigma_pmr_err_hi"])) / np.sqrt(2)
    mod_h = _sigma_1d_kms(jeans, D, rh) / (KMS_PER_MASYR_KPC * D)

    rg = np.asarray(mix["r_median"])
    gaia = np.sqrt(0.5 * (np.asarray(mix["sigma_pmr"]) ** 2 + np.asarray(mix["mean_pmr"]) ** 2
                          + np.asarray(mix["sigma_pmt"]) ** 2 + np.asarray(mix["mean_pmt"]) ** 2))
    gaia_err = np.asarray(mix["sigma_pm_err"])
    mod_g = _sigma_1d_kms(jeans, D, rg) / (KMS_PER_MASYR_KPC * D)

    fig, (ax, axq) = plt.subplots(2, 1, figsize=(9.5, 7.6), sharex=True,
                                  gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.08})
    ax.axvspan(346, 500, color=style.SERIES[1], alpha=0.12, lw=0)
    ax.text(352, -13, "no reliable data\nfrom either instrument", fontsize=8, color=style.SERIES[1])
    ax.errorbar(rh, 100 * (hst / mod_h - 1), yerr=100 * hst_err / mod_h, fmt="o", ms=5, color=style.SERIES[0],
                ecolor=style.SERIES[0], elinewidth=1.2, lw=0, label="HST oMEGACat (the model's anchor)")
    ax.errorbar(rg, 100 * (gaia / mod_g - 1), yerr=100 * gaia_err / mod_g, fmt="D", ms=6, color=style.INK,
                ecolor=style.INK, elinewidth=1.4, lw=0, label="Gaia EDR3, our measurement")
    ax.axhline(0, color=style.INK, lw=1)
    ax.set_xscale("log"); ax.set_xlim(140, 2600)
    ax.set_ylabel("deviation from the same unscaled model  [%]")
    ax.set_title("HST ends at 346 arcsec, Gaia becomes usable near 500: do they join?", fontsize=11)
    ax.legend(fontsize=8.5, loc="upper left"); add_pc_axis(ax, distance_kpc)

    frac = [((q.r_arcsec >= lo) & (q.r_arcsec < hi)).sum() / max(((s.r_arcsec >= lo) & (s.r_arcsec < hi)).sum(), 1)
            for lo, hi in zip(edges[:-1], edges[1:])]
    axq.step(rg, 100 * np.array(frac), where="mid", color=style.SERIES[2], lw=2, label="Gaia stars passing the quality flag")
    axq.step(rh, 100 * np.ones_like(rh), where="mid", color=style.SERIES[0], lw=2, label="HST coverage")
    axq.axvspan(346, 500, color=style.SERIES[1], alpha=0.12, lw=0)
    axq.axhline(15, color=style.INK_SECONDARY, lw=0.8, ls="--")
    axq.set_xscale("log"); axq.set_xlabel("R  [arcsec]"); axq.set_ylabel("usable stars  [%]")
    axq.legend(fontsize=8)
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_residual_significance(path: Path | str = "plots/outer_residual_significance.png",
                               k1: str = "K1_noDM_composite", edges: np.ndarray | None = None,
                               distance_kpc: float = 5.43) -> Path:
    """Is the apparent wiggle in the Gaia residuals real?

    Top: the deviation of our Gaia measurement from the no-dark-matter model, in total
    second moments so that no rotation curve enters the comparison, with the constant-offset
    fit over the range where Gaia is reliable. Bottom: the fraction of catalogue stars that
    survive the astrometric quality flag, which is what makes the innermost annuli useless.
    """
    import json

    from scipy.stats import chi2 as chi2_dist

    from ..kinematics.report import _family_for
    style.apply()
    edges = np.array([300., 380, 460, 540, 630, 730, 850, 1000, 1200, 1500, 1900, 2400.]) if edges is None else edges
    mix = our_mixture_profile(edges, distance_kpc=distance_kpc)
    s = load_members(exact=True, distance_kpc=distance_kpc)
    q = s.select((s.quality_flag & QUALITY_BIT) > 0)

    summary = json.loads((results_dir() / "fits" / k1 / "summary.json").read_text())
    fam = _family_for(summary)
    x = np.array([summary["parameters"][n]["ml"] for n in fam.names])
    jeans, D, scales = fam.build(fam.to_dict(x))
    r = np.asarray(mix["r_median"])
    total = np.sqrt(0.5 * (np.asarray(mix["sigma_pmr"]) ** 2 + np.asarray(mix["mean_pmr"]) ** 2
                           + np.asarray(mix["sigma_pmt"]) ** 2 + np.asarray(mix["mean_pmt"]) ** 2))
    model = _sigma_1d_kms(jeans, D, r) / (KMS_PER_MASYR_KPC * D) * scales.get("GaiaEDR3", 1.0)
    res = 100 * (total / model - 1); rerr = 100 * np.asarray(mix["sigma_pm_err"]) / model

    fig, (ax, axq) = plt.subplots(2, 1, figsize=(9.5, 7.4), sharex=True,
                                  gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.08})
    good = (r > 460) & (r < 1500)
    w = 1 / rerr[good] ** 2
    c = float(np.sum(w * res[good]) / np.sum(w))
    chi2 = float(np.sum(w * (res[good] - c) ** 2)); dof = int(good.sum() - 1)
    ax.axhspan(c - 1, c + 1, color=style.SERIES[0], alpha=0.15, lw=0)
    ax.axhline(c, color=style.SERIES[0], lw=2,
               label=r"constant %+.1f %% over 460$-$1500 arcsec: $\chi^2$ = %.1f / %d, p = %.2f"
                     % (c, chi2, dof, 1 - chi2_dist.cdf(chi2, dof)))
    ax.axhline(0, color=style.INK, lw=1)
    ax.errorbar(r, res, yerr=rerr, fmt="D", ms=6, color=style.INK, ecolor=style.INK, elinewidth=1.4, lw=0,
                label="our Gaia measurement (total second moment)")
    for i in np.flatnonzero(~good):
        ax.annotate("%.1f$\sigma$" % (res[i] / rerr[i]), (r[i], res[i]), textcoords="offset points",
                    xytext=(0, 12 if res[i] > 0 else -18), ha="center", fontsize=8, color=style.SERIES[1])
    ax.axvspan(edges[0], 460, color=style.COLOR_FIELD, alpha=0.3, lw=0)
    ax.axvspan(1500, edges[-1], color=style.SERIES[1], alpha=0.08, lw=0)
    ax.text(330, ax.get_ylim()[1] * 0.75, "Gaia unusable here\n(< 15 % of stars pass\nthe quality flag)",
            fontsize=8, color=style.INK_SECONDARY)
    ax.text(1560, 12, "genuine rise\n(4-5 sigma)", fontsize=8, color=style.SERIES[1])
    ax.set_ylabel("deviation from the no-DM model  [%]")
    ax.set_xscale("log"); ax.legend(fontsize=8.5, loc="lower right")
    ax.set_title("Is the wiggle real? Residuals with their errors", fontsize=11)
    add_pc_axis(ax, distance_kpc)

    frac = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        a = ((s.r_arcsec >= lo) & (s.r_arcsec < hi)).sum()
        b = ((q.r_arcsec >= lo) & (q.r_arcsec < hi)).sum()
        frac.append(b / max(a, 1))
    axq.step(r, 100 * np.array(frac), where="mid", color=style.SERIES[2], lw=2)
    axq.axhline(15, color=style.INK_SECONDARY, lw=0.8, ls="--")
    axq.axvspan(edges[0], 460, color=style.COLOR_FIELD, alpha=0.3, lw=0)
    axq.set_xscale("log"); axq.set_xlabel("R  [arcsec]")
    axq.set_ylabel("stars passing the\nquality flag  [%]")
    for i, (rr, ff) in enumerate(zip(r, frac)):
        if ff < 0.2:
            axq.annotate("%.1f %%" % (100 * ff), (rr, 100 * ff), textcoords="offset points",
                         xytext=(0, 8), ha="center", fontsize=8, color=style.SERIES[2])
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_offset_explained(path: Path | str = "plots/outer_offset_explained.png",
                          k1: str = "K1_noDM_composite", distance_kpc: float = 5.43) -> Path:
    """Why the Gaia points sit above the model beyond 500 arcsec.

    Left: the deviation split into what each step contributes -- the fitted instrument scale,
    our measurement against the published profile, and what is left over. Right: the model's
    projected anisotropy against the measured one, which is where the leftover comes from at
    large radius.
    """
    import json

    from ..kinematics.report import _family_for
    style.apply()
    edges = np.array([460., 540, 630, 730, 850, 1000, 1200, 1500, 1900, 2400.])
    mix = our_mixture_profile(edges, distance_kpc=distance_kpc)
    pub = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    summary = json.loads((results_dir() / "fits" / k1 / "summary.json").read_text())
    fam = _family_for(summary)
    x = np.array([summary["parameters"][n]["ml"] for n in fam.names])
    jeans, D, scales = fam.build(fam.to_dict(x))
    s_g = scales.get("GaiaEDR3", 1.0)
    pc = distance_kpc * 1e3 / 206264.806
    r = np.asarray(mix["r_median"])
    ours = np.asarray(mix["sigma_pm"])
    total = np.sqrt(0.5 * (np.asarray(mix["sigma_pmr"]) ** 2 + np.asarray(mix["mean_pmr"]) ** 2
                           + np.asarray(mix["sigma_pmt"]) ** 2 + np.asarray(mix["mean_pmt"]) ** 2))
    err = np.asarray(mix["sigma_pm_err"])
    published = np.interp(r, np.asarray(pub["r"]), np.asarray(pub["sigma_pm"]))
    model = _sigma_1d_kms(jeans, D, r) / (KMS_PER_MASYR_KPC * D)

    fig, (ax, ab) = plt.subplots(1, 2, figsize=(13, 5))
    ax.axhline(0, color=style.INK, lw=1)
    ax.errorbar(r, 100 * (total / model - 1), yerr=100 * err / model, fmt="D", ms=6, color=style.INK,
                ecolor=style.INK, elinewidth=1.3, lw=0, label="total: our measurement vs the unscaled model")
    ax.errorbar(r, 100 * (total / (model * s_g) - 1), yerr=100 * err / (model * s_g), fmt="s", ms=5,
                color=style.SERIES[0], ecolor=style.SERIES[0], elinewidth=1.2, lw=0,
                label="after the fitted Gaia scale (%.3f)" % s_g)
    ax.plot(r, 100 * (ours / published - 1), "^--", ms=5, color=style.SERIES[1], lw=1.2,
            label="our measurement vs the published profile")
    ax.plot(r, 100 * (total / ours - 1), "v:", ms=5, color=style.SERIES[2], lw=1.2,
            label="what the mean motions add (total vs dispersion)")
    ax.axvspan(1500, 2400, color=style.SERIES[1], alpha=0.08, lw=0)
    ax.text(1540, 21, "shape mismatch\ngrows here", fontsize=8, color=style.SERIES[1])
    ax.set_xscale("log"); ax.set_xlabel("R  [arcsec]"); ax.set_ylabel("contribution to the offset  [%]")
    ax.set_title("What makes up the Gaia offset", fontsize=10.5); ax.legend(fontsize=7.5, loc="upper left")
    add_pc_axis(ax, distance_kpc)

    model_ratio = np.array([jeans.dispersions_kms(np.array([rr * pc]))["pmt"][0]
                            / jeans.dispersions_kms(np.array([rr * pc]))["pmr"][0] for rr in r])
    meas = np.asarray(mix["sigma_pmt"]) / np.asarray(mix["sigma_pmr"])
    merr = meas * np.hypot(np.asarray(mix["sigma_pmt_err"]) / np.asarray(mix["sigma_pmt"]),
                           np.asarray(mix["sigma_pmr_err"]) / np.asarray(mix["sigma_pmr"]))
    ab.plot(r, model_ratio, "-", color=style.SERIES[0], lw=2, label=r"K1 model ($\beta_\infty$ = %.2f, radial)"
            % summary["parameters"]["beta_inf"]["ml"])
    ab.errorbar(r, meas, yerr=merr, fmt="D", ms=6, color=style.INK, ecolor=style.INK, elinewidth=1.3, lw=0,
                label="measured")
    ab.axhline(1.0, color=style.INK_SECONDARY, lw=0.9, ls="--")
    ab.set_xscale("log"); ab.set_xlabel("R  [arcsec]"); ab.set_ylabel(r"$\sigma_T/\sigma_R$")
    ab.set_title("The model is radial where the data are tangential", fontsize=10.5)
    ab.legend(fontsize=8); add_pc_axis(ab, distance_kpc)
    ab.text(0.03, 0.08, "radial orbits", transform=ab.transAxes, fontsize=8, color=style.INK_SECONDARY)
    ab.text(0.03, 0.9, "tangential orbits", transform=ab.transAxes, fontsize=8, color=style.INK_SECONDARY)

    fig.suptitle("Why the Gaia dispersions sit above the model beyond 500 arcsec", y=0.99)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_datasets_unscaled(path: Path | str = "plots/datasets_unscaled.png", model_x: str | None = None,
                           distance_kpc: float = 5.43) -> Path:
    """Every dataset against one model, with no per-instrument rescaling anywhere.

    The reference is a K1 (no dark matter) model fitted with the instrument nuisances
    removed, so each dataset's disagreement is visible instead of being absorbed. Upper
    panel: the dispersions themselves; lower: the deviation of each dataset from the same
    model, in per cent.
    """
    import numpy as np

    from ..kinematics import FitProblem, KinematicData, NoDarkMatterModel
    style.apply()
    fam = NoDarkMatterModel(tracer="composite", instruments=())
    x = np.load(model_x or (results_dir() / "fits" / "K1_noscale_ml_x.npy"))
    jeans, D, _ = fam.build(fam.to_dict(x))
    data = _datasets_in_kms(D, contamination_modelled=True)
    # also show the published EDR3 profile -- that is the version the fits were given
    from ..kinematics.likelihood import load_profile
    pub = load_profile("gaia_edr3_pm", r_min_arcsec=300.0, n_max=None)
    keep = np.linspace(0, pub.n - 1, 14).astype(int)
    k = KMS_PER_MASYR_KPC * D
    data.append(dict(name="Gaia EDR3, published profile (what the fits were given)", r=pub.r[keep],
                     sigma=pub.value[keep] * k, err=0.5 * (pub.err_lo + pub.err_hi)[keep] * k, n=None,
                     kind="pm", color=style.SERIES_EXTRA, marker="v", streaming2=pub.streaming2[keep] * k**2))
    R = np.geomspace(2.0, 2400.0, 220)

    fig, (ax, axr, axg) = plt.subplots(3, 1, figsize=(9.5, 10.2), sharex=True,
                                       gridspec_kw={"height_ratios": [2.2, 1.6, 1.1], "hspace": 0.07})
    for kind, ls, lab in (("pm", "-", "proper motion"), ("los", "--", "line of sight")):
        ax.plot(R, _sigma_1d_kms(jeans, D, R, kind), color=style.SERIES[0], lw=2, ls=ls, alpha=0.9,
                label="K1 no dark matter, no instrument scales, %s" % lab)
    for d in data:
        ax.errorbar(d["r"], d["sigma"], yerr=d["err"], fmt=d["marker"], ms=4.5, color=d["color"],
                    ecolor=d["color"], elinewidth=1, capsize=0, lw=0, label=d["name"], zorder=5)
    ax.set_yscale("log"); ax.set_ylabel("1-D velocity dispersion  [km/s]")
    ax.legend(fontsize=7.5, ncol=2, loc="lower left"); ax.tick_params(labelbottom=False)
    ax.set_title("Every dataset against the same model, no rescaling anywhere", fontsize=11)
    add_pc_axis(ax, D)

    axr.axhline(0, color=style.INK, lw=1)
    for d in data:
        model = _sigma_1d_kms(jeans, D, d["r"], d["kind"])
        if d.get("streaming2") is not None:
            model = np.sqrt(np.maximum(model**2 - d["streaming2"], 1e-6))
        elif "Gaia" in d["name"]:
            from ..kinematics.likelihood import pm_rotation_curve
            rot = pm_rotation_curve(d["r"]) * KMS_PER_MASYR_KPC * D
            model = np.sqrt(np.maximum(model**2 - 0.5 * rot**2, 1e-6))
        axr.errorbar(d["r"], 100 * (d["sigma"] / model - 1), yerr=100 * d["err"] / model, fmt=d["marker"],
                     ms=4.5, color=d["color"], ecolor=d["color"], elinewidth=1, capsize=0, lw=0, label=d["name"])
    axr.axvspan(346, 500, color=style.SERIES[1], alpha=0.10, lw=0)
    axr.text(352, 17, "no reliable data", fontsize=7.5, color=style.SERIES[1])
    axr.set_ylabel("deviation from the model  [%]")
    axr.set_xscale("log"); axr.set_ylim(-25, 30); axr.legend(fontsize=7.5, ncol=2, loc="upper left")
    axr.tick_params(labelbottom=False)

    # do the two Gaia releases agree with each other?
    dr2 = Table.read(processed_dir() / "kinematics" / "baumgardt2019_ocen_pm_dispersion.ecsv")
    pubp = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    rr = np.asarray(dr2["r"]); d = np.asarray(dr2["sigma_pm"])
    de = 0.5 * (np.asarray(dr2["sigma_pm_err_lo"]) + np.asarray(dr2["sigma_pm_err_hi"]))
    e3 = np.interp(rr, np.asarray(pubp["r"]), np.asarray(pubp["sigma_pm"]))
    axg.errorbar(rr, 100 * (e3 / d - 1), yerr=100 * e3 / d * de / d, fmt="o", ms=5, color=style.SERIES[2],
                 ecolor=style.SERIES[2], elinewidth=1.2, lw=0, label="EDR3 (published) / DR2 (Baumgardt+ 2019)")
    axg.axhline(0, color=style.INK, lw=1)
    axg.set_xscale("log"); axg.set_xlabel("R  [arcsec]"); axg.set_ylabel("EDR3 vs DR2  [%]")
    axg.legend(fontsize=8, loc="upper left"); axg.set_ylim(-18, 18)
    axg.text(0.55, 0.08, "the two Gaia releases disagree with a radial trend, not a constant",
             transform=axg.transAxes, fontsize=8, color=style.INK_SECONDARY, ha="center")
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def hst_gaia_excess_table(edges=(0., 100., 200., 300., 460.), g_max: float = 17.0):
    """Per-star Gaia-minus-HST scatter, and the dispersion each reports from the same stars.

    HST proper motions are *relative* (locally corrected), so the difference is taken about
    its own median in every annulus and only widths are compared, never zero points. The
    matched sample contains no Gaia star that passes the astrometric quality flag, so this
    measures the **unflagged** Gaia catalogue in the crowded region.
    """
    from ..kinematics.outer_profile import dispersion_ml
    from ..selection.hst_gaia_match import load_match
    t = load_match()
    r = np.asarray(t["r_arcsec"]); g = np.asarray(t["g_mag"])
    sel = (np.asarray(t["membership_prob"]) > 0.9) & (np.asarray(t["hst_quality"]) > 0) & (g < g_max)
    cols = {k: np.asarray(t[k], float) for k in
            ("hst_pmra", "hst_pmdec", "gaia_pmra", "gaia_pmdec", "hst_pmra_error",
             "hst_pmdec_error", "gaia_pmra_error", "gaia_pmdec_error")}
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = sel & (r >= lo) & (r < hi)
        if m.sum() < 30:
            continue
        sh, sg, ex2, q2 = [], [], [], []
        for ax in ("pmra", "pmdec"):
            h, gg = cols["hst_" + ax][m], cols["gaia_" + ax][m]
            he, ge = cols["hst_%s_error" % ax][m], cols["gaia_%s_error" % ax][m]
            sh.append(dispersion_ml(h, he)[0] ** 2); sg.append(dispersion_ml(gg, ge)[0] ** 2)
            d = (gg - np.median(gg)) - (h - np.median(h))
            ex2.append(np.var(d) - np.median(he ** 2 + ge ** 2)); q2.append(np.median(he ** 2 + ge ** 2))
        n = int(m.sum())
        sh1, sg1 = np.sqrt(np.mean(sh)), np.sqrt(np.mean(sg))
        rows.append((lo, hi, float(np.median(r[m])), n, sh1, sg1, sg1 / sh1,
                     (sg1 / sh1) * np.sqrt(1.0 / n), np.sqrt(max(np.mean(ex2), 0.0)), np.sqrt(np.mean(q2))))
    return Table(rows=rows, names=("r_lower", "r_upper", "r_median", "n_stars", "sigma_hst",
                                   "sigma_gaia", "ratio", "ratio_err", "excess_noise", "quoted"))


def hst_gaia_overlap_table(lo: float = 300.0, hi: float = 380.0, distance_kpc: float = 5.43):
    """HST against **quality-flagged** Gaia in the one annulus where both are usable.

    HST's high-quality astrometry ends at 380 arcsec and Gaia's quality flag passes no star
    inside 300 arcsec, so 300-380 arcsec is the entire overlap. The two samples are not the
    same stars here (no flagged Gaia star has an HST counterpart), so this compares the
    dispersion each instrument measures in the same annulus.
    """
    from astropy.table import Table as _T
    from ..kinematics.outer_profile import dispersion_ml, load_members
    from ..paths import raw_dir

    def comb(v1, e1, v2, e2):
        a = dispersion_ml(v1, e1); b = dispersion_ml(v2, e2)
        s = np.sqrt(0.5 * (a[0] ** 2 + b[0] ** 2))
        return s, 0.5 * np.hypot(a[1], b[1])

    hst = _T.read(raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits")
    hr = np.hypot((np.asarray(hst["RA"], float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC)),
                  np.asarray(hst["DEC"], float) - OCEN_DEC) * 3600.0
    hm = (hr >= lo) & (hr < hi) & (np.asarray(hst["selection_hq_astrometry_and_membership"], float) > 0) \
        & np.isfinite(np.asarray(hst["pmra_corrected"], float))
    sh, esh = comb(np.asarray(hst["pmra_corrected"], float)[hm], np.asarray(hst["pmra_corrected_err"], float)[hm],
                   np.asarray(hst["pmdec_corrected"], float)[hm], np.asarray(hst["pmdec_corrected_err"], float)[hm])
    s = load_members(exact=True, distance_kpc=distance_kpc)
    base = (s.r_arcsec >= lo) & (s.r_arcsec < hi) & (s.prob > 0.9)
    rows = [("HST, high-quality astrometry", int(hm.sum()), sh, esh, 1.0, 0.0)]
    for label, m in (("Gaia, quality flag", base & ((s.quality_flag & QUALITY_BIT) > 0)),
                     ("Gaia, P > 0.99 and G < 17", base & (s.prob > 0.99) & (s.g_mag < 17)),
                     ("Gaia, G < 16", base & (s.g_mag < 16)),
                     ("Gaia, G < 17", base & (s.g_mag < 17)),
                     ("Gaia, no quality cut", base)):
        if m.sum() < 5:
            continue
        sg, esg = comb(s.mu_r[m], s.err_r[m], s.mu_t[m], s.err_t[m])
        rows.append((label, int(m.sum()), sg, esg, sg / sh, (sg / sh) * np.hypot(esg / sg, esh / sh)))
    return Table(rows=rows, names=("sample", "n_stars", "sigma", "sigma_err", "ratio", "ratio_err"))


def plot_hst_gaia_star_by_star(path: Path | str = "plots/hst_gaia_star_by_star.png",
                               distance_kpc: float = 5.43) -> Path:
    """Do the two proper-motion catalogues agree?

    Top left: the per-star Gaia-minus-HST difference against what the quoted errors allow.
    Top right: the excess scatter this implies, as a function of radius. Bottom left: the
    dispersion each instrument reports from the identical stars, with the quality-flagged
    comparison in the one annulus where both catalogues are usable. Bottom right: the
    fraction of Gaia stars that survive the astrometric quality flag, which is what empties
    the inner region.
    """
    from ..kinematics.outer_profile import load_members
    from ..selection.hst_gaia_match import load_match
    style.apply()
    t = load_match()
    r = np.asarray(t["r_arcsec"]); g = np.asarray(t["g_mag"])
    sel0 = (np.asarray(t["membership_prob"]) > 0.9) & (np.asarray(t["hst_quality"]) > 0) & (g < 17)
    he = 0.5 * (np.asarray(t["hst_pmra_error"]) + np.asarray(t["hst_pmdec_error"]))
    ge = 0.5 * (np.asarray(t["gaia_pmra_error"]) + np.asarray(t["gaia_pmdec_error"]))
    da = np.asarray(t["gaia_pmra"]) - np.median(np.asarray(t["gaia_pmra"])[sel0]) \
        - (np.asarray(t["hst_pmra"]) - np.median(np.asarray(t["hst_pmra"])[sel0]))

    fig, ((a1, a2), (a3, a4)) = plt.subplots(2, 2, figsize=(11.6, 8.6))
    bins = np.linspace(-4, 4, 65); mid = 0.5 * (bins[1:] + bins[:-1])
    a1.hist(da[sel0], bins=bins, color=style.COLOR_FIELD, alpha=0.8, label="Gaia $-$ HST, same stars")
    rms = float(np.std(da[sel0]))
    exp = float(np.hypot(np.median(he[sel0]), np.median(ge[sel0])))
    n, dx = int(sel0.sum()), bins[1] - bins[0]
    a1.plot(mid, n * dx * np.exp(-0.5 * (mid / exp) ** 2) / (exp * np.sqrt(2 * np.pi)), color=style.SERIES[1],
            lw=2, label="allowed by the quoted errors (%.2f mas/yr)" % exp)
    a1.plot(mid, n * dx * np.exp(-0.5 * (mid / rms) ** 2) / (rms * np.sqrt(2 * np.pi)), color=style.SERIES[0],
            lw=2, ls="--", label="observed (%.2f mas/yr, %.1f$\\times$ wider)" % (rms, rms / exp))
    a1.set_yscale("log"); a1.set_ylim(0.5, 4.0 * n * dx / (exp * np.sqrt(2 * np.pi)))
    a1.set_xlabel(r"$\mu_{\alpha*}$ difference  [mas/yr]"); a1.set_ylabel("stars per bin")
    a1.set_title("%d matched stars, G < 17, median separation %.3f arcsec"
                 % (n, np.median(np.asarray(t["separation_arcsec"])[sel0])), fontsize=9.5)
    a1.legend(fontsize=7.5)

    ex = hst_gaia_excess_table()
    a2.plot(ex["r_median"], ex["excess_noise"], "o-", color=style.SERIES[0], lw=2, ms=7,
            label="Gaia scatter beyond its quoted errors")
    a2.plot(ex["r_median"], ex["quoted"], "s--", color=style.INK_SECONDARY, lw=1.5, ms=5,
            label="quoted (HST and Gaia combined)")
    for x_, y_ in zip(ex["r_median"], ex["excess_noise"]):
        a2.annotate("%.2f" % y_, (x_, y_), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=8)
    a2.set_xscale("log"); a2.set_xlabel("R  [arcsec]"); a2.set_ylabel("PM scatter  [mas/yr]")
    a2.set_ylim(0, None); a2.legend(fontsize=8)
    a2.set_title("Undeclared Gaia scatter falls outwards: crowding", fontsize=9.5)
    add_pc_axis(a2, distance_kpc)

    a3.errorbar(ex["r_median"], ex["sigma_hst"], fmt="o-", color=style.SERIES[0], lw=1.8, ms=6, label="HST")
    a3.errorbar(ex["r_median"], ex["sigma_gaia"], fmt="D-", color=style.INK, lw=1.8, ms=6,
                label="Gaia, no quality cut")
    for x_, a_, b_ in zip(ex["r_median"], ex["sigma_hst"], ex["sigma_gaia"]):
        a3.annotate("%+.0f %%" % (100 * (b_ / a_ - 1)), (x_, b_), textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8, color=style.SERIES[1])
    ov = hst_gaia_overlap_table()
    row = ov[ov["sample"] == "Gaia, quality flag"]
    if len(row):
        a3.errorbar([340.0], [row["sigma"][0]], yerr=[row["sigma_err"][0]], fmt="*", color=style.SERIES[2],
                    ms=17, lw=1.6, capsize=3, label="Gaia, quality flag (300-380 arcsec)")
        a3.annotate("ratio %.2f $\\pm$ %.2f\nfrom %d stars" % (row["ratio"][0], row["ratio_err"][0], row["n_stars"][0]),
                    (340.0, row["sigma"][0]), textcoords="offset points", xytext=(-6, -30), ha="right",
                    fontsize=8, color=style.SERIES[2])
    a3.set_xscale("log"); a3.set_xlabel("R  [arcsec]"); a3.set_ylabel("1-D PM dispersion  [mas/yr]")
    a3.set_title("Dispersion from the identical stars", fontsize=9.5)
    a3.legend(fontsize=8, loc="upper right")
    a3.margins(y=0.18)
    add_pc_axis(a3, distance_kpc)

    s = load_members(); qq = s.select((s.quality_flag & QUALITY_BIT) > 0)
    e2 = np.array([0., 100, 200, 300, 380, 460, 700, 1000, 1500, 2400.])
    frac = [((qq.r_arcsec >= lo) & (qq.r_arcsec < hi)).sum() / max(((s.r_arcsec >= lo) & (s.r_arcsec < hi)).sum(), 1)
            for lo, hi in zip(e2[:-1], e2[1:])]
    mid2 = np.sqrt(np.maximum(e2[:-1], 1) * e2[1:])
    a4.step(mid2, 100 * np.array(frac), where="mid", color=style.SERIES[2], lw=2)
    a4.axvspan(300, 380, color=style.SERIES[2], alpha=0.12)
    a4.set_xscale("log"); a4.set_xlabel("R  [arcsec]"); a4.set_ylabel("Gaia stars passing the quality flag  [%]")
    a4.set_title("Why EDR3 stops: crowding", fontsize=9.5)
    a4.text(115, 40, "no star passes\ninside 200 arcsec", fontsize=8, color=style.SERIES[2])
    a4.text(330, 8, "only overlap\nwith usable HST", fontsize=7.5, color=style.SERIES[2], ha="center")
    add_pc_axis(a4, distance_kpc)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_pm_datasets(path: Path | str = "plots/pm_datasets_edr3_rebuild.png",
                     distance_kpc: float = 5.43) -> Path:
    """Every proper-motion dispersion measurement on one axis, and what changed.

    Top: HST and Gaia DR2 as published, the published Gaia EDR3 spline with its confidence
    band, our first (old) EDR3 measurement which deconvolved the raw catalogue errors, and
    our new one built from stars whose errors are small next to the signal. The starred pair
    is the 300-380 arcsec annulus, the only place where HST's high-quality astrometry and
    Gaia's quality flag both have usable stars, so it is the one direct instrument test.
    Bottom: the same, divided by the published EDR3 spline.
    """
    from ..kinematics.likelihood import load_profile
    from ..kinematics.outer_gaia import load_edr3_profile
    from ..kinematics.vb2021_replication import published_profile
    style.apply()
    rp, sp = published_profile()
    t = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    lo_b, hi_b = np.asarray(t["sigma_pm_p16"], float), np.asarray(t["sigma_pm_p84"], float)
    old = our_mixture_profile()
    new = load_edr3_profile()
    ov = hst_gaia_overlap_table()
    hst, dr2 = load_profile("hst_pm_combined"), load_profile("gaia_dr2_pm")

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.6, 8.4), sharex=True,
                                 gridspec_kw={"height_ratios": [2.3, 1]})
    a1.fill_between(rp[rp > 0], lo_b[rp > 0], hi_b[rp > 0], color=style.INK_SECONDARY, alpha=0.18, lw=0)
    a1.plot(rp[rp > 0], sp[rp > 0], color=style.INK_SECONDARY, lw=2,
            label="Gaia EDR3, published spline (Vasiliev & Baumgardt 2021)")
    a1.axvspan(0.5, 300, color=style.SERIES[1], alpha=0.06, lw=0)
    a1.errorbar(hst.r, hst.value, yerr=[hst.err_lo, hst.err_hi], fmt="o", ms=4,
                color=style.SERIES[1], lw=1, label="HST (oMEGACat)")
    a1.errorbar(dr2.r, dr2.value, yerr=[dr2.err_lo, dr2.err_hi], fmt="s", ms=5,
                color=style.SERIES_EXTRA, lw=1, label="Gaia DR2 (Baumgardt+ 2019)")
    a1.errorbar(old["r_median"], old["sigma_pm"], yerr=old["sigma_pm_err"], fmt="v--", ms=6,
                color=style.COLOR_FIELD, lw=1.3, label="ours, old: all quality stars, raw catalogue errors")
    a1.errorbar(new["r_median"], new["sigma_pm"], yerr=new["sigma_pm_err"], fmt="D-", ms=7,
                color=style.SERIES[0], lw=2, capsize=3, zorder=5,
                label="ours, new: err < 0.4$\\sigma$, error-model independent, R > 300\"")
    row_h = ov[ov["sample"] == "HST, high-quality astrometry"][0]
    row_g = ov[ov["sample"] == "Gaia, quality flag"][0]
    a1.errorbar([318.0], [row_h["sigma"]], yerr=[row_h["sigma_err"]], fmt="*", ms=19,
                color=style.SERIES[1], zorder=6, label="HST, quality cut, 300-380\"  (N = %d)" % row_h["n_stars"])
    a1.errorbar([338.0], [row_g["sigma"]], yerr=[row_g["sigma_err"]], fmt="*", ms=19,
                color=style.SERIES[2], zorder=6, capsize=3,
                label="Gaia EDR3, quality cut, same annulus  (N = %d)" % row_g["n_stars"])
    a1.annotate("same stars' quality standard,\nGaia/HST = %.2f $\\pm$ %.2f" % (row_g["ratio"], row_g["ratio_err"]),
                (338.0, row_g["sigma"]), textcoords="offset points", xytext=(10, -56), fontsize=8,
                color=style.SERIES[2])
    a1.text(165, 0.80, "the quality flag passes\n5 Gaia stars in total\ninside 300 arcsec",
            fontsize=8.5, color=style.SERIES[1], ha="center")
    a1.annotate("53 and 441 flagged stars,\nerrors only 5-9 % of the signal",
                (415.0, 0.78), fontsize=8, color=style.SERIES[0], ha="center")
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_ylabel("1-D PM dispersion  [mas/yr]")
    a1.set_ylim(0.18, 0.95)
    a1.legend(fontsize=7.6, loc="lower left", framealpha=0.95, borderpad=0.6)
    a1.set_title("Proper-motion dispersion of $\\omega$ Cen: what enters the fits", fontsize=11)
    add_pc_axis(a1, distance_kpc)

    def ratio(r, v, e):
        ref = np.interp(r, rp, sp)
        return np.asarray(v) / ref, np.asarray(e) / ref

    a2.axhline(1.0, color=style.INK_SECONDARY, lw=2)
    a2.fill_between(rp[rp > 0], lo_b[rp > 0] / sp[rp > 0], hi_b[rp > 0] / sp[rp > 0],
                    color=style.INK_SECONDARY, alpha=0.18, lw=0)
    for r_, v_, e_, fmt, col, lab in (
            (hst.r, hst.value, hst.err_lo, "o", style.SERIES[1], "HST"),
            (dr2.r, dr2.value, dr2.err_lo, "s", style.SERIES_EXTRA, "Gaia DR2"),
            (old["r_median"], old["sigma_pm"], old["sigma_pm_err"], "v--", style.COLOR_FIELD, "ours, old"),
            (new["r_median"], new["sigma_pm"], new["sigma_pm_err"], "D-", style.SERIES[0], "ours, new")):
        y, ye = ratio(np.asarray(r_, float), v_, e_)
        a2.errorbar(r_, y, yerr=ye, fmt=fmt, ms=5 if fmt[0] != "D" else 7, color=col,
                    lw=1.8 if fmt[0] == "D" else 1.0, label=lab, zorder=5 if fmt[0] == "D" else 3)
    a2.axvspan(0.5, 300, color=style.SERIES[1], alpha=0.06, lw=0)
    a2.set_xscale("log"); a2.set_xlabel("R  [arcsec]")
    a2.set_ylabel("ratio to published\nEDR3 spline")
    a2.set_ylim(0.80, 1.30); a2.set_xlim(100, 2700)
    a2.legend(fontsize=8, ncol=4, loc="lower left", framealpha=0.95)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_periphery(path: Path | str = "plots/periphery_where_the_cluster_ends.png",
                   distance_kpc: float = 5.43) -> Path:
    """Where does omega Cen stop? The profile from the centre out to 400 pc.

    The Gaia member catalogue is truncated at 63 pc, which is a retrieval radius and not a
    boundary, and the Jacobi radius is larger. Beyond that edge the Pristine periphery sample
    carries proper motions whose membership is itself proper-motion based, so it is shown as
    a band between two fixed selection windows rather than as points; the spectroscopic
    line-of-sight velocities have no such selection and are the honest probe.
    """
    from ..kinematics.likelihood import load_profile
    from ..kinematics.outer_gaia import load_edr3_profile
    from ..kinematics.periphery import (GAIA_EDGE_DEG, KMS_PER_MASYR_KPC, R_JACOBI_NOW_PC,
                                        R_JACOBI_PERI_PC, periphery_los_profile,
                                        periphery_pm_profile)
    style.apply()
    k = KMS_PER_MASYR_KPC * distance_kpc
    pc = lambda arcsec: np.asarray(arcsec, float) * distance_kpc * 1e3 / 206264.806

    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    from ..kinematics.hst_profile import hst_profile
    from .style import dataset_label, dataset_style
    hst = load_profile("hst_pm_combined")
    ax.errorbar(pc(hst.r), np.asarray(hst.value) * k, yerr=np.asarray(hst.err_lo) * k,
                linestyle="none", lw=1, label="HST, published oMEGACat profile",
                **dataset_style("hst", size=3.5, fitted=False))
    # our own HST measurement carries the profile from 300 to 340 arcsec, where the
    # published one stops, so that it genuinely overlaps Gaia rather than stopping short
    ours = hst_profile(edges_arcsec=tuple(np.concatenate(
        [np.geomspace(2.0, 150.0, 12), [200., 250., 300., 340.]])))
    ax.errorbar(pc(ours["r_median"]), np.asarray(ours["sigma_pm"]) * k,
                yerr=np.asarray(ours["sigma_pm_err"]) * k, linestyle="-", lw=2, capsize=3,
                zorder=6, label=dataset_label("hst") + ", ours (flagged stars)",
                **dataset_style("hst", size=7.5))
    g = load_edr3_profile()
    ax.errorbar(pc(g["r_median"]), np.asarray(g["sigma_pm"]) * k,
                yerr=np.asarray(g["sigma_pm_err"]) * k, linestyle="-", lw=2, capsize=3,
                label=dataset_label("gaia_edr3"), **dataset_style("gaia_edr3", size=7))
    # the one annulus where both surveys have quality-selected stars
    tests = hst_gaia_overlap_tests(distance_kpc)
    clean = tests[tests["flagged_only"]][0]
    ax.errorbar([pc(clean["r_gaia"])], [clean["sigma_gaia"] * k],
                yerr=[clean["sigma_gaia_err"] * k], linestyle="none", lw=1.8, capsize=4,
                zorder=7, label="Gaia, same annulus as HST's last point (%d stars)" % clean["n_gaia"],
                **dataset_style("gaia_edr3", size=11))
    lo_pc, hi_pc = pc(clean["r_lower"]), pc(clean["r_upper"])
    ax.axvspan(lo_pc, hi_pc, color=style.SERIES[2], alpha=0.22, lw=0)
    ax.annotate("both surveys\nmeasure here", (np.sqrt(lo_pc * hi_pc), 21.0), fontsize=8,
                ha="center", va="top", color=style.SERIES[2])
    lo_w, hi_w = periphery_pm_profile(0.8), periphery_pm_profile(1.2)
    ax.fill_between(lo_w["r_pc"], lo_w["sigma_kms"], hi_w["sigma_kms"],
                    color=style.SERIES[2], alpha=0.20, lw=0,
                    label="Pristine, between 0.8 and 1.2 mas/yr windows")
    ax.plot(lo_w["r_pc"], lo_w["sigma_kms"], color=style.INK_SECONDARY, lw=1.2, ls=":")
    ax.plot(hi_w["r_pc"], hi_w["sigma_kms"], color=style.INK_SECONDARY, lw=1.2, ls=":")
    los = periphery_los_profile()
    ax.errorbar(los["r_pc"], los["sigma_kms"], yerr=los["sigma_kms_err"], linestyle="none",
                lw=1.8, capsize=4, zorder=6, label=dataset_label("spectro", fitted=False),
                **dataset_style("spectro", size=13, fitted=False))
    for x_, lab, ls in ((pc(GAIA_EDGE_DEG * 3600), "Gaia catalogue edge", "-"),
                        (R_JACOBI_PERI_PC, "$r_J$ at pericentre", "--"),
                        (R_JACOBI_NOW_PC, "$r_J$ now", ":")):
        ax.axvline(x_, color=style.INK_SECONDARY, lw=1.2, ls=ls, alpha=0.8)
        ax.text(x_ * 1.03, 21.5, lab, rotation=90, fontsize=8, color=style.INK_SECONDARY, va="top")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("r  [pc]"); ax.set_ylabel("velocity dispersion  [km/s]")
    ax.set_xlim(0.2, 500); ax.set_ylim(3.5, 24)
    ax.set_yticks([4, 6, 8, 10, 15, 20]); ax.set_yticklabels(["4", "6", "8", "10", "15", "20"])
    ax.legend(fontsize=8.2, loc="lower left")
    ax.set_title("Where does $\\omega$ Cen stop? The dispersion flattens near the Jacobi radius",
                 fontsize=11)
    add_arcsec_axis(ax, distance_kpc)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_periphery_density(path: Path | str = "plots/periphery_overdensity.png",
                           distance_kpc: float = 5.43) -> Path:
    """Are the stars measured in the periphery actually there? A star-count test.

    Left: the sky distribution of cluster-like stars, selected on proper motion and CaHK
    metallicity, with the spectroscopic members overplotted and the preferred axis of the
    outermost candidates drawn. Right: their background-subtracted surface density against
    radius, with the background measured from eight control windows in proper-motion space.
    """
    from ..kinematics.periphery import (GAIA_EDGE_DEG, R_JACOBI_PERI_PC, periphery_los_profile)
    from ..selection.periphery_density import (TAIL_AXIS_PA_DEG, along_across, axial_rayleigh,
                                               candidate_masks, density_profile, load_periphery)
    style.apply()
    t, r, th = load_periphery()
    sig, _ = candidate_masks(t)
    prof = density_profile()
    rj_deg = R_JACOBI_PERI_PC / (distance_kpc * 1e3) * 180 / np.pi

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.2, 5.9))
    x = (np.asarray(t["ra"], float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC))
    y = np.asarray(t["dec"], float) - OCEN_DEC
    far = r > 5.2
    a1.scatter(x[sig & ~far], y[sig & ~far], s=14, color=style.SERIES[0], alpha=0.85,
               label="cluster-like: PM within 0.8 mas/yr and [Fe/H] < $-$1.2")
    spec = Table.read(processed_dir() / "tails" / "kuzma2026_spectroscopy.ecsv")
    sx = (np.asarray(spec["ra"], float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC))
    sy = np.asarray(spec["dec"], float) - OCEN_DEC
    sm = np.asarray(spec["member_flag"]) == "True"
    a1.scatter(sx[sm], sy[sm], s=34, facecolors="none", edgecolors=style.SERIES[2], lw=1.4,
               label="spectroscopic members")
    for rad, lab, ls in ((GAIA_EDGE_DEG, "Gaia edge", "-"), (rj_deg, "$r_J$ pericentre", "--")):
        a1.add_patch(plt.Circle((0, 0), rad, fill=False, color=style.INK_SECONDARY, lw=1.3, ls=ls))
        a1.annotate(lab, (0, rad), fontsize=8, color=style.INK_SECONDARY, ha="center", va="bottom")
    ang = np.radians(TAIL_AXIS_PA_DEG)
    a1.plot([-5.2 * np.cos(ang), 5.2 * np.cos(ang)], [-5.2 * np.sin(ang), 5.2 * np.sin(ang)],
            color=style.SERIES[1], lw=1.6, ls="-.", alpha=0.8,
            label="preferred axis of the outer candidates (PA %.0f$^\\circ$)" % TAIL_AXIS_PA_DEG)
    a1.set_xlim(5.4, -5.4); a1.set_ylim(-5.4, 5.4); a1.set_aspect("equal")
    a1.set_xlabel(r"$\Delta\alpha\cos\delta$  [deg]"); a1.set_ylabel(r"$\Delta\delta$  [deg]")
    a1.legend(fontsize=7.6, loc="upper left"); a1.set_title("Where the candidates are", fontsize=10)

    ok = np.asarray(prof["excess"]) > 0
    a2.errorbar(np.asarray(prof["r_pc"])[ok], np.asarray(prof["excess"])[ok],
                yerr=np.asarray(prof["excess_err"])[ok], fmt="o-", ms=7, lw=1.8,
                color=style.SERIES[0], capsize=3, label="excess over the control windows")
    for row in prof:
        if row["significance"] > 1.5:
            a2.annotate("%.1f$\\sigma$" % row["significance"], (row["r_pc"], row["excess"]),
                        textcoords="offset points", xytext=(6, 8), fontsize=8, color=style.SERIES[0])
    a2.axhline(0, color=style.INK_SECONDARY, lw=1)
    for i, row in enumerate(periphery_los_profile()):
        a2.axvline(row["r_pc"], color=style.SERIES[2], lw=1.1, ls=":", alpha=0.9)
        a2.annotate("%d spec.\nstars" % row["n_stars"], (row["r_pc"] * 0.97, (150, 40, 150)[i % 3]),
                    fontsize=7.5, color=style.SERIES[2], ha="right", rotation=0)
    a2.axvline(R_JACOBI_PERI_PC, color=style.INK_SECONDARY, lw=1.3, ls="--")
    a2.annotate("$r_J$ pericentre", (R_JACOBI_PERI_PC * 1.04, 0.8), rotation=90, fontsize=8,
                color=style.INK_SECONDARY)
    aa = along_across()
    a2.annotate("beyond 133 pc the circular average is empty,\nbut along PA %.0f$^\\circ$ the excess is "
                "%.2f $\\pm$ %.2f /deg$^2$ (%.1f$\\sigma$)\nagainst %.2f across it"
                % (TAIL_AXIS_PA_DEG, aa["excess"][0], aa["excess_err"][0], aa["significance"][0],
                   aa["excess"][1]), (0.97, 0.62), xycoords="axes fraction", ha="right",
                fontsize=8, color=style.SERIES[1])
    a2.set_xscale("log"); a2.set_yscale("log"); a2.set_ylim(0.05, 300); a2.set_xlim(45, 520)
    a2.set_xlabel("r  [pc]"); a2.set_ylabel("excess surface density  [stars deg$^{-2}$]")
    a2.legend(fontsize=8, loc="upper right")
    a2.set_title("The cluster truncates near the Jacobi radius", fontsize=10)
    add_arcsec_axis(a2, distance_kpc, unit="arcmin")
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_extended_profile(path: Path | str = "plots/profile_extended_pristine.png",
                          distance_kpc: float = 5.43) -> Path:
    """The dispersion profile extended outwards with Pristine, under our own mixture model.

    Left: HST, our Gaia EDR3 measurement, the Pristine mixture measurement with no
    proper-motion selection, the same catalogue read through a fixed proper-motion window for
    contrast, and the spectroscopic line-of-sight points. Right: Pristine over Gaia on the
    identical annuli, which is a like-for-like test of two independent catalogues, selections
    and field models.
    """
    from ..kinematics.outer_gaia import load_edr3_profile
    from ..kinematics.periphery import (KMS_PER_MASYR_KPC, R_JACOBI_PERI_PC,
                                        periphery_los_profile, periphery_pm_profile)
    from ..kinematics.likelihood import load_profile
    from ..kinematics.pristine_profile import pristine_profile
    style.apply()
    k = KMS_PER_MASYR_KPC * distance_kpc
    pc = lambda a: np.asarray(a, float) * distance_kpc * 1e3 / 206264.806
    g = load_edr3_profile()
    edges = np.array([300., 380., 460., 582., 737., 934., 1182., 1497., 1895., 2413.])
    matched = pristine_profile(edges_deg=tuple(edges / 3600.0), min_stars=25)
    wide = pristine_profile()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.4, 5.8),
                                 gridspec_kw={"width_ratios": [1.5, 1]})
    from .style import dataset_label, dataset_style
    hst = load_profile("hst_pm_combined")
    a1.errorbar(pc(hst.r), np.asarray(hst.value) * k, yerr=np.asarray(hst.err_lo) * k,
                linestyle="none", lw=0.9, label=dataset_label("hst"),
                **dataset_style("hst", size=4.0))
    a1.errorbar(pc(g["r_median"]), np.asarray(g["sigma_pm"]) * k,
                yerr=np.asarray(g["sigma_pm_err"]) * k, linestyle="-", lw=2, capsize=3,
                label=dataset_label("gaia_edr3"), **dataset_style("gaia_edr3", size=6.5))
    ok = np.asarray(wide["reliable"], bool)
    a1.errorbar(np.asarray(wide["r_pc"])[ok], np.asarray(wide["sigma_kms"])[ok],
                yerr=np.asarray(wide["sigma_kms_err"])[ok], linestyle="-", lw=2, capsize=3,
                zorder=6, label=dataset_label("pristine") + ", same mixture, no PM cut",
                **dataset_style("pristine", size=9))
    a1.errorbar(np.asarray(wide["r_pc"])[~ok], np.asarray(wide["sigma_kms"])[~ok],
                yerr=np.asarray(wide["sigma_kms_err"])[~ok], linestyle="none", lw=1.1,
                label="Pristine, fewer than 10 cluster stars: not a measurement",
                **dataset_style("pristine", size=8, fitted=False))
    w = periphery_pm_profile(0.8)
    a1.plot(w["r_pc"], w["sigma_kms"], ":", color=style.SERIES[2], lw=1.6, alpha=0.8,
            label="Pristine read through a fixed 0.8 mas/yr window")
    los = periphery_los_profile()
    a1.errorbar(los["r_pc"], los["sigma_kms"], yerr=los["sigma_kms_err"], linestyle="none",
                lw=1.6, capsize=4, zorder=7, label=dataset_label("spectro", fitted=False),
                **dataset_style("spectro", size=13, fitted=False))
    a1.axvline(R_JACOBI_PERI_PC, color=style.INK_SECONDARY, lw=1.2, ls="--")
    a1.annotate("$r_J$ pericentre", (R_JACOBI_PERI_PC * 1.05, 4.2), rotation=90, fontsize=8,
                color=style.INK_SECONDARY)
    a1.axvspan(100, 260, color=style.COLOR_FIELD, alpha=0.22, lw=0)
    a1.annotate("cluster fraction below 5 %:\nthe mixture has nothing to fit",
                (160, 17), fontsize=8, color=style.INK_SECONDARY, ha="center")
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_xlim(0.5, 300); a1.set_ylim(3, 24)
    a1.set_yticks([4, 6, 8, 10, 15, 20]); a1.set_yticklabels(["4", "6", "8", "10", "15", "20"])
    a1.set_xlabel("r  [pc]"); a1.set_ylabel("velocity dispersion  [km/s]")
    a1.legend(fontsize=7.6, loc="lower left")
    a1.set_title("Profile extended outwards with Pristine", fontsize=10.5)
    add_arcsec_axis(a1, distance_kpc, unit="arcmin")

    gr = pc(g["r_median"]); ratio, rerr = [], []
    for row in matched:
        j = int(np.argmin(np.abs(gr - row["r_pc"])))
        gs, ge = float(g["sigma_pm"][j]), float(g["sigma_pm_err"][j])
        ps, pe = float(row["sigma_pm"]), float(row["sigma_pm_err"])
        ratio.append(ps / gs); rerr.append((ps / gs) * np.hypot(pe / ps, ge / gs))
    ratio, rerr = np.array(ratio), np.array(rerr)
    wmean = np.sum(ratio / rerr ** 2) / np.sum(1 / rerr ** 2)
    wsig = 1 / np.sqrt(np.sum(1 / rerr ** 2))
    a2.axhline(1.0, color=style.INK_SECONDARY, lw=2)
    a2.axhspan(wmean - wsig, wmean + wsig, color=style.SERIES[2], alpha=0.18, lw=0)
    a2.axhline(wmean, color=style.SERIES[2], lw=1.6, ls="--",
               label="weighted mean %.3f $\\pm$ %.3f" % (wmean, wsig))
    a2.errorbar(matched["r_pc"], ratio, yerr=rerr, linestyle="none", lw=1.5, capsize=3,
                **dataset_style("pristine", size=8))
    a2.set_xscale("log"); a2.set_xlabel("r  [pc]")
    a2.set_ylabel("Pristine / Gaia, identical annuli")
    a2.set_ylim(0.78, 1.22); a2.legend(fontsize=8.5, loc="upper left")
    a2.set_title("Two catalogues, two field models, one answer", fontsize=10.5)
    add_arcsec_axis(a2, distance_kpc, unit="arcmin")
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def plot_master_datasets(path: Path | str = "plots/master_datasets.png",
                         datasets: str | None = None, distance_kpc: float = 5.43) -> Path:
    """Everything the likelihood is shown, in one figure, in km/s.

    One symbol and one colour per dataset throughout (:data:`ocen_dm.plotting.style.
    DATASET_STYLE`); the radial and tangential components of the same stars share that
    symbol and are told apart by the marker fill, left for radial and right for tangential,
    so they read as one measurement split rather than as two instruments. Data shown for
    context but not fitted keep their symbol and are lightened.

    Top: every dataset in the default likelihood, converted to km/s at the fitted distance
    so proper motions and line-of-sight velocities share an axis. Middle: the anisotropy
    implied by the components. Bottom: what each dataset covers and how many points it
    contributes.
    """
    from ..cli import DEFAULT_DATASETS
    from ..kinematics.likelihood import load_profile
    from ..kinematics.outer_gaia import load_edr3_profile
    from ..kinematics.periphery import (KMS_PER_MASYR_KPC, R_JACOBI_PERI_PC,
                                        periphery_los_profile)
    from ..kinematics.pristine_profile import pristine_profile
    from ..kinematics.vb2021_replication import published_profile
    from .style import dataset_label, dataset_style
    style.apply()
    keys = (datasets or DEFAULT_DATASETS).split(",")
    k = KMS_PER_MASYR_KPC * distance_kpc
    pc = lambda a: np.asarray(a, float) * distance_kpc * 1e3 / 206264.806
    g = load_edr3_profile()
    pr = pristine_profile()
    ok = np.asarray(pr["reliable"], bool)

    def draw(ax, x, y, yerr, key, component="combined", fitted=True, size=6.5, ls="none", **kw):
        st = dataset_style(key, component, size=size, fitted=fitted)
        return ax.errorbar(x, y, yerr=yerr, linestyle=ls, lw=1.4, capsize=2,
                           label=dataset_label(key, component, fitted), zorder=5 if fitted else 3,
                           **st, **kw)

    fig, (a1, a3, a2) = plt.subplots(3, 1, figsize=(10.6, 11.6), sharex=True,
                                     gridspec_kw={"height_ratios": [2.5, 1.15, 1.0]})

    # --- context: real data that exists but is not fitted --------------------------
    rp, sp = published_profile()
    a1.plot(pc(rp[rp > 0]), sp[rp > 0] * k, color=style.SERIES[0], lw=1.4, ls="--", alpha=0.5,
            label="Gaia EDR3 published spline [not fitted]")
    draw(a1, np.asarray(pr["r_pc"])[ok], np.asarray(pr["sigma_kms"])[ok],
         np.asarray(pr["sigma_kms_err"])[ok], "pristine", fitted=False, size=8)
    los = periphery_los_profile()
    draw(a1, los["r_pc"], los["sigma_kms"], los["sigma_kms_err"], "spectro", fitted=False, size=13)

    # --- the likelihood's own data -------------------------------------------------
    spans = []
    for key in keys:
        p = load_profile(key)
        scale = 1.0 if p.kind == "los" else k
        draw(a1, pc(p.r), np.asarray(p.value) * scale,
             [np.asarray(p.err_lo) * scale, np.asarray(p.err_hi) * scale], key,
             size=4.5 if p.instrument in ("HST", "MUSE") else 7.0)
        st = dataset_style(key)
        spans.append((dataset_label(key), pc(p.r).min(), pc(p.r).max(), len(p.r), st["color"]))

    # our Gaia measurement, split into the components the likelihood does not receive
    if "gaia_edr3_ours" in keys:
        for col, comp in (("sigma_pmr", "radial"), ("sigma_pmt", "tangential")):
            st = dataset_style("gaia_edr3", comp, size=5.5)
            a1.plot(pc(g["r_median"]), np.asarray(g[col]) * k, ls=":", lw=1.0, **st,
                    label="Gaia EDR3 (our measurement), %s [measured, not fitted apart]" % comp)

    a1.axvspan(0.03, 200 * distance_kpc * 1e3 / 206264.806, color=style.SERIES[0],
               alpha=0.05, lw=0)
    a1.annotate("no Gaia EDR3 star passes\nthe quality flag here", (0.55, 24.0), fontsize=8.5,
                color=style.SERIES[0], ha="center", va="top")
    a1.axvline(R_JACOBI_PERI_PC, color=style.INK_SECONDARY, lw=1.2, ls=":")
    a1.annotate("$r_J$ pericentre", (R_JACOBI_PERI_PC * 1.05, 4.3), rotation=90, fontsize=8,
                color=style.INK_SECONDARY)
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_ylim(3.2, 26); a1.set_xlim(0.03, 300)
    a1.set_yticks([4, 5, 6, 8, 10, 15, 20]); a1.set_yticklabels(["4", "5", "6", "8", "10", "15", "20"])
    a1.set_ylabel("velocity dispersion  [km/s]")
    a1.legend(fontsize=7.6, loc="lower left", framealpha=0.93)
    a1.set_title("$\\omega$ Cen: every dataset in the likelihood", fontsize=12)
    add_arcsec_axis(a1, distance_kpc)

    # --- anisotropy ----------------------------------------------------------------
    hr, ht = load_profile("hst_pm_radial"), load_profile("hst_pm_tangential")
    rh = np.asarray(ht.value) / np.asarray(hr.value)
    eh = rh * np.hypot(np.asarray(ht.err_lo) / np.asarray(ht.value),
                       np.asarray(hr.err_lo) / np.asarray(hr.value))
    st = dataset_style("hst", "combined", size=4.0)
    a3.errorbar(pc(hr.r), rh, yerr=eh, linestyle="none", lw=0.9, label="HST", **st)
    st = dataset_style("gaia_edr3", "combined", size=6.5)
    a3.plot(pc(g["r_median"]), np.asarray(g["sigma_pmt"]) / np.asarray(g["sigma_pmr"]),
            ls="-", lw=1.8, label="Gaia EDR3, ours", **st)
    st = dataset_style("pristine", "combined", size=7.0, fitted=False)
    a3.plot(np.asarray(pr["r_pc"])[ok],
            np.asarray(pr["sigma_pmt"])[ok] / np.asarray(pr["sigma_pmr"])[ok],
            ls="none", label="Pristine [not fitted]", **st)
    a3.axhline(1.0, color=style.INK_SECONDARY, lw=1.5)
    a3.annotate("tangentially biased", (75, 1.22), fontsize=8, color=style.INK_SECONDARY, ha="right")
    a3.annotate("radially biased", (75, 0.755), fontsize=8, color=style.INK_SECONDARY, ha="right")
    a3.axvline(R_JACOBI_PERI_PC, color=style.INK_SECONDARY, lw=1.2, ls=":")
    a3.set_xscale("log"); a3.set_ylim(0.72, 1.28)
    a3.set_ylabel(r"$\sigma_T / \sigma_R$")
    a3.legend(fontsize=8, loc="upper left", ncol=3, framealpha=0.93)

    # --- coverage ------------------------------------------------------------------
    for i, (label, lo, hi, n, colour) in enumerate(spans):
        y = len(spans) - i
        a2.plot([lo, hi], [y, y], lw=7, color=colour, solid_capstyle="round", alpha=0.9)
        a2.annotate("%d points" % n, (hi, y), textcoords="offset points", xytext=(9, -3),
                    fontsize=8.5, color=colour, va="center")
        a2.annotate(label, (lo, y), textcoords="offset points", xytext=(-9, -3), fontsize=8.5,
                    color=colour, ha="right", va="center")
    a2.axvline(R_JACOBI_PERI_PC, color=style.INK_SECONDARY, lw=1.2, ls=":")
    a2.set_ylim(0.3, len(spans) + 0.7); a2.set_yticks([])
    a2.set_xscale("log"); a2.set_xlabel("r  [pc]"); a2.set_xlim(0.03, 300)
    a2.set_title("radial coverage and weight: %d points in total" % sum(s[3] for s in spans),
                 fontsize=10)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path


def hst_gaia_overlap_profile(edges_arcsec=(300.0, 380.0, 460.0)) -> Table:
    """HST and Gaia EDR3 on identical annuli, corrected to a common effective radius.

    Inside one annulus the two samples do not sit at the same radius: HST's coverage falls
    outwards while Gaia's quality flag passes more stars outwards, so their median radii
    differ by 30-40 arcsec. Both are moved to the annulus's geometric midpoint using the
    local logarithmic slope of the profile before the ratio is formed.
    """
    from ..kinematics.hst_profile import hst_profile
    from ..kinematics.outer_gaia import load_edr3_profile
    from ..kinematics.vb2021_replication import published_profile
    h = hst_profile(edges_arcsec=edges_arcsec)
    h_raw = hst_profile(edges_arcsec=edges_arcsec, correct_unflagged=False)
    g = load_edr3_profile()
    rp, sp = published_profile()
    keep = rp > 0
    slope_grid = np.gradient(np.log(sp[keep]), np.log(rp[keep]))
    rows = []
    for hr in h:
        j = int(np.argmin(np.abs(np.asarray(g["r_lower"]) - hr["r_lower"])))
        gr = g[j]
        mid = float(np.sqrt(hr["r_lower"] * hr["r_upper"]))
        sl = float(np.interp(mid, rp[keep], slope_grid))
        hc = float(hr["sigma_pm"]) * (mid / float(hr["r_median"])) ** sl
        gc = float(gr["sigma_pm"]) * (mid / float(gr["r_median"])) ** sl
        ratio = gc / hc
        err = ratio * float(np.hypot(gr["sigma_pm_err"] / gr["sigma_pm"],
                                     hr["sigma_pm_err"] / hr["sigma_pm"]))
        raw_h = float(h_raw["sigma_pm"][list(h["r_lower"]).index(hr["r_lower"])])
        rows.append((hr["r_lower"], hr["r_upper"], mid, hr["r_median"], gr["r_median"], sl,
                     hr["sigma_pm"], hr["sigma_pm_err"], gr["sigma_pm"], gr["sigma_pm_err"],
                     hc, gc, float(gr["sigma_pm"] / raw_h), ratio, err,
                     int(hr["n_stars"]), int(gr["n_stars"]), raw_h, float(hr["f_unflagged"])))
    return Table(rows=rows, names=("r_lower", "r_upper", "r_mid", "r_hst", "r_gaia", "slope",
                                   "sigma_hst", "sigma_hst_err", "sigma_gaia", "sigma_gaia_err",
                                   "sigma_hst_mid", "sigma_gaia_mid", "ratio_raw", "ratio",
                                   "ratio_err", "n_hst", "n_gaia", "sigma_hst_raw",
                                   "f_unflagged"))


def hst_gaia_overlap_tests(distance_kpc: float = 5.43) -> Table:
    """Three routes to the same comparison, differing in how HST's unflagged stars are used.

    1. ``300-340", both flagged`` -- HST's own quality flag and Gaia's, nothing corrected.
       The two samples land at 311 and 318 arcsec, so the radius correction is negligible.
    2. ``300-380", both flagged`` -- triples the Gaia sample; HST still flagged-only, so its
       effective radius stays at 311 arcsec and a 45-arcsec slope correction does part of
       the work.
    3. ``300-460", HST corrected`` -- uses HST's unflagged stars with the
       ``1 + f_unflagged (k - 1)`` correction, the only route that reaches 460 arcsec.
    """
    from ..kinematics.hst_profile import hst_profile
    from ..kinematics.outer_gaia import load_edr3_profile
    from ..kinematics.outer_profile import dispersion_2d, load_members
    from ..kinematics.vb2021_replication import published_profile
    from ..selection.field_template import field_density_2d
    rp, sp = published_profile(); keep = rp > 0
    grid = np.gradient(np.log(sp[keep]), np.log(rp[keep]))
    slope = lambda r: float(np.interp(r, rp[keep], grid))

    cat = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    s = load_members(exact=True, distance_kpc=distance_kpc)
    qf = np.asarray(cat["quality_flag"], int)
    err = 0.5 * (s.err_r + s.err_t)
    usable = ((qf & 2) > 0) & (err < 0.4 * np.interp(s.r_arcsec, rp, sp))
    dens = field_density_2d()

    def gaia(lo, hi):
        m = usable & (s.r_arcsec >= lo) & (s.r_arcsec < hi)
        o = dispersion_2d(s, m, dens, depth_var=0.0, field_at=s.absolute_pm)
        return (int(m.sum()), float(np.median(s.r_arcsec[m])),
                float(np.sqrt(0.5 * (o["sigma_r"] ** 2 + o["sigma_t"] ** 2))),
                float(0.5 * np.hypot(o["sigma_r_err"], o["sigma_t_err"])))

    rows = []
    for label, lo, hi, flag_only in (("300-340\"\nboth flagged", 300., 340., True),
                                     ("300-380\"\nboth flagged", 300., 380., True),
                                     ("300-460\"\nHST corrected", 300., 460., False)):
        h = hst_profile(edges_arcsec=(lo, hi), require_flag=flag_only,
                        correct_unflagged=not flag_only, min_stars=50)[0]
        ng, rg, sg, eg = gaia(lo, hi)
        mid = float(np.sqrt(lo * hi)); sl = slope(mid)
        hc = float(h["sigma_pm"]) * (mid / float(h["r_median"])) ** sl
        gc = sg * (mid / rg) ** sl
        ratio = gc / hc
        rerr = ratio * float(np.hypot(eg / sg, h["sigma_pm_err"] / h["sigma_pm"]))
        rows.append((label, lo, hi, int(h["n_stars"]), float(h["r_median"]),
                     float(h["sigma_pm"]), float(h["sigma_pm_err"]), ng, rg, sg, eg,
                     ratio, rerr, bool(flag_only)))
    return Table(rows=rows, names=("label", "r_lower", "r_upper", "n_hst", "r_hst", "sigma_hst",
                                   "sigma_hst_err", "n_gaia", "r_gaia", "sigma_gaia",
                                   "sigma_gaia_err", "ratio", "ratio_err", "flagged_only"))


def plot_hst_gaia_overlap(path: Path | str = "plots/hst_gaia_overlap.png",
                          distance_kpc: float = 5.43) -> Path:
    """The HST/Gaia overlap, built from each survey's own quality selection.

    Left: HST measured from flagged stars alone, which reaches 340 arcsec, against Gaia,
    which starts at 300. The extension to 466 arcsec using HST's rejected stars is drawn
    faintly as a cross-check, not as a measurement. Right: the comparison by three routes
    that treat those rejected stars differently.
    """
    from ..kinematics.hst_profile import hst_profile
    from ..kinematics.likelihood import load_profile
    from ..kinematics.outer_gaia import load_edr3_profile
    from .style import dataset_label, dataset_style
    style.apply()
    flagged = hst_profile(edges_arcsec=tuple(np.concatenate(
        [np.geomspace(40.0, 150.0, 5), [200., 250., 300., 340.]])))
    extended = hst_profile(edges_arcsec=(340., 380., 420., 466.), correct_unflagged=True)
    g = load_edr3_profile()
    pub = load_profile("hst_pm_combined")
    tests = hst_gaia_overlap_tests(distance_kpc)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.4, 5.7),
                                 gridspec_kw={"width_ratios": [1.55, 1]})
    a1.axvspan(300, 340, color=style.SERIES[2], alpha=0.18, lw=0)
    a1.annotate("the overlap\n300-340\"", (319, 0.745), fontsize=9.5, ha="center",
                va="top", color=style.SERIES[2])
    a1.errorbar(pub.r, np.asarray(pub.value), yerr=np.asarray(pub.err_lo), linestyle="none",
                lw=0.9, label="HST, published oMEGACat profile (ends at 300\")",
                **dataset_style("hst", size=4.0, fitted=False))
    a1.errorbar(flagged["r_median"], flagged["sigma_pm"], yerr=flagged["sigma_pm_err"],
                linestyle="-", lw=2, capsize=3, zorder=6,
                label="HST, ours, flagged stars only (to 340\")",
                **dataset_style("hst", size=8.5))
    a1.errorbar(extended["r_median"], extended["sigma_pm"], yerr=extended["sigma_pm_err"],
                linestyle=":", lw=1.4, capsize=3,
                label="HST, rejected stars with the 8 % correction (cross-check only)",
                **dataset_style("hst", size=7, fitted=False))
    a1.errorbar(g["r_median"][:4], g["sigma_pm"][:4], yerr=g["sigma_pm_err"][:4],
                linestyle="-", lw=2, capsize=3, zorder=6, label=dataset_label("gaia_edr3"),
                **dataset_style("gaia_edr3", size=8.5))
    a1.axvline(340, color=style.INK_SECONDARY, lw=1, ls=":")
    a1.annotate("HST quality flag ends", (347, 0.70), rotation=90, fontsize=7.5,
                color=style.INK_SECONDARY, va="top")
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_xlim(140, 700); a1.set_ylim(0.30, 0.78)
    a1.set_yticks([0.3, 0.4, 0.5, 0.6, 0.7]); a1.set_yticklabels(["0.3", "0.4", "0.5", "0.6", "0.7"])
    a1.set_xlabel("R  [arcsec]"); a1.set_ylabel("1-D PM dispersion  [mas/yr]")
    a1.legend(fontsize=8.0, loc="lower left"); add_pc_axis(a1, distance_kpc)
    a1.set_title("Each survey's own quality selection", fontsize=10.5)

    x = np.arange(len(tests))
    a2.axhline(1.0, color=style.INK_SECONDARY, lw=2)
    for i, row in enumerate(tests):
        st = dataset_style("gaia_edr3", size=12, fitted=bool(row["flagged_only"]))
        a2.errorbar([i], [row["ratio"]], yerr=[row["ratio_err"]], linestyle="none", lw=1.8,
                    capsize=5, zorder=6, **st)
        a2.annotate("%.3f $\\pm$ %.3f" % (row["ratio"], row["ratio_err"]), (i, row["ratio"]),
                    textcoords="offset points", xytext=(0, 16), ha="center", fontsize=8.5,
                    color=style.SERIES[0])
        a2.annotate("%d HST\n%d Gaia\nr = %.0f\" / %.0f\"" % (row["n_hst"], row["n_gaia"],
                                                                row["r_hst"], row["r_gaia"]),
                    (i, 0.775), ha="center", fontsize=7.8, color=style.INK_SECONDARY)
    a2.set_xticks(x); a2.set_xticklabels(list(tests["label"]), fontsize=8.5)
    a2.set_xlim(-0.55, len(tests) - 0.45); a2.set_ylim(0.74, 1.24)
    a2.set_ylabel("Gaia EDR3 / HST")
    a2.set_title("Three routes, three uses of the rejected stars", fontsize=10.5)
    a2.annotate("filled: no correction on either side", (0.5, 0.955), xycoords="axes fraction",
                ha="center", fontsize=8, color=style.INK_SECONDARY)
    fig.tight_layout()
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path

