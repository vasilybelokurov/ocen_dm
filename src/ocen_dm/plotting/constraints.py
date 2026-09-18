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

__all__ = ["plot_constraint_map", "plot_outer_tracer_audit", "plot_contamination_model", "plot_annulus_fits",
           "plot_method_comparison",
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
