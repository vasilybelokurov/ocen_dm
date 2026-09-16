"""Overview figures of every dataset ingested in Milestone 1.

Each function reads a processed product from ``data/processed`` and writes one
PNG to ``plots/``. Nothing here changes any data. Coordinates are plotted as
delivered; the cluster centre used for annotation is the SIMBAD/Harris (2010)
position of NGC 5139, RA 201.697, Dec -47.480 deg.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

from ..paths import processed_dir, project_root, raw_dir
from . import style

__all__ = ["plot_all", "PLOTS"]

#: NGC 5139 centre, SIMBAD (J2000): 13h26m47.28s -47d28m46.1s
OCEN_RA, OCEN_DEC = 201.697, -47.4795
#: oMEGACat VI kinematic distance, used to label a secondary axis in pc and to
#: convert proper motions to km/s for display: v = 4.74047 * mu * D
OCEN_DISTANCE_KPC = 5.43
ARCSEC_TO_PC = OCEN_DISTANCE_KPC * 1e3 / 206265.0
MASYR_TO_KMS = 4.74047 * OCEN_DISTANCE_KPC
#: systemic heliocentric velocity, taken as the median of the high-quality
#: oMEGACat sample (232.6 km/s), used only to centre displays
OCEN_VSYS = 232.6
MEMBER_THRESHOLD = 0.5


def plots_dir() -> Path:
    """Return ``<root>/plots``, creating it."""
    path = project_root() / "plots"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _load(name: str, subdir: str) -> Table:
    return Table.read(processed_dir() / subdir / f"{name}.ecsv")


def _save(fig, name: str) -> Path:
    path = plots_dir() / f"{name}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def _tangent(ra, dec):
    """Offsets from the cluster centre in degrees (cos-dec corrected)."""
    dx = (np.asarray(ra) - OCEN_RA) * np.cos(np.radians(OCEN_DEC))
    dy = np.asarray(dec) - OCEN_DEC
    return dx, dy


# ---------------------------------------------------------------- sky ---------
def plot_sky_overview() -> Path:
    """All tail datasets on the sky: a wide field, and a zoom on the periphery.

    Fimbulthul lies ~20 degrees north of the cluster, so one panel cannot show
    both it and the 5-degree periphery survey legibly.
    """
    per = _load("kuzma2025_periphery", "tails")
    spec = _load("kuzma2026_spectroscopy", "tails")
    fim = _load("fimbulthul_members", "tails")

    members = np.asarray(per["membership_prob"]) > MEMBER_THRESHOLD
    spec_member = np.asarray(spec["member_flag"]).astype(str) == "True"
    dx, dy = _tangent(per["ra"], per["dec"])
    sx, sy = _tangent(spec["ra"], spec["dec"])
    fx, fy = _tangent(fim["ra"], fim["dec"])

    fig, (ax_wide, ax_zoom) = plt.subplots(1, 2, figsize=(13, 6.2),
                                           gridspec_kw={"width_ratios": [1, 1.15]})

    # -- wide field: where everything is relative to the cluster
    footprint = plt.Circle((0, 0), 5.0, fc=style.COLOR_FIELD, ec="none", alpha=0.45)
    ax_wide.add_patch(footprint)
    ax_wide.scatter(dx[members], dy[members], s=3, c=style.COLOR_KUZMA2025, linewidths=0,
                    label="Kuzma & Ishigaki 2025 members")
    ax_wide.scatter(fx, fy, s=10, marker="^", c=style.COLOR_FIMBULTHUL, linewidths=0,
                    label=f"Fimbulthul candidates, Ibata+ 2019 (n={len(fim)})")
    ax_wide.plot(0, 0, marker="+", ms=12, mew=1.5, color=style.INK, ls="none", label="NGC 5139 centre")
    ax_wide.text(0, -6.3, "Gaia + Pristine footprint (5°)", ha="center", va="top",
                 color=style.INK_SECONDARY, fontsize=9)
    ax_wide.set_xlim(12, -12)
    ax_wide.set_ylim(-9, 28)
    ax_wide.set_aspect("equal")
    ax_wide.set_xlabel("ΔRA cos(Dec)  [deg]  (east to the left)")
    ax_wide.set_ylabel("ΔDec  [deg]")
    ax_wide.set_title("Wide field: cluster, periphery survey and Fimbulthul")
    ax_wide.legend(loc="upper right", markerscale=1.8)

    # -- zoom: the periphery survey and the spectroscopic fields
    ax_zoom.scatter(dx[~members], dy[~members], s=1.2, c=style.COLOR_FIELD, alpha=0.3,
                    linewidths=0, rasterized=True,
                    label=f"field, P_mem ≤ {MEMBER_THRESHOLD} (n={(~members).sum()})")
    ax_zoom.scatter(dx[members], dy[members], s=7, c=style.COLOR_KUZMA2025, linewidths=0,
                    label=f"members, P_mem > {MEMBER_THRESHOLD} (n={members.sum()})")
    ax_zoom.scatter(sx[~spec_member], sy[~spec_member], s=12, marker="x", c=style.INK_SECONDARY,
                    linewidths=0.6, label=f"Kuzma+ 2026 targets, non-members (n={(~spec_member).sum()})")
    ax_zoom.scatter(sx[spec_member], sy[spec_member], s=26, marker="o", facecolors="none",
                    edgecolors=style.COLOR_KUZMA2026, linewidths=1.0,
                    label=f"Kuzma+ 2026 members (n={spec_member.sum()})")
    ax_zoom.plot(0, 0, marker="+", ms=14, mew=1.5, color=style.INK, ls="none")
    ax_zoom.set_xlim(5.6, -5.6)
    ax_zoom.set_ylim(-5.6, 5.6)
    ax_zoom.set_aspect("equal")
    ax_zoom.set_xlabel("ΔRA cos(Dec)  [deg]  (east to the left)")
    ax_zoom.set_title("Periphery: Gaia + Pristine members and FLAMES fields")
    ax_zoom.legend(loc="upper left", markerscale=1.6, fontsize=8.5)

    fig.suptitle("Omega Centauri: every ingested tail dataset on the sky", y=0.99)
    fig.tight_layout()
    return _save(fig, "sky_overview")


# ------------------------------------------------------ Kuzma & Ishigaki 2025 -
def plot_kuzma2025_proper_motions() -> Path:
    """Proper-motion diagram: field density versus members."""
    per = _load("kuzma2025_periphery", "tails")
    pmra, pmdec = np.asarray(per["pmra"]), np.asarray(per["pmdec"])
    members = np.asarray(per["membership_prob"]) > MEMBER_THRESHOLD

    fig, ax = plt.subplots(figsize=(7, 6.5))
    hb = ax.hexbin(pmra[~members], pmdec[~members], gridsize=80, extent=(-15, 5, -15, 5),
                   cmap=style.SEQUENTIAL, mincnt=1, bins="log", linewidths=0.1)
    ax.scatter(pmra[members], pmdec[members], s=5, c=style.COLOR_KUZMA2026, linewidths=0,
               label=f"members, P_mem > {MEMBER_THRESHOLD} (n={members.sum()})")
    cb = fig.colorbar(hb, ax=ax, pad=0.02)
    cb.set_label("field stars per cell (log)", color=style.INK_SECONDARY)
    cb.outline.set_visible(False)
    ax.set_xlabel("μ_α cos δ  [mas / yr]")
    ax.set_ylabel("μ_δ  [mas / yr]")
    ax.set_title("Kuzma & Ishigaki 2025: proper motions — field density and members")
    ax.legend(loc="lower left", markerscale=2.0)
    return _save(fig, "kuzma2025_proper_motions")


def plot_kuzma2025_cmd_and_metallicity() -> Path:
    """Dereddened CMD and photometric [Fe/H] for members against the field."""
    per = _load("kuzma2025_periphery", "tails")
    members = np.asarray(per["membership_prob"]) > MEMBER_THRESHOLD
    g0 = np.asarray(per["g0_mag"])
    color = np.asarray(per["bp0_mag"]) - np.asarray(per["rp0_mag"])
    feh = np.asarray(per["feh"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5))
    ax1.scatter(color[~members], g0[~members], s=1.5, c=style.COLOR_FIELD, alpha=0.4,
                linewidths=0, rasterized=True, label="field")
    ax1.scatter(color[members], g0[members], s=6, c=style.COLOR_KUZMA2025, linewidths=0,
                label=f"members (n={members.sum()})")
    ax1.invert_yaxis()
    ax1.set_xlabel("(BP − RP)₀  [mag]")
    ax1.set_ylabel("G₀  [mag]")
    ax1.set_xlim(0.3, 2.2)
    ax1.set_title("Dereddened colour–magnitude diagram")
    ax1.legend(loc="upper right", markerscale=2.0)

    bins = np.linspace(-4, 0.5, 46)
    ax2.hist(feh[~members & np.isfinite(feh)], bins=bins, density=True, histtype="stepfilled",
             color=style.COLOR_FIELD, alpha=0.7, label="field")
    ax2.hist(feh[members & np.isfinite(feh)], bins=bins, density=True, histtype="step",
             color=style.COLOR_KUZMA2025, lw=1.8, label="members")
    ax2.set_xlabel("[Fe/H] from synthetic CaHK  [dex]")
    ax2.set_ylabel("normalised density")
    ax2.set_title("Photometric metallicity")
    ax2.legend(loc="upper left")
    fig.suptitle("Kuzma & Ishigaki 2025 periphery catalogue", y=1.02)
    return _save(fig, "kuzma2025_cmd_metallicity")


# ------------------------------------------------------------ Kuzma+ 2026 -----
def plot_kuzma2026_spectroscopy() -> Path:
    """Velocity–metallicity plane and velocity histogram of the spectroscopic targets."""
    spec = _load("kuzma2026_spectroscopy", "tails")
    v = np.asarray(spec["vlos"])
    ev = np.asarray(spec["vlos_error"])
    feh = np.asarray(spec["feh"])
    efeh = np.asarray(spec["feh_error"])
    member = np.asarray(spec["member_flag"]).astype(str) == "True"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5.5), gridspec_kw={"width_ratios": [1.5, 1]})
    ax1.errorbar(v[~member], feh[~member], xerr=ev[~member], yerr=efeh[~member], fmt="o",
                 ms=3, mfc=style.COLOR_FIELD, mec="none", ecolor=style.COLOR_FIELD,
                 elinewidth=style.ERR_PT, alpha=0.7, label=f"non-members (n={(~member).sum()})")
    ax1.errorbar(v[member], feh[member], xerr=ev[member], yerr=efeh[member], fmt="o",
                 ms=4, mfc=style.COLOR_KUZMA2026, mec="none", ecolor=style.COLOR_KUZMA2026,
                 elinewidth=style.ERR_PT, label=f"members (n={member.sum()})")
    vm = v[member]
    ax1.axvline(np.mean(vm), color=style.INK_SECONDARY, lw=0.8, ls="--")
    ax1.text(np.mean(vm) + 4, -3.7, f"⟨v⟩ = {np.mean(vm):.1f} km/s\nσ = {np.std(vm):.1f} km/s",
             color=style.INK_SECONDARY, fontsize=9, va="bottom")
    ax1.set_xlabel("heliocentric velocity  [km / s]")
    ax1.set_ylabel("[Fe/H]  [dex]")
    ax1.set_title("Velocity–metallicity plane")
    ax1.legend(loc="lower left")

    bins = np.arange(-150, 400, 10)
    ax2.hist(v[~member], bins=bins, color=style.COLOR_FIELD, histtype="stepfilled", alpha=0.7,
             label="non-members")
    ax2.hist(v[member], bins=bins, color=style.COLOR_KUZMA2026, histtype="stepfilled", alpha=0.9,
             label="members")
    ax2.set_xlabel("heliocentric velocity  [km / s]")
    ax2.set_ylabel("stars per 10 km/s")
    ax2.set_title("Velocity distribution")
    ax2.legend(loc="upper left")
    fig.suptitle("Kuzma et al. 2026: FLAMES spectroscopy of the periphery and tails", y=1.02)
    return _save(fig, "kuzma2026_spectroscopy")


# -------------------------------------------------------------- Fimbulthul ----
def plot_fimbulthul() -> Path:
    """Fimbulthul candidates: positions with proper-motion vectors, and their CMD."""
    fim = _load("fimbulthul_members", "tails")
    ra, dec = np.asarray(fim["ra"]), np.asarray(fim["dec"])
    pmra, pmdec = np.asarray(fim["pmra"]), np.asarray(fim["pmdec"])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), gridspec_kw={"width_ratios": [1.6, 1]})
    # 1 mas/yr drawn as 0.25 deg; the reference arrow states the scale.
    q = ax1.quiver(ra, dec, pmra, pmdec, color=style.COLOR_FIMBULTHUL, width=0.002,
                   scale=4.0, scale_units="xy", angles="xy", alpha=0.8)
    ax1.quiverkey(q, 0.12, 0.08, 4.0, "4 mas/yr", labelpos="E", coordinates="axes",
                  color=style.INK_SECONDARY, labelcolor=style.INK_SECONDARY, fontproperties={"size": 9})
    ax1.scatter(ra, dec, s=8, c=style.COLOR_FIMBULTHUL, linewidths=0, label=f"candidates (n={len(fim)})")
    ax1.plot(OCEN_RA, OCEN_DEC, marker="+", ms=14, mew=1.5, color=style.INK, ls="none",
             label="NGC 5139")
    ax1.set_xlim(217.5, 194.5)
    ax1.set_ylim(-49.5, -19.5)
    ax1.set_aspect("equal")
    ax1.set_xlabel("RA  [deg]")
    ax1.set_ylabel("Dec  [deg]")
    ax1.set_title("Sky positions with Gaia proper-motion vectors")
    ax1.legend(loc="center left")

    ax2.scatter(np.asarray(fim["bp_rp0"]), np.asarray(fim["g0_mag"]), s=10,
                c=style.COLOR_FIMBULTHUL, linewidths=0)
    ax2.invert_yaxis()
    ax2.set_xlabel("(BP − RP)₀  [mag]")
    ax2.set_ylabel("G₀  [mag]")
    ax2.set_title("Dereddened CMD")
    fig.suptitle("Fimbulthul stream candidates (Ibata et al. 2019, via VizieR J/other/NatAs/3.667)", y=1.02)
    return _save(fig, "fimbulthul")


# ------------------------------------------------------------- oMEGACat VI ----
def _asym(table: Table, name: str) -> np.ndarray:
    return np.vstack([np.asarray(table[f"{name}_err_lo"]), np.asarray(table[f"{name}_err_hi"])])


def _add_pc_axis(ax) -> None:
    secax = ax.secondary_xaxis("top", functions=(lambda a: a * ARCSEC_TO_PC, lambda p: p / ARCSEC_TO_PC))
    secax.set_xlabel(f"r  [pc]  (D = {OCEN_DISTANCE_KPC} kpc)", color=style.INK_SECONDARY)
    secax.tick_params(colors=style.INK_SECONDARY)


def plot_omegacat_profiles() -> Path:
    """Proper-motion and line-of-sight kinematic profiles with asymmetric errors."""
    rad = _load("omegacat_vi_pm_radial", "kinematics")
    tan = _load("omegacat_vi_pm_tangential", "kinematics")
    com = _load("omegacat_vi_pm_combined", "kinematics")
    los = _load("omegacat_vi_los_dispersion", "kinematics")
    rot = _load("omegacat_vi_los_rotation", "kinematics")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    (ax_pm, ax_an), (ax_los, ax_rot) = axes
    r_pm = np.asarray(rad["r_median"])
    r_los = np.asarray(los["r_median"])

    for table, key, color, label in (
        (com, "sigma_pmc", style.SERIES[0], "combined σ_PM"),
        (rad, "sigma_pmr", style.SERIES[1], "radial σ_R"),
        (tan, "sigma_pmt", style.SERIES[2], "tangential σ_T"),
    ):
        ax_pm.errorbar(r_pm, table[key], yerr=_asym(table, key), fmt="o-", ms=style.MARKER_PT,
                       color=color, ecolor=color, elinewidth=style.ERR_PT, capsize=2, label=label)
    ax_pm.set_ylabel("proper-motion dispersion  [mas / yr]")
    ax_pm.set_title("PM dispersion, adaptive log bins (n = 40)")
    ax_pm.legend(loc="upper right")
    _add_pc_axis(ax_pm)

    ratio = np.asarray(rad["sigma_pmr"]) / np.asarray(tan["sigma_pmt"])
    err = ratio * np.sqrt(
        (0.5 * (rad["sigma_pmr_err_lo"] + rad["sigma_pmr_err_hi"]) / rad["sigma_pmr"]) ** 2
        + (0.5 * (tan["sigma_pmt_err_lo"] + tan["sigma_pmt_err_hi"]) / tan["sigma_pmt"]) ** 2
    )
    ax_an.errorbar(r_pm, ratio, yerr=err, fmt="o", ms=style.MARKER_PT, color=style.SERIES[0],
                   ecolor=style.SERIES[0], elinewidth=style.ERR_PT, capsize=2)
    ax_an.axhline(1.0, color=style.INK_SECONDARY, lw=0.8, ls="--")
    ax_an.text(r_pm[1], 1.015, "isotropic", color=style.INK_SECONDARY, fontsize=9)
    ax_an.set_ylabel("σ_R / σ_T")
    ax_an.set_title("Anisotropy diagnostic (symmetrised errors, propagated)")
    _add_pc_axis(ax_an)

    ax_los.errorbar(r_los, los["sigma_los"], yerr=_asym(los, "sigma_los"), fmt="o-",
                    ms=style.MARKER_PT, color=style.SERIES[0], ecolor=style.SERIES[0],
                    elinewidth=style.ERR_PT, capsize=2)
    ax_los.set_ylabel("σ_LOS  [km / s]")
    ax_los.set_xlabel("r  [arcsec]")
    ax_los.set_title("Line-of-sight dispersion, MUSE (n = 29)")

    ax_rot.errorbar(r_los, rot["v_rot"], yerr=_asym(rot, "v_rot"), fmt="o-", ms=style.MARKER_PT,
                    color=style.SERIES[0], ecolor=style.SERIES[0], elinewidth=style.ERR_PT, capsize=2)
    ax_rot.set_ylabel("rotation amplitude v_rot  [km / s]")
    ax_rot.set_xlabel("r  [arcsec]")
    ax_rot.set_title("Line-of-sight rotation")

    for ax in axes.ravel():
        ax.set_xscale("log")
    fig.suptitle("oMEGACat VI kinematic profiles (Zenodo 10.5281/zenodo.14978551)", y=0.995)
    fig.tight_layout()
    return _save(fig, "omegacat_vi_profiles")


def plot_omegacat_rotation_axis() -> Path:
    """Position angle of the rotation axis versus radius."""
    rot = _load("omegacat_vi_los_rotation", "kinematics")
    r = np.asarray(rot["r_median"])
    theta = np.asarray(rot["theta_0"])
    # The angle is periodic in 360 deg; the file stores it in (-180, 180], which
    # splits one direction across the two ends of the axis. Wrap into a window
    # centred on the well-determined outer value so the trend reads continuously.
    centre = np.median(theta[-10:])
    wrapped = (theta - centre + 180.0) % 360.0 + centre - 180.0
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.errorbar(r, wrapped, yerr=_asym(rot, "theta_0"), fmt="o", ms=style.MARKER_PT,
                color=style.SERIES[0], ecolor=style.SERIES[0], elinewidth=style.ERR_PT, capsize=2)
    ax.axhline(centre, color=style.INK_SECONDARY, lw=0.8, ls="--")
    ax.text(r[0], centre + 6, f"outer median {centre:.0f}°", color=style.INK_SECONDARY, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("r  [arcsec]")
    ax.set_ylabel(f"θ₀  [deg], wrapped to ({centre - 180:.0f}°, {centre + 180:.0f}°]")
    ax.set_title("oMEGACat VI: position angle of the rotation axis")
    _add_pc_axis(ax)
    return _save(fig, "omegacat_vi_rotation_axis")


# --------------------------------------------------- velocities along the tails ---
def plot_kuzma2026_vlos_along_tails() -> Path:
    """Line-of-sight velocity against signed distance along the tidal tails.

    Sign comes from the target names: fields I/M/O (inner, middle, outer) and
    L/T (leading, trailing). Leading is drawn positive.
    """
    spec = _load("kuzma2026_spectroscopy", "tails")
    ra, dec = np.asarray(spec["ra"]), np.asarray(spec["dec"])
    v, ev = np.asarray(spec["vlos"]), np.asarray(spec["vlos_error"])
    member = np.asarray(spec["member_flag"]).astype(str) == "True"
    ids = np.asarray(spec["star_id"]).astype(str)
    field = np.array([i[:2] for i in ids])              # OT, MT, IT, IL, ML, OL
    leading = np.array([f[1] == "L" for f in field])
    dx, dy = _tangent(ra, dec)
    sep = np.hypot(dx, dy) * np.where(leading, 1.0, -1.0)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.axhline(OCEN_VSYS, color=style.INK_SECONDARY, lw=0.8, ls="--")
    ax.text(3.35, OCEN_VSYS + 6, f"systemic {OCEN_VSYS} km/s", color=style.INK_SECONDARY,
            fontsize=9, ha="right")
    ax.errorbar(sep[~member], v[~member], yerr=ev[~member], fmt="o", ms=3, mfc=style.COLOR_FIELD,
                mec="none", ecolor=style.COLOR_FIELD, elinewidth=style.ERR_PT, alpha=0.7,
                label=f"non-members (n={(~member).sum()})")
    ax.errorbar(sep[member], v[member], yerr=ev[member], fmt="o", ms=4.5, mfc=style.COLOR_KUZMA2026,
                mec="none", ecolor=style.COLOR_KUZMA2026, elinewidth=style.ERR_PT,
                label=f"members (n={member.sum()})")
    # per-field mean of the members, with the error of the mean, labelled directly
    for name in ("OT", "MT", "IT", "IL", "ML", "OL"):
        pick = member & (field == name)
        if pick.sum() < 2:
            ax.text(np.mean(sep[field == name]), 290, f"{name}\n({pick.sum()} member)",
                    ha="center", va="top", color=style.INK_SECONDARY, fontsize=8.5)
            continue
        m, e = np.mean(v[pick]), np.std(v[pick]) / np.sqrt(pick.sum())
        x = np.mean(sep[pick])
        ax.errorbar(x, m, yerr=e, fmt="s", ms=8, mfc=style.SURFACE, mec=style.INK, mew=1.2,
                    ecolor=style.INK, elinewidth=1.4, capsize=4, zorder=5)
        ax.text(x, 290, f"{name}\n{m:.1f} ± {e:.1f}\n(n={pick.sum()})", ha="center", va="top",
                color=style.INK, fontsize=8.5)
    ax.set_ylim(-160, 300)
    ax.set_title("Kuzma et al. 2026: line-of-sight velocities along the tidal tails", pad=14)
    ax.set_xlim(-3.5, 3.5)
    ax.set_xlabel("signed angular distance from the centre  [deg]   (trailing ← 0 → leading)")
    ax.set_ylabel("heliocentric velocity  [km / s]")
    ax.legend(loc="lower left")
    return _save(fig, "kuzma2026_vlos_along_tails")


# ------------------------------------------------------ oMEGACat 3D velocities ----
def _omegacat_hq_3d() -> dict[str, np.ndarray]:
    """High-quality stars with both proper motions and v_los, in tangent-plane arcsec."""
    cat = Table.read(raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits")
    hq = np.asarray(cat["selection_hq_pm_and_los"]) == 1
    ra, dec = np.asarray(cat["RA"])[hq], np.asarray(cat["DEC"])[hq]
    dx, dy = _tangent(ra, dec)
    x, y = dx * 3600.0, dy * 3600.0
    pmra = np.asarray(cat["pmra_corrected"])[hq]
    pmde = np.asarray(cat["pmdec_corrected"])[hq]
    r = np.hypot(x, y)
    ux, uy = x / r, y / r
    return {
        "x": x, "y": y, "r": r,
        "v_los": np.asarray(cat["vlos"])[hq],
        "v_rad": (pmra * ux + pmde * uy) * MASYR_TO_KMS,      # plane-of-sky radial
        "v_tan": (-pmra * uy + pmde * ux) * MASYR_TO_KMS,     # plane-of-sky tangential
        "v_ra": pmra * MASYR_TO_KMS, "v_de": pmde * MASYR_TO_KMS,
    }


def _fit_kinematic_axis(x, y, v, r_min=30.0):
    """Fit ``v = v0 + A sin(phi - phi0)`` with ``phi`` the position angle from North
    through East, using stars outside ``r_min`` arcsec. Returns ``(v0, A, phi0)`` in
    km/s, km/s and degrees; ``phi0`` is the zero-velocity (rotation-axis) direction."""
    r = np.hypot(x, y)
    keep = r > r_min
    phi = np.arctan2(x[keep], y[keep])            # PA: 0 at North (+y), +90 at East (+x)
    design = np.column_stack([np.ones(keep.sum()), np.sin(phi), np.cos(phi)])
    coef, *_ = np.linalg.lstsq(design, v[keep], rcond=None)
    v0, a, b = coef                               # v = v0 + a sin(phi) + b cos(phi)
    amp = np.hypot(a, b)
    phi0 = np.degrees(np.arctan2(-b, a)) % 180.0  # sin(phi - phi0) expansion
    return float(v0), float(amp), float(phi0)


def plot_omegacat_3d_velocities() -> Path:
    """The 3D velocity data of the core: LOS rotation map with the fitted kinematic
    axis, rotation curves from PM and LOS, the three velocity components, and the
    PM-versus-LOS dispersion comparison that sets the kinematic distance."""
    d = _omegacat_hq_3d()
    n = len(d["x"])
    dv = d["v_los"] - OCEN_VSYS
    v0, amp, phi0 = _fit_kinematic_axis(d["x"], d["y"], d["v_los"])

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 11))
    (ax_los, ax_rot), (ax_hist, ax_disp) = axes
    extent = (-400, 400, -300, 300)

    hb = ax_los.hexbin(d["x"], d["y"], C=dv, vmin=-8, vmax=8, gridsize=28, extent=extent,
                       cmap=style.DIVERGING, reduce_C_function=np.mean, mincnt=8,
                       linewidths=0.2, edgecolors=style.SURFACE)
    cb = fig.colorbar(hb, ax=ax_los, pad=0.02)
    cb.set_label(f"⟨v_LOS⟩ − {OCEN_VSYS} km/s  [km / s]", color=style.INK_SECONDARY)
    cb.outline.set_visible(False)
    # the fitted zero-velocity line, i.e. the projected rotation axis
    t = np.linspace(-320, 320, 2)
    ax_los.plot(t * np.sin(np.radians(phi0)), t * np.cos(np.radians(phi0)), color=style.INK,
                lw=1.2, ls="--")
    ax_los.text(-390, -280, f"fitted rotation axis: PA = {phi0:.0f}° (N through E), r > 30″\n"
                f"amplitude {amp:.1f} km/s   |   paper θ₀ = 104.3 ± 1.4°",
                fontsize=8.5, color=style.INK, ha="right", va="bottom",
                bbox=dict(fc=style.SURFACE, ec="none", alpha=0.85))
    ax_los.set_title(f"Line-of-sight velocity field (mean per cell, n = {n:,} stars)")
    ax_los.set_xlabel("ΔRA cos(Dec)  [arcsec]  (east to the left)")
    ax_los.set_ylabel("ΔDec  [arcsec]")
    ax_los.set_xlim(400, -400)
    ax_los.set_ylim(-300, 300)
    ax_los.set_aspect("equal")
    ax_los.grid(False)

    # rotation curves: plane-of-sky from PM (these PMs are relative to the bulk
    # motion, so any signal here is differential rotation) and LOS from the paper
    edges = np.geomspace(5.0, 330.0, 14)
    mids, mean_t, err_t = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        pick = (d["r"] >= lo) & (d["r"] < hi)
        if pick.sum() < 30:
            continue
        mids.append(np.median(d["r"][pick]))
        mean_t.append(np.mean(d["v_tan"][pick]))
        err_t.append(np.std(d["v_tan"][pick]) / np.sqrt(pick.sum()))
    ax_rot.errorbar(mids, mean_t, yerr=err_t, fmt="o-", ms=style.MARKER_PT, color=style.SERIES[1],
                    ecolor=style.SERIES[1], elinewidth=style.ERR_PT, capsize=2,
                    label="⟨v_tan⟩ from proper motions (this sample, plane of sky)")
    rot = _load("omegacat_vi_los_rotation", "kinematics")
    ax_rot.errorbar(rot["r_median"], rot["v_rot"], yerr=_asym(rot, "v_rot"), fmt="o-",
                    ms=style.MARKER_PT, color=style.SERIES[0], ecolor=style.SERIES[0],
                    elinewidth=style.ERR_PT, capsize=2, label="v_rot from line-of-sight velocities (paper)")
    ax_rot.axhline(0, color=style.INK_SECONDARY, lw=0.8, ls="--")
    ax_rot.set_xscale("log")
    ax_rot.set_xlabel("r  [arcsec]")
    ax_rot.set_ylabel("rotation velocity  [km / s]")
    ax_rot.set_title("Rotation: plane of sky versus line of sight")
    ax_rot.legend(loc="upper left", fontsize=8.5)
    _add_pc_axis(ax_rot)

    bins = np.linspace(-80, 80, 81)
    for values, color, label in (
        (d["v_rad"], style.SERIES[1], f"v_R (PM, plane-of-sky radial)   σ = {np.std(d['v_rad']):.1f}"),
        (d["v_tan"], style.SERIES[2], f"v_T (PM, plane-of-sky tangential)   σ = {np.std(d['v_tan']):.1f}"),
        (dv, style.SERIES[0], f"v_LOS − systemic   σ = {np.std(dv):.1f}"),
    ):
        ax_hist.hist(values, bins=bins, histtype="step", lw=1.6, color=color, label=label)
    ax_hist.set_xlabel(f"velocity  [km / s]   (proper motions converted at D = {OCEN_DISTANCE_KPC} kpc)")
    ax_hist.set_ylabel("stars per 2 km/s")
    ax_hist.set_title("The three velocity components of the same stars")
    ax_hist.legend(loc="upper left", fontsize=8.5)

    com = _load("omegacat_vi_pm_combined", "kinematics")
    los = _load("omegacat_vi_los_dispersion", "kinematics")
    spm = np.asarray(com["sigma_pmc"]) * MASYR_TO_KMS
    epm = _asym(com, "sigma_pmc") * MASYR_TO_KMS
    ax_disp.errorbar(com["r_median"], spm, yerr=epm, fmt="o-", ms=style.MARKER_PT,
                     color=style.SERIES[1], ecolor=style.SERIES[1], elinewidth=style.ERR_PT,
                     capsize=2, label=f"σ_PM,combined × 4.74 D  (D = {OCEN_DISTANCE_KPC} kpc)")
    ax_disp.errorbar(los["r_median"], los["sigma_los"], yerr=_asym(los, "sigma_los"), fmt="o-",
                     ms=style.MARKER_PT, color=style.SERIES[0], ecolor=style.SERIES[0],
                     elinewidth=style.ERR_PT, capsize=2, label="σ_LOS")
    ax_disp.set_xscale("log")
    ax_disp.set_xlabel("r  [arcsec]")
    ax_disp.set_ylabel("velocity dispersion  [km / s]")
    ax_disp.set_title("PM and LOS dispersions agree when D is right: the kinematic distance")
    ax_disp.legend(loc="lower left")
    _add_pc_axis(ax_disp)

    fig.suptitle("oMEGACat VI: three-dimensional kinematics of the core "
                 "(HST proper motions + MUSE line-of-sight velocities)", y=0.995)
    fig.tight_layout()
    return _save(fig, "omegacat_vi_3d_velocities")

PLOTS = (
    plot_sky_overview,
    plot_kuzma2025_proper_motions,
    plot_kuzma2025_cmd_and_metallicity,
    plot_kuzma2026_spectroscopy,
    plot_kuzma2026_vlos_along_tails,
    plot_fimbulthul,
    plot_omegacat_profiles,
    plot_omegacat_rotation_axis,
    plot_omegacat_3d_velocities,
)


def plot_all(only: Iterable[str] | None = None) -> list[Path]:
    """Write every overview figure and return the paths."""
    style.apply()
    wanted = set(only) if only else None
    out = []
    for fn in PLOTS:
        if wanted and fn.__name__ not in wanted:
            continue
        out.append(fn())
    return out
