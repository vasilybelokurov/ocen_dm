#!/usr/bin/env python3
"""Detailed figures of the measurements, read through the same loaders the likelihood uses.

Figures (plots/data_<name>.png):
  kinematics  the five kinematic profiles of the baseline DF fit (ours: HST PM radial and
              tangential, MUSE line-of-sight, Gaia EDR3 PM radial and tangential) and the
              published profiles that are not fitted (oMEGACat VI PMs, Vasiliev & Baumgardt
              2021 and Baumgardt+2019 Gaia PMs), with the rotation term each fitted profile
              carries (sqrt(streaming2)), against projected radius (arcsec; pc on top);
  rotation    the rotation inputs of the streaming correction: Vasiliev & Baumgardt 2021 PM
              rotation curve (median, 16-84 and 2-98 per cent) and the MUSE LOS rotation;
  counts      the Poisson star-count products (HST F625W<19, Gaia G<17): surface number
              density with Poisson errors, stars per bin and footprint coverage;
  stream54    STREAMFINDER (Ibata+2024) stream 54 members: sky, proper motions and radial
              velocities along the stream (chord great-circle frame phi1).
Needs $OCEN_CATALOGUES for the stream figure.

Usage: python bin/plot_data_products.py [--only kinematics,rotation,counts,stream54]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from ocen_dm.cluster import OCEN_DISTANCE_KPC
from ocen_dm.kinematics.likelihood import load_profile

D_KPC = OCEN_DISTANCE_KPC                                   # distance used for the pc axis
PC_PER_ARCSEC = D_KPC*1e3/206264.806
KMS_PER_MASYR = 4.740470*D_KPC
PROC = ROOT/"data/processed/kinematics"
BASE = dict(hst="black", muse="C2", gaia="C0")


def pc_axis(ax):
    sec = ax.secondary_xaxis("top", functions=(lambda x: x*PC_PER_ARCSEC, lambda x: x/PC_PER_ARCSEC))
    sec.set_xlabel(f"R [pc] (D = {D_KPC} kpc)")


def errbar(ax, p, **kw):
    ax.errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt=kw.pop("fmt", "o"), ms=kw.pop("ms", 3.5),
                capsize=0, lw=0.8, **kw)


def kinematics():
    fig, axs = plt.subplots(2, 2, figsize=(14, 9.5), constrained_layout=True)
    for ax, comp, title in ((axs[0, 0], "radial", "Proper motion, radial"), (axs[0, 1], "tangential", "Proper motion, tangential")):
        pub = load_profile(f"hst_pm_{comp}")
        errbar(ax, pub, color="0.65", fmt="s", ms=3, label=f"oMEGACat VI published (not fitted, {pub.n} bins)")
        ours = load_profile(f"hst_pm_{comp}_ours")
        errbar(ax, ours, color=BASE["hst"], label=f"HST, ours (fitted, {ours.n} bins)")
        g = load_profile(f"gaia_edr3_ours_{comp}")
        errbar(ax, g, color=BASE["gaia"], label=f"Gaia EDR3, ours (fitted, {g.n} bins)")
        for p, c in ((ours, BASE["hst"]), (g, BASE["gaia"])):
            if p.streaming2 is not None and np.any(p.streaming2 > 0):
                ax.plot(p.r, np.sqrt(p.streaming2), color=c, ls="--", lw=1, label=f"rotation term in the fit ({p.instrument})")
        ax.set(xscale="log", ylabel=r"$\sigma$ [mas/yr]", title=title)
        sec = ax.secondary_yaxis("right", functions=(lambda y: y*KMS_PER_MASYR, lambda y: y/KMS_PER_MASYR))
        sec.set_ylabel("km/s")
    ax = axs[1, 0]
    los = load_profile("muse_los_dispersion")
    errbar(ax, los, color=BASE["muse"], label=f"MUSE LOS dispersion (fitted, {los.n} bins)")
    ax.plot(los.r, np.sqrt(los.streaming2), color=BASE["muse"], ls="--", lw=1, label=r"rotation term in the fit ($v_{\rm rot}/\sqrt{2}$)")
    ax.set(xscale="log", ylabel=r"$\sigma_{\rm LOS}$ [km/s]", title="Line of sight (MUSE, oMEGACat VI)")
    ax = axs[1, 1]
    for key, c, m, lab in (("gaia_edr3_pm", "C1", "D", "Vasiliev & Baumgardt 2021 EDR3 (not fitted)"),
                           ("gaia_dr2_pm", "C3", "^", "Baumgardt+2019 DR2 (not fitted)"),
                           ("hst_pm_combined", "0.65", "s", "oMEGACat VI combined (not fitted)")):
        p = load_profile(key)
        errbar(ax, p, color=c, fmt=m, ms=3, label=lab)
    for comp, mk in (("radial", "o"), ("tangential", "v")):
        g = load_profile(f"gaia_edr3_ours_{comp}")
        errbar(ax, g, color=BASE["gaia"], fmt=mk, label=f"Gaia EDR3 ours, {comp} (fitted)")
    ax.set(xscale="log", ylabel=r"$\sigma$ [mas/yr]", title="Published 1D proper-motion profiles vs ours")
    sec = ax.secondary_yaxis("right", functions=(lambda y: y*KMS_PER_MASYR, lambda y: y/KMS_PER_MASYR))
    sec.set_ylabel("km/s")
    for ax in axs.flat:
        ax.set_xlabel("R [arcsec]")
        pc_axis(ax)
        ax.legend(fontsize=8)
    fig.suptitle("Kinematic profiles: the 89 bins of the baseline fit (colour) and published profiles not fitted (grey/others)")
    return fig


def rotation():
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), constrained_layout=True)
    t = Table.read(PROC/"vasiliev2021_ocen_pm_profiles.ecsv")
    r = np.asarray(t["r"], float)
    ax = axs[0]
    ax.fill_between(r, t["vrot_pm_p02"], t["vrot_pm_p98"], color="C1", alpha=0.15, lw=0, label="2-98%")
    ax.fill_between(r, t["vrot_pm_p16"], t["vrot_pm_p84"], color="C1", alpha=0.35, lw=0, label="16-84%")
    ax.plot(r, t["vrot_pm"], color="C1", label="median")
    ax.set(xscale="log", xlim=(10, 2400), ylabel=r"$v_{\rm rot}$ [mas/yr]", xlabel="R [arcsec]",
           title="PM rotation curve, Vasiliev & Baumgardt 2021 (used for the PM rotation term)")
    sec = ax.secondary_yaxis("right", functions=(lambda y: y*KMS_PER_MASYR, lambda y: y/KMS_PER_MASYR))
    sec.set_ylabel("km/s")
    ax = axs[1]
    m = Table.read(PROC/"omegacat_vi_los_rotation.ecsv")
    ax.errorbar(m["r_median"], m["v_rot"], yerr=[m["v_rot_err_lo"], m["v_rot_err_hi"]], fmt="o", ms=3.5, color="C2", lw=0.8)
    ax.set(xscale="log", ylabel=r"$v_{\rm rot}$ [km/s]", xlabel="R [arcsec]",
           title="LOS rotation amplitude, MUSE (oMEGACat VI; LOS rotation term)")
    for ax in axs:
        pc_axis(ax)
        ax.legend(fontsize=8) if ax.get_legend_handles_labels()[0] else None
    return fig


def counts():
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    for name, c, lab in (("ocen_counts_hst_f625w19", BASE["hst"], "HST oMEGACat, F625W < 19"),
                         ("ocen_counts_gaia_g17", BASE["gaia"], "Gaia EDR3 (V&B 2021 field), G < 17")):
        t = Table.read(PROC/f"{name}.ecsv")
        r = np.asarray(t["r_median"], float)
        xerr = [r-t["r_lower"], t["r_upper"]-r]
        axs[0].errorbar(r, t["density_per_arcmin2"], xerr=xerr, yerr=t["poisson_err_per_arcmin2"], fmt="o", ms=3.5,
                        color=c, lw=0.8, label=f"{lab}: {int(np.sum(t['n_stars']))} stars, {len(t)} bins")
        axs[1].step(np.r_[t["r_lower"], t["r_upper"][-1]], np.r_[t["n_stars"], t["n_stars"][-1]], where="post", color=c, label=lab)
        axs[2].step(np.r_[t["r_lower"], t["r_upper"][-1]], np.r_[t["coverage"], t["coverage"][-1]], where="post", color=c, label=lab)
    axs[0].set(xscale="log", yscale="log", ylabel=r"stars arcmin$^{-2}$", title="Surface number density (Poisson errors)")
    axs[1].set(xscale="log", yscale="log", ylabel="stars per bin", title="Counts per bin")
    axs[2].set(xscale="log", ylim=(0.9, 1.01), ylabel="covered fraction of the annulus", title="Footprint coverage")
    for ax in axs:
        ax.set_xlabel("R [arcsec]")
        pc_axis(ax)
        ax.legend(fontsize=8)
    fig.suptitle("Star-count products (tracer density in the likelihood; the field level is fitted for Gaia)")
    return fig


def stream54():
    from ocen_dm.streams.path import project_gc
    import astropy.coordinates as coord
    t = Table.read(os.path.expandvars("$OCEN_CATALOGUES/streamfinder_ibata2024_dr3.fits"))
    t = t[np.asarray(t["Stream"]) == 54]
    g = coord.SkyCoord(np.asarray(t["RAdeg"], float), np.asarray(t["DEdeg"], float), unit="deg").galactic
    l = np.where(g.l.deg > 180, g.l.deg-360, g.l.deg)
    b = g.b.deg
    GC = json.loads((ROOT/"results/plot_data/stream54_gc_frame_chord.json").read_text())
    phi1, phi2 = project_gc(l, b, GC)
    rv = np.asarray(t["VHel"], float)
    has_rv = np.asarray(t["e_VHel"], float) < 300
    fig, axs = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    ax = axs[0, 0]
    sc = ax.scatter(l, b, c=np.asarray(t["Gmag"], float), s=3, cmap="viridis_r")
    ax.plot(309.10-360, 14.97, "r*", ms=12, label=r"$\omega$ Cen")
    ax.set(xlabel="l [deg]", ylabel="b [deg]", title=f"Stream 54 members on the sky ({len(t)}; colour: G)")
    ax.invert_xaxis()
    fig.colorbar(sc, ax=ax, label="G [mag]")
    ax.legend(fontsize=8)
    for ax, col, lab in ((axs[0, 1], "pmRA", r"$\mu_{\alpha*}$"), (axs[1, 0], "pmDE", r"$\mu_\delta$")):
        ax.scatter(phi1, np.asarray(t[col], float), s=3, color="k", alpha=0.5)
        ax.set(xlabel=r"$\phi_1$ [deg] (chord great-circle frame)", ylabel=f"{lab} [mas/yr]",
               title=f"{lab} along the stream")
    ax = axs[1, 1]
    ax.scatter(phi1[has_rv], rv[has_rv], s=14, color="C3")
    ax.set(xlabel=r"$\phi_1$ [deg]", ylabel=r"$v_{\rm hel}$ [km/s]", title=f"Radial velocities ({int(has_rv.sum())} members)")
    fig.suptitle("STREAMFINDER (Ibata+2024) stream 54 (Fimbulthul): the data the stream models are fitted to")
    return fig


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--only", default="kinematics,rotation,counts,stream54")
    a = p.parse_args()
    for name in a.only.split(","):
        fig = globals()[name]()
        out = ROOT/"plots"/f"data_{name}.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(out)


if __name__ == "__main__":
    main()
