"""Predicted DM density inside omega Cen if the dark matter within the nucleus's Jacobi radius is
retained: rho(r) = rho_NFW(r; M200, c) exp(-r / r_J)  (EN21-style exponential truncation at r_J;
the sharp-cut alternative is rho_NFW for r < r_J, 0 outside). Compared with the stellar density
of the rung-0 K1 model and the 0.1-3 Msun/pc^3 range discussed in docs/dm_capture_constraints.tex
for the deep white-dwarf field at 20 pc.

Usage: PYTHONPATH=src python -m ocen_dm.tails.dm_density_plots
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
from .tidal_tracks import NFW

MSUN_PC3_TO_GEV_CM3 = 37.97
M200S = (1e9, 1e10, 1e11)
CS = (5.0, 10.0, 15.0)
RJ_CLASS3 = (33.0, 36.0)      # pc, from tidal_tracks budget (0.4 kpc pericentre)
RJ_CLASS12 = (40.0, 70.0)     # pc (1.25-2.0 kpc pericentres)
RADII_PC = (3.0, 10.0, 20.0)  # inner, intermediate, deep WD field (Scalco et al. 2024)
WD_BAND = (0.1, 3.0)          # Msun/pc^3, plausible range in dm_capture_constraints.tex


def rho_retained(halo: NFW, r_pc, rj_pc, sharp=False):
    """DM density [Msun/pc^3] at r inside the cluster for the initial NFW cusp retained inside r_J."""
    r = np.asarray(r_pc, float) / 1e3
    rho = halo.rho(r) / 1e9                                  # Msun/kpc^3 -> Msun/pc^3
    if sharp:
        return np.where(np.asarray(r_pc) < rj_pc, rho, 0.0)
    return rho * np.exp(-np.asarray(r_pc) / rj_pc)


def stellar_density(r_pc, m_star=2.88e6, distance_kpc=5.43):
    """Stellar mass density [Msun/pc^3] from the project's MGE light model scaled to the rung-0 K1 M_star."""
    from ..light_model import build_stellar_mge, fit_mge_projected, load_tracer_profile
    fit = fit_mge_projected(load_tracer_profile("composite"), sigma_range_arcsec=(7.0, 3000.0))
    stars = build_stellar_mge(fit, distance_kpc, m_star)
    return np.asarray(stars.density(np.asarray(r_pc, float)))


def _gev_axis(ax):
    sec = ax.secondary_yaxis("right", functions=(lambda y: y * MSUN_PC3_TO_GEV_CM3, lambda y: y / MSUN_PC3_TO_GEV_CM3))
    sec.set_ylabel(r"$\rho_{\rm DM}$ [GeV cm$^{-3}$]")
    return sec


def figure_profiles(path="plots/dm_density_in_rj_profiles.png"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    r = np.logspace(0, np.log10(300), 400)
    rho_star = stellar_density(r)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    cols = {5.0: "C0", 10.0: "C1", 15.0: "C3"}
    for ax, m200 in zip(axes, M200S):
        for c in CS:
            h = NFW(m200, c)
            ax.plot(r, rho_retained(h, r, 70.0), color=cols[c], lw=1.4, label="c = %.0f, $r_J$ = 70 pc (classes 1-2)" % c)
            ax.plot(r, rho_retained(h, r, 35.0), color=cols[c], lw=1.4, ls="--", label="c = %.0f, $r_J$ = 35 pc (class 3)" % c)
        ax.plot(r, rho_star, "k", lw=2, label=r"stars (K1 fit, $M_\star$ = 2.9e6)")
        ax.axhspan(*WD_BAND, color="grey", alpha=0.15)
        for x, lab in zip(RADII_PC, ("3 pc", "10 pc", "20 pc: deep WD field")):
            ax.axvline(x, color="k", ls=":", lw=0.7); ax.text(x * 1.05, 2e3, lab, rotation=90, fontsize=7, va="top")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1, 300); ax.set_ylim(1e-2, 5e3)
        ax.set_xlabel("r [pc]"); ax.set_title(r"initial NFW $M_{200}$ = %.0e M$_\odot$ (z = 2)" % m200, fontsize=10); ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"$\rho$ [M$_\odot$ pc$^{-3}$]"); axes[0].legend(fontsize=6.5, loc="lower left", ncol=1)
    _gev_axis(axes[2])
    fig.suptitle(r"DM retained inside the nucleus's Jacobi radius: $\rho_{\rm NFW}(r)\,e^{-r/r_J}$; grey band 0.1-3 M$_\odot$ pc$^{-3}$ (WD-heating note)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def figure_vs_rj(path="plots/dm_density_in_rj_vs_rj.png"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rj = np.logspace(1, np.log10(200), 200)
    rho_star = stellar_density(np.array(RADII_PC))
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    cols = {1e9: "C0", 1e10: "C1", 1e11: "C3"}; ls = {5.0: ":", 10.0: "-", 15.0: "--"}
    for ax, r0, rs in zip(axes, RADII_PC, rho_star):
        for m200 in M200S:
            for c in CS:
                h = NFW(m200, c)
                ax.plot(rj, rho_retained(h, r0, rj), color=cols[m200], ls=ls[c], lw=1.3,
                        label=(r"$M_{200}$ = %.0e, c = %.0f" % (m200, c)) if r0 == RADII_PC[0] else None)
        ax.axvspan(*RJ_CLASS3, color="C6", alpha=0.3); ax.text(RJ_CLASS3[0], 2e2, "class 3", color="C6", fontsize=8, rotation=90)
        ax.axvspan(*RJ_CLASS12, color="C0", alpha=0.12); ax.text(RJ_CLASS12[1] * 1.02, 2e2, "classes 1-2", color="C0", fontsize=8, rotation=90)
        ax.axhspan(*WD_BAND, color="grey", alpha=0.15)
        ax.axhline(rs, color="k", lw=2); ax.text(11, rs * 1.15, r"stars at this radius (%.0f M$_\odot$/pc$^3$)" % rs, fontsize=7)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(10, 200); ax.set_ylim(1e-2, 5e3)
        ax.set_xlabel(r"Jacobi radius of the nucleus $r_J$ [pc]"); ax.grid(alpha=0.3)
        ax.set_title("DM density at r = %.0f pc%s" % (r0, " (deep WD field)" if r0 == 20 else ""), fontsize=10)
    axes[0].set_ylabel(r"$\rho_{\rm DM}(r)$ [M$_\odot$ pc$^{-3}$]"); axes[0].legend(fontsize=6.5, loc="lower right", ncol=1)
    _gev_axis(axes[2])
    fig.suptitle(r"$\rho_{\rm DM}(r) = \rho_{\rm NFW}(r)\,e^{-r/r_J}$ for DM retained inside the nucleus's Jacobi radius (colour: $M_{200}$; line style: c = 5 dotted, 10 solid, 15 dashed)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def grid_table(path="results/tails/dm_density_grid.ecsv"):
    from astropy.table import Table
    rows = []
    rs = stellar_density(np.array(RADII_PC))
    for m200 in M200S:
        for c in CS:
            h = NFW(m200, c)
            for rj in (35.0, 50.0, 70.0, 100.0):
                row = dict(m200=m200, c=c, rj_pc=rj, m_dm_in_rj=h.mass(rj / 1e3))
                for r0, s in zip(RADII_PC, rs):
                    row["rho_dm_%dpc" % r0] = float(rho_retained(h, r0, rj)); row["rho_star_%dpc" % r0] = float(s)
                rows.append(row)
    t = Table(rows=rows); Path(path).parent.mkdir(parents=True, exist_ok=True); t.write(path, overwrite=True); return t


if __name__ == "__main__":
    Path("plots").mkdir(exist_ok=True)
    figure_profiles(); figure_vs_rj(); t = grid_table()
    for r in t:
        if r["rj_pc"] in (35.0, 70.0):
            print("M200 %.0e c %4.1f rJ %3.0f: M_DM(<rJ) %.1e | rho_DM(3,10,20 pc) = %6.1f %5.2f %5.2f Msun/pc3 | stars %6.0f %5.1f %5.2f" % (
                r["m200"], r["c"], r["rj_pc"], r["m_dm_in_rj"], r["rho_dm_3pc"], r["rho_dm_10pc"], r["rho_dm_20pc"], r["rho_star_3pc"], r["rho_star_10pc"], r["rho_star_20pc"]))
