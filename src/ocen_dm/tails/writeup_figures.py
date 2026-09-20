"""Figures for docs/progenitor_orbits.tex. Usage: PYTHONPATH=src python -m ocen_dm.tails.writeup_figures
Writes plots/po_*.png and results/tails/po_comparison.ecsv."""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table

from .progenitor_orbits import (plausible_orbit_set, class2_friction_backwards_fast, exponential_stripping,
                                _agama_potential, M_NUCLEUS_MSUN, CLASS2_SET)
from . import bar_migration as bm

PLOTS = Path("plots"); OUT = Path("results/tails")
COL = {"class1": ["C0", "C1", "C2"], "class2": ["C3", "C4", "C5"], "class3": ["C6", "C7", "C8"]}
_G = 4.30091e-6


def per_gyr_elements(o, edges=None):
    """peri, apo, and peak tidal strength G M(<r_peri)/r_peri^3 per look-back-time bin."""
    t = -o.t_gyr
    if edges is None:
        edges = np.arange(0.0, t.max() + 1e-9, 1.0)
    pot = _agama_potential("McMillan17")
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (t >= lo) & (t < hi)
        if m.sum() < 10:
            continue
        r = o.r[m]; ip = np.argmin(r); rp = r[ip]
        Menc = float(pot.enclosedMass(rp))
        rows.append(dict(t_mid=0.5 * (lo + hi), peri=rp, apo=r.max(), tidal=_G * Menc / rp**3,
                         rj_nucleus=rp * (M_NUCLEUS_MSUN / (3 * Menc)) ** (1 / 3),
                         rj_host1e9=rp * (1e9 / (3 * Menc)) ** (1 / 3)))
    return Table(rows=rows)


def main():
    PLOTS.mkdir(exist_ok=True); OUT.mkdir(parents=True, exist_ok=True)
    S = plausible_orbit_set()

    # ---------------- Fig: class 1 (zoom 3 Gyr + meridional plane)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for o, c in zip(S["class1"], COL["class1"]):
        m = -o.t_gyr < 3
        ax[0].plot(-o.t_gyr[m], o.r[m], c=c, lw=0.8, label="%s: peri %.2f, apo %.2f, e %.2f" % (o.potential, o.r_peri, o.r_apo, o.eccentricity))
        ax[1].plot(o.R[m], o.z[m], c=c, lw=0.3, alpha=0.8)
    ax[0].set_xlabel("look-back time [Gyr]"); ax[0].set_ylabel("r [kpc]"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].set_xlabel("R [kpc]"); ax[1].set_ylabel("z [kpc]"); ax[1].set_aspect("equal"); ax[1].grid(alpha=0.3)
    ax[1].set_title("meridional plane, last 3 Gyr", fontsize=10)
    fig.tight_layout(); fig.savefig(PLOTS / "po_class1.png", dpi=150); plt.close(fig)

    # ---------------- Fig: class 2 mass histories and r(t) with constant-mass comparison
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for (label, m_inf, tau, t_inf), o, c in zip(CLASS2_SET, S["class2"], COL["class2"]):
        h = exponential_stripping(m_inf, tau, t_inf)
        tl = np.linspace(0, 10, 500)
        ax[0].plot(tl, [h(x) for x in tl], c=c, label=label)
        ax[1].plot(-o.t_gyr, o.r, c=c, lw=0.5, label=label)
    for m_inf, c in ((1e10, "grey"), (1e11, "k")):
        oc = class2_friction_backwards_fast(m_inf, 1.0 * (m_inf / 1e10) ** (1 / 3), t_gyr=10.0)
        ax[2].plot(-oc.t_gyr, oc.r, c=c, lw=0.5, label="constant M = %.0e" % m_inf)
    ax[2].plot(-S["class2"][2].t_gyr, S["class2"][2].r, c=COL["class2"][2], lw=0.5, label="stripped 1e11, tau 0.75")
    ax[0].set_yscale("log"); ax[0].set_ylabel(r"bound mass $M(t)$ [M$_\odot$]"); ax[0].set_xlabel("look-back time [Gyr]")
    ax[0].axhline(M_NUCLEUS_MSUN, color="k", ls=":", lw=0.7); ax[0].text(0.2, 4.5e6, "nucleus", fontsize=8)
    for a in ax[1:]:
        a.set_yscale("log"); a.set_ylabel("r [kpc]"); a.set_xlabel("look-back time [Gyr]")
        a.axhline(100, color="grey", ls="--", lw=0.7)
    ax[1].set_title("stripped histories (the picks)", fontsize=10); ax[2].set_title("why constant mass fails", fontsize=10)
    for a in ax: a.legend(fontsize=7); a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(PLOTS / "po_class2.png", dpi=150); plt.close(fig)

    # ---------------- Fig: class 2 scan heat map (apocentre in the last Gyr before infall)
    taus = np.arange(0.5, 3.01, 0.25); masses = (1e10, 3e10, 1e11, 2e11)
    grid = np.zeros((len(masses), len(taus)))
    for i, m_inf in enumerate(masses):
        rh = 1.0 * (m_inf / 1e10) ** (1 / 3)
        for j, tau in enumerate(taus):
            o = class2_friction_backwards_fast(m_inf, rh, t_gyr=10.0, mass_history=exponential_stripping(m_inf, tau, 10.0))
            grid[i, j] = o.r[o.t_gyr < -9.0].max()
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    im = ax.imshow(np.log10(grid), aspect="auto", origin="lower", cmap="viridis",
                   extent=[taus[0] - 0.125, taus[-1] + 0.125, -0.5, len(masses) - 0.5])
    for i in range(len(masses)):
        for j in range(len(taus)):
            ok = 50 <= grid[i, j] <= 150
            ax.text(taus[j], i, "%.0f" % grid[i, j], ha="center", va="center", fontsize=7,
                    color="w" if not ok else "yellow", fontweight="bold" if ok else "normal")
    ax.set_yticks(range(len(masses))); ax.set_yticklabels(["%.0e" % m for m in masses])
    ax.set_xlabel(r"stripping e-folding time $\tau$ [Gyr]"); ax.set_ylabel(r"$M_{\rm inf}$ [M$_\odot$]")
    fig.colorbar(im, label=r"log$_{10}$ apocentre 9-10 Gyr ago [kpc]"); ax.set_title("bold: 50-150 kpc (virial radius of the young MW)", fontsize=9)
    fig.tight_layout(); fig.savefig(PLOTS / "po_class2_scan.png", dpi=150); plt.close(fig)
    Table(grid, names=["tau_%.2f" % t for t in taus]).write(OUT / "po_class2_scan.ecsv", overwrite=True)

    # ---------------- Fig: class 3 bar history + meridional planes before/after migration
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for w, ls in ((22, ":"), (24, "-"), (26, "--"), (30, "-."), (35, (0, (1, 1)))):
        h = bm.BarHistory(float(w)); ax[0].plot(h.tf - h.t, h.omega, "k", ls=ls, lw=1, label=r"$\Omega_{b,0}$ = %d  ($\Omega_1$ = %.0f)" % (w, h.omega_1))
    h = bm.BarHistory(24.0); a2 = ax[0].twinx(); a2.plot(h.tf - h.t, h.frac, c="C1", lw=1); a2.set_ylabel("bar amplitude fraction", color="C1")
    ax[0].set_xlabel("look-back time [Gyr]"); ax[0].set_ylabel(r"$\Omega_b$ [km/s/kpc]"); ax[0].legend(fontsize=7); ax[0].set_ylim(15, 120); ax[0].grid(alpha=0.3)
    o = S["class3"][1]; t = -o.t_gyr
    for a, (lo, hi, ttl) in zip(ax[1:], ((5, 8, "pre-migration, 5-8 Gyr ago"), (0, 2, "last 2 Gyr"))):
        m = (t >= lo) & (t < hi)
        a.plot(o.R[m], o.z[m], c=COL["class3"][1], lw=0.25, alpha=0.8); a.set_aspect("equal")
        a.set_xlabel("R [kpc]"); a.set_ylabel("z [kpc]"); a.set_title(ttl + " (sample %d)" % o.notes["sample"], fontsize=10); a.grid(alpha=0.3)
        a.set_xlim(0, 13); a.set_ylim(-6, 6)
    fig.tight_layout(); fig.savefig(PLOTS / "po_class3_bar.png", dpi=150); plt.close(fig)

    # ---------------- Fig: cross-class comparison (per-Gyr peri/apo, tidal strength, Jacobi radius)
    fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    rows = []
    for cls in ("class1", "class2", "class3"):
        for o, c in zip(S[cls], COL[cls]):
            tab = per_gyr_elements(o)
            lab = o.label if cls != "class1" else "class1 " + o.potential
            ax[0, 0].plot(tab["t_mid"], tab["peri"], "o-", c=c, ms=3, lw=0.8, label=lab)
            ax[0, 1].plot(tab["t_mid"], tab["apo"], "o-", c=c, ms=3, lw=0.8)
            ax[1, 0].plot(tab["t_mid"], tab["tidal"], "o-", c=c, ms=3, lw=0.8)
            ax[1, 1].plot(tab["t_mid"], tab["rj_nucleus"] * 1e3, "o-", c=c, ms=3, lw=0.8)
            for r in tab:
                rows.append(dict(cls=cls, label=lab, **{k: float(r[k]) for k in tab.colnames}))
    ax[0, 0].set_ylabel("pericentre per Gyr [kpc]"); ax[0, 0].set_yscale("log"); ax[0, 0].legend(fontsize=6.5, ncol=2)
    ax[0, 1].set_ylabel("apocentre per Gyr [kpc]"); ax[0, 1].set_yscale("log")
    ax[1, 0].set_ylabel(r"$G M(<r_{\rm peri})/r_{\rm peri}^3$ [(km/s/kpc)$^2$]"); ax[1, 0].set_yscale("log")
    ax[1, 1].set_ylabel(r"Jacobi radius of the $3.55\times10^6$ M$_\odot$ nucleus at pericentre [pc]"); ax[1, 1].set_yscale("log")
    for a in ax.ravel(): a.grid(alpha=0.3)
    for a in ax[1]: a.set_xlabel("look-back time [Gyr]")
    fig.suptitle("Per-Gyr orbital elements of the nine orbits (tides evaluated in McMillan 2017 for all)", fontsize=11)
    fig.tight_layout(); fig.savefig(PLOTS / "po_comparison.png", dpi=150); plt.close(fig)
    Table(rows=rows).write(OUT / "po_comparison.ecsv", overwrite=True)

    # ---------------- Fig: E and Lz histories (relative to today) for all nine
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
    pot = _agama_potential("McMillan17")
    for cls in ("class1", "class2"):
        for o, c in zip(S[cls], COL[cls]):
            if cls == "class1" and o.potential != "McMillan17":
                continue          # E would not be conserved when evaluated in a different potential
            xyz = np.column_stack([o.R, np.zeros_like(o.R), o.z])
            E = 0.5 * (o.vR**2 + o.vT**2 + o.vz**2) + pot.potential(xyz)
            Lz = o.R * o.vT
            sl = slice(None, None, max(1, len(E) // 4000))
            ax[0].plot(-o.t_gyr[sl], (E[sl] - E[0]) / 1e5, c=c, lw=0.6, label=o.label if cls == "class2" else "class1 " + o.potential)
            ax[1].plot(-o.t_gyr[sl], Lz[sl], c=c, lw=0.6)
    fid = np.load(OUT / "bar_migration_fiducial.npz")
    lb = 8.0 - fid["t"]
    for i, c in zip(fid["picks"], COL["class3"]):
        ax[0].plot(lb, (fid["E"][:, i] - fid["E"][-1, i]) / 1e5, c=c, lw=1.0, label="class3 sample %d" % i)
        ax[1].plot(lb, fid["Lz"][:, i], c=c, lw=1.0)       # oCen_bar frame is already prograde-positive (omega Cen -533)
    ax[0].set_ylabel(r"$E(t) - E_{\rm today}$ [$10^5$ km$^2$ s$^{-2}$]")
    ax[1].set_ylabel(r"$L_z$ [kpc km/s], prograde positive")
    ax[1].axhline(0, color="k", lw=0.5)
    for a in ax: a.set_xlabel("look-back time [Gyr]"); a.grid(alpha=0.3)
    ax[0].legend(fontsize=7); ax[0].set_xlim(0, 10)
    fig.tight_layout(); fig.savefig(PLOTS / "po_E_Lz.png", dpi=150); plt.close(fig)
    print("figures written")


if __name__ == "__main__":
    main()
