#!/usr/bin/env python3
"""Ratio of the two star-count products, N(F625W < 19) / N(G < 17), measured where HST and Gaia overlap.

The DF fits give the HST and Gaia count profiles separate free amplitudes (a_HST, a_Gaia). For one cluster with one
stellar population their ratio is fixed by the luminosity function, so it can be measured and compared with the fit.

Method (read-only):
- Both catalogues around the same centre (cluster.OCEN_RA/DEC, catalogue sky coordinates) and counted in the same area:
  the 3-arcsec cells occupied by any oMEGACat star. The Gaia uniform field (1.60 arcmin^-2, fitted) is subtracted.
- Completeness: Gaia stars are matched to HST after removing the epoch offset (median over a 0.4" first pass), 0.06"
  tolerance (as selection/hst_gaia_match.py). Gaia completeness = fraction of HST stars with F625W < F_eq - 0.3 that
  have a Gaia counterpart; HST completeness = fraction of oMEGACat stars with a measured F625W.
- F_eq, the F625W equivalent of G = 17, is the median F625W - G of matched stars at 16.5 < G < 17.5 beyond 150".
- HST-only check: N(F625W < 19) / N(F625W < F_eq), the ratio a complete G < 17 sample would give.
- Fitted ratios a_HST/a_Gaia are read from the report.json of the base fit and the two fixed-BH fits.

Output: plots/count_ratio_overlap.png and results/plot_data/count_ratio_overlap.json.
Usage: PYTHONPATH=src python bin/check_count_ratio_overlap.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.spatial import cKDTree

from ocen_dm.cluster import OCEN_DEC, OCEN_RA
from ocen_dm.paths import processed_dir, raw_dir

CELL = 3.          # arcsec, footprint cells
B_FIELD = 1.60     # arcmin^-2, G < 17 uniform field of the fits
FITS = [("no BH", "dftwo_counts_area_20261013", "observed_no_halo_start3"),
        (r"BH 8,200 M$_\odot$", "dfbh8200_restart_20261015", None),
        (r"BH $4\times10^4$ M$_\odot$", "dfbh40000_restart_20261015", None)]


def fitted_ratio(batch, job):
    rep = json.loads((ROOT/"results/df"/batch/"report.json").read_text())
    jobs = rep["jobs"] if job is None else [j for j in rep["jobs"] if j["job"] == job]
    a = jobs[0]["gates"]["count_amplitudes"]
    return a["ocen_counts_hst_f625w19"]/a["ocen_counts_gaia_g17"]


def cell_key(x, y):
    return (np.floor(x/CELL).astype(np.int64)+10**5)*10**6+np.floor(y/CELL).astype(np.int64)+10**5


def main():
    c = np.cos(np.radians(OCEN_DEC))
    h = fits.open(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")[1].data
    hx = (np.asarray(h["RA"], float)-OCEN_RA)*c*3600; hy = (np.asarray(h["DEC"], float)-OCEN_DEC)*3600
    f = np.asarray(h["f625w"], float)
    hok = np.isfinite(hx) & np.isfinite(hy)
    hx, hy, f = hx[hok], hy[hok], f[hok]; hR = np.hypot(hx, hy)
    g = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    gx = (np.asarray(g["ra"], float)-OCEN_RA)*c*3600; gy = (np.asarray(g["dec"], float)-OCEN_DEC)*3600
    G = np.asarray(g["g_mag"], float); gR = np.hypot(gx, gy)

    # shared footprint and its area per annulus
    occ = np.unique(cell_key(hx, hy))
    gin = np.isin(cell_key(gx, gy), occ)
    off = ((np.arange(10)+.5)/10-.5)*CELL
    cx = (occ//10**6-10**5+.5)*CELL; cy = (occ % 10**6-10**5+.5)*CELL
    Rp = np.sort(np.hypot((cx[:, None, None]+off[None, :, None]+0*off[None, None, :]).ravel(),
                          (cy[:, None, None]+0*off[None, :, None]+off[None, None, :]).ravel()))
    dA = (CELL/10)**2/3600.

    # epoch offset, then a tight match
    htree = cKDTree(np.column_stack([hx, hy]))
    d0, j0 = htree.query(np.column_stack([gx, gy]), distance_upper_bound=0.4); o0 = np.isfinite(d0)
    ox = np.median(gx[o0]-hx[j0[o0]]); oy = np.median(gy[o0]-hy[j0[o0]])
    d1, j1 = htree.query(np.column_stack([gx-ox, gy-oy]), distance_upper_bound=0.06)
    m = np.isfinite(d1) & (G > 16.5) & (G < 17.5) & (gR > 150)
    dm = f[j1[m]]-G[m]; dm = dm[np.isfinite(dm)]
    colour = float(np.median(dm)); f_eq = 17.+colour
    gtree = cKDTree(np.column_stack([gx-ox, gy-oy]))

    edges = np.arange(50., 351., 25.)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        A = (np.searchsorted(Rp, hi)-np.searchsorted(Rp, lo))*dA
        k = (hR >= lo) & (hR < hi)
        nh = int(np.sum(k & (f < 19))); nheq = int(np.sum(k & (f < f_eq)))
        kg = gin & (gR >= lo) & (gR < hi)
        ng = int(np.sum(kg & (G < 17))); ng16 = int(np.sum(kg & (G < 16)))
        kb = k & (f > 12) & (f < f_eq-0.3)
        dd, _ = gtree.query(np.column_stack([hx[kb], hy[kb]]), distance_upper_bound=0.06)
        fld = B_FIELD*A; r = nh/(ng-fld)
        rows.append(dict(r_lo=lo, r_hi=hi, area_arcmin2=A, n_hst19=nh, n_gaia17=ng, n_gaia16=ng16, field=fld,
                         ratio=r, ratio_err=r*np.sqrt(1/nh+ng/(ng-fld)**2), gaia_complete=float(np.isfinite(dd).mean()),
                         n_bright=int(kb.sum()), hst_complete=float(np.isfinite(f[k]).mean()), lf_ratio=nh/nheq,
                         lf_ratio_err=np.sqrt(nh*(nh-nheq)/nheq**3)))
    # the single clean overlap annulus 250-300"
    lo, hi = 250., 300.
    A = (np.searchsorted(Rp, hi)-np.searchsorted(Rp, lo))*dA
    nh = int(np.sum((hR >= lo) & (hR < hi) & (f < 19))); ng = int(np.sum(gin & (gR >= lo) & (gR < hi) & (G < 17)))
    r250 = nh/(ng-B_FIELD*A); e250 = r250*np.sqrt(1/nh+ng/(ng-B_FIELD*A)**2)
    fitted = {lab: fitted_ratio(b, j) for lab, b, j in FITS}

    out = dict(colour_f625w_minus_g=colour, f_eq=f_eq, epoch_offset_arcsec=[ox, oy], overlap_250_300=dict(ratio=r250, err=e250,
               n_hst19=nh, n_gaia17=ng, area_arcmin2=A), fitted_ratio={k.replace("$", ""): v for k, v in fitted.items()}, bins=rows)
    (ROOT/"results/plot_data/count_ratio_overlap.json").write_text(json.dumps(out, indent=1))

    R = np.array([0.5*(r["r_lo"]+r["r_hi"]) for r in rows])
    gc = np.array([r["gaia_complete"] for r in rows]); hc = np.array([r["hst_complete"] for r in rows])
    rr = np.array([r["ratio"] for r in rows]); re = np.array([r["ratio_err"] for r in rows])
    lf = np.array([r["lf_ratio"] for r in rows]); lfe = np.array([r["lf_ratio_err"] for r in rows])
    good = (gc > 0.9) & (hc > 0.9)

    fig, ax = plt.subplots(2, 1, figsize=(6.0, 7.4), sharex=True, gridspec_kw=dict(hspace=0.06, height_ratios=[1, 1.5]))
    ax[0].plot(R, gc, "o-", color="C0", label=rf"Gaia: bright HST stars ($F625W<{f_eq-0.3:.1f}$) found in Gaia")
    ax[0].plot(R, hc, "s-", color="C3", label="HST: oMEGACat stars with a measured F625W")
    ax[0].axhline(1, color="k", lw=0.5); ax[0].set_ylim(0, 1.08); ax[0].set_ylabel("completeness")
    ax[0].legend(loc="lower left", fontsize=7.5)
    ax[1].errorbar(R[good], rr[good], re[good], fmt="o", color="k", ms=5, label="direct, both catalogues > 0.9 complete")
    ax[1].errorbar(R[~good], rr[~good], re[~good], fmt="o", mfc="none", color="0.5", ms=5, label="direct, one catalogue < 0.9 complete")
    ax[1].errorbar(R, lf, lfe, fmt="s", mfc="none", color="C2", ms=4,
                   label=rf"HST only: $N(F625W<19)/N(F625W<{f_eq:.2f})$")
    ax[1].errorbar(275, r250, e250, xerr=25, fmt="D", color="C1", ms=6, zorder=6, label=rf"250–300″: {r250:.2f} $\pm$ {e250:.2f}")
    for (lab, v), ls in zip(fitted.items(), ("-", "--", ":")):
        ax[1].axhline(v, color="C4", ls=ls, lw=1.1, label=rf"fitted $a_{{\rm HST}}/a_{{\rm Gaia}}$, {lab}: {v:.2f}")
    ax[1].set_ylim(5, 25); ax[1].set_ylabel(r"$N(F625W<19)\,/\,N(G<17)$")
    ax[1].set_xlabel("R [arcsec]"); ax[1].set_xlim(45, 355)
    ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=7, ncol=2, frameon=False)
    for j, a in enumerate(ax):
        a.text(0.015, 0.95, "ab"[j], transform=a.transAxes, va="top", fontweight="bold")
    fig.savefig(ROOT/"plots/count_ratio_overlap.png", dpi=150, bbox_inches="tight")
    print(f"F625W-G = {colour:.3f}, F_eq = {f_eq:.2f}; 250-300: {r250:.2f} +- {e250:.2f}; fitted", {k: round(v, 2) for k, v in fitted.items()})
    for r in rows:
        print(f"{r['r_lo']:4.0f}-{r['r_hi']:<4.0f} ratio {r['ratio']:6.2f}+-{r['ratio_err']:.2f}  LF {r['lf_ratio']:5.2f}  "
              f"Gaia {r['gaia_complete']:.2f}  HST {r['hst_complete']:.2f}  G17/G16 {r['n_gaia17']/max(r['n_gaia16'], 1):.2f}")


if __name__ == "__main__":
    main()
