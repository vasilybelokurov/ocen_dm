#!/usr/bin/env python3
"""Completeness of the Gaia G < 17 count sample against radius, measured with HST fields.

Gaia completeness c_G = fraction of HST stars brighter than the HST equivalent of G = 17 that are present in the
Vasiliev & Baumgardt (2021) EDR3 catalogue used for the Gaia count product (ocen_counts_gaia_g17).

HST references:
- oMEGACat (Haeberle et al. 2025), F625W, to ~350" (as bin/check_count_ratio_corrected.py).
- HST Large Programme on omega Cen (GO-14118+14662, PI Bedin), $OCEN_CATALOGUES/omegacen_hst_lp:
  F2 + F3 (Scalco et al. 2021, MNRAS 505, 3549), F1, F4, F5 (Scalco, Bedin & Vesperini 2024, A&A 688, A180); F606W.
For each field: epoch offset = median Gaia - HST position offset (0.5" first pass), then a 0.1" match. The HST
equivalent of G = 17 is 17 + median(m_HST - G) of matched stars at 16 < G < 17.5 (F606W) or 16.5 < G < 17.5 (F625W).
Binomial errors; for k = n the 68% lower limit 0.32^(1/n) is quoted. Bins: 250-300" plus the Gaia count-product bins.

Output: plots/gaia_completeness_lp.png, results/plot_data/gaia_completeness_lp.json.
Usage: PYTHONPATH=src python bin/check_gaia_completeness_lp.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.spatial import cKDTree

from ocen_dm.cluster import OCEN_DEC, OCEN_RA
from ocen_dm.paths import processed_dir, raw_dir

LP = Path(os.path.expandvars("$OCEN_CATALOGUES/omegacen_hst_lp")).expanduser()
LP_FIELDS = {"F5": ("P06/CATALOG.F5/CATALOG.F5", 3), "F4": ("P06/CATALOG.F4/CATALOG.F4", 3),
             "F2+F3": ("P04/CATALOG/CATALOG", 4), "F1": ("P06/CATALOG.F1/CATALOG.F1", 3)}
TOL = 0.1


def gaia():
    g = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    c = np.cos(np.radians(OCEN_DEC))
    return (np.asarray(g["ra"], float)-OCEN_RA)*c*3600, (np.asarray(g["dec"], float)-OCEN_DEC)*3600, np.asarray(g["g_mag"], float)


def load_lp(sub, ic):
    d = LP/sub
    a = np.genfromtxt(d/"ID_XY_RD_PM.dat", comments="#", usecols=(ic, ic+1))
    m = np.genfromtxt(d/"F606W.m1.dat", comments="#", usecols=(0,))
    m[m < -90] = np.nan
    c = np.cos(np.radians(OCEN_DEC))
    return (a[:, 0]-OCEN_RA)*c*3600, (a[:, 1]-OCEN_DEC)*3600, m


def load_omegacat():
    h = fits.open(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")[1].data
    c = np.cos(np.radians(OCEN_DEC))
    x = (np.asarray(h["RA"], float)-OCEN_RA)*c*3600; y = (np.asarray(h["DEC"], float)-OCEN_DEC)*3600
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok], np.asarray(h["f625w"], float)[ok]


def found_flags(x, y, m, gx, gy, G, colour_range):
    """Match HST stars to Gaia; return (R, found, m_eq, colour, offset) for HST stars brighter than the G = 17 equivalent."""
    t = cKDTree(np.column_stack([x, y]))
    near = np.hypot(gx-np.median(x), gy-np.median(y)) < max(400., 1.2*np.ptp(x))
    d0, j0 = t.query(np.column_stack([gx[near], gy[near]]), distance_upper_bound=0.5); o = np.isfinite(d0)
    ox = float(np.median(gx[near][o]-x[j0[o]])); oy = float(np.median(gy[near][o]-y[j0[o]]))
    d1, j1 = t.query(np.column_stack([gx[near]-ox, gy[near]-oy]), distance_upper_bound=TOL); o1 = np.isfinite(d1)
    Gm, mm = G[near][o1], m[j1[o1]]
    k = (Gm > colour_range[0]) & (Gm < colour_range[1]) & np.isfinite(mm)
    colour = float(np.median(mm[k]-Gm[k])); m_eq = 17.+colour
    sel = np.isfinite(m) & (m < m_eq)
    gt = cKDTree(np.column_stack([gx-ox, gy-oy]))
    dd, _ = gt.query(np.column_stack([x[sel], y[sel]]), distance_upper_bound=TOL)
    return np.hypot(x[sel], y[sel]), np.isfinite(dd), m_eq, colour, (ox, oy), int(k.sum())


def frac(k, n):
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k/n
    if k == n:
        return p, 1.-0.32**(1./n), 0.
    s = np.sqrt(p*(1-p)/n)
    return p, s, s


def main():
    gx, gy, G = gaia()
    prod = json.loads((processed_dir()/"kinematics"/"ocen_counts_gaia_g17.json").read_text())
    edges = [250.]+list(prod["r_lower"])+[prod["r_upper"][-1]]
    sources = {"oMEGACat (F625W)": found_flags(*load_omegacat(), gx, gy, G, (16.5, 17.5))}
    for name, (sub, ic) in LP_FIELDS.items():
        sources[f"LP {name} (F606W)"] = found_flags(*load_lp(sub, ic), gx, gy, G, (16., 17.5))
    out = dict(match_tolerance_arcsec=TOL, bin_edges=edges, sources={}, combined_lp=[])
    for name, (R, f, m_eq, col, off, nc) in sources.items():
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (R >= lo) & (R < hi)
            p, el, eh = frac(int(f[k].sum()), int(k.sum()))
            rows.append(dict(r_lo=lo, r_hi=hi, n=int(k.sum()), found=int(f[k].sum()), c=p, err_lo=el, err_hi=eh))
        out["sources"][name] = dict(m_eq=m_eq, colour=col, n_colour=nc, epoch_offset_arcsec=off, r_range=[float(R.min()), float(R.max())], bins=rows)
    # combined over the Large Programme fields (independent of oMEGACat)
    Rl = np.concatenate([v[0] for k, v in sources.items() if k.startswith("LP")])
    fl = np.concatenate([v[1] for k, v in sources.items() if k.startswith("LP")])
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (Rl >= lo) & (Rl < hi)
        p, el, eh = frac(int(fl[k].sum()), int(k.sum()))
        out["combined_lp"].append(dict(r_lo=lo, r_hi=hi, n=int(k.sum()), found=int(fl[k].sum()), c=p, err_lo=el, err_hi=eh))
    k = Rl >= 500
    out["combined_lp_beyond_500"] = dict(n=int(k.sum()), found=int(fl[k].sum()), c=float(fl[k].mean()),
                                         err=float(np.sqrt(fl[k].mean()*(1-fl[k].mean())/k.sum())))
    (ROOT/"results/plot_data/gaia_completeness_lp.json").write_text(json.dumps(out, indent=1))

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for e in edges[1:]:
        ax.axvline(e, color="0.9", lw=0.6, zorder=0)
    for i, (name, s) in enumerate(out["sources"].items()):
        b = [r for r in s["bins"] if r["n"] >= 10]
        Rm = np.array([np.sqrt(r["r_lo"]*r["r_hi"]) for r in b])*(1+0.012*(i-2))
        ax.errorbar(Rm, [r["c"] for r in b], yerr=[[r["err_lo"] for r in b], [r["err_hi"] for r in b]], fmt="os^vD"[i],
                    ms=4.5, capsize=0, lw=0.9, color=f"C{i}", mfc="none" if i == 0 else None,
                    label=f"{name}, {s['r_range'][0]:.0f}–{s['r_range'][1]:.0f}″")
    b = [r for r in out["combined_lp"] if r["n"] >= 10]
    for r in b:
        ax.plot([r["r_lo"], r["r_hi"]], [r["c"]]*2, color="k", lw=1.2, label="LP fields combined, per Gaia bin" if r is b[0] else None)
    ax.axvspan(300, edges[-1], color="C7", alpha=0.06, lw=0)
    ax.text(310, 0.715, "Gaia count bins (grey lines)", fontsize=7, color="0.35")
    ax.set_xscale("log"); ax.set_xlim(240, 1300); ax.set_ylim(0.7, 1.04); ax.axhline(1, color="k", lw=0.5)
    ax.set_xticks([250, 300, 400, 500, 700, 1000]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xticklabels(["250", "300", "400", "500", "700", "1000"])
    ax.set_xlabel("R [arcsec]"); ax.set_ylabel(r"Gaia completeness, $G<17$")
    ax.legend(fontsize=6.8, loc="lower right")
    fig.savefig(ROOT/"plots/gaia_completeness_lp.png", dpi=150, bbox_inches="tight")
    for name, s in out["sources"].items():
        print(f"{name}: m_eq {s['m_eq']:.2f} (colour {s['colour']:+.3f}, n {s['n_colour']}), R {s['r_range'][0]:.0f}-{s['r_range'][1]:.0f}")
        for r in s["bins"]:
            if r["n"]:
                print(f"   {r['r_lo']:5.0f}-{r['r_hi']:<5.0f} {r['found']:4d}/{r['n']:<4d} {r['c']:.3f} -{r['err_lo']:.3f} +{r['err_hi']:.3f}")
    print("LP combined:")
    for r in out["combined_lp"]:
        if r["n"]:
            print(f"   {r['r_lo']:5.0f}-{r['r_hi']:<5.0f} {r['found']:4d}/{r['n']:<4d} {r['c']:.3f} -{r['err_lo']:.3f} +{r['err_hi']:.3f}")
    print("LP beyond 500:", out["combined_lp_beyond_500"])


if __name__ == "__main__":
    main()
