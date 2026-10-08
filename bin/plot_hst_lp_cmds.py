#!/usr/bin/env python3
"""CMDs of the HST Large Programme on omega Cen outer fields (Bedin's public e-material), members only.

Data (~/data/catalogues/omegacen_hst_lp, see README_provenance.txt):
  F2, F3 (Paper IV, Scalco+2021, arXiv:2107.08726): P04/CATALOG/CATALOG/*.dat
  F1, F4, F5 (Scalco, Bedin & Vesperini 2024, A&A 688, A180): P06/CATALOG.F?/CATALOG.F?/*.dat
All files of one catalogue are row-aligned. Photometry: KS2 method 2 (m2) in WFC3/UVIS F606W and F814W
(Vega). Selection: membership probability >= 90 per cent (the catalogues' PM-based value, given in per cent), found in both filters,
QFIT >= 0.8 in both, not saturated in the long exposures (flag from the m1 files; m2/m3 carry none). Radius from the centre (201.696991, -47.479472)
(Anderson & van der Marel 2010).
Absolute magnitudes: D = 5.43 kpc, E(B-V) = 0.12, A_V = 3.1 E(B-V); A_F606W, A_F814W from the MIST
bolometric-correction tables for WFC3/UVIS (median over the main sequence), as in bin/build_mass_function.py.
Top row: Hess diagram of members with the MIST main sequence (12.5 Gyr, [Fe/H] = -1.53); bottom row: zoom
on the white-dwarf region as points.

Usage: python bin/plot_hst_lp_cmds.py [--pmin 90] [--qfit 0.8] [--method m2]
Output: plots/cmd_hst_lp_fields.png, results/plot_data/cmd_hst_lp_fields.json
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
from astropy.coordinates import SkyCoord

BASE = Path(os.path.expanduser("~/data/catalogues/omegacen_hst_lp"))
CENTRE = (201.696991, -47.479472)
BANDS = ["WFC3_UVIS_F606W", "WFC3_UVIS_F814W"]
RAMP = LinearSegmentedColormap.from_list("blue", ["#e8f0fb", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
INK, MUTED = "#1f1f1e", "#6b6a64"
CB = np.arange(-1.0, 2.6, 0.02)
MB = np.arange(2.0, 14.0, 0.05)


def extinction_and_isochrone(age=12.5, feh=-1.53, ebv=0.12):
    import minimint
    ii = minimint.Interpolator(BANDS, data_prefix=os.path.expanduser("~/data/isochrones/minimint_hst")+"/")
    mass = np.geomspace(.09, 1.2, 1200)
    res = ii(mass, np.log10(age*1e9), feh)
    good = np.isfinite(res[BANDS[0]]) & (res["phase"] <= 0)
    p0 = np.column_stack((res["logteff"][good], res["logg"][good], np.full(good.sum(), feh), np.zeros(good.sum())))
    p1 = p0.copy(); p1[:, 3] = 3.1*ebv
    bc0, bc1 = ii.bolomInt(p0), ii.bolomInt(p1)
    a6 = float(np.median(bc0[BANDS[0]]-bc1[BANDS[0]])); a8 = float(np.median(bc0[BANDS[1]]-bc1[BANDS[1]]))
    return a6, a8, res[BANDS[0]][good], res[BANDS[0]][good]-res[BANDS[1]][good]


def load_field(kind, field, method):
    if kind == "P04":
        d = BASE/"P04/CATALOG/CATALOG"
        ast = np.loadtxt(d/"ID_XY_RD_PM.dat", comments="#")
        sel = ast[:, 1] == int(field[1])
        ra, de, prob = ast[:, 4], ast[:, 5], ast[:, 26]     # 27 columns; membership probability is the last
    else:
        d = BASE/f"P06/CATALOG.{field}/CATALOG.{field}"
        ast = np.loadtxt(d/"ID_XY_RD_PM.dat", comments="#")
        sel = np.ones(len(ast), bool)
        ra, de, prob = ast[:, 3], ast[:, 4], ast[:, 5]
    ph = {b: np.loadtxt(d/f"{b}.{method}.dat", comments="#") for b in ("F606W", "F814W")}
    sat = {b: np.loadtxt(d/f"{b}.m1.dat", comments="#", usecols=9) for b in ("F606W", "F814W")}   # flag only in m1
    if any(len(v) != len(ast) for v in list(ph.values())+list(sat.values())):
        raise ValueError(f"{field}: photometry and astrometry files are not row-aligned")
    out = dict(ra=ra[sel], dec=de[sel], prob=prob[sel])
    for b in ph:
        out[f"{b}_mag"], out[f"{b}_qfit"], out[f"{b}_sat"] = ph[b][sel][:, 0], ph[b][sel][:, 2], sat[b][sel]
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pmin", type=float, default=90., help="membership probability threshold [per cent; the catalogues give 0-100]")
    p.add_argument("--qfit", type=float, default=0.8)
    p.add_argument("--method", default="m2")
    p.add_argument("--distance", type=float, default=5.43)
    args = p.parse_args()
    a6, a8, iso_M, iso_c = extinction_and_isochrone()
    dm0 = 5*np.log10(args.distance*1e3/10)
    fields = [("P06", "F5"), ("P06", "F4"), ("P04", "F2"), ("P04", "F3"), ("P06", "F1")]
    fig, axes = plt.subplots(2, len(fields), figsize=(3.5*len(fields), 11))
    out = dict(distance_kpc=args.distance, dm0=dm0, A_F606W=a6, A_F814W=a8, pmin=args.pmin, qfit=args.qfit,
               method=args.method, fields=[])
    c0 = SkyCoord(*CENTRE, unit="deg")
    for i, (kind, field) in enumerate(fields):
        f = load_field(kind, field, args.method)
        m6, m8 = f["F606W_mag"], f["F814W_mag"]
        ok = (m6 > -90) & (m8 > -90) & (f["F606W_qfit"] >= args.qfit) & (f["F814W_qfit"] >= args.qfit) \
            & (f["F606W_sat"] == 0) & (f["F814W_sat"] == 0)
        mem = ok & (f["prob"] >= args.pmin)
        r = SkyCoord(f["ra"], f["dec"], unit="deg").separation(c0).arcmin
        M = m6 - dm0 - a6
        col = (m6 - m8) - (a6 - a8)
        wd = mem & (col < 0.6) & (M > 8)
        out["fields"].append(dict(field=field, n=int(len(m6)), n_quality=int(ok.sum()), n_members=int(mem.sum()),
                                  r_arcmin_median=float(np.median(r)), n_wd_box=int(wd.sum()),
                                  F606W_member_99pct=float(np.percentile(m6[mem], 99)) if mem.any() else None))
        ax = axes[0, i]
        H, _, _ = np.histogram2d(col[mem], M[mem], bins=(CB, MB))
        if H.sum():
            ax.pcolormesh(CB, MB, H.T, cmap=RAMP, norm=LogNorm(vmin=1, vmax=max(H.max(), 2)), rasterized=True)
        ax.plot(iso_c-(a6-a8)*0, iso_M, "-", color="#eb6834", lw=1.2, label="MIST MS")
        ax.add_patch(plt.Rectangle((-1.0, 8), 1.6, 5.5, fill=False, ec=MUTED, lw=0.8, ls="--"))
        ax.set_xlim(CB[0], CB[-1]); ax.set_ylim(MB[-1], MB[0])
        ax.set_title(f"{field} (r ~ {np.median(r):.1f}'): N = {mem.sum():,}\nmembers (P >= {args.pmin:.0f}%), QFIT >= {args.qfit}",
                     fontsize=9.5, color=INK, loc="left")
        if i == 0:
            ax.legend(fontsize=8, frameon=False, loc="lower left")
        ax = axes[1, i]
        ax.plot(col[wd], M[wd], ".", color="#2a78d6", ms=2.5, alpha=0.6, mew=0, rasterized=True)
        ax.set_xlim(-1.0, 0.6); ax.set_ylim(13.5, 8)
        ax.set_title(f"{field}: white-dwarf region, N = {wd.sum():,}", fontsize=9.5, color=INK, loc="left")
        for row in (0, 1):
            a = axes[row, i]
            a.tick_params(colors=MUTED, labelsize=9)
            for s in ("top", "right"):
                a.spines[s].set_visible(False)
            a.set_xlabel("(F606W - F814W)$_0$", color=INK)
            if i == 0:
                a.set_ylabel("$M_{F606W}$", color=INK)
    fig.suptitle(f"omega Cen HST Large Programme outer fields (WFC3/UVIS, KS2 {args.method}), PM members; "
                 f"D = {args.distance} kpc, E(B-V) = 0.12 (A_F606W = {a6:.3f}, A_F814W = {a8:.3f})\n"
                 "Top: log counts (dashed box = zoom). Bottom: white-dwarf region", fontsize=10.5, color=INK)
    fig.tight_layout()
    fig.savefig(ROOT/"plots/cmd_hst_lp_fields.png", dpi=130)
    (ROOT/"results/plot_data/cmd_hst_lp_fields.json").write_text(json.dumps(out, indent=1))
    for r_ in out["fields"]:
        print(r_)


if __name__ == "__main__":
    main()
