#!/usr/bin/env python3
"""Completeness-corrected ratio of the two star-count products, N(F625W < 19) / N(G < 17), from 2 to 350 arcsec.

Companion to bin/check_count_ratio_overlap.py (same centre, shared footprint, field and epoch offset). Completeness is
measured in the magnitude range of each sample:
- Gaia, c_G(R): fraction of oMEGACat stars with F625W < F_eq (the HST equivalent of G < 17) that have a Gaia counterpart
  within 0.06" (also 0.10", as a tolerance check), so it is weighted by the sample's own luminosity function. Also given
  in three magnitude slices.
- HST, c_H(R): fraction of stars with F814W < F814_eq (the F814W equivalent of F625W = 19, median colour at
  18.5 < F625W < 19) that have a measured F625W. Stars missing in both filters are not seen by this proxy.
Corrected ratio: [N_H19 / c_H] / [(N_G17 - field) / c_G]. If the completeness is right it is flat with radius and equals
the HST-only ratio N(F625W < 19) / N(F625W < F_eq).
The fitted a_HST/a_Gaia contains the completeness of each product (HST 2-250", Gaia 300-2400"); it is converted to a
complete-sample ratio as fitted x c_G,out / c_H,in, with c_H,in the count-weighted c_H over 2-250" and c_G,out between
c_G(250-300") and 1 (Gaia beyond 300" cannot be checked against HST).

Output: plots/count_ratio_corrected.png, results/plot_data/count_ratio_corrected.json.
Usage: PYTHONPATH=src python bin/check_count_ratio_corrected.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.spatial import cKDTree

from check_count_ratio_overlap import B_FIELD, CELL, FITS, cell_key, fitted_ratio
from ocen_dm.cluster import OCEN_DEC, OCEN_RA
from ocen_dm.paths import processed_dir, raw_dir

EDGES = np.array([2., 5., 10., 20., 35., 50., 75., 100., 125., 150., 175., 200., 225., 250., 275., 300., 325., 350.])


def main():
    c = np.cos(np.radians(OCEN_DEC))
    h = fits.open(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")[1].data
    hx = (np.asarray(h["RA"], float)-OCEN_RA)*c*3600; hy = (np.asarray(h["DEC"], float)-OCEN_DEC)*3600
    ok = np.isfinite(hx) & np.isfinite(hy)
    hx, hy = hx[ok], hy[ok]; f = np.asarray(h["f625w"], float)[ok]; f8 = np.asarray(h["f814w"], float)[ok]
    hR = np.hypot(hx, hy)
    g = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    gx = (np.asarray(g["ra"], float)-OCEN_RA)*c*3600; gy = (np.asarray(g["dec"], float)-OCEN_DEC)*3600
    G = np.asarray(g["g_mag"], float); gR = np.hypot(gx, gy)

    occ = np.unique(cell_key(hx, hy)); gin = np.isin(cell_key(gx, gy), occ)
    off = ((np.arange(10)+.5)/10-.5)*CELL
    cx = (occ//10**6-10**5+.5)*CELL; cy = (occ % 10**6-10**5+.5)*CELL
    Rp = np.sort(np.hypot((cx[:, None, None]+off[None, :, None]+0*off[None, None, :]).ravel(),
                          (cy[:, None, None]+0*off[None, :, None]+off[None, None, :]).ravel()))
    dA = (CELL/10)**2/3600.

    htree = cKDTree(np.column_stack([hx, hy]))
    d0, j0 = htree.query(np.column_stack([gx, gy]), distance_upper_bound=0.4); o0 = np.isfinite(d0)
    ox = np.median(gx[o0]-hx[j0[o0]]); oy = np.median(gy[o0]-hy[j0[o0]])
    d1, j1 = htree.query(np.column_stack([gx-ox, gy-oy]), distance_upper_bound=0.06)
    m = np.isfinite(d1) & (G > 16.5) & (G < 17.5) & (gR > 150)
    dm = f[j1[m]]-G[m]; f_eq = 17.+float(np.median(dm[np.isfinite(dm)]))
    k = np.isfinite(f) & np.isfinite(f8) & (f > 18.5) & (f < 19.)
    f8_eq = 19.-float(np.median(f[k]-f8[k]))
    gtree = cKDTree(np.column_stack([gx-ox, gy-oy]))
    slices = [(-np.inf, 15.5), (15.5, 16.5), (16.5, f_eq)]

    def found(sel, tol):
        d, _ = gtree.query(np.column_stack([hx[sel], hy[sel]]), distance_upper_bound=tol)
        return np.isfinite(d)

    rows = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        A = (np.searchsorted(Rp, hi)-np.searchsorted(Rp, lo))*dA
        kr = (hR >= lo) & (hR < hi)
        nh = int(np.sum(kr & (f < 19))); nheq = int(np.sum(kr & (f < f_eq)))
        kg = gin & (gR >= lo) & (gR < hi); ng = int(np.sum(kg & (G < 17))); fld = B_FIELD*A
        ks = kr & (f < f_eq)
        cg = float(found(ks, 0.06).mean()) if ks.sum() else np.nan
        cg10 = float(found(ks, 0.10).mean()) if ks.sum() else np.nan
        cgs = [float(found(kr & (f >= a) & (f < b), 0.06).mean()) if np.sum(kr & (f >= a) & (f < b)) else np.nan for a, b in slices]
        k8 = kr & np.isfinite(f8) & (f8 < f8_eq)
        ch = float(np.isfinite(f[k8]).mean()) if k8.sum() else np.nan
        raw = nh/(ng-fld) if ng > fld else np.nan
        err = raw*np.sqrt(1/nh+ng/(ng-fld)**2) if ng > fld else np.nan
        corr = raw*cg/ch
        rows.append(dict(r_lo=lo, r_hi=hi, n_hst19=nh, n_gaia17=ng, field=fld, raw=raw, raw_err=err, c_gaia=cg,
                         c_gaia_tol010=cg10, c_gaia_slices=cgs, n_gaia_sample=int(ks.sum()), c_hst=ch, n_hst_proxy=int(k8.sum()),
                         corrected=corr, corrected_err=corr*np.sqrt(1/nh+ng/(ng-fld)**2+(1-cg)/max(cg*ks.sum(), 1)) if ng > fld else np.nan,
                         lf_ratio=nh/nheq, lf_ratio_err=np.sqrt(nh*(nh-nheq)/nheq**3)))
    R = np.array([0.5*(r["r_lo"]+r["r_hi"]) for r in rows])
    # count-weighted HST completeness inside the HST product (2-250") and Gaia completeness at 250-300"
    inn = [r for r in rows if r["r_hi"] <= 250]
    ch_in = sum(r["n_hst19"] for r in inn)/sum(r["n_hst19"]/r["c_hst"] for r in inn)
    o = [r for r in rows if 250 <= r["r_lo"] < 300]
    cg_out = sum(r["c_gaia"]*r["n_gaia_sample"] for r in o)/sum(r["n_gaia_sample"] for r in o)
    fitted = {lab: fitted_ratio(b, j) for lab, b, j in FITS}
    band = {lab: (v*cg_out/ch_in, v/ch_in) for lab, v in fitted.items()}
    good = np.array([(r["c_gaia"] > 0.5) & (r["c_hst"] > 0.5) for r in rows])
    cr = np.array([r["corrected"] for r in rows]); cre = np.array([r["corrected_err"] for r in rows])
    w = good & (R > 150) & (R < 275)
    mean_corr = float(np.sum(cr[w]/cre[w]**2)/np.sum(1/cre[w]**2)); mean_err = float(np.sum(1/cre[w]**2)**-0.5)
    lf = np.array([r["lf_ratio"] for r in rows]); lfe = np.array([r["lf_ratio_err"] for r in rows])
    out = dict(f_eq=f_eq, f814_eq=f8_eq, c_hst_inside_250=ch_in, c_gaia_250_300=cg_out, corrected_mean_150_275=[mean_corr, mean_err],
               fitted_ratio={k.replace("$", ""): v for k, v in fitted.items()},
               fitted_complete_sample_band={k.replace("$", ""): v for k, v in band.items()}, bins=rows)
    (ROOT/"results/plot_data/count_ratio_corrected.json").write_text(json.dumps(out, indent=1))

    fig, ax = plt.subplots(2, 1, figsize=(6.0, 7.4), sharex=True, gridspec_kw=dict(hspace=0.06, height_ratios=[1, 1.5]))
    cgv = np.array([r["c_gaia"] for r in rows]); chv = np.array([r["c_hst"] for r in rows])
    ax[0].plot(R, cgv, "o-", color="C0", label=rf"Gaia, sample $F625W<{f_eq:.2f}$ (0.06″)")
    ax[0].plot(R, [r["c_gaia_tol010"] for r in rows], "o", mfc="none", color="C0", ms=4, label="Gaia, 0.10″ match")
    for i, (a, b) in enumerate(slices):
        lab = rf"Gaia, ${'' if np.isinf(a) else f'{a:.1f}<'}F625W<{b:.1f}$"
        ax[0].plot(R, [r["c_gaia_slices"][i] for r in rows], ":", color=f"C{i+5}", lw=1, label=lab)
    ax[0].plot(R, chv, "s-", color="C3", label=rf"HST: F625W measured, of $F814W<{f8_eq:.2f}$")
    ax[0].axhline(1, color="k", lw=0.5); ax[0].set_ylim(0, 1.08); ax[0].set_ylabel("completeness")
    ax[0].legend(loc="upper left", bbox_to_anchor=(0.03, 0.9), fontsize=6.5, ncol=1, framealpha=0.9)
    strong = good & (cgv > 0.8) & (R < 300); weak = good & (cgv <= 0.8) & (R < 300); edge = good & (R > 300)
    ax[1].errorbar(R[strong], cr[strong], cre[strong], fmt="o", color="k", ms=5, label=r"direct, corrected, $c_{\rm G}>0.8$")
    ax[1].errorbar(R[weak], cr[weak], cre[weak], fmt="o", mfc="none", color="k", ms=5, label=r"direct, corrected, $0.5<c_{\rm G}<0.8$")
    ax[1].errorbar(R[edge], cr[edge], cre[edge], fmt="x", color="0.5", ms=5, label="direct, HST mosaic edge (>300″)")
    ax[1].errorbar(R, lf, lfe, fmt="s", mfc="none", color="C2", ms=4, label=rf"HST only: $N(F625W<19)/N(F625W<{f_eq:.2f})$")
    lab0 = next(iter(band)); lo_, hi_ = band[lab0]
    ax[1].axhspan(lo_, hi_, color="C4", alpha=0.2, lw=0,
                  label=rf"fitted $a_{{\rm HST}}/a_{{\rm Gaia}}$ as complete samples: {lo_:.1f}–{hi_:.1f} (no BH; BH fits within 0.2)")
    ax[1].fill_between([150, 275], mean_corr-mean_err, mean_corr+mean_err, color="k", alpha=0.25, lw=0,
                       label=rf"corrected mean, 150–275″: {mean_corr:.2f} $\pm$ {mean_err:.2f}")
    ax[1].set_ylim(5, 20); ax[1].set_ylabel(r"$N(F625W<19)\,/\,N(G<17)$, complete samples")
    ax[1].set_xscale("log"); ax[1].set_xlabel("R [arcsec]"); ax[1].set_xlim(3, 360)
    ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=7, ncol=1, frameon=False)
    for j, a in enumerate(ax):
        a.text(0.015, 0.95, "ab"[j], transform=a.transAxes, va="top", fontweight="bold")
    fig.savefig(ROOT/"plots/count_ratio_corrected.png", dpi=150, bbox_inches="tight")
    print(f"F_eq {f_eq:.2f}, F814_eq {f8_eq:.2f}, c_H(2-250) {ch_in:.3f}, c_G(250-300) {cg_out:.3f}, "
          f"corrected mean 150-275: {mean_corr:.2f} +- {mean_err:.2f}")
    print("fitted complete-sample bands", {k: tuple(round(x, 2) for x in v) for k, v in band.items()})
    for r in rows:
        print(f"{r['r_lo']:5.0f}-{r['r_hi']:<4.0f} raw {r['raw']:7.2f}  corr {r['corrected']:6.2f}+-{r['corrected_err']:.2f}  "
              f"LF {r['lf_ratio']:5.2f}  cG {r['c_gaia']:.3f} (0.1\" {r['c_gaia_tol010']:.3f}; slices "
              + " ".join(f"{x:.2f}" for x in r["c_gaia_slices"]) + f")  cH {r['c_hst']:.3f}")


if __name__ == "__main__":
    main()
