#!/usr/bin/env python3
"""Luminosity functions of the stars detected by HST (oMEGACat) and JWST (NIRCam GO 8322 F200W) in the omega Cen core,
on a common magnitude scale (2026-10-20).

The bands differ, so each telescope's magnitudes are mapped onto the other's with the median relation of the
cross-matched stars (bin/jwst_hst_depth.py catalogues and matching: 0.08" after the median offset), fitted in
0.25-mag bins where both are unsaturated (JWST F200W(3-px AB) > 18; brighter JWST magnitudes are corrupted by
saturation). Stars outside the fitted range are not mapped (left out of the mapped histograms; the range is shaded).
Rows: top in F625W (HST measured; JWST-only stars mapped), bottom in F200W (JWST measured; HST-only stars mapped).
Columns: annuli 2-4.1", 6.7-17.6", 20-60". Curves: all HST, all JWST, matched, HST-only, JWST-only, per arcsec^2 per mag.
Only the area covered by the JWST mosaic is used for both catalogues, so the comparison is like for like.
Output: plots/jwst_hst_lf_compare.png, results/plot_data/jwst_hst_lf_compare.json
Usage: PYTHONPATH=src python bin/jwst_hst_lf_compare.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"bin"), str(ROOT/"src")]
from jwst_hst_depth import ANN, crossmatch, hst_catalogue, jwst_catalogue      # noqa: E402

FN = "jw08322-o001_t001_nircam_clear-f200w_i2d.fits"
SAT = 18.0          # JWST F200W(3-px AB) brighter than this: saturation-affected
BINS = np.arange(12., 27.01, 0.25)


def median_relation(x, y, lo, hi, step=0.25, nmin=30):
    """Median of y in bins of x over [lo, hi); returns (centres, medians) where each bin has >= nmin stars, made monotone."""
    e = np.arange(lo, hi+step, step); c, m = [], []
    for a, b in zip(e[:-1], e[1:]):
        k = (x >= a) & (x < b)
        if k.sum() >= nmin:
            c.append(0.5*(a+b)); m.append(np.median(y[k]))
    c, m = np.array(c), np.maximum.accumulate(np.array(m))
    return c, m


def footprint_mask(J, pts, cell=1.0):
    """True for points in a 1"-cell that holds at least one JWST detection (the mosaic footprint, at this resolution)."""
    occ = set(zip(np.floor(J["dra"]/cell).astype(int), np.floor(J["ddec"]/cell).astype(int)))
    return np.array([(i, j) in occ for i, j in zip(np.floor(pts[0]/cell).astype(int), np.floor(pts[1]/cell).astype(int))])


def main():
    D = Path(os.path.expandvars(os.environ["OCEN_JWST"]))
    H = hst_catalogue(); J = jwst_catalogue(D/FN)
    H = {k: v[footprint_mask(J, (H["dra"], H["ddec"]))] for k, v in H.items()}
    H["R"] = np.hypot(H["dra"], H["ddec"]); J["R"] = np.hypot(J["dra"], J["ddec"])
    jh, off = crossmatch(H, J); hj, _ = crossmatch(J, H)
    g = jh >= 0; hm, jm = H["mag"][g], J["mag"][jh[g]]
    ok = jm > SAT
    cJ, mJ = median_relation(jm[ok], hm[ok], SAT, 25.)             # F200W -> F625W
    cH, mH = median_relation(hm[ok], jm[ok], 15., 26.)             # F625W -> F200W (pairs with unsaturated JWST mags)
    to625 = lambda m: np.where((m >= cJ[0]) & (m <= cJ[-1]), np.interp(m, cJ, mJ), np.nan)
    to200 = lambda m: np.where((m >= cH[0]) & (m <= cH[-1]), np.interp(m, cH, mH), np.nan)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 4, figsize=(20, 9), constrained_layout=True, gridspec_kw=dict(width_ratios=(1, 1, 1, 0.8)))
    out = dict(offset_arcsec=off.tolist(), saturation_f200w=SAT, relation_f200w_to_f625w=dict(x=cJ.tolist(), y=mJ.tolist()),
               relation_f625w_to_f200w=dict(x=cH.tolist(), y=mH.tolist()), annuli={})
    matched_h = np.zeros(len(H["mag"]), bool); matched_h[g] = True
    matched_j = hj >= 0
    for i, (lo, hi) in enumerate(ANN):
        kh = (H["R"] >= lo) & (H["R"] < hi); kj = (J["R"] >= lo) & (J["R"] < hi)
        A = np.pi*(hi*hi-lo*lo)            # annulus area; both catalogues cut to the same footprint
        lf = lambda m: np.histogram(m[np.isfinite(m)], BINS)[0]/A/0.25
        rows = {
            "F625W": dict(hst_all=lf(H["mag"][kh]), jwst_all=lf(to625(J["mag"][kj])), matched=lf(H["mag"][kh & matched_h]),
                          hst_only=lf(H["mag"][kh & ~matched_h]), jwst_only=lf(to625(J["mag"][kj & ~matched_j]))),
            "F200W": dict(hst_all=lf(to200(H["mag"][kh])), jwst_all=lf(J["mag"][kj]), matched=lf(J["mag"][kj & matched_j]),
                          hst_only=lf(to200(H["mag"][kh & ~matched_h])), jwst_only=lf(J["mag"][kj & ~matched_j]))}
        out["annuli"][f"{lo:g}-{hi:g}"] = {b: {k: v.tolist() for k, v in r.items()} | dict(n_hst=int(kh.sum()), n_jwst=int(kj.sum()))
                                           for b, r in rows.items()}
        for row, (band, r) in enumerate(rows.items()):
            x = ax[row, i]; c = BINS[:-1]
            x.step(c, r["hst_all"], where="post", color="k", lw=1.6, label=f"HST all (N={kh.sum()})")
            x.step(c, r["jwst_all"], where="post", color="C3", lw=1.6, label=f"JWST all (N={kj.sum()})")
            x.step(c, r["matched"], where="post", color="C0", lw=1, label="matched")
            x.step(c, r["hst_only"], where="post", color="0.5", lw=1, ls="--", label="HST only")
            x.step(c, r["jwst_only"], where="post", color="C1", lw=1, ls="--", label="JWST only")
            if band == "F625W":
                x.axvline(19, color="k", ls=":", lw=1); x.axvspan(12, cJ[0], color="0.92", zorder=0)
            else:
                x.axvspan(12, SAT, color="0.92", zorder=0)
            x.set(yscale="log", ylim=(3e-3, 8), xlim=(12, 27), xlabel=("HST F625W" if band == "F625W" else "JWST F200W (3-px AB)")
                  + (" [JWST mapped]" if band == "F625W" else " [HST mapped]"), ylabel="stars / arcsec^2 / mag",
                  title=f"{lo:g}-{hi:g}\" ({band} scale)")
            if i == 0:
                x.legend(fontsize=7, loc="upper left")
    for row, (x, y, xl, yl) in enumerate(((jm, hm, "JWST F200W (3-px AB)", "HST F625W"), (hm, jm, "HST F625W", "JWST F200W (3-px AB)"))):
        a = ax[row, 3]; a.scatter(x, y, s=0.5, c="0.6", rasterized=True)
        cc, mm = (cJ, mJ) if row == 0 else (cH, mH); a.plot(cc, mm, "r-", lw=2, label="median relation (fitted, JWST > 18)")
        a.set(xlabel=xl, ylabel=yl, xlim=(12, 26), ylim=(12, 26), title="mapping from matched stars"); a.legend(fontsize=7)
    fig.suptitle("omega Cen core: luminosity functions of stars detected by HST (oMEGACat) and JWST (NIRCam 8322 F200W), same footprint; "
                 "grey band: saturated / unmapped", fontsize=10)
    fig.savefig(ROOT/"plots/jwst_hst_lf_compare.png", dpi=110); print(ROOT/"plots/jwst_hst_lf_compare.png")
    (ROOT/"results/plot_data/jwst_hst_lf_compare.json").write_text(json.dumps(out, indent=1, default=float))
    for a_, r in out["annuli"].items():
        f = r["F625W"]; c = BINS[:-1]; tot = lambda y, lo_, hi_: float(np.sum(np.array(y)[(c >= lo_) & (c < hi_)])*0.25)
        print(f"{a_:9s} N_HST {r['F625W']['n_hst']:6d} N_JWST {r['F625W']['n_jwst']:6d} | F625W 19-23 per arcsec^2: HST {tot(f['hst_all'],19,23):.2f} "
              f"JWST {tot(f['jwst_all'],19,23):.2f}; 23-26: HST {tot(f['hst_all'],23,26):.2f} JWST {tot(f['jwst_all'],23,26):.2f}")


if __name__ == "__main__":
    main()
