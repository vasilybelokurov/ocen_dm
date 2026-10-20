#!/usr/bin/env python3
"""Depth of the HST (oMEGACat) and JWST (NIRCam) catalogues of the omega Cen core, by radius (2026-10-20).

Before trusting a central count drop seen in both, compare how deep each catalogue reaches where.
HST: oMEGACat VI catalogue (catalog_and_selections.fits), all stars with measured F625W (the C1 count selection has no
quality flag). JWST: detections of bin/jwst_star_counts.detect (DAOStarFinder 5 sigma, 3-px aperture AB magnitudes),
GO 8322 F200W and GO 4343 F200W, with sub-pixel sky positions from the mosaic WCS.
Panels, in three annuli (2-4.1", 6.7-17.6" = the deficit test regions, and 20-60" = the normalisation region):
  (a) HST F625W luminosity function per arcsec^2;  (b) JWST F200W luminosity function per arcsec^2;
  (c) fraction of HST stars with a JWST counterpart (within 0.08", after removing the median offset) vs F625W;
  (d) fraction of JWST stars with an HST counterpart vs F200W;
  (e) F625W vs F200W of matched stars, with the F625W < 19 cut and the JWST matched cut.
Output: plots/jwst_hst_depth.png, results/plot_data/jwst_hst_depth.json
Usage: PYTHONPATH=src python bin/jwst_hst_depth.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"bin"), str(ROOT/"src")]
from ocen_dm.cluster import OCEN_DEC, OCEN_RA          # noqa: E402

ANN = ((2.0, 4.1), (6.7, 17.6), (20., 60.))
JWST = {"8322 F200W": "jw08322-o001_t001_nircam_clear-f200w_i2d.fits", "4343 F200W": "jw04343-o002_t001_nircam_clear-f200w_i2d.fits"}
MATCH = 0.08                                            # arcsec
RMAX = 62.


def jwst_catalogue(path, fwhm=2.1):
    """Detections with sky offsets (dRA*, dDec) [arcsec] from the adopted centre and 3-px AB magnitudes."""
    from astropy.io import fits
    from astropy.wcs import WCS
    from photutils.aperture import CircularAperture, aperture_photometry
    from photutils.background import Background2D, MedianBackground
    from photutils.detection import DAOStarFinder
    with fits.open(path, memmap=True) as h:
        w = WCS(h["SCI"].header); pixsr = h["SCI"].header["PIXAR_SR"]
        x0, y0 = w.world_to_pixel_values(OCEN_RA, OCEN_DEC)
        sc = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600.); n = int(RMAX/sc)+2
        ny, nx = h["SCI"].shape
        xa, xb, ya, yb = max(int(x0)-n, 0), min(int(x0)+n, nx), max(int(y0)-n, 0), min(int(y0)+n, ny)   # clipped at the mosaic edge
        img = np.array(h["SCI"].section[ya:yb, xa:xb], float); wht = np.array(h["WHT"].section[ya:yb, xa:xb])
    img[(wht <= 0) | (img == 0)] = np.nan; mask = ~np.isfinite(img)
    bkg = Background2D(np.where(mask, 0., img), 64, filter_size=3, mask=mask, bkg_estimator=MedianBackground(), exclude_percentile=80.)
    sub = np.where(mask, 0., img-bkg.background)
    src = DAOStarFinder(threshold=5., fwhm=fwhm, exclude_border=True)(sub/bkg.background_rms, mask=mask)
    xy = np.column_stack((src["xcentroid"], src["ycentroid"]))
    flux = aperture_photometry(sub, CircularAperture(xy, r=3.))["aperture_sum"].data*pixsr*1e6
    ok = flux > 0
    ra, dec = w.pixel_to_world_values(xy[ok, 0]+xa, xy[ok, 1]+ya)
    return dict(dra=(ra-OCEN_RA)*np.cos(np.deg2rad(OCEN_DEC))*3600., ddec=(dec-OCEN_DEC)*3600., mag=-2.5*np.log10(flux[ok]/3631.))


def hst_catalogue():
    from astropy.table import Table
    from ocen_dm.paths import raw_dir
    t = Table.read(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")
    m = np.isfinite(np.asarray(t["f625w"], float))
    ra, dec = np.asarray(t["RA"], float)[m], np.asarray(t["DEC"], float)[m]
    dra = (ra-OCEN_RA)*np.cos(np.deg2rad(OCEN_DEC))*3600.; ddec = (dec-OCEN_DEC)*3600.
    k = np.hypot(dra, ddec) < RMAX
    return dict(dra=dra[k], ddec=ddec[k], mag=np.asarray(t["f625w"], float)[m][k])


def crossmatch(a, b):
    """Nearest-neighbour match a -> b after removing the median offset of bright pairs; returns index (or -1) and offset."""
    from scipy.spatial import cKDTree
    tb = cKDTree(np.column_stack((b["dra"], b["ddec"])))
    pa = np.column_stack((a["dra"], a["ddec"]))
    d, j = tb.query(pa, distance_upper_bound=0.3)
    good = np.isfinite(d)
    off = np.median(np.column_stack((b["dra"][j[good]], b["ddec"][j[good]]))-pa[good], axis=0)
    d, j = tb.query(pa+off, distance_upper_bound=MATCH)
    j[~np.isfinite(d)] = -1
    return j, off


def main():
    D = Path(os.path.expandvars(os.environ["OCEN_JWST"]))
    H = hst_catalogue(); H["R"] = np.hypot(H["dra"], H["ddec"])
    cut = json.loads((ROOT/"results/plot_data/jwst_star_counts.json").read_text())["jwst"]
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5), constrained_layout=True); ax = ax.ravel()
    out = dict(annuli=ANN, match_arcsec=MATCH, jwst={})
    hb = np.arange(12., 27.01, 0.25); jb = np.arange(12., 27.01, 0.25)
    for i, (lo, hi) in enumerate(ANN):
        A = np.pi*(hi*hi-lo*lo); k = (H["R"] >= lo) & (H["R"] < hi)
        n, _ = np.histogram(H["mag"][k], hb); ax[0].step(hb[:-1], n/A/0.25, where="post", color=f"C{i}", label=f"{lo:g}-{hi:g}\" (N={k.sum()})")
    ax[0].axvline(19, color="k", ls=":", lw=1); ax[0].set(yscale="log", xlabel="HST F625W", ylabel="stars / arcsec^2 / mag",
                                                          title="(a) HST oMEGACat luminosity function", xlim=(12, 27)); ax[0].legend(fontsize=8)
    for s, (lab, fn) in enumerate(JWST.items()):
        J = jwst_catalogue(D/fn); J["R"] = np.hypot(J["dra"], J["ddec"]); mc = cut[lab]["matched_cut"]
        jh, off = crossmatch(H, J); hj, _ = crossmatch(J, H)
        rec = dict(n_detected=len(J["mag"]), offset_arcsec=off.tolist(), matched_cut=mc, annuli={})
        ls = "-" if s == 0 else "--"
        for i, (lo, hi) in enumerate(ANN):
            A = np.pi*(hi*hi-lo*lo); kj = (J["R"] >= lo) & (J["R"] < hi); kh = (H["R"] >= lo) & (H["R"] < hi)
            n, _ = np.histogram(J["mag"][kj], jb)
            ax[1].step(jb[:-1], n/A/0.25, where="post", color=f"C{i}", ls=ls, label=f"{lab} {lo:g}-{hi:g}\" (N={kj.sum()})")
            # matched fractions
            fh = [np.mean(jh[kh & (H["mag"] >= a) & (H["mag"] < a+0.5)] >= 0) if (kh & (H["mag"] >= a) & (H["mag"] < a+0.5)).sum() >= 5 else np.nan
                  for a in hb[:-1:2]]
            fj = [np.mean(hj[kj & (J["mag"] >= a) & (J["mag"] < a+0.5)] >= 0) if (kj & (J["mag"] >= a) & (J["mag"] < a+0.5)).sum() >= 5 else np.nan
                  for a in jb[:-1:2]]
            ax[2].plot(hb[:-1:2]+0.25, fh, "o"+ls, ms=3, color=f"C{i}", label=f"{lab}: {lo:g}-{hi:g}\"")
            ax[3].plot(jb[:-1:2]+0.25, fj, "o"+ls, ms=3, color=f"C{i}", label=f"{lab}: {lo:g}-{hi:g}\"")
            nh19 = int((kh & (H["mag"] < 19)).sum()); nj = int((kj & (J["mag"] < mc)).sum())
            rec["annuli"][f"{lo:g}-{hi:g}"] = dict(n_hst_f625w_lt19=nh19, n_jwst_lt_matched=nj,
                                                   frac_hst19_in_jwst=float(np.mean(jh[kh & (H["mag"] < 19)] >= 0)) if nh19 else None,
                                                   frac_jwst_cut_in_hst=float(np.mean(hj[kj & (J["mag"] < mc)] >= 0)) if nj else None)
        if s == 0:
            g = jh >= 0
            ax[4].scatter(H["mag"][g], J["mag"][jh[g]], s=1, c=H["R"][g], cmap="viridis", vmin=0, vmax=60, rasterized=True)
            ax[4].axvline(19, color="k", ls=":"); ax[4].axhline(mc, color="r", ls=":")
            ax[4].set(xlabel="HST F625W", ylabel=f"JWST {lab} (3-px AB)", xlim=(12, 26), ylim=(12, 25),
                      title="(e) matched stars (colour: R 0-60\"); dotted: F625W<19, JWST matched cut")
        ax[1].axvline(mc, color="k" if s == 0 else "0.5", ls=":", lw=1)
        out["jwst"][lab] = rec
        print(lab, "offset", np.round(off, 3), "| annuli:", rec["annuli"])
    ax[1].set(yscale="log", xlabel="JWST F200W (3-px AB)", ylabel="stars / arcsec^2 / mag", title="(b) JWST luminosity function (solid 8322, dashed 4343)",
              xlim=(12, 27)); ax[1].legend(fontsize=6)
    for x, t, xl in ((ax[2], "(c) fraction of HST stars with a JWST match", "HST F625W"), (ax[3], "(d) fraction of JWST stars with an HST match", "JWST F200W (3-px AB)")):
        x.set(xlabel=xl, ylabel="matched fraction", ylim=(0, 1.05), xlim=(12, 27), title=t); x.legend(fontsize=6); x.grid(alpha=0.3)
    ax[2].axvline(19, color="k", ls=":")
    # (f) the two-catalogue counts at the cut vs radius
    e = np.array(json.loads((ROOT/"results/plot_data/jwst_star_counts.json").read_text())["edges_arcsec"]); r = 0.5*(e[1:]+e[:-1])
    A = np.pi*(e[1:]**2-e[:-1]**2)
    nh, _ = np.histogram(H["R"][H["mag"] < 19], e); ax[5].errorbar(r, nh/A, np.sqrt(nh)/A, fmt="ks-", ms=3, label="HST F625W<19 (raw, all)")
    for lab in JWST:
        p = cut[lab]["profiles"]["matched/adopted"]; n = np.array(p["counts"])
        ax[5].errorbar(r, n/A, np.sqrt(n)/A, fmt="o-", ms=3, label=f"JWST {lab} < {cut[lab]['matched_cut']:.2f}")
    ax[5].set(xscale="log", xlabel="R [arcsec]", ylabel="stars / arcsec^2", title="(f) raw densities at the cuts (no completeness)")
    ax[5].axvspan(2.0, 3.2, color="0.9", zorder=0); ax[5].legend(fontsize=7)
    fig.suptitle("omega Cen core: depth of HST oMEGACat vs JWST NIRCam detections, by radius", fontsize=11)
    fig.savefig(ROOT/"plots/jwst_hst_depth.png", dpi=110); print(ROOT/"plots/jwst_hst_depth.png")
    (ROOT/"results/plot_data/jwst_hst_depth.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
