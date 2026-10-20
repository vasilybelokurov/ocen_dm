#!/usr/bin/env python3
"""Central star-count profile of omega Cen from JWST NIRCam, in the bins of the HST counts (2026-10-20).

Pixel statistics of the core are dominated by a few bright giants (bin/jwst_sb_profile.py), so the dip question is
answered with resolved counts instead. Stars are detected with DAOStarFinder on the background-subtracted significance
image (Background2D, 64-px boxes, median; threshold 5 sigma; FWHM 2.1 px in F200W, 1.4 px in F115W, SW pixel 0.031"),
measured with a 3-px aperture (AB magnitudes from the MJy/sr image and PIXAR_SR; no aperture correction, so the
magnitudes are offset by a constant), and counted in the HST count bins (ocen_counts_hst_f625w19_c: 2-250").
Magnitude cuts: the "matched" cut gives the same density as the completeness-corrected HST F625W < 19 counts over
20-60"; brighter and fainter cuts (+-1 mag) show whether the inner shape depends on depth (crowding incompleteness
grows faintwards and inwards). No artificial-star completeness: a central drop that grows with depth and is absent for
the bright cut points to crowding, not to the stellar density.
Centres: adopted (ocen_dm.cluster) and the oMEGACat VI kinematic centre.
Output: results/plot_data/jwst_star_counts.json, plots/jwst_star_counts.png
Usage: PYTHONPATH=src python bin/jwst_star_counts.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"bin"), str(ROOT/"src")]
from jwst_sb_profile import KINEMATIC, cutout, hst_counts      # noqa: E402

MOSAICS = {"8322 F200W": ("jw08322-o001_t001_nircam_clear-f200w_i2d.fits", 2.1),
           "8322 F115W": ("jw08322-o001_t001_nircam_clear-f115w_i2d.fits", 1.4),
           "4343 F200W": ("jw04343-o002_t001_nircam_clear-f200w_i2d.fits", 2.1)}
RMAX = 62.


def detect(path, fwhm):
    from astropy.io import fits
    from photutils.aperture import CircularAperture, aperture_photometry
    from photutils.background import Background2D, MedianBackground
    from photutils.detection import DAOStarFinder
    img, dra, ddec, scale, _ = cutout(path, RMAX)
    pixsr = fits.getheader(path, "SCI")["PIXAR_SR"]
    mask = ~np.isfinite(img)
    bkg = Background2D(np.where(mask, 0., img), 64, filter_size=3, mask=mask, bkg_estimator=MedianBackground(), exclude_percentile=80.)
    sub = np.where(mask, 0., img-bkg.background)
    src = DAOStarFinder(threshold=5., fwhm=fwhm, exclude_border=True)(sub/bkg.background_rms, mask=mask)
    xy = np.column_stack((src["xcentroid"], src["ycentroid"]))
    flux = aperture_photometry(sub, CircularAperture(xy, r=3.))["aperture_sum"].data*pixsr*1e6      # Jy
    ok = flux > 0
    ix, iy = np.round(xy[ok, 0]).astype(int), np.round(xy[ok, 1]).astype(int)
    return dict(dra=dra[iy, ix], ddec=ddec[iy, ix], mag=-2.5*np.log10(flux[ok]/3631.), n=int(ok.sum()))


def density(s, edges, mcut, centre=(0., 0.)):
    R = np.hypot(s["dra"]-centre[0], s["ddec"]-centre[1]); k = s["mag"] < mcut
    n, _ = np.histogram(R[k], edges); area = np.pi*(edges[1:]**2-edges[:-1]**2)
    return n, n/area


def main():
    D = Path(os.path.expandvars(os.environ["OCEN_JWST"]))
    h = hst_counts(); k = (h["r_hi"] <= RMAX); e = np.r_[h["r_lo"][k], h["r_hi"][k][-1]]
    hr = 0.5*(e[1:]+e[:-1]); hd, he = h["dens"][k], h["err"][k]; ref = (hr >= 20) & (hr < 60)
    out = dict(edges_arcsec=e.tolist(), hst=dict(dens=hd.tolist(), err=he.tolist(), counts=h["counts"][k].tolist()), jwst={})
    for lab, (fn, fwhm) in MOSAICS.items():
        s = detect(D/fn, fwhm)
        grid = np.arange(14., 24., 0.05)                                   # matched cut: equal density over 20-60"
        dens20 = np.array([density(s, e, m)[1][ref].mean() for m in grid])
        mm = float(np.interp(hd[ref].mean(), dens20, grid))
        cuts = {"matched-1": mm-1, "matched": mm, "matched+1": mm+1}
        r = dict(file=fn, n_detected=s["n"], matched_cut=mm, profiles={})
        for cl, mc in cuts.items():
            for cen, c0 in (("adopted", (0., 0.)), ("kinematic", KINEMATIC)):
                n, dd = density(s, e, mc, c0)
                r["profiles"][f"{cl}/{cen}"] = dict(mcut=mc, counts=n.tolist(), dens=dd.tolist())
        out["jwst"][lab] = r
        p = r["profiles"]["matched/adopted"]
        print(f"{lab}: {s['n']} detections; matched cut m_AB(3px) < {mm:.2f}; inner bins (R, N_JWST, N_HST):",
              [(round(a, 1), b, c) for a, b, c in zip(e[:6], p["counts"][:6], h["counts"][k][:6])])
    (ROOT/"results/plot_data/jwst_star_counts.json").write_text(json.dumps(out, indent=1, default=float))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    nrm = lambda y: np.asarray(y)/np.mean(np.asarray(y)[ref])
    for i, (lab, r) in enumerate(out["jwst"].items()):
        x = ax[i]
        x.errorbar(hr, nrm(hd), he/hd*nrm(hd), fmt="ks", ms=4, capsize=2, label="HST F625W<19 (completeness-corrected)", zorder=5)
        for j, cl in enumerate(("matched-1", "matched", "matched+1")):
            p = r["profiles"][f"{cl}/adopted"]; n = np.array(p["counts"]); y = nrm(p["dens"])
            x.errorbar(hr*(1+0.02*(j-1)), y, y/np.sqrt(np.maximum(n, 1)), fmt="o-", ms=3, capsize=2, color=f"C{j}",
                       label=f"JWST m<{p['mcut']:.2f} ({cl})")
            q = r["profiles"][f"{cl}/kinematic"]; x.plot(hr, nrm(q["dens"]), "--", color=f"C{j}", lw=0.8,
                                                         label="  kinematic centre" if j == 1 else None)
        x.axvspan(2.0, 3.2, color="0.9", zorder=0); x.axhline(1, color="0.6", lw=0.5)
        x.set(xscale="log", xlabel="R [arcsec]", ylabel="density / mean(20-60\")", ylim=(0.4, 1.6), title=f"JWST {lab}: star counts")
        x.legend(fontsize=7, loc="lower left")
    fig.suptitle("omega Cen central star counts: JWST NIRCam detections (no artificial-star completeness) vs HST; grey = 2-3.2\"", fontsize=10)
    fig.savefig(ROOT/"plots/jwst_star_counts.png", dpi=110); print(ROOT/"plots/jwst_star_counts.png")


if __name__ == "__main__":
    main()
