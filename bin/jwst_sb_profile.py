#!/usr/bin/env python3
"""Central surface-brightness profile of omega Cen from public JWST NIRCam mosaics (2026-10-20).

Question: does JWST confirm the drop of the HST star counts in the innermost bins (2-3.2", 6 and 12 stars against
10.8 and 17.5 expected for a flat core)? The counts' drop moved with the adopted centre (bin/check_count_centre.py).

Data: level-3 NIRCam i2d mosaics (SCI in MJy/sr) in $OCEN_JWST: GO 8322 (Haeberle; 2025-07-30; F115W, F200W) and
GO 4343 (Kargaltsev; 2024-07-10; F200W). Pixels with no data (NaN or zero weight; saturated star cores drop out of the
resampled mosaic this way) are excluded; the excluded fraction per annulus is reported.

Profile: annuli around the adopted centre (ocen_dm.cluster: Anderson & van der Marel 2010 / oMEGACat), 0.5" wide inside
10" and log-spaced outside. Three statistics per annulus, each computed in 16 azimuthal sectors:
  mean    - mean of the valid pixels: integrated light, dominated by the few bright giants (shot noise);
  clipped - mean after iterative 3-sigma (MAD) clipping of bright pixels: the light of the faint, unresolved stars;
  median  - median pixel: the same, most robust.
The value is the mean over sectors and the error the sector scatter / sqrt(16) (azimuthal stochasticity). The same
profiles at the oMEGACat VI kinematic centre (dRA, dDec = -1.31", +1.44"; Haeberle et al. 2025) test the centre choice.
Comparison: completeness-corrected HST star counts (ocen_counts_hst_f625w19_c) and the J-counts MGE tracer, all
normalised to their mean over 10-20".
Output: results/plot_data/jwst_sb_profile.json, plots/jwst_sb_profile.png
Usage: PYTHONPATH=src python bin/jwst_sb_profile.py [--rmax 60]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from ocen_dm.cluster import OCEN_DEC, OCEN_RA          # noqa: E402

MOSAICS = {"8322 F115W": "jw08322-o001_t001_nircam_clear-f115w_i2d.fits",
           "8322 F200W": "jw08322-o001_t001_nircam_clear-f200w_i2d.fits",
           "4343 F200W": "jw04343-o002_t001_nircam_clear-f200w_i2d.fits"}
KINEMATIC = (-1.31, 1.44)            # oMEGACat VI kinematic centre offset (dRA*, dDec), arcsec
NSECT = 16
MJYSR_TO_JY_ARCSEC2 = 1e6/((180/np.pi*3600)**2)


def edges(rmax):
    inner = np.arange(0., 10.01, 0.5)
    outer = np.geomspace(10., rmax, 15)[1:]
    return np.r_[inner, outer]


def cutout(path, rmax):
    """SCI cutout [MJy/sr] (NaN where no data), and the pixel offsets (dRA*, dDec) [arcsec] from the adopted centre."""
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=True) as h:
        w = WCS(h["SCI"].header); sci = h["SCI"].section
        wht = h["WHT"].section if "WHT" in h else None
        x0, y0 = w.world_to_pixel_values(OCEN_RA, OCEN_DEC)
        scale = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600.)
        n = int(np.ceil(rmax/scale))+2; ny, nx = h["SCI"].shape
        xa, xb = max(int(x0)-n, 0), min(int(x0)+n, nx); ya, yb = max(int(y0)-n, 0), min(int(y0)+n, ny)
        img = np.array(sci[ya:yb, xa:xb], float)
        if wht is not None:
            img[np.array(wht[ya:yb, xa:xb]) <= 0] = np.nan
        img[img == 0] = np.nan
        yy, xx = np.mgrid[ya:yb, xa:xb]
        ra, dec = w.pixel_to_world_values(xx, yy)
    dra = (ra-OCEN_RA)*np.cos(np.deg2rad(OCEN_DEC))*3600.; ddec = (dec-OCEN_DEC)*3600.
    return img, dra, ddec, scale, dict(x0=float(x0), y0=float(y0), date=None)


def clipped_mean(v, n_iter=5):
    v = v[np.isfinite(v)]
    for _ in range(n_iter):
        if v.size < 5:
            break
        m = np.median(v); s = 1.4826*np.median(np.abs(v-m))
        k = v < m+3*s
        if k.all():
            break
        v = v[k]
    return float(np.mean(v)) if v.size else np.nan


def profile(img, dra, ddec, e, centre=(0., 0.)):
    x, y = dra-centre[0], ddec-centre[1]; R = np.hypot(x, y); th = np.arctan2(y, x)
    sect = np.floor((th+np.pi)/(2*np.pi)*NSECT).astype(int) % NSECT
    out = {k: [] for k in ("r_lo", "r_hi", "r_mid", "n_pix", "valid_frac")}
    for s in ("mean", "clipped", "median"):
        out[s] = []; out[s+"_err"] = []
    for a, b in zip(e[:-1], e[1:]):
        k = (R >= a) & (R < b); v = img[k]; ok = np.isfinite(v)
        out["r_lo"].append(a); out["r_hi"].append(b); out["r_mid"].append(float(np.sqrt(0.5*(a*a+b*b))))
        out["n_pix"].append(int(k.sum())); out["valid_frac"].append(float(ok.mean()) if k.any() else 0.)
        for s, f in (("mean", lambda q: float(np.nanmean(q)) if np.isfinite(q).any() else np.nan),
                     ("clipped", clipped_mean), ("median", lambda q: float(np.nanmedian(q)) if np.isfinite(q).any() else np.nan)):
            per = np.array([f(img[k & (sect == j)]) for j in range(NSECT)])
            per = per[np.isfinite(per)]
            if per.size >= NSECT//2 and ok.mean() > 0.5:
                out[s].append(float(per.mean())); out[s+"_err"].append(float(per.std(ddof=1)/np.sqrt(per.size)))
            else:
                out[s].append(np.nan); out[s+"_err"].append(np.nan)
    return out


def hst_counts():
    from ocen_dm.kinematics.counts import CountProfile
    from ocen_dm.paths import processed_dir
    c = CountProfile.from_dict(json.loads((processed_dir()/"kinematics/ocen_counts_hst_f625w19_c.json").read_text()))
    comp = np.ones(c.n) if c.completeness is None else np.asarray(c.completeness)
    dens = np.asarray(c.counts)/(np.asarray(c.area_arcsec2)*comp); err = np.sqrt(np.maximum(c.counts, 1))/(np.asarray(c.area_arcsec2)*comp)
    return dict(r_lo=np.asarray(c.r_lower), r_hi=np.asarray(c.r_upper), counts=np.asarray(c.counts), dens=dens, err=err)


def mge_counts(R):
    j = json.loads((ROOT/"data/processed/kinematics/tracer_mge_counts_c.json").read_text())
    s, f = np.array(j["sigmas_arcsec"]), np.array(j["fractions"])
    return (f/(2*np.pi*s**2)*np.exp(-np.asarray(R)[:, None]**2/(2*s**2))).sum(1)


def norm(r, y, lo=10., hi=20.):
    k = (np.asarray(r) >= lo) & (np.asarray(r) < hi) & np.isfinite(y)
    return np.asarray(y)/np.mean(np.asarray(y)[k])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--rmax", type=float, default=60.); a = ap.parse_args()
    D = Path(os.path.expandvars(os.environ.get("OCEN_JWST", "")))
    if not D.is_dir():
        raise SystemExit("OCEN_JWST not set or missing (see .env)")
    e = edges(a.rmax); res = {}
    for lab, fn in MOSAICS.items():
        if not (D/fn).exists():
            print(f"missing {fn}; skipped"); continue
        img, dra, ddec, scale, info = cutout(D/fn, a.rmax)
        res[lab] = dict(file=fn, pixel_arcsec=scale, adopted=profile(img, dra, ddec, e), kinematic=profile(img, dra, ddec, e, KINEMATIC))
        p = res[lab]["adopted"]
        print(f"{lab}: pixel {scale:.4f}\" | valid fraction inside 3\": {np.mean(p['valid_frac'][:6]):.2f}")
    h = hst_counts()
    out = dict(centre_deg=(OCEN_RA, OCEN_DEC), kinematic_offset_arcsec=KINEMATIC, edges_arcsec=e.tolist(), profiles=res,
               hst_counts={k: np.asarray(v).tolist() for k, v in h.items()},
               units="SCI MJy/sr; 1 MJy/sr = %.4e Jy/arcsec^2" % MJYSR_TO_JY_ARCSEC2)
    (ROOT/"results/plot_data/jwst_sb_profile.json").write_text(json.dumps(out, indent=1, default=float))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    hr = 0.5*(h["r_lo"]+h["r_hi"]); hn = norm(hr, h["dens"]); he = h["err"]/h["dens"]*hn
    Rm = np.geomspace(0.3, a.rmax, 200); mg = norm(Rm, mge_counts(Rm))
    for col, stat in enumerate(("clipped", "mean")):
        for row, (lo, hi, xs) in enumerate(((0, a.rmax, "log"), (0, 12, "linear"))):
            x = ax[row, col]
            x.errorbar(hr, hn, he, fmt="ks", ms=4, capsize=2, label="HST counts F625W<19 (completeness-corr.)", zorder=5)
            x.plot(Rm, mg, "k:", lw=1, label="J-counts MGE")
            for i, (lab, r) in enumerate(res.items()):
                p = r["adopted"]; y = np.array(p[stat]); yn = norm(p["r_mid"], y); en = np.array(p[stat+"_err"])/y*yn
                x.errorbar(np.array(p["r_mid"])*(1+0.01*i), yn, en, fmt="o-", ms=3, capsize=2, color=f"C{i}", label=f"JWST {lab}")
                q = r["kinematic"]; yk = norm(q["r_mid"], np.array(q[stat]))
                x.plot(q["r_mid"], yk, "--", color=f"C{i}", lw=0.8, label=f"  kinematic centre" if row == 0 and col == 0 else None)
            x.axvspan(2.0, 3.2, color="0.9", zorder=0); x.set_xscale(xs)
            x.set(xlim=(0.25 if xs == "log" else 0, hi), ylim=(0.5, 1.6), xlabel="R [arcsec]",
                  ylabel="Sigma / mean(10-20\")", title=f"{stat} pixel statistic" + (" (inner 12\")" if row else ""))
            x.axhline(1, color="0.6", lw=0.5)
            if row == 0 and col == 0:
                x.legend(fontsize=7, loc="lower left")
    fig.suptitle("omega Cen central profile: JWST NIRCam (16-sector errors) vs HST star counts; grey band = HST counts' 2-3.2\" bins", fontsize=10)
    fig.savefig(ROOT/"plots/jwst_sb_profile.png", dpi=110); print(ROOT/"plots/jwst_sb_profile.png")
    for lab, r in res.items():
        p = r["adopted"]; y = norm(p["r_mid"], np.array(p["clipped"])); m = norm(p["r_mid"], np.array(p["mean"]))
        print(f"{lab}:  R   clipped   mean   valid")
        for i in range(10):
            print(f"   {p['r_lo'][i]:4.1f}-{p['r_hi'][i]:4.1f}  {y[i]:5.3f}  {m[i]:5.3f}  {p['valid_frac'][i]:.2f}")


if __name__ == "__main__":
    main()
