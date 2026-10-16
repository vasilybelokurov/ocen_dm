#!/usr/bin/env python3
"""Count products with per-bin detection completeness, and the HST/Gaia amplitude tie (2026-10-16).

Inputs: the schema-2 count products (bin/build_count_profiles.py), the Gaia completeness against HST
(bin/check_gaia_completeness_lp.py functions), and the complete-sample count ratio measured by
bin/check_count_ratio_corrected.py (results/plot_data/count_ratio_corrected.json; run that script first).

- HST (ocen_counts_hst_f625w19 -> ocen_counts_hst_f625w19_c): c_H per bin = fraction of oMEGACat stars with
  F814W < F814_eq (the F814W equivalent of F625W = 19) that have a measured F625W, in the product's own bins
  and pixel centre; inner bins are pooled until each group holds >= 500 such stars (Poisson noise).
- Gaia (ocen_counts_gaia_g17 -> ocen_counts_gaia_g17_c): c_G = fraction of HST stars brighter than the HST
  equivalent of G = 17 present in the Vasiliev & Baumgardt catalogue. Bin 1 (300-350"): oMEGACat and the HST
  Large Programme fields combined. Every later bin: one value, the Large Programme fields combined over 350-1112"
  (the measurements per bin are consistent with a constant; beyond 1112" no HST field exists, the same value is
  assumed).
- Tie (ocen_counts_tie_c.json): a_Gaia = a_HST / r for complete samples, r ~ N(ratio, sigma) from the corrected
  ratio at 150-275".

Usage: PYTHONPATH=src python bin/build_count_completeness.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np
from astropy.io import fits

from check_gaia_completeness_lp import LP_FIELDS, found_flags, gaia, load_lp, load_omegacat
from ocen_dm.kinematics.counts import CountProfile, CountTie
from ocen_dm.paths import processed_dir, raw_dir

KIN = processed_dir()/"kinematics"
MIN_PROXY = 500


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hst_completeness(profile):
    h = fits.open(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")[1].data
    x = -0.04*(np.asarray(h["x"], float)-15000.); y = 0.04*(np.asarray(h["y"], float)-15000.)
    f = np.asarray(h["f625w"], float); f8 = np.asarray(h["f814w"], float)
    R = np.hypot(x, y)
    k = np.isfinite(f) & np.isfinite(f8) & (f > 18.5) & (f < 19.)
    f8_eq = 19.-float(np.median(f[k]-f8[k]))
    proxy = np.isfinite(f8) & (f8 < f8_eq)
    n = [int(np.sum(proxy & (R >= lo) & (R < hi))) for lo, hi in zip(profile.r_lower, profile.r_upper)]
    # pool bins from the centre outwards until each group holds >= MIN_PROXY proxy stars (inner bins hold 6-500)
    groups, cur = [], []
    for i in range(profile.n):
        cur.append(i)
        if sum(n[j] for j in cur) >= MIN_PROXY or i == profile.n-1:
            groups.append(cur); cur = []
    c = np.empty(profile.n)
    for g in groups:
        s = proxy & (R >= profile.r_lower[g[0]]) & (R < profile.r_upper[g[-1]])
        c[g] = float(np.isfinite(f[s]).mean())
    return c, dict(method=f"fraction of stars with F814W < {f8_eq:.3f} that have a measured F625W; bins pooled from the "
                          f"centre until each group holds >= {MIN_PROXY} such stars", n_proxy=n,
                   groups=[[int(g[0]), int(g[-1])] for g in groups], f814w_equivalent=f8_eq)


def gaia_completeness(profile):
    gx, gy, G = gaia()
    om = found_flags(*load_omegacat(), gx, gy, G, (16.5, 17.5))
    lp = [found_flags(*load_lp(sub, ic), gx, gy, G, (16., 17.5)) for sub, ic in LP_FIELDS.values()]
    R = np.concatenate([v[0] for v in lp]); F = np.concatenate([v[1] for v in lp])
    lo1, hi1 = profile.r_lower[0], profile.r_upper[0]
    k1o = (om[0] >= lo1) & (om[0] < hi1); k1l = (R >= lo1) & (R < hi1)
    n1, f1 = int(k1o.sum()+k1l.sum()), int(om[1][k1o].sum()+F[k1l].sum())
    ko = (R >= hi1)
    n2, f2 = int(ko.sum()), int(F[ko].sum())
    c = np.full(profile.n, f2/n2); c[0] = f1/n1
    meta = dict(method="fraction of HST stars brighter than the HST equivalent of G = 17 present in the VB21 catalogue "
                       "(0.1 arcsec match after the epoch offset)",
                bin1=dict(r=[float(lo1), float(hi1)], found=f1, n=n1, c=f1/n1, err=float(np.sqrt(f1/n1*(1-f1/n1)/n1)),
                          sources="oMEGACat + HST Large Programme fields"),
                outer=dict(r=[float(hi1), float(R.max())], found=f2, n=n2, c=f2/n2, err=float(np.sqrt(f2/n2*(1-f2/n2)/n2)),
                           sources="HST Large Programme fields F1-F5", note="assumed constant beyond the last HST field"))
    return c, meta


def write(profile, completeness, meta, name):
    row = profile.to_dict()
    row.update(name=name, completeness=completeness.tolist(),
               note=(row.get("note") or "")+" Per-bin detection completeness added 2026-10-16 (bin/build_count_completeness.py).",
               meta=dict(row["meta"], parent_product=profile.name, completeness=meta,
                         built_utc=datetime.now(timezone.utc).isoformat()))
    CountProfile.from_dict(row)    # validate
    path = KIN/f"{name}.json"
    path.write_text(json.dumps(row, indent=1)+"\n")
    return path


def main():
    hst = CountProfile.from_dict(json.loads((KIN/"ocen_counts_hst_f625w19.json").read_text()))
    gai = CountProfile.from_dict(json.loads((KIN/"ocen_counts_gaia_g17.json").read_text()))
    ch, mh = hst_completeness(hst)
    cg, mg = gaia_completeness(gai)
    for k, v in (("parent_sha256", sha256(KIN/"ocen_counts_hst_f625w19.json")),):
        mh[k] = v
    mg["parent_sha256"] = sha256(KIN/"ocen_counts_gaia_g17.json")
    ph = write(hst, ch, mh, "ocen_counts_hst_f625w19_c")
    pg = write(gai, cg, mg, "ocen_counts_gaia_g17_c")
    src = ROOT/"results/plot_data/count_ratio_corrected.json"
    ratio, sigma = json.loads(src.read_text())["corrected_mean_150_275"]
    tie = CountTie("ocen_counts_hst_f625w19_c", "ocen_counts_gaia_g17_c", round(ratio, 4), round(sigma, 4),
                   note="N(F625W<19)/N(G<17) for complete samples, inverse-variance mean at 150-275 arcsec of the "
                        "completeness-corrected ratio (bin/check_count_ratio_corrected.py, 2026-10-16)")
    (KIN/"ocen_counts_tie_c.json").write_text(json.dumps(tie.to_dict(), indent=1)+"\n")
    print(ph.name, "c_H:", np.round(ch, 3).tolist())
    print(pg.name, "c_G:", np.round(cg, 3).tolist(), mg["bin1"], mg["outer"])
    print("tie:", tie.to_dict())


if __name__ == "__main__":
    main()
