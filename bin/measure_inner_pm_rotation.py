#!/usr/bin/env python3
"""Step 1 of the rotation plan (2026-11-13): inner plane-of-sky rotation from Gaia-minus-HST proper motions.

HST (oMEGACat) proper motions are locally corrected: the local mean motion, including rotation, is removed. Gaia PMs
are absolute. For the same star, d = PM_Gaia - PM_HST = systemic motion + local mean field (rotation) + noise, with the
star's own random motion cancelled. This is the method of oMEGACat II (Haeberle et al. 2024, arXiv:2404.03722,
Sec. VIII, Table 10), which used Gaia FPR; here Gaia EDR3 (Vasiliev & Baumgardt 2021 members), which is sparse and
noisy inside ~200".

Steps: positional match (same procedure as ocen_dm.selection.hst_gaia_match: epoch offset from a 0.4" pass, then
0.06"); d per star; systemic = weighted median of d over all matches; tangential/radial components of d about AvdM10
(tangential positive counter-clockwise on the sky, i.e. from North towards East); per 30" ring a Gaussian ML mean with
free intrinsic scatter (the EDR3 crowding scatter is not described by its errors) and the quoted errors; bootstrap
errors. Compared with Table 10 (median tangential PM per ring).
Usage: PYTHONPATH=src python bin/measure_inner_pm_rotation.py [--all-hst]
  --all-hst: drop the HST selection_hq_astrometry_and_membership cut (keeps Gaia membership > 0.9); more inner stars.
Output: plots/inner_pm_rotation_20261113[_allhst].png, results/plot_data/inner_pm_rotation_20261113[_allhst].json
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
from astropy.table import Table
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from ocen_dm.cluster import OCEN_DEC, OCEN_DISTANCE_KPC, OCEN_RA
from ocen_dm.paths import processed_dir, raw_dir

KMS = 4.74047*OCEN_DISTANCE_KPC
# oMEGACat II Table 10 (Haeberle et al. 2024): ring edges ["], median radius ["], <mu_T> [mas/yr] and error
T10 = np.array([[0, 30, 20.66, .0548, .0324], [30, 60, 47.6, .086, .0181], [60, 90, 76.65, .1478, .0131],
                [90, 120, 105.91, .2088, .0105], [120, 150, 135.99, .2556, .0093], [150, 180, 165.6, .2562, .0083],
                [180, 210, 195.11, .2678, .008], [210, 240, 225.38, .2648, .0078], [240, 270, 254.06, .2797, .0083],
                [270, 300, 282.86, .2619, .0103], [300, 330, 309.92, .2658, .0176], [330, 360, 340.91, .2455, .0421]])


def match():
    hst = Table.read(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")
    cosd = np.cos(np.radians(OCEN_DEC))
    hx = (np.asarray(hst["RA"], float)-OCEN_RA)*cosd*3600; hy = (np.asarray(hst["DEC"], float)-OCEN_DEC)*3600
    ha, hd = np.asarray(hst["pmra_corrected"], float), np.asarray(hst["pmdec_corrected"], float)
    hea, hed = np.asarray(hst["pmra_corrected_err"], float), np.asarray(hst["pmdec_corrected_err"], float)
    hq = np.asarray(hst["selection_hq_astrometry_and_membership"], bool)
    fin = np.isfinite(hx) & np.isfinite(hy) & np.isfinite(ha) & np.isfinite(hd)
    idx = np.flatnonzero(fin)
    tree = cKDTree(np.column_stack([hx[fin], hy[fin]]))
    g = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    gx = (np.asarray(g["ra"], float)-OCEN_RA)*cosd*3600; gy = (np.asarray(g["dec"], float)-OCEN_DEC)*3600
    d0, j0 = tree.query(np.column_stack([gx, gy]), distance_upper_bound=0.4)
    ok0 = np.isfinite(d0)
    ox, oy = np.median(gx[ok0]-hx[fin][j0[ok0]]), np.median(gy[ok0]-hy[fin][j0[ok0]])
    dist, j = tree.query(np.column_stack([gx-ox, gy-oy]), distance_upper_bound=0.06)
    ok = np.isfinite(dist)
    gi, hi = np.flatnonzero(ok), idx[j[ok]]
    return dict(x=hx[hi], y=hy[hi], hst_q=hq[hi], memb=np.asarray(g["membership_prob"], float)[gi],
                g_mag=np.asarray(g["g_mag"], float)[gi],
                da=np.asarray(g["pmra"], float)[gi]-ha[hi], dd=np.asarray(g["pmdec"], float)[gi]-hd[hi],
                ea=np.hypot(np.asarray(g["pmra_error"], float)[gi], hea[hi]),
                ed=np.hypot(np.asarray(g["pmdec_error"], float)[gi], hed[hi]))


def ml_mean(v, e):
    """Gaussian ML mean with free intrinsic scatter s: v_i ~ N(m, e_i^2 + s^2)."""
    def nll(p):
        var = e**2+np.exp(2*p[1])
        return 0.5*np.sum((v-p[0])**2/var+np.log(var))
    p0 = [np.median(v), np.log(max(np.std(v), 1e-3))]
    return minimize(nll, p0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-8)).x


def main():
    all_hst = "--all-hst" in sys.argv
    tag = "_allhst" if all_hst else ""
    m = match()
    sel = np.isfinite(m["da"]) & np.isfinite(m["dd"]) & (m["memb"] > 0.9) & (m["hst_q"] | all_hst)
    for k in m:
        m[k] = m[k][sel]
    r = np.hypot(m["x"], m["y"])
    sys_a, sys_d = ml_mean(m["da"], m["ea"])[0], ml_mean(m["dd"], m["ed"])[0]
    ua, ud = m["da"]-sys_a, m["dd"]-sys_d
    # unit vectors: x = East (RA*), y = North; radial = (x, y)/r, tangential = (-y, x)/r (N -> E is clockwise
    # on the sky as seen with East left; sign is irrelevant for the amplitude comparison)
    vr = (ua*m["x"]+ud*m["y"])/r
    vt = (-ua*m["y"]+ud*m["x"])/r
    et = np.hypot(m["ea"]*m["y"], m["ed"]*m["x"])/r
    er = np.hypot(m["ea"]*m["x"], m["ed"]*m["y"])/r
    rng = np.random.default_rng(1)
    rows = []
    for lo, hi, rmed, t10, t10e in T10:
        k = (r >= lo) & (r < hi)
        n = int(k.sum())
        if n < 15:
            rows.append(dict(lo=lo, hi=hi, n=n)); continue
        mt, st = ml_mean(vt[k], et[k]); mr, _ = ml_mean(vr[k], er[k])
        boot = np.array([ml_mean(vt[k][b], et[k][b])[0] for b in rng.integers(0, n, (200, n))])
        bootr = np.array([ml_mean(vr[k][b], er[k][b])[0] for b in rng.integers(0, n, (200, n))])
        rows.append(dict(lo=lo, hi=hi, n=n, r_med=float(np.median(r[k])), mu_t=float(mt), mu_t_err=float(boot.std()),
                         mu_r=float(mr), mu_r_err=float(bootr.std()), scatter_t=float(np.exp(st)),
                         t10=float(t10), t10_err=float(t10e)))
    print(f"matched members (memb>0.9{'' if all_hst else ', HST hq'}): {sel.sum()}; systemic d = ({sys_a:.3f}, {sys_d:.3f}) mas/yr")
    sign = np.sign(np.nansum([x.get("mu_t", 0)*x.get("n", 0) for x in rows]))
    print(" ring    n   r_med   mu_T(ours)        Table10          diff/sig   mu_R(ours)")
    passed = 0; tested = 0
    for x in rows:
        if "mu_t" not in x:
            print(f"{x['lo']:4.0f}-{x['hi']:3.0f} {x['n']:4d}  too few"); continue
        d = (sign*x["mu_t"]-x["t10"])/np.hypot(x["mu_t_err"], x["t10_err"])
        tested += 1; passed += abs(d) < 2
        x["diff_sigma"] = float(d)
        print(f"{x['lo']:4.0f}-{x['hi']:3.0f} {x['n']:4d} {x['r_med']:6.1f}  {sign*x['mu_t']:.3f}+-{x['mu_t_err']:.3f}  "
              f"{x['t10']:.3f}+-{x['t10_err']:.3f}  {d:+5.1f}   {x['mu_r']:+.3f}+-{x['mu_r_err']:.3f}")
    print(f"within 2 sigma of Table 10: {passed}/{tested}")

    vb = Table.read(processed_dir()/"kinematics"/"vasiliev2021_ocen_pm_profiles.ecsv")
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.errorbar(T10[:, 2], T10[:, 3]*KMS, T10[:, 4]*KMS, fmt="s", color="k", ms=4, label="oMEGACat II Table 10 (Gaia FPR − HST)")
    ok = [x for x in rows if "mu_t" in x]
    ax.errorbar([x["r_med"] for x in ok], [sign*x["mu_t"]*KMS for x in ok], [x["mu_t_err"]*KMS for x in ok],
                fmt="o", color="C3", ms=4, label="ours (Gaia EDR3 − HST), tangential")
    ax.errorbar([x["r_med"]+4 for x in ok], [x["mu_r"]*KMS for x in ok], [x["mu_r_err"]*KMS for x in ok],
                fmt="o", mfc="none", color="C3", ms=4, alpha=.6, label="ours, radial (expect 0)")
    ax.plot(vb["r"], vb["vrot_pm"]*KMS, color="C0", lw=1.5, label="Vasiliev & Baumgardt 2021 curve (current HST term)")
    ax.axhline(0, color="k", lw=.5)
    ax.set_xlim(0, 370); ax.set_xlabel("projected radius R [arcsec]"); ax.set_ylabel("mean PM rotation [km/s] (D = 5.43 kpc)")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.savefig(ROOT/f"plots/inner_pm_rotation_20261113{tag}.png", dpi=140, bbox_inches="tight")
    (ROOT/f"results/plot_data/inner_pm_rotation_20261113{tag}.json").write_text(json.dumps(dict(
        n=int(sel.sum()), systemic=[float(sys_a), float(sys_d)], sign=float(sign), rings=rows), indent=1))


if __name__ == "__main__":
    main()
