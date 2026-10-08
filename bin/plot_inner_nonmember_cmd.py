#!/usr/bin/env python3
"""CMD of the proper-motion non-member candidates within 10'' of the omega Cen centre (oMEGACat).

Selection as bin/plot_field_counts_outside.py: PM errors < 0.5 mas/yr, chi2 > 25 against the cluster
(sigma = 0.81 mas/yr inside 30''), |mu_rel| > 4 mas/yr. Candidates are split by the half of the PM plane they
fall in: towards the field centroid (-2.58, +5.37) mas/yr ("field side") or away from it ("far side", which holds
no field population and traces cluster leakage). Stars in the 3-mas/yr field window and the mirrored control
window (bin/plot_field_counts_inner.py) are ringed.
Left: R < 10''; background = all other stars there with PMs (cluster). Right: same selection at 150-270'',
where the field dominates, as the reference for where real field stars sit in the CMD.
Frame: M_F625W and (F625W - F814W)_0 with D = 5.43 kpc, A_F625W = 0.316, A_F814W = 0.219 (MIST, as in
bin/plot_cmd_radial.py).

Usage: python bin/plot_inner_nonmember_cmd.py
Output: plots/omegacat_inner10_nonmember_cmd.png, results/plot_data/omegacat_inner10_nonmember_cmd.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

DM0, A6, A8 = 13.674, 0.316, 0.219
FC = np.array([-2.58, 5.37])
INK, MUTED = "#1f1f1e", "#6b6a64"


def main():
    t = fits.open(ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits", memmap=True)[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.); y = 0.04*(np.asarray(t["y"], float)-15000.)
    R = np.hypot(x, y)
    ma, md = np.asarray(t["pmra_corrected"], float), np.asarray(t["pmdec_corrected"], float)
    ea, ed = np.asarray(t["pmra_corrected_err"], float), np.asarray(t["pmdec_corrected_err"], float)
    f6, f8 = np.asarray(t["f625w"], float), np.asarray(t["f814w"], float)
    M = f6 - DM0 - A6; C = (f6 - f8) - (A6 - A8)
    phot = np.isfinite(f6) & np.isfinite(f8)
    s = np.interp(R, [15, 45, 90, 160, 270], [0.81, 0.78, 0.73, 0.66, 0.58])
    par = np.isfinite(ma) & (ea < 0.5) & (ed < 0.5)
    chi2 = ma**2/(s**2+ea**2) + md**2/(s**2+ed**2)
    mu = np.hypot(ma, md)
    cand = par & (chi2 > 25) & (mu > 4)
    side = (ma*FC[0]+md*FC[1]) > 0
    win_f = cand & (np.hypot(ma-FC[0], md-FC[1]) < 3)
    win_c = cand & (np.hypot(ma+FC[0], md+FC[1]) < 3)
    fig, ax = plt.subplots(1, 2, figsize=(13, 7.2), sharey=True)
    out = {}
    for a, (lo, hi, title) in zip(ax, ((0, 10, "R < 10''"), (150, 270, "150-270'' (reference: field dominates)"))):
        ring = (R >= lo) & (R < hi)
        bg = ring & np.isfinite(ma) & ~cand & phot
        a.plot(C[bg], M[bg], ".", color="#c9c8c0", ms=2 if lo == 0 else 0.6, alpha=0.7, rasterized=True,
               label=f"other stars with PMs (cluster), N = {bg.sum():,}")
        fs = ring & cand & side & phot; fa = ring & cand & ~side & phot
        a.plot(C[fs], M[fs], "o", color="#2a78d6", ms=7 if lo == 0 else 3, mec="white", mew=0.6, alpha=0.9,
               label=f"|mu| > 4, field side, N = {fs.sum()}")
        a.plot(C[fa], M[fa], "s", color="#eb6834", ms=6 if lo == 0 else 3, mec="white", mew=0.6, alpha=0.9,
               label=f"|mu| > 4, far side (leakage), N = {fa.sum()}")
        if lo == 0:
            wf = ring & win_f & phot; wc = ring & win_c & phot
            a.plot(C[wf], M[wf], "o", mfc="none", mec="#0d366b", ms=15, mew=2, label=f"in field window, N = {wf.sum()}")
            a.plot(C[wc], M[wc], "o", mfc="none", mec="#6b6a64", ms=15, mew=2, label=f"in control window, N = {wc.sum()}")
            nophot = ring & cand & ~phot
            out["inner10"] = dict(n_candidates=int((ring & cand).sum()), n_field_side=int((ring & cand & side).sum()),
                                  n_far_side=int((ring & cand & ~side).sum()), n_without_photometry=int(nophot.sum()),
                                  stars=[dict(id=int(t["ID"][k]), R=float(R[k]), pmra=float(ma[k]), pmdec=float(md[k]),
                                              mu=float(mu[k]), side="field" if side[k] else "far", M_F625W=float(M[k]),
                                              colour0=float(C[k]), field_window=bool(win_f[k]), control_window=bool(win_c[k]))
                                         for k in np.where(ring & cand)[0]])
        a.set_xlim(-0.6, 2.6); a.set_ylim(11.5, -1.5)
        a.set_xlabel("(F625W - F814W)$_0$", color=INK)
        a.set_title(title, fontsize=10.5, color=INK, loc="left")
        a.legend(fontsize=8, frameon=False, loc="lower right" if lo == 0 else "upper right", markerscale=1.0)
        a.tick_params(colors=MUTED, labelsize=9)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    ax[0].set_ylabel("$M_{F625W}$", color=INK)
    fig.suptitle("oMEGACat: CMD of proper-motion non-member candidates (|mu| > 4 mas/yr, chi2 > 25, PM err < 0.5); "
                 "D = 5.43 kpc, E(B-V) = 0.12", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(ROOT/"plots/omegacat_inner10_nonmember_cmd.png", dpi=130)
    (ROOT/"results/plot_data/omegacat_inner10_nonmember_cmd.json").write_text(json.dumps(out, indent=1))
    print({k: v for k, v in out["inner10"].items() if k != "stars"})


if __name__ == "__main__":
    main()
