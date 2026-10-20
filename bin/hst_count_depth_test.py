#!/usr/bin/env python3
"""Does the central HST count deficit survive deeper magnitude cuts? (2026-10-20)

oMEGACat VI catalogue, all stars with measured F625W (no quality flag, as the C1 count product), centre pixel (15000,
15000) at 0.04"/px (the adopted centre) and the oMEGACat VI kinematic centre (dRA, dDec = -1.31", +1.44").
Cuts F625W < 19, 20, 21, 22, 23. For each cut:
  - density profile in the C1 bins inside 60" (coverage 1 inside 30"; beyond, the same for every cut, so it cancels
    in ratios), normalised to its mean over 20-60";
  - N(<cut)/N(<19) per bin: flat if the fainter stars are as complete as the bright ones at every radius AND have the
    same spatial distribution; a central drop means either crowding incompleteness or mass segregation (the two are not
    separable without artificial stars);
  - the inner deficit: N(2.0-4.1") against a flat core at the 6.7-17.6" density (binomial given the total), as in
    bin/jwst_star_counts.py.
Output: plots/hst_count_depth_test.png, results/plot_data/hst_count_depth_test.json
Usage: PYTHONPATH=src python bin/hst_count_depth_test.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
CUTS = (19, 20, 21, 22, 23)
KINEMATIC = (-1.31, 1.44)
INNER, REF = (2.0, 4.15), (6.6, 17.7)


def main():
    from astropy.table import Table
    from scipy.stats import binom
    from ocen_dm.kinematics.counts import CountProfile
    from ocen_dm.paths import processed_dir, raw_dir
    t = Table.read(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")
    m = np.asarray(t["f625w"], float); ok = np.isfinite(m)
    x = -0.04*(np.asarray(t["x"], float)[ok]-15000.); y = 0.04*(np.asarray(t["y"], float)[ok]-15000.); m = m[ok]
    c = CountProfile.from_dict(json.loads((processed_dir()/"kinematics/ocen_counts_hst_f625w19_c.json").read_text()))
    e = np.r_[c.r_lower, c.r_upper[-1]]; e = e[e <= 61.]; r = np.sqrt(e[1:]*e[:-1]); A = np.pi*(e[1:]**2-e[:-1]**2)
    lo, hi = e[:-1], e[1:]; inn = (lo >= INNER[0]-1e-6) & (hi <= INNER[1]); ref = (lo >= REF[0]) & (hi <= REF[1]); nrm = (lo >= 20.) & (hi <= 61.)
    f = A[inn].sum()/(A[inn].sum()+A[ref].sum())
    out = dict(edges_arcsec=e.tolist(), inner=[float(lo[inn][0]), float(hi[inn][-1])], reference=[float(lo[ref][0]), float(hi[ref][-1])],
               cuts={})
    for cen, c0 in (("adopted", (0., 0.)), ("kinematic", KINEMATIC)):
        R = np.hypot(x-c0[0], y-c0[1])
        n19 = np.histogram(R[m < 19], e)[0]
        for cut in CUTS:
            n = np.histogram(R[m < cut], e)[0]; Ni, Nr = int(n[inn].sum()), int(n[ref].sum())
            out["cuts"].setdefault(str(cut), {})[cen] = dict(
                counts=n.tolist(), dens_norm=(n/A/np.mean((n/A)[nrm])).tolist(), ratio_to_19=(n/np.maximum(n19, 1)).tolist(),
                n_inner=Ni, expected_flat=f*(Ni+Nr), deficit=Ni/(f*(Ni+Nr)), p=float(binom.cdf(Ni, Ni+Nr, f)))
    (ROOT/"results/plot_data/hst_count_depth_test.json").write_text(json.dumps(out, indent=1, default=float))
    print(f'inner {out["inner"]} vs flat at {out["reference"]} arcsec')
    for cut in CUTS:
        a, k = out["cuts"][str(cut)]["adopted"], out["cuts"][str(cut)]["kinematic"]
        print(f"F625W < {cut}: adopted N_in {a['n_inner']:5d} exp {a['expected_flat']:7.1f} ratio {a['deficit']:.2f} p {a['p']:.2g} | "
              f"kinematic N_in {k['n_inner']:5d} exp {k['expected_flat']:7.1f} ratio {k['deficit']:.2f} p {k['p']:.2g}")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    for i, cut in enumerate(CUTS):
        for cen, ls in (("adopted", "o-"), ("kinematic", "--")):
            d = out["cuts"][str(cut)][cen]; n = np.array(d["counts"]); y_ = np.array(d["dens_norm"])
            if cen == "adopted":
                ax[0].errorbar(r*(1+0.015*(i-2)), y_, y_/np.sqrt(np.maximum(n, 1)), fmt=ls, ms=3, capsize=2, color=f"C{i}", label=f"F625W < {cut}")
                ax[1].plot(r, d["ratio_to_19"], ls, ms=3, color=f"C{i}", label=f"N(<{cut})/N(<19)")
            else:
                ax[0].plot(r, y_, ls, lw=0.8, color=f"C{i}", label="kinematic centre" if i == 2 else None)
    for cen, mk in (("adopted", "o"), ("kinematic", "s")):
        ax[2].plot(CUTS, [out["cuts"][str(c_)][cen]["deficit"] for c_ in CUTS], mk+"-", label=f"{cen} centre")
        for c_ in CUTS:
            d = out["cuts"][str(c_)][cen]
            ax[2].annotate(f"p={d['p']:.1g}", (c_, d["deficit"]), fontsize=7, xytext=(3, 4 if cen == "adopted" else -10), textcoords="offset points")
    ax[0].axvspan(*INNER, color="0.9", zorder=0); ax[0].axvspan(*REF, color="0.95", zorder=0); ax[0].axhline(1, color="0.6", lw=0.5)
    ax[0].set(xscale="log", xlabel="R [arcsec]", ylabel="density / mean(20-60\")", title="HST oMEGACat counts by depth (shaded: inner, reference)")
    ax[0].legend(fontsize=7)
    ax[1].set(xscale="log", xlabel="R [arcsec]", ylabel="N(<cut) / N(<19)", title="faint-to-bright ratio per bin (adopted centre)"); ax[1].legend(fontsize=7)
    ax[2].axhline(1, color="0.6", lw=0.5)
    ax[2].set(xlabel="F625W cut", ylabel="N(2-4.1\") / flat expectation from 6.7-17.6\"", title="inner deficit vs depth (binomial p)", ylim=(0.5, 1.2))
    ax[2].legend(fontsize=8)
    fig.savefig(ROOT/"plots/hst_count_depth_test.png", dpi=110); print(ROOT/"plots/hst_count_depth_test.png")


if __name__ == "__main__":
    main()
