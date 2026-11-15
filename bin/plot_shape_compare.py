#!/usr/bin/env python3
"""Compare the projected-shape measurements of omega Cen (2026-11-15): isophotal q'(R) of the unbinned fits
(bin/fit_shape_counts.py variants), the ring-table small-eps estimates (bin/shape_counts.py), and the literature
(van de Ven et al. 2006, A&A 445, 513, Table 2: q'_j of their MGE, set to the observed flattening of Geyer et al. 1983 at
R = sigma'_j; https://ui.adsabs.harvard.edu/abs/2006A&A...445..513V; values read from the paper, arXiv astro-ph/0509228).

Usage: python bin/plot_shape_compare.py   -> plots/shape_compare_20261115.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VDV_S = np.array([0.153, 1.477, 2.525, 3.691, 5.219, 7.534, 11.069, 17.547])*60      # major-axis sigma' (arcsec)
VDV_Q = np.array([1.000, 0.934, 0.877, 0.848, 0.850, 0.836, 0.866, 0.926])
LABEL = {"pafree_h19_g17members_r2400": "fit: PA free, HST<19, Gaia G<17 PM members",
         "pafixed_h19_g17members_r2400": "fit: PA fixed 104.2",
         "pafree_h18_g16members_r2400": "fit: bright (HST<18, G<16 members)",
         "pafree_h19_g17all_r2400": "fit: Gaia all stars (field term)",
         "pafree_h19_g17all_r2400_quad": "fit: Gaia all stars, curved field"}


def main():
    fig, ax = plt.subplots(2, 1, figsize=(9, 9), sharex=True, gridspec_kw=dict(height_ratios=(2, 1)))
    rows = []
    for j, (tag, lab) in enumerate(LABEL.items()):
        f = ROOT/f"results/shape/shape_fit_{tag}.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        R, q, e = (np.array(d[k]) for k in ("R_iso", "q_iso", "q_iso_err"))
        ax[0].plot(R, q, "-", lw=2 if j == 0 else 1, label=lab, c=f"C{j}")
        if j == 0:
            ax[0].fill_between(R, q-e, q+e, color="C0", alpha=.2)
        if d.get("pa_err", 0) > 0:                              # global PA of the fit, drawn as a band (free-PA fits)
            ax[1].axhspan(d["pa"]-d["pa_err"], d["pa"]+d["pa_err"], color=f"C{j}", alpha=.15, label=f"{lab}: global PA")
        rows.append((lab, d))
    tab = json.loads((ROOT/"results/shape/ring_table_20261115.json").read_text())
    for name, mk in (("HST F625W<19", "s"), ("Gaia G<17 PM members (s_int 0.5)", "D")):
        t = tab[name]
        r = np.array([np.sqrt(x["r_lower"]*x["r_upper"]) for x in t])
        ax[0].errorbar(r, 1-np.array([x["eps"] for x in t]), [x["eps_err"] for x in t], fmt=mk, c="gray", ms=4, lw=.7,
                       mfc="none", label=f"ring table (small eps): {name}")
        ax[1].errorbar(r, [x["m2"]["pa"] for x in t], [x["m2"]["pa_err"] for x in t], fmt=mk, c="gray", ms=3, lw=.6, mfc="none")
    ax[0].plot(VDV_S, VDV_Q, "k*", ms=11, label="van de Ven+06 MGE q'_j at R = sigma'_j (from Geyer+83)")
    ax[0].axhline(np.cos(np.radians(44)), c="r", ls=":", lw=1, label="cos 44 deg (oblate limit at i = 44)")
    ax[0].axhline(1, c="k", lw=.5); ax[0].set_ylim(0.6, 1.2); ax[0].set_ylabel("projected axis ratio q'")
    ax[0].legend(fontsize=7, loc="lower left"); ax[0].grid(alpha=.3)
    ax[1].axhline(104.2, c="k", ls="--", lw=.8, label="kinematic axis 104.2"); ax[1].set_ylim(0, 180)
    ax[1].legend(fontsize=6, loc="lower left"); ax[1].set_ylabel("PA (deg E of N)")
    ax[1].set_xscale("log"); ax[1].set_xlabel("R (arcsec; fits: major-axis radius of the isophote)"); ax[1].grid(alpha=.3)
    ax[1].set_xlim(50, 3000)
    plt.tight_layout(); p = ROOT/"plots/shape_compare_20261115.png"; plt.savefig(p, dpi=100); print(p)
    for lab, d in rows:
        qi = np.interp(VDV_S[1:], d["R_iso"], d["q_iso"])
        print(f"{lab}: PA {d['pa']:.1f}+-{d.get('pa_err', 0):.1f}; q_iso at vdV radii {np.round(qi, 3).tolist()}; "
              f"field lnF {d['x'][d['names'].index('lnF')]:.2f}")
    print("van de Ven+06:", VDV_Q[1:].tolist(), "at", np.round(VDV_S[1:]).tolist())


if __name__ == "__main__":
    main()
