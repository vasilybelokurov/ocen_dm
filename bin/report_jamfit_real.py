#!/usr/bin/env python3
"""Data / model / residual plot of the real-data JAM least-squares fits, K1 vs NFW (2026-11-14).

Reads results/rotation/jamfit_real_K1_lsq.json and jamfit_real_NFW_lsq.json (bin/fit_jam_halo.py --data real
--method lsq); radii and errors come from tune_jam_to_data.Setup (the corrected real data the fits used; ~3 GB).

Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/report_jamfit_real.py
Output: plots/jamfit_real_lsq_20261114.png
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

import jam_halo_models as H
from tune_jam_to_data import Setup

FITS = {"K1 (no DM)": "jamfit_real_K1_lsq.json", "NFW": "jamfit_real_NFW_lsq.json"}
COL = {"K1 (no DM)": "C0", "NFW": "C3"}


def main():
    S = Setup(run="jeans_K1_turnover_counts_c", models=H.jam_models)
    fits = {k: json.loads((ROOT/"results/rotation"/v).read_text()) for k, v in FITS.items()}
    profs = S.data.profiles
    fig, ax = plt.subplots(2, len(profs), figsize=(4.2*len(profs), 7), sharex="col",
                           gridspec_kw=dict(height_ratios=(2.2, 1)))
    for j, p in enumerate(profs):
        a, b = ax[0, j], ax[1, j]
        a.errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt="ko", ms=3, lw=.8, label="data")
        for k, d in fits.items():
            m = np.array(d["pred"][p.name])
            a.plot(p.r, m, "-", c=COL[k], label=f"{k}: chi2 {d['stats'][p.name]['chi2']:.1f}")
            res = (m-p.value)/np.where(m > p.value, p.err_hi, p.err_lo)
            b.plot(p.r, res, "o-", c=COL[k], ms=3, lw=.8)
        a.set_title(p.name, fontsize=10); a.legend(fontsize=7); a.set_xscale("log")
        b.axhline(0, c="k", lw=.5); b.set_xlabel("R (arcsec)"); b.set_ylim(-4, 4)
    ax[0, 0].set_ylabel("dispersion (mas/yr or km/s)"); ax[1, 0].set_ylabel("(model-data)/err")
    lines = []
    for k, d in fits.items():
        t = d["theta"]
        fs = sorted(s["f"] for s in d["starts"])
        lines.append(f"{k}: -lnpost {d['neg_lnpost']:.2f}, chi2 {d['chi2_total']:.1f}, q' {d['qobs']:.3f}, "
                     f"D {d['D']:.3f} kpc, M* {t['M_star']:.3g}, Mrem {t['M_rem']:.3g} (a {t['a_rem']:.2f} pc), "
                     f"Mbh {t['M_bh']:.3g}" + (f", M_dm(<100pc) {t['M_dm_100']:.3g}, r_s {t['r_s']:.0f} pc"
                                                if "M_dm_100" in t else "")
                     + f"; starts within 0.5 of best: {sum(f < fs[0]+0.5 for f in fs)}/{len(fs)}")
    fig.suptitle("\n".join(lines), fontsize=9)
    plt.tight_layout(rect=(0, 0, 1, 0.9))
    out = ROOT/"plots/jamfit_real_lsq_20261114.png"
    plt.savefig(out, dpi=100)
    print("\n".join(lines)); print(out)


if __name__ == "__main__":
    main()
