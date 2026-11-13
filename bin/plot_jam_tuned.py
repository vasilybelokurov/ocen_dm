#!/usr/bin/env python3
"""Plot the tuned flattened JAM models (bin/tune_jam_to_data.py) against omega Cen's 89 binned data.

Data: the snapshot with the corrected HST rotation terms (Setup._corrected_data). Models: saved best-fit
predictions in results/rotation/jam_tuned_<tag>.json. Output: plots/jam_tuned_redo_20261113.png.
Usage: PYTHONPATH=src:bin python bin/plot_jam_tuned.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from tune_jam_to_data import Setup

TAGS = [("q1.00_r1.00_i90_from_q1.00_r1.00_i90_v1", "round (q'=1)", "k", "-"),
        ("q0.95_r0.95_i44_from_q0.90_r0.90_i44_v1", "q'=0.95, i=44 (unstable)", "C0", "--"),
        ("q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1", "q'=0.90, i=44", "C3", "-"),
        ("q0.85_r0.85_i44_from_q0.90_r0.90_i44_v1", "q'=0.85, i=44", "C2", ":")]     # best of the 2026-11-13 redo
S = Setup(per_bin=10)          # only the data (corrected terms and errors) are used here
prof = {p.name: p for p in S.data.profiles}
res = {t: json.loads((ROOT/f"results/rotation/jam_tuned_{t}.json").read_text()) for t, *_ in TAGS}
fig = plt.figure(figsize=(13, 7.2))
gs = fig.add_gridspec(2, 3, height_ratios=[2.2, 1], hspace=0.06, wspace=0.3)
order = [("muse_los_dispersion", 0), ("hst_pm_radial_ours", 1), ("gaia_edr3_ours_radial", 1),
         ("hst_pm_tangential_ours", 2), ("gaia_edr3_ours_tangential", 2)]
axes = {c: (fig.add_subplot(gs[0, c]), fig.add_subplot(gs[1, c])) for c in range(3)}
for name, c in order:
    p = prof[name]; a, r = axes[c]; mk = "s" if "gaia" in name else "o"
    a.errorbar(p.r, p.value, [p.err_lo, p.err_hi], fmt=mk, color="0.4", ms=3, lw=0.7,
               mfc="none" if "gaia" in name else "0.4")
    for t, lab, col, ls in TAGS:
        pred = np.array(res[t]["pred"][name])
        a.plot(p.r, pred, ls, color=col, lw=1.2)
        e = np.where(pred > p.value, p.err_hi, p.err_lo)
        r.plot(p.r, (p.value-pred)/e, mk+ls, color=col, ms=2.5, lw=0.7)
    for q in (a, r): q.set_xscale("log")
    a.set_xticklabels([]); r.axhline(0, color="k", lw=0.5); r.set_ylim(-4, 4); r.set_xlabel("projected radius [arcsec]")
axes[0][0].set_ylabel(r"$\sigma_{\rm los}$ [km/s]")
axes[1][0].set_ylabel(r"$\sigma_{\rm R}$ [mas/yr] (HST ●, Gaia □)")
axes[2][0].set_ylabel(r"$\sigma_{\rm T}$ [mas/yr] (HST ●, Gaia □)")
axes[0][1].set_ylabel("(data − model)/err")
for t, lab, col, ls in TAGS:
    chi2 = sum(v["chi2"] for v in res[t]["stats"].values())
    axes[0][0].plot([], [], ls, color=col, label=f"{lab}: χ²={chi2:.0f}")
axes[0][0].legend(fontsize=8, title="tuned JAM, no DM (89 bins)", title_fontsize=8)
out = ROOT/"plots/jam_tuned_redo_20261113.png"; fig.savefig(out, dpi=130, bbox_inches="tight"); print(out)
