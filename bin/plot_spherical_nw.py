#!/usr/bin/env python3
"""Spherical Jeans refits after the Gaia node-weighting fix (2026-11-16): data / ML model / residuals for turnover and
spline beta, K1 vs NFW (results/fits/recover_real_{K1,NFW}{,_spline}_nw). Each run is replayed from its own saved
config and data snapshot. Usage: PYTHONPATH=src:bin python bin/plot_spherical_nw.py -> plots/sph_nw_fits_20261116.png"""
import json, sys, warnings
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src"))
from ocen_dm.kinematics.report import load_run_metadata, problem_for
warnings.filterwarnings("ignore")
RUNS = {"turnover K1": "recover_real_K1_nw", "turnover NFW": "recover_real_NFW_nw",
        "spline K1": "recover_real_K1_spline_nw", "spline NFW": "recover_real_NFW_spline_nw"}


def main():
    fits = {}
    for lab, d in RUNS.items():
        P = problem_for(load_run_metadata(ROOT/"results/fits"/d)); x = np.load(ROOT/"results/fits"/d/"ml_x.npy")
        s = json.loads((ROOT/"results/fits"/d/"summary.json").read_text())
        fits[lab] = (P, P.predict(x), s, P.loglike_vector(x))
    profs = next(iter(fits.values()))[0].data.profiles
    fig, ax = plt.subplots(2, len(profs), figsize=(4.2*len(profs), 7.5), sharex="col", gridspec_kw=dict(height_ratios=(2.2, 1)))
    for j, p in enumerate(profs):
        ax[0, j].errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt="ko", ms=3, lw=.8)
        for i, (lab, (P, pred, s, lnl)) in enumerate(fits.items()):
            m = np.asarray(pred[p.name]); res = (m-p.value)/np.where(m > p.value, p.err_hi, p.err_lo)
            ls = "--" if "turnover" in lab else "-"
            ax[0, j].plot(p.r, m, ls, c=f"C{i}", lw=1.3, label=f"{lab}: chi2 {np.sum(res**2):.0f}")
            ax[1, j].plot(p.r, res, ls, c=f"C{i}", lw=1)
        ax[0, j].set_title(p.name, fontsize=9); ax[0, j].legend(fontsize=6); ax[0, j].set_xscale("log")
        ax[1, j].axhline(0, c="k", lw=.5); ax[1, j].set_ylim(-6, 6); ax[1, j].set_xlabel("R (arcsec)")
    ax[0, 0].set_ylabel("dispersion (mas/yr or km/s)"); ax[1, 0].set_ylabel("(model-data)/err")
    fig.suptitle("spherical Jeans (UltraNest ML), Gaia nodes membership-weighted: turnover (dashed) vs spline beta (solid)", fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, .96)); f = ROOT/"plots/sph_nw_fits_20261116.png"; plt.savefig(f, dpi=90); plt.close()
    for lab, (P, pred, s, lnl) in fits.items():
        print(f"{lab:14s} lnL {lnl:7.2f} (saved {s['lnL_max']:.2f}) logZ {s['logz']:.2f}+-{s['logzerr']:.2f}")
    print(f)


if __name__ == "__main__":
    main()
