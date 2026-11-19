#!/usr/bin/env python3
"""Best spherical Jeans fit vs best JAM fit (2026-11-19): data / model / residuals, plus ML parameters.

Spherical: spline beta, no halo (K1), UltraNest ML (results/fits/recover_real_K1_spline_nw; highest evidence).
JAM: shape A (global q'), turnover beta, no halo (K1), 18-start lsq (results/rotation/jamfit_real_K1_lsq_nw.json;
the NFW version gains only 0.78 lnL). Both after the Gaia node-weighting fix (commit 26399d8).
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/plot_best_models.py -> plots/best_models_20261119.png
"""
import json, sys, warnings
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src"))
from ocen_dm.kinematics.report import load_run_metadata, problem_for
warnings.filterwarnings("ignore")
SPH, JAM = "recover_real_K1_spline_nw", "jamfit_real_K1_lsq_nw.json"


def main():
    d = ROOT/"results/fits"/SPH
    P = problem_for(load_run_metadata(d)); x = np.load(d/"ml_x.npy"); s = json.loads((d/"summary.json").read_text())
    sph = {k: np.asarray(v) for k, v in P.predict(x).items()}
    J = json.loads((ROOT/"results/rotation"/JAM).read_text())
    profs = P.data.profiles
    for p in profs:                                         # same data in both fits?
        assert np.allclose(J["data_used"][p.name]["value"], p.value) and np.allclose(J["data_used"][p.name]["err_lo"], p.err_lo), p.name
    fits = {"spherical Jeans, spline beta, no DM": (sph, "C0"), "JAM (global q', i=44), turnover beta, no DM": (J["pred"], "C3")}
    fig, ax = plt.subplots(2, len(profs), figsize=(4.2*len(profs), 7.5), sharex="col", gridspec_kw=dict(height_ratios=(2.2, 1)))
    for j, p in enumerate(profs):
        ax[0, j].errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt="ko", ms=3, lw=.8)
        for lab, (pred, c) in fits.items():
            m = np.asarray(pred[p.name]); res = (m-p.value)/np.where(m > p.value, p.err_hi, p.err_lo)
            ax[0, j].plot(p.r, m, c=c, lw=1.4, label=f"{lab.split(',')[0]}: chi2 {np.sum(res**2):.0f} / {len(m)} bins")
            ax[1, j].plot(p.r, res, "o-", c=c, ms=2.5, lw=1)
        ax[0, j].set_title(p.name, fontsize=9); ax[0, j].legend(fontsize=7); ax[0, j].set_xscale("log")
        ax[1, j].axhline(0, c="k", lw=.5); ax[1, j].axhspan(-1, 1, color="0.9", zorder=0)
        ax[1, j].set_ylim(-5, 5); ax[1, j].set_xlabel("R (arcsec)")
    ax[0, 0].set_ylabel("dispersion (mas/yr; MUSE km/s)"); ax[1, 0].set_ylabel("(model - data) / err")
    fig.suptitle("Best models, no dark matter: blue = spherical Jeans (spline beta), red = JAM (global q', turnover beta)", fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, .96)); f = ROOT/"plots/best_models_20261119.png"; plt.savefig(f, dpi=90); plt.close()
    print("spherical ML:", dict(zip([q["name"] for q in s["model"]["parameters"]], np.round(x, 4))), f"lnL {s['lnL_max']:.2f} logZ {s['logz']:.2f}")
    print("JAM ML:", {k: round(v, 4) for k, v in J["theta"].items()}, f"q' {J['qobs']:.3f} lnL {-J['neg_lnpost']:.2f} chi2 {J['chi2_total']:.1f}")
    print(f)


if __name__ == "__main__":
    main()
