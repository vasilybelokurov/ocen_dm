#!/usr/bin/env python3
"""Inner-node fits (beta(2") = -0.5, b10 <= 0; 2026-11-19): data / model / residuals for the spherical spline (local MAP,
results/limits/sph_local_{K1,NFW}_inner.json) and JAM C6 (results/rotation/jamfit_real_{K1,NFW}_lsq_qhyb600_remstars_
bspline6_in2_nw.json), each without (K1) and with an NFW halo.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/plot_inner_fits.py -> plots/inner_fits_20261119.png
"""
import json, sys, warnings
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
from sph_local import Local
warnings.filterwarnings("ignore")


def main():
    fits = {}
    for i, k in enumerate(("K1", "NFW")):
        L = Local(k); r = json.loads((ROOT/f"results/limits/sph_local_{k}_inner.json").read_text())["best"]
        fits[f"spherical {k}"] = ({n: np.asarray(v) for n, v in L.P.predict(np.array(r["x"])).items()}, f"C{i}", "-")
        J = json.loads((ROOT/f"results/rotation/jamfit_real_{k}_lsq_qhyb600_remstars_bspline6_in2_nw.json").read_text())
        fits[f"JAM C6 {k}"] = (J["pred"], f"C{i+2}", "--")
    profs = L.P.data.profiles
    for p in profs:
        assert np.allclose(J["data_used"][p.name]["value"], p.value) and np.allclose(J["data_used"][p.name]["err_lo"], p.err_lo), p.name
    fig, ax = plt.subplots(2, len(profs), figsize=(4.2*len(profs), 7.5), sharex="col", gridspec_kw=dict(height_ratios=(2.2, 1)))
    for j, p in enumerate(profs):
        ax[0, j].errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt="ko", ms=3, lw=.8)
        for lab, (pred, c, ls) in fits.items():
            m = np.asarray(pred[p.name]); res = (m-p.value)/np.where(m > p.value, p.err_hi, p.err_lo)
            ax[0, j].plot(p.r, m, c=c, ls=ls, lw=1.3, label=f"{lab}: chi2 {np.sum(res**2):.1f}/{len(m)}")
            ax[1, j].plot(p.r, res, "o", c=c, ls=ls, ms=2.5, lw=.9)
        ax[0, j].set_title(p.name, fontsize=9); ax[0, j].legend(fontsize=6.5); ax[0, j].set_xscale("log")
        ax[1, j].axhline(0, c="k", lw=.5); ax[1, j].axhspan(-1, 1, color="0.9", zorder=0)
        ax[1, j].set_ylim(-5, 5); ax[1, j].set_xlabel("R (arcsec)")
    ax[0, 0].set_ylabel("dispersion (mas/yr; MUSE km/s)"); ax[1, 0].set_ylabel("(model - data) / err")
    fig.suptitle('Inner beta node fits: spherical spline (solid) and JAM C6 (dashed), no DM vs NFW', fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, .96)); f = ROOT/"plots/inner_fits_20261119.png"; plt.savefig(f, dpi=90); plt.close(); print(f)


if __name__ == "__main__":
    main()
