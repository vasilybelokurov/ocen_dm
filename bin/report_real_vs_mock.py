#!/usr/bin/env python3
"""Real-data spherical fit (corrected HST rotation terms) against the flattened-mock fingerprint (2026-11-14).

1. Residuals (data - ML model)/err of the spherical K1 and NFW fits to the REAL data (results/fits/recover_real_*),
   overlaid on the K1 residuals of the spherical fits to the noiseless tuned JAM mocks (round, q'0.90, q'0.85).
   Same data errors everywhere, so the curves are directly comparable. Correlation coefficient of the real K1
   residuals with each mock's K1 residuals (all 89 bins) is printed and shown.
2. Halo gain: dlnL and dlnZ (NFW - K1), real vs mocks.
3. Real spherical K1 posterior vs the flattened JAM fit to the same data (tuned q'0.90 and round): M(<1,5,20,50 pc),
   M_bh, D. For the JAM fit the mass is the oblate model inside spheres (recover_tuned_mock.oblate_enclosed).

Usage: PYTHONPATH=src:bin python bin/report_real_vs_mock.py
Output: plots/real_vs_mock_residuals_20261114.png, results/rotation/real_vs_mock_20261114.json
"""
from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

import jam_mock_rotation_test as J
from ocen_dm.kinematics import FitProblem, KinematicData
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config
from recover_tuned_mock import R_MASS, oblate_enclosed, summarise
from tune_jam_to_data import Setup

MOCKS = (("round", "q1.00_r1.00_i90_from_q1.00_r1.00_i90_v1", "0.6"),
         ("q'0.90", "q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1", "C1"),
         ("q'0.85", "q0.85_r0.85_i44_from_q0.90_r0.90_i44_v1", "C2"))
FAMS = {"K1": "jeans_K1_turnover_counts_c", "NFW": "jeans_K2_nfw_turnover_counts_c"}


def resid(fam, data, run_dir):
    pred = FitProblem(fam, data).predict(np.load(run_dir/"ml_x.npy"))
    out = {}
    for p in data.profiles:
        e = np.where(pred[p.name] > p.value, p.err_hi, p.err_lo)
        out[p.name] = (p.value-pred[p.name])/e
    return out


def main():
    S = Setup(per_bin=10)
    fams = {k: family_from_config(family_config(problem_for(load_run_metadata(ROOT/"results/fits"/v)).family))
            for k, v in FAMS.items()}
    names = [p.name for p in S.data.profiles]
    real = {k: resid(fams[k], S.data, ROOT/f"results/fits/recover_real_{k}") for k in FAMS}
    summ = {k: json.loads((ROOT/f"results/fits/recover_real_{k}/summary.json").read_text()) for k in FAMS}
    out = dict(real=dict(chi2_K1=summ["K1"]["chi2_ml"], chi2_NFW=summ["NFW"]["chi2_ml"],
                         dlnL=summ["NFW"]["lnL_max"]-summ["K1"]["lnL_max"], dlnZ=summ["NFW"]["logz"]-summ["K1"]["logz"],
                         dlnZ_err=float(np.hypot(summ["NFW"]["logzerr"], summ["K1"]["logzerr"]))), mocks={})
    rv = np.concatenate([real["K1"][n] for n in names])
    mres = {}
    for lab, tag, _ in MOCKS:
        T = json.loads((ROOT/f"results/rotation/jam_tuned_{tag}.json").read_text())
        md = KinematicData(tuple(replace(p, value=np.array(T["pred"][p.name])) for p in S.data.profiles))
        mres[lab] = resid(fams["K1"], md, ROOT/f"results/fits/recover_{tag}_K1")
        mv = np.concatenate([mres[lab][n] for n in names])
        ms = {k: json.loads((ROOT/f"results/fits/recover_{tag}_{k}/summary.json").read_text()) for k in FAMS}
        out["mocks"][lab] = dict(corr_with_real=float(np.corrcoef(rv, mv)[0, 1]),
                                 dlnL=ms["NFW"]["lnL_max"]-ms["K1"]["lnL_max"], dlnZ=ms["NFW"]["logz"]-ms["K1"]["logz"])

    fig, axes = plt.subplots(1, 5, figsize=(21, 4.2))
    for ax, n, p in zip(axes, names, S.data.profiles):
        for lab, _, c in MOCKS:
            ax.plot(p.r, mres[lab][n], "-", c=c, lw=2, alpha=0.8,
                    label=f"mock {lab} (r={out['mocks'][lab]['corr_with_real']:+.2f})")
        ax.plot(p.r, real["K1"][n], "o-", c="k", ms=4, lw=1, label="REAL, K1")
        ax.plot(p.r, real["NFW"][n], "s--", c="C3", ms=3, lw=1, label="REAL, NFW")
        ax.axhline(0, c="0.5", lw=0.5); ax.set_xscale("log"); ax.set_title(n, fontsize=9)
        ax.set_xlabel("R [arcsec]")
    axes[0].set_ylabel("(data - spherical ML model)/err")
    axes[0].legend(fontsize=7)
    fig.suptitle(f"Spherical Jeans residuals: real data (corrected HST rotation) vs noiseless flattened no-DM JAM mocks; "
                 f"real dlnZ(NFW-K1) {out['real']['dlnZ']:+.2f}")
    fig.tight_layout()
    fig.savefig(ROOT/"plots/real_vs_mock_residuals_20261114.png", dpi=110)

    # spherical posterior on the real data vs the flattened JAM fits to the same data
    post = Table.read(ROOT/"results/fits/recover_real_K1/posterior.ecsv", format="ascii.ecsv")
    X = np.array([post[n] for n in fams["K1"].names]).T
    ss = [summarise(fams["K1"], X[j]) for j in np.random.default_rng(1).choice(len(X), 400, replace=False)]
    sph = {k: np.percentile(np.array([s[k] for s in ss], float), [16, 50, 84], axis=0).tolist()
           for k in ("M_r", "M_bh", "D")}
    out["spherical_real_K1"] = sph
    print(f"REAL: K1 chi2 {summ['K1']['chi2_ml_total']}, NFW chi2 {summ['NFW']['chi2_ml_total']}; "
          f"dlnL {out['real']['dlnL']:+.2f}; dlnZ {out['real']['dlnZ']:+.2f} +- {out['real']['dlnZ_err']:.2f}")
    for lab in out["mocks"]:
        m = out["mocks"][lab]
        print(f"mock {lab}: corr(real K1 resid, mock K1 resid) {m['corr_with_real']:+.2f}; dlnL {m['dlnL']:+.2f}; "
              f"dlnZ {m['dlnZ']:+.2f}")
    print(f"spherical K1 (real): M(<1,5,20,50) {np.round(np.array(sph['M_r'][1])/1e6, 4)} e6 "
          f"(+- {np.round((np.array(sph['M_r'][2])-np.array(sph['M_r'][0]))/2e6, 4)}); "
          f"M_bh {sph['M_bh'][1]:.0f} +- {(sph['M_bh'][2]-sph['M_bh'][0])/2:.0f}; D {sph['D'][1]:.3f} +- "
          f"{(sph['D'][2]-sph['D'][0])/2:.3f}")
    out["jam_fits"] = {}
    for lab, tag, _ in MOCKS:
        T = json.loads((ROOT/f"results/rotation/jam_tuned_{tag}.json").read_text())
        jeans, D, _ = S.fam.build(dict(zip(T["names"], T["x"])))
        lum, pot, mbh = J.jam_models(jeans, D, T["qobs"], T["qrem"])
        Mr = oblate_enclosed(pot, mbh, D, T["inc"], R_MASS)
        out["jam_fits"][lab] = dict(M_r=Mr.tolist(), M_bh=float(mbh), D=float(D))
        print(f"JAM {lab} (best fit, same data): M(<1,5,20,50) {np.round(Mr/1e6, 4)} e6; M_bh {mbh:.0f}; D {D:.3f}")
    (ROOT/"results/rotation/real_vs_mock_20261114.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
