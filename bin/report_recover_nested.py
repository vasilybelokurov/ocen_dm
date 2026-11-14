#!/usr/bin/env python3
"""Report the UltraNest recoveries of the tuned flattened JAM mocks (bin/recover_tuned_mock_nested.py).

For each truth (round control, q'0.90, q'0.85) and each spherical family (K1 no DM, NFW):
posterior median and 16/84 of M(<1,5,20,50 pc), M_tot, M_bh, D, beta(4,15,25,40 pc) from 400 posterior draws, against
the truth (oblate mass inside spheres + BH; shell-averaged spherical-equivalent beta; recover_tuned_mock.py), the bias
in units of that posterior's own s.d., lnL_max and lnZ for NFW - K1, and the NFW halo mass posterior.
Plot: mock data, ML model of each family, normalised residuals; one row per truth.

Usage: PYTHONPATH=src:bin python bin/report_recover_nested.py
Output: results/rotation/recover_nested_summary.json, plots/recover_nested_20261114.png
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
from recover_tuned_mock import R_BETA, R_MASS, beta_sph_truth, oblate_enclosed, summarise
from tune_jam_to_data import Setup

TRUTHS = (("round", "q1.00_r1.00_i90_from_q1.00_r1.00_i90_v1"),
          ("q'0.90", "q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1"),
          ("q'0.85", "q0.85_r0.85_i44_from_q0.90_r0.90_i44_v1"))
FAMS = {"K1": "jeans_K1_turnover_counts_c", "NFW": "jeans_K2_nfw_turnover_counts_c"}
KEYS = ("M_r", "M_tot", "M_bh", "D", "beta")


def truth_of(S, T):
    theta = dict(zip(T["names"], T["x"]))
    jeans, D, _ = S.fam.build(theta)
    lum, pot, mbh = J.jam_models(jeans, D, T["qobs"], T["qrem"])
    return dict(M_r=oblate_enclosed(pot, mbh, D, T["inc"], R_MASS).tolist(),
                M_tot=float(J.mge_total(*pot, D*1e3*np.pi/648000.)+mbh), M_bh=float(mbh), D=float(D),
                beta=beta_sph_truth(lum, pot, mbh, D, T["inc"], jeans.anisotropy, R_BETA).tolist())


def main():
    S = Setup(per_bin=10)
    fams = {k: family_from_config(family_config(problem_for(load_run_metadata(ROOT/"results/fits"/v)).family))
            for k, v in FAMS.items()}
    rng = np.random.default_rng(1)
    out = {}
    fig, axes = plt.subplots(2*len(TRUTHS), 5, figsize=(20, 4.2*len(TRUTHS)), sharex="col",
                             gridspec_kw=dict(height_ratios=[2, 1]*len(TRUTHS)))
    for it, (lab, tag) in enumerate(TRUTHS):
        T = json.loads((ROOT/f"results/rotation/jam_tuned_{tag}.json").read_text())
        mdata = KinematicData(tuple(replace(p, value=np.array(T["pred"][p.name])) for p in S.data.profiles))
        truth = truth_of(S, T)
        res = dict(truth=truth)
        for fk, fam in fams.items():
            d = ROOT/f"results/fits/recover_{tag}_{fk}"
            summ = json.loads((d/"summary.json").read_text())
            post = Table.read(d/"posterior.ecsv", format="ascii.ecsv")
            X = np.array([post[n] for n in fam.names]).T
            ss = [summarise(fam, X[j]) for j in rng.choice(len(X), min(400, len(X)), replace=False)]
            st = {}
            for k in KEYS+("M_dm_100",):
                a = np.array([s[k] for s in ss], float)
                p16, p50, p84 = np.percentile(a, [16, 50, 84], axis=0)
                st[k] = dict(p16=np.atleast_1d(p16).tolist(), p50=np.atleast_1d(p50).tolist(),
                             p84=np.atleast_1d(p84).tolist(), sd=np.atleast_1d(a.std(0)).tolist())
                if k in truth:
                    st[k]["bias_sd"] = ((np.atleast_1d(p50)-np.atleast_1d(truth[k]))/np.atleast_1d(a.std(0))).tolist()
            x_ml = np.load(d/"ml_x.npy")
            res[fk] = dict(lnL_max=summ["lnL_max"], logz=summ["logz"], logzerr=summ["logzerr"],
                           chi2_ml=summ["chi2_ml"], chi2_ml_total=summ["chi2_ml_total"], stats=st,
                           M_dm_100_95=float(np.percentile(np.asarray(post["M_dm_100"]), 95)) if fk == "NFW" else None)
            pred = FitProblem(fam, mdata).predict(x_ml)
            for ip, p in enumerate(mdata.profiles):
                ax, rx = axes[2*it, ip], axes[2*it+1, ip]
                if fk == "K1":
                    ax.errorbar(p.r, p.value, [p.err_lo, p.err_hi], fmt="o", ms=3, c="k", label="mock")
                c = "C0" if fk == "K1" else "C3"
                ax.plot(p.r, pred[p.name], c=c, label=f"{fk} ML")
                e = np.where(pred[p.name] > p.value, p.err_hi, p.err_lo)
                rx.plot(p.r, (p.value-pred[p.name])/e, "o-" if fk == "K1" else "s--", ms=3, c=c)
                rx.axhline(0, c="0.5", lw=0.5)
                ax.set_xscale("log")
                if it == 0:
                    ax.set_title(p.name, fontsize=9)
                if ip == 0:
                    ax.set_ylabel(f"{lab}\nsigma [mas/yr; MUSE km/s]"); rx.set_ylabel("(mock-model)/err")
                if it == len(TRUTHS)-1:
                    rx.set_xlabel("R [arcsec]")
        axes[2*it, 0].legend(fontsize=8)
        res["dlnL_NFW_K1"] = res["NFW"]["lnL_max"]-res["K1"]["lnL_max"]
        res["dlnZ_NFW_K1"] = res["NFW"]["logz"]-res["K1"]["logz"]
        out[lab] = res
        k = res["K1"]["stats"]
        print(f"== {lab}: K1 chi2 {res['K1']['chi2_ml_total']:.1f}, NFW chi2 {res['NFW']['chi2_ml_total']:.1f}; "
              f"dlnL {res['dlnL_NFW_K1']:+.2f}; dlnZ {res['dlnZ_NFW_K1']:+.2f}; NFW M_dm_100 median "
              f"{res['NFW']['stats']['M_dm_100']['p50'][0]:.3g}, 95% {res['NFW']['M_dm_100_95']:.3g}")
        for fk in ("K1", "NFW"):
            s = res[fk]["stats"]
            print(f"  {fk}: M(<1,5,20,50)/truth-1 {np.round(np.array(s['M_r']['p50'])/truth['M_r']-1, 3)} "
                  f"[sd {np.round(s['M_r']['bias_sd'], 1)}]; M_tot {s['M_tot']['p50'][0]/truth['M_tot']-1:+.3f} "
                  f"[{s['M_tot']['bias_sd'][0]:+.1f}]; M_bh {s['M_bh']['p50'][0]/truth['M_bh']-1:+.3f} "
                  f"[{s['M_bh']['bias_sd'][0]:+.1f}]; D {s['D']['p50'][0]-truth['D']:+.3f} [{s['D']['bias_sd'][0]:+.1f}]")
            print(f"       beta(4,15,25,40) {np.round(s['beta']['p50'], 2)} truth {np.round(truth['beta'], 2)} "
                  f"[sd {np.round(s['beta']['bias_sd'], 1)}]")
    fig.suptitle("Spherical Jeans (UltraNest) fits to noiseless tuned no-DM JAM mocks: K1 (blue), NFW (red)")
    fig.tight_layout()
    fig.savefig(ROOT/"plots/recover_nested_20261114.png", dpi=110)
    (ROOT/"results/rotation/recover_nested_summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
