#!/usr/bin/env python3
"""Report of the least-squares JAM mock tests (bin/fit_jam_halo.py --method lsq), 2026-11-14.

For each fit JSON: best -lnpost vs the truth's, the spread of the multi-start results ("n at or below truth + 0.1":
starts scoring <= truth + 0.1, which includes fits that beat the truth, as noisy-free fits with free parameters can),
parameter offsets, and the enclosed mass M(<r) of fit vs truth (total and per component), within spheres, of the
flattened JAM model (jam_halo_models.enclosed_flat: each Gaussian deprojected at the inclination; 2026-11-15, before:
a spherical proxy), plus residuals (model-data)/data.
Truth -lnpost values are computed by fit_jam_halo.Problem.objective at the truth (see JOURNAL 2026-11-14).

Usage: PYTHONPATH=src:bin python bin/report_jamfit_lsq.py
Output: plots/jamfit_lsq_mocks_20261114.png, results/rotation/jamfit_lsq_mocks_summary.json
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
from ocen_dm.kinematics.report import load_run_metadata, problem_for

Q85, Q90 = "q0.85_r0.85_i44_from_q0.90_r0.90_i44_v1", "q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1_injM1.35e+06_rs300"
CASES = [  # (label, json, truth -lnpost)
    ("q'0.85 no DM, K1 fit", f"jamfit_{Q85}_K1_lsq.json", -218.594),
    ("q'0.85 no DM, NFW fit", f"jamfit_{Q85}_NFW_lsq.json", -218.594),
    ("q'0.90 + NFW, NFW fit", f"jamfit_{Q90}_NFW_lsq.json", -217.282),
    ("q'0.90 + NFW, K1 fit", f"jamfit_{Q90}_K1_lsq.json", -217.282),
]
FAM = {k: problem_for(load_run_metadata(ROOT/"results/fits"/r)).family
       for k, r in dict(K1="jeans_K1_turnover_counts_c", NFW="jeans_K2_nfw_turnover_counts_c").items()}
R = np.geomspace(0.5, 70., 60)          # pc; data reach ~67 pc


def enclosed(theta, qobs, qrem, inc):
    """Per-component M(<r) [Msun] of the flattened model within spheres of radius R."""
    fam = FAM["NFW"] if "M_dm_100" in theta else FAM["K1"]
    jeans, _, _ = fam.build({n: theta[n] for n in fam.names})
    return H.enclosed_flat(jeans, qobs, qrem, inc, R)


def main():
    fig, ax = plt.subplots(len(CASES), 4, figsize=(19, 3.4*len(CASES)))
    summary = {}
    for i, (lab, fn, ftrue) in enumerate(CASES):
        d = json.loads((ROOT/"results/rotation"/fn).read_text())
        th, tr = dict(d["theta"], qobs=d["qobs"]), d["truth"]
        if "M_dm_100" not in tr:
            tr = dict(tr)
        fs = np.sort([s["f"] for s in d["starts"]])
        keys = [k for k in th if k in tr and tr[k]]
        off = {k: th[k]/tr[k]-1 for k in keys}
        a = ax[i, 0]
        a.barh(range(len(keys)), [100*off[k] for k in keys]); a.set_yticks(range(len(keys)))
        a.set_yticklabels(["q'" if k == "qobs" else k for k in keys], fontsize=8); a.axvline(0, c="k", lw=.5)
        a.set_xlim(-60, 60); a.set_xlabel("fit/truth - 1 (%)"); a.set_title(lab, fontsize=10)
        ef = enclosed(th, d["qobs"], d["qobs"], d["inc"])
        et = enclosed(tr, tr["qobs"], tr["qrem"], tr["inc"])
        tf, tt = sum(ef.values()), sum(et.values())
        a = ax[i, 1]
        a.plot(R, 100*(tf/tt-1), "k-", lw=2, label="total")
        for k in sorted(set(ef) | set(et)):
            if k in ef and k in et:
                a.plot(R, 100*(ef[k]/et[k]-1), lw=1, label=k)
            elif k in et:
                a.plot(R, -100*et[k]/tt, "--", lw=1, label=f"{k} (truth only, share of total)")
            else:
                a.plot(R, 100*ef[k]/tt, "--", lw=1, label=f"{k} (fit only, share of total)")
        a.axhline(0, c="k", lw=.5); a.set_xscale("log"); a.set_ylim(-30, 30); a.set_xlabel("r (pc)")
        a.set_ylabel("M(<r) fit/truth - 1 (%)"); a.legend(fontsize=6)
        a = ax[i, 2]
        a.plot(fs, "o"); a.axhline(ftrue, c="r", label=f"truth {ftrue:.2f}"); a.set_ylim(ftrue-3, ftrue+15)
        a.set_xlabel("start (sorted)"); a.set_ylabel("-lnpost"); a.legend(fontsize=8)
        a.set_title(f"best {fs[0]:.3f}; {np.sum(fs <= ftrue+0.1)}/{len(fs)} at or below truth + 0.1", fontsize=9)
        a = ax[i, 3]
        x0 = 0
        for name, p in d["pred"].items():
            v, s = np.array(d["data_values"][name]), d["stats"][name]
            res = np.array(p)-v
            a.plot(np.arange(x0, x0+len(v)), res/np.maximum(np.abs(v), 1e-9)*100, ".-", ms=3,
                   label=f"{name} chi2 {s['chi2']:.2f}")
            x0 += len(v)
        a.axhline(0, c="k", lw=.5); a.set_xlabel("bin (all datasets)"); a.set_ylabel("(model-data)/data (%)")
        a.legend(fontsize=6)
        m5, m20, m50 = (float(np.interp(r, R, tf/tt-1)) for r in (5, 20, 50))
        summary[lab] = dict(json=fn, best=float(fs[0]), truth=ftrue, starts=fs.tolist(),
                            n_at_or_below_truth_plus_0p1=int(np.sum(fs <= ftrue+0.1)), chi2=d["chi2_total"], qobs=d["qobs"],
                            q_truth=tr.get("qobs"), offsets=off, M_lt_5_20_50_frac=[m5, m20, m50],
                            M_dm_100=th.get("M_dm_100"), r_s=th.get("r_s"))
        print(f"{lab}: best {fs[0]:.3f} (truth {ftrue:.3f}), {np.sum(fs <= ftrue+0.1)}/{len(fs)} starts at or below truth + 0.1; "
              f"chi2 {d['chi2_total']:.2f}; q' {d['qobs']:.3f}; M(<5,20,50 pc) {100*m5:+.1f} {100*m20:+.1f} "
              f"{100*m50:+.1f}%; M_dm_100 {th.get('M_dm_100', 0):.3g}; r_s {th.get('r_s', 0):.0f}")
    plt.tight_layout()
    out = ROOT/"plots/jamfit_lsq_mocks_20261114.png"
    plt.savefig(out, dpi=100)
    (ROOT/"results/rotation/jamfit_lsq_mocks_summary.json").write_text(json.dumps(summary, indent=1))
    print(out)


if __name__ == "__main__":
    main()
