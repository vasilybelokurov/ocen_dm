#!/usr/bin/env python3
"""Summarise the shape-fit mock tests (bin/mock_shape_counts.py; docs/JAM_FIT_PLAN.md step 4d; 2026-11-15).

Per test and radius: mean bias <q_fit - q_true>, its standard error, s.d. of the pull (q_fit - q_true)/err and the
fraction of mocks within 1 sigma; PA bias. Pass (plan): |bias| < 0.01 and pull s.d. ~1 (we flag pull s.d. > 1.3).

Usage: python bin/report_shape_mocks.py  -> plots/shape_mocks_20261115.png, results/shape/mocks/summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TESTS = ("fit", "round", "step", "steep", "offset", "m1", "field0", "field2x", "crowdflat")
RS = (100, 205, 420, 600, 860, 1230, 1760)


def main():
    fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    summ = {}
    for j, t in enumerate(TESTS):
        f = ROOT/f"results/shape/mocks/mock_{t}.json"
        if not f.exists():
            continue
        recs = json.loads(f.read_text()); tr = json.loads(f.with_suffix(".truth.json").read_text())
        R = np.array(tr["R"]); qt = np.array(tr["q_true"])
        q = np.array([r["q_iso"] for r in recs]); e = np.array([r["q_iso_err"] for r in recs])
        d = q-qt; pull = d/e; n = len(recs)
        b, se, ps = d.mean(0), d.std(0, ddof=1)/np.sqrt(n), pull.std(0, ddof=1)
        ax[0].errorbar(R*(1+0.01*j), b, se, fmt="-o", ms=3, c=f"C{j}", label=f"{t} (n={n})")
        ax[1].plot(R, ps, "-", c=f"C{j}")
        pa = np.array([r["pa"]-tr["pa_true"] for r in recs])
        idx = [int(np.argmin(abs(R-r))) for r in RS]
        summ[t] = dict(n=n, R=[float(R[i]) for i in idx], q_true=[round(float(qt[i]), 3) for i in idx],
                       bias=[round(float(b[i]), 4) for i in idx], bias_se=[round(float(se[i]), 4) for i in idx],
                       pull_sd=[round(float(ps[i]), 2) for i in idx],
                       frac_1sig=[round(float(np.mean(abs(pull[:, i]) < 1)), 2) for i in idx],
                       median_err=[round(float(np.median(e[:, i])), 4) for i in idx],
                       pa_bias=round(float(pa.mean()), 2), pa_bias_se=round(float(pa.std(ddof=1)/np.sqrt(n)), 2))
        print(f"{t:7s} n={n:2d} PA bias {summ[t]['pa_bias']:+.2f}+-{summ[t]['pa_bias_se']:.2f}")
        for i, r in enumerate(summ[t]["R"]):
            print(f"   R {r:5.0f}\" q_true {summ[t]['q_true'][i]:.3f} bias {summ[t]['bias'][i]:+.4f}+-{summ[t]['bias_se'][i]:.4f}"
                  f" err {summ[t]['median_err'][i]:.3f} pull s.d. {summ[t]['pull_sd'][i]:.2f} in 1sig {summ[t]['frac_1sig'][i]:.2f}")
    for a in ax:
        a.set_xscale("log"); a.grid(alpha=.3)
    ax[0].axhspan(-0.01, 0.01, color="g", alpha=.1); ax[0].axhline(0, c="k", lw=.6)
    ax[0].set_ylabel("mean q'_fit - q'_true"); ax[0].legend(fontsize=8, ncol=2); ax[0].set_ylim(-0.06, 0.06)
    ax[1].axhline(1, c="k", lw=.6); ax[1].axhline(1.3, c="r", ls=":", lw=.8); ax[1].set_ylim(0, 3)
    ax[1].set_ylabel("s.d. of pull"); ax[1].set_xlabel("R (arcsec, isophote major axis)")
    plt.tight_layout(); p = ROOT/"plots/shape_mocks_20261115.png"; plt.savefig(p, dpi=100); print(p)
    (ROOT/"results/shape/mocks/summary.json").write_text(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
