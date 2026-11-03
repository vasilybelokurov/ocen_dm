#!/usr/bin/env python3
"""Mock test of the model-side centre prediction (src/ocen_dm/centre_model.py).

Question: if the cluster centre is really at c (not AvdM10), does our binned estimator, run on
the real stars in the FIXED AvdM10 bins, measure what CentreGeometry predicts?

Mock: the real flagged HST stars (positions, per-star PM errors) and the real MUSE stars
(positions, LOS errors). True velocities drawn from an analytic, steeply rising velocity
profile centred at c (radial/tangential about c), plus Gaussian measurement errors.
Common random numbers: the same normals are used for the mock centred at AvdM10 (c = 0), so

    Delta_est  = est(mock_c) - est(mock_0)          (our mixture / ML estimator, fixed bins)
    Delta_pred = sqrt(<q(c)>) - sqrt(<q(0)>)        (per-star prediction, data axes about AvdM10)

differ only by estimator response, not by the stars' noise. Pass: |<Delta_est - Delta_pred>|
over seeds < 0.2 of the bin error in every bin (and consistent with zero given the seed scatter).
Also reported: est(mock_0) - sqrt(<q(0)>) (baseline estimator bias) and Delta_pred/err (the signal).

Usage: python bin/centre_model_mock.py [--seeds 4] [--workers 4] [--rmax 36.3]
Output: results/plot_data/centre_model_mock_<label>.json, plots/centre_model_mock_<label>.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import replace
from multiprocessing import get_context
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np

CENTRES = [(2., 0.), (-1.5, 1.5), (0., -4.)]


# analytic velocity profile: steeper than any of our DF models inside 30" (makes the test sensitive)
def s2r(r):
    return 0.62**2+0.30**2*(8./np.maximum(r, 1e-3))          # (mas/yr)^2


def s2t(r):
    return 0.93*s2r(r)


def s2los(r):
    return 18.**2+8.**2*(8./np.maximum(r, 1e-3))             # (km/s)^2


def _moments():
    from ocen_dm.centre_model import ModelMoments
    return ModelMoments.from_functions(s2r=s2r, s2t=s2t, s2los=s2los, Sigma=lambda r: 1./np.maximum(r, 1e-3))


def _setup(rmax):
    from ocen_dm.centre_marg import _hst_catalogue
    from ocen_dm.kinematics.hst_profile import load_hst_product, load_hst_sample
    t = load_hst_product()
    keep = np.asarray(t["r_upper"]) <= rmax
    edges = np.r_[np.asarray(t["r_lower"])[keep], np.asarray(t["r_upper"])[keep][-1]]
    s = load_hst_sample(require_flag=True)
    x, y = s.r_arcsec*np.sin(s.phi), s.r_arcsec*np.cos(s.phi)
    m = s.r_arcsec < edges[-1]+1
    cat = _hst_catalogue()
    return edges, s, x, y, m, cat


def _hst_mock(s, x, y, c, z):
    """Sample with PMs drawn about centre c (true) and projected on AvdM10 axes (data)."""
    dx, dy = x-c[0], y-c[1]
    rc = np.hypot(dx, dy)
    ux, uy = dx/rc, dy/rc
    vR, vT = z[0]*np.sqrt(s2r(rc)), z[1]*np.sqrt(s2t(rc))
    a = vR*ux-vT*uy+s.err_a*z[2]
    d = vR*uy+vT*ux+s.err_d*z[3]
    r0 = np.hypot(x, y)
    cp, sp = x/r0, y/r0
    return replace(s, mu_a=a, mu_d=d, mu_r=a*cp+d*sp, mu_t=-a*sp+d*cp)


def _muse_est(lx, ly, ev, edges_lo, edges_hi, c, z):
    from ocen_dm.centre_marg import _ml_sigma
    dx, dy = lx-c[0], ly-c[1]
    rc = np.hypot(dx, dy)
    v = z[0]*np.sqrt(s2los(rc))+ev*z[1]
    R0, phi = np.hypot(lx, ly), np.arctan2(lx, ly)
    out, n = [], []
    for lo, hi in zip(edges_lo, edges_hi):
        k = (R0 >= lo) & (R0 < hi)
        out.append(_ml_sigma(v[k], ev[k], phi[k])); n.append(int(k.sum()))
    return np.array(out), np.array(n)


def run_seed(args):
    seed, rmax = args
    from ocen_dm.kinematics.hst_profile import hst_profile
    edges, s, x, y, m, cat = _setup(rmax)
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((4, len(x)))
    zl = rng.standard_normal((2, len(cat["lx"])))
    muse_p = _muse_profile()
    out = {}
    for c in [(0., 0.)]+CENTRES:
        t0 = time.time()
        tab = hst_profile(edges_arcsec=edges, require_flag=True, min_stars=0, sample=_hst_mock(s, x, y, c, z))
        sl, nl = _muse_est(cat["lx"], cat["ly"], cat["ev"], muse_p[0], muse_p[1], c, zl)
        out[str(c)] = dict(pmr=np.asarray(tab["sigma_pmr"]).tolist(), pmr_err=np.asarray(tab["sigma_pmr_err"]).tolist(),
                           pmt=np.asarray(tab["sigma_pmt"]).tolist(), pmt_err=np.asarray(tab["sigma_pmt_err"]).tolist(),
                           n=np.asarray(tab["n_stars"]).tolist(), los=sl.tolist(), los_n=nl.tolist(),
                           seconds=time.time()-t0)
        print(f"seed {seed} c {c}: {time.time()-t0:.0f}s", flush=True)
    return seed, out


def _muse_profile():
    from ocen_dm.kinematics.likelihood import load_profile
    p = load_profile("muse_los_dispersion")
    return np.asarray(p.r_lower), np.asarray(p.r_upper)


def predictions(rmax):
    """sqrt(<q(c)>) per bin for HST radial/tangential (fixed bins, stars) and MUSE."""
    from ocen_dm.centre_model import PointSet, projected_variance
    edges, s, x, y, m, cat = _setup(rmax)
    mom = _moments()
    hst = PointSet.from_positions(x, y, edges[:-1], edges[1:])
    lo, hi = _muse_profile()
    mus = PointSet.from_positions(cat["lx"], cat["ly"], lo, hi)
    out = {}
    for c in [(0., 0.)]+CENTRES:
        out[str(c)] = dict(
            pmr=np.sqrt(hst.mean(projected_variance(mom, hst.x, hst.y, c, "r"))).tolist(),
            pmt=np.sqrt(hst.mean(projected_variance(mom, hst.x, hst.y, c, "t"))).tolist(),
            los=np.sqrt(mus.mean(projected_variance(mom, mus.x, mus.y, c, "los"))).tolist())
    r_hst = np.sqrt(edges[:-1]*edges[1:])
    r_muse = np.sqrt(np.maximum(lo, 0.3)*hi)
    return out, r_hst.tolist(), r_muse.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--rmax", type=float, default=36.3)
    ap.add_argument("--label", default="20261103")
    a = ap.parse_args()
    t0 = time.time()
    pred, r_hst, r_muse = predictions(a.rmax)
    with get_context("spawn").Pool(a.workers) as pool:
        est = dict(pool.map(run_seed, [(s, a.rmax) for s in range(a.seeds)]))
    summary = {}
    for c in CENTRES:
        row = {}
        for comp in ("pmr", "pmt", "los"):
            d_est = np.array([np.subtract(est[s][str(c)][comp], est[s][str((0., 0.))][comp]) for s in est])
            d_pred = np.subtract(pred[str(c)][comp], pred[str((0., 0.))][comp])
            if comp == "los":
                sig = np.array(est[0][str((0., 0.))]["los"]); n = np.array(est[0][str((0., 0.))]["los_n"])
                err = sig/np.sqrt(2*n)                                   # approximate bin error
            else:
                err = np.mean([est[s][str((0., 0.))][comp+"_err"] for s in est], axis=0)
            res = (d_est-d_pred)/err
            base = (np.array([est[s][str((0., 0.))][comp] for s in est])-np.array(pred[str((0., 0.))][comp]))/err
            row[comp] = dict(signal=(d_pred/err).tolist(), resid_mean=res.mean(0).tolist(),
                             resid_sem=(res.std(0, ddof=1)/np.sqrt(len(res))).tolist() if len(res) > 1 else None,
                             baseline_bias_mean=base.mean(0).tolist(), err=err.tolist())
        row["pass"] = bool(all(np.all(np.abs(row[k]["resid_mean"]) < 0.2) for k in ("pmr", "pmt", "los")))
        summary[str(c)] = row
        print(c, "pass" if row["pass"] else "FAIL",
              {k: (np.round(np.max(np.abs(row[k]["signal"])), 2), np.round(np.max(np.abs(row[k]["resid_mean"])), 3))
               for k in ("pmr", "pmt", "los")}, flush=True)
    out = dict(centres=CENTRES, r_hst=r_hst, r_muse=r_muse, predictions=pred, estimates=est, summary=summary,
               seeds=a.seeds, rmax=a.rmax, elapsed_s=time.time()-t0)
    path = ROOT/"results"/"plot_data"/f"centre_model_mock_{a.label}.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"wrote {path.relative_to(ROOT)} in {time.time()-t0:.0f}s")
    plot(out, ROOT/"plots"/f"centre_model_mock_{a.label}.png")


def plot(out, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(13, 7), sharex="col")
    cols = ["C0", "C1", "C2"]
    for j, (comp, r, title) in enumerate((("pmr", out["r_hst"], "HST PM radial"), ("pmt", out["r_hst"], "HST PM tangential"),
                                          ("los", out["r_muse"], "MUSE LOS"))):
        for c, col in zip(out["centres"], cols):
            s = out["summary"][str(tuple(c))][comp]
            ax[0, j].plot(r, s["signal"], "o-", color=col, ms=3, label=f"c = ({c[0]:+.1f}, {c[1]:+.1f})\"")
            e = s["resid_sem"] if s["resid_sem"] is not None else None
            ax[1, j].errorbar(r, s["resid_mean"], yerr=e, fmt="o", color=col, ms=3, capsize=2)
        ax[0, j].set_title(title)
        ax[0, j].axhline(0, color="k", lw=.5)
        ax[1, j].axhspan(-.2, .2, color="0.85")
        ax[1, j].axhline(0, color="k", lw=.5)
        ax[1, j].set_xscale("log")
        ax[1, j].set_xlabel("R about AvdM10 [arcsec]")
    ax[0, 0].set_ylabel("signal: predicted change / bin error")
    ax[1, 0].set_ylabel("(estimated - predicted change) / bin error")
    ax[0, 0].legend(fontsize=8)
    fig.suptitle("Mock: true centre at c, data binned about AvdM10 (mean over seeds; band = 0.2 error pass)")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
