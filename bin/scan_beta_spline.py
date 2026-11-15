#!/usr/bin/env python3
"""Can a flexible anisotropy profile reconcile the MEASURED shape with the kinematics? (2026-11-15; adversarial
review step 3; diagnostic, not a fit of the model).

At FIXED masses/distance (base K1 fit, results/rotation/jamfit_real_K1_lsq_base.json) only the anisotropy beta(r) is
changed: piecewise linear in log r through nodes at NODES arcsec (constant outside). The score uses the ratios
sigma_T/sigma_R of HST and Gaia per bin, which cancel the mass and distance scale (PM ~ sqrt(M)/D for both components;
only the subtracted streaming term breaks this slightly). The 5 node values are optimised (scipy least_squares, trf,
bounds [-1, 0.95], 3 starts) on that ratio chi2 for each shape x inclination:
  shapes: matched q'_k (results/shape/q_matched_r1300.json; remnants interpolated, or round), global q' 0.891 (control).
Afterwards the full data chi2 is evaluated after an analytic rescale: v^2 -> lam v^2 (mass) and PM^2 -> mu PM^2
(mu = (D0/D)^2, distance), on a (lam, mu) grid; calibration: the base model must return its own chi2 (~75).
Necessary DF conditions (spherical, An & Evans 2006, https://ui.adsabs.harvard.edu/abs/2006ApJ...642..752A):
gamma(r) >= 2 beta(r) (central cusp-anisotropy, used here at all r as a warning only) and
(1-beta)(2-beta) - (d beta/d ln r)/2 >= 0 (separable augmented density; An, Van Hese & Baes 2012).

Usage: OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=src:bin python bin/scan_beta_spline.py --workers 3
Output: results/rotation/beta_spline_scan.json, plots/beta_spline_scan_20261115.png
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np

NODES = np.array([10., 100., 250., 700., 2000.])               # arcsec
STARTS = ([0., 0., 0., 0., 0.], [0., 0., 0.3, 0.5, 0.5], [0., 0.3, 0.8, 0.3, -0.5])
BASE = "results/rotation/jamfit_real_K1_lsq_base.json"
QM = "results/shape/q_matched_r1300.json"
PAIRS = (("hst_pm_tangential_ours", "hst_pm_radial_ours"), ("gaia_edr3_ours_tangential", "gaia_edr3_ours_radial"))
_S = None


class SplineBeta:
    """beta(r) linear in log r through (nodes [pc], values), constant outside; .beta(r_pc) as the JAM code expects."""

    def __init__(self, nodes_pc, vals):
        self.lx, self.v = np.log(nodes_pc), np.asarray(vals, float)

    def beta(self, r):
        return np.interp(np.log(np.maximum(np.asarray(r, float), 1e-6)), self.lx, self.v)

    def dlog(self, r):
        lr = np.log(np.asarray(r, float)); h = 1e-4
        return (np.interp(lr+h, self.lx, self.v)-np.interp(lr-h, self.lx, self.v))/(2*h)


def setup():
    global _S
    if _S is None:
        import jam_halo_models as H
        from tune_jam_to_data import Setup
        _S = Setup(run="jeans_K1_turnover_counts_c", models=H.jam_models)
    return _S


def shape(name):
    import jam_halo_models as H
    if name == "global":
        q = json.loads((ROOT/BASE).read_text())["qobs"]; return q, q
    Q = H.QMatched.from_json(ROOT/QM)
    return Q, (Q if name == "matched" else 1.0)


def predict(theta, qs, qr, inc, anis=None):
    """Prediction with the base theta; anis (SplineBeta) replaces the fitted anisotropy if given."""
    S = setup()
    if anis is None:
        return S.predict(theta, qs, qr, inc)[0]
    orig = S.grid_moments
    S.grid_moments = lambda lum, pot, mbh, D, _a, inc_: orig(lum, pot, mbh, D, anis, inc_)
    try:
        return S.predict(theta, qs, qr, inc)[0]
    finally:
        S.grid_moments = orig


def ratio_resid(pred, P):
    out = []
    for t, r in PAIRS:
        pt, pr = P[t], P[r]
        rho = pt.value/pr.value
        e = rho*np.hypot(0.5*(pt.err_lo+pt.err_hi)/pt.value, 0.5*(pr.err_lo+pr.err_hi)/pr.value)
        out.append((pred[t]/pred[r]-rho)/e)
    return np.concatenate(out)


def full_chi2(pred, P):
    """Per-dataset chi2 after the best (lam, mu) rescale of the second moments (grid search)."""
    tot = {n: pred[n]**2+(0 if p.streaming2 is None else p.streaming2) for n, p in P.items()}
    best = (np.inf, None)
    for lam in np.linspace(0.7, 1.4, 141):
        for mu in np.linspace(0.85, 1.15, 61):
            c = {}
            for n, p in P.items():
                x = lam*(1. if p.kind == "los" else mu)*tot[n]-(0 if p.streaming2 is None else p.streaming2)
                if np.any(x <= 0):
                    c = None; break
                v = np.sqrt(x); e = np.where(v > p.value, p.err_hi, p.err_lo); c[n] = float(np.sum(((v-p.value)/e)**2))
            if c and sum(c.values()) < best[0]:
                best = (sum(c.values()), dict(lam=float(lam), mu=float(mu), chi2=c))
    return best[1]


def df_check(theta, anis_beta, dlog, D):
    """Necessary DF conditions on 2-2600" for the stellar tracer: min of gamma - 2 beta and of R2."""
    from ocen_dm.kinematics.backends import _split_mass_model
    from ocen_dm.kinematics.df_consistency import mge_log_slopes
    S = setup(); jeans, D, _ = S.fam.build(theta); pc = D*1e3/206264.806
    tracer = _split_mass_model(jeans.mass)[0][0]
    r = np.geomspace(2., 2600., 200)*pc
    g, _ = mge_log_slopes(tracer, r)
    b = anis_beta(r); db = dlog(r)
    r2 = (1-b)*(2-b)-db/2
    return dict(min_gamma_minus_2beta=float(np.min(g-2*b)), r_at_min=float(r[np.argmin(g-2*b)]/pc),
                min_R2=float(np.min(r2)))


def run(job):
    from scipy.optimize import least_squares
    shape_name, inc = job
    S = setup(); P = {p.name: p for p in S.data.profiles}
    A = json.loads((ROOT/BASE).read_text()); th = A["theta"]
    jeans, D, _ = S.fam.build(th); pc = D*1e3/206264.806
    qs, qr = shape(shape_name)
    out = dict(shape=shape_name, inc=inc)
    try:                                                      # fitted (turnover) anisotropy, no change
        p0 = predict(th, qs, qr, inc)
        out["turnover"] = dict(ratio_chi2=float(np.sum(ratio_resid(p0, P)**2)), full=full_chi2(p0, P),
                               TR={t: (p0[t]/p0[r]).tolist() for t, r in PAIRS})
    except ValueError as e:
        out["turnover"] = dict(error=str(e))

    def res(v):
        try:
            return ratio_resid(predict(th, qs, qr, inc, SplineBeta(NODES*pc, v)), P)
        except ValueError:
            return np.full(30, 30.)
    runs = []
    for s0 in STARTS:
        sol = least_squares(res, s0, bounds=(-1., 0.95), method="trf", diff_step=1e-2, x_scale=0.3, max_nfev=200)
        runs.append(dict(x=sol.x.tolist(), chi2=float(np.sum(sol.fun**2)), nfev=int(sol.nfev), status=int(sol.status)))
    best = min(runs, key=lambda r: r["chi2"])
    sb = SplineBeta(NODES*pc, best["x"])
    pb = predict(th, qs, qr, inc, sb)
    out["spline"] = dict(best=best, runs=runs, full=full_chi2(pb, P), TR={t: (pb[t]/pb[r]).tolist() for t, r in PAIRS},
                         df=df_check(th, sb.beta, sb.dlog, D))
    print(shape_name, inc, "turnover ratio chi2", round(out["turnover"].get("ratio_chi2", np.nan), 1),
          "| spline ratio chi2", [round(r["chi2"], 1) for r in runs], "best beta", np.round(best["x"], 2).tolist(), flush=True)
    return out


def plot(results, P):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    incs = sorted({r["inc"] for r in results}); shapes = ["global", "matched", "matched_remround"]
    fig, ax = plt.subplots(3, len(incs), figsize=(5*len(incs), 11), sharey="row")
    col = dict(zip(shapes, ("C0", "C3", "C1")))
    for j, inc in enumerate(incs):
        for k, (t, r) in enumerate(PAIRS):
            pt, pr = P[t], P[r]; rho = pt.value/pr.value
            e = rho*np.hypot(0.5*(pt.err_lo+pt.err_hi)/pt.value, 0.5*(pr.err_lo+pr.err_hi)/pr.value)
            ax[k, j].errorbar(pr.r, rho, e, fmt="ko", ms=3, lw=.8, label="data")
            for res in results:
                if res["inc"] != inc:
                    continue
                c = col[res["shape"]]
                if "TR" in res["turnover"]:
                    ax[k, j].plot(pr.r, res["turnover"]["TR"][t], ":", c=c, lw=1.2,
                                  label=f"{res['shape']} fitted beta ({res['turnover']['ratio_chi2']:.0f})")
                ax[k, j].plot(pr.r, res["spline"]["TR"][t], "-", c=c, lw=1.5,
                              label=f"{res['shape']} spline beta ({res['spline']['best']['chi2']:.0f})")
            ax[k, j].set_xscale("log"); ax[k, j].axhline(1, c="k", lw=.4)
            ax[k, j].set_title(f"i = {inc:g}: {'HST' if k == 0 else 'Gaia'} sigma_T/sigma_R", fontsize=9)
        ax[0, j].legend(fontsize=6)
        rr = np.geomspace(3, 2500, 200)
        for res in results:
            if res["inc"] == inc:
                ax[2, j].plot(rr, np.interp(np.log(rr), np.log(NODES), res["spline"]["best"]["x"]), c=col[res["shape"]],
                              label=res["shape"])
        ax[2, j].plot(NODES, 0*NODES, "k|", ms=10); ax[2, j].set_xscale("log"); ax[2, j].axhline(0, c="k", lw=.4)
        ax[2, j].set_xlabel("r (arcsec; 3D radius)"); ax[2, j].set_title("best spline beta(r) (fixed masses)", fontsize=9)
        ax[2, j].legend(fontsize=7)
    ax[0, 0].set_ylabel("sigma_T/sigma_R (HST)"); ax[1, 0].set_ylabel("sigma_T/sigma_R (Gaia)"); ax[2, 0].set_ylabel("beta")
    fig.suptitle("Fixed masses (base K1): shape vs anisotropy; numbers = T/R chi2 (HST 21 + Gaia 9 bins)", fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, .97)); p = ROOT/"plots/beta_spline_scan_20261115.png"; plt.savefig(p, dpi=90)
    return p


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--inc", type=float, nargs="+", default=[44., 50., 55.]); a = ap.parse_args()
    jobs = [(s, i) for i in a.inc for s in ("global", "matched", "matched_remround")]
    with ProcessPoolExecutor(a.workers, mp_context=get_context("spawn")) as ex:
        results = list(ex.map(run, jobs))
    (ROOT/"results/rotation/beta_spline_scan.json").write_text(json.dumps(dict(nodes_arcsec=NODES.tolist(), base=BASE,
                                                                                qmatched=QM, results=results), indent=1))
    P = {p.name: p for p in setup().data.profiles}
    print(plot(results, P))
    for r in results:
        f = r["spline"]["full"]
        print(f"{r['shape']:17s} i{r['inc']:g}: T/R chi2 turnover {r['turnover'].get('ratio_chi2', np.nan):6.1f} -> spline "
              f"{r['spline']['best']['chi2']:6.1f}; full chi2 (lam {f['lam']:.2f}, mu {f['mu']:.2f}) "
              + " ".join(f"{k.split('_')[0]}{k.split('_')[-1][:1]} {v:.0f}" for k, v in f["chi2"].items())
              + f"; DF min(gamma-2b) {r['spline']['df']['min_gamma_minus_2beta']:.2f} min R2 {r['spline']['df']['min_R2']:.2f}")
        t = r["turnover"].get("full")
        if t:
            print(f"{'':17s}      turnover full chi2: " + " ".join(f"{v:.0f}" for v in t["chi2"].values())
                  + f" (sum {sum(t['chi2'].values()):.0f}, lam {t['lam']:.2f}, mu {t['mu']:.2f})")


if __name__ == "__main__":
    main()
