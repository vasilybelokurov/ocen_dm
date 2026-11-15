#!/usr/bin/env python3
"""Match per-Gaussian q'_k of the JAM stellar MGE to the measured isophotal q'(R) (2026-11-15).

Why: giving each Gaussian the measured q' at R = sigma_major' (jam_halo_models.QProfile) does not reproduce the measured
isophotes once broad Gaussians with different q' blend (model flattest at ~860" instead of ~420"; isophote chi2 510 vs 2).
Here q'_k are fitted directly: minimise sum_R [(q_iso,model(R) - q_iso,meas(R))/err(R)]^2 over rmin..rmax, with
q'_k in [qlo, 1]. The relative stellar MGE weights and widths are fixed by the count profile (theta only rescales the
mass; checked: identical to 1e-16 across fits), so one matched set serves every fit and inclination.
Ties: Gaussians narrower than --tie-in (default 80", < 1% of the light, unconstrained) take the q' of the next one; the
outermost (1500", ~0 mass) takes the q' of the 772" one.
Errors: diagonal q_iso_err (neighbouring points of the smooth 4-node profile are correlated, so the chi2 is not a
goodness-of-fit and the parameter errors are indicative only).

Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/match_q_isophotes.py [--rmax 1300] [--qlo 0.729]
Output: results/shape/q_matched_r<rmax>.json, plots/q_matched_r<rmax>_20261115.png
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, least_squares

import jam_halo_models as H
from ocen_dm.kinematics.backends import _split_mass_model
from ocen_dm.kinematics.report import load_run_metadata, problem_for

SHAPE = "results/shape/shape_fit_pafree_h19_g17members_r2400.json"
BASE = "results/rotation/jamfit_real_K1_lsq_base.json"


def stellar_mge():
    """(fractions, circular widths in arcsec) of the JAM stellar MGE (at the base-fit theta; theta-independent)."""
    pb = problem_for(load_run_metadata(ROOT/"results/fits/jeans_K1_turnover_counts_c"))
    th = json.loads((ROOT/BASE).read_text())["theta"]
    jeans, D, _ = pb.family.build(th); pc = D*1e3/H.J.ARCSEC_PER_RAD
    g = _split_mass_model(jeans.mass)[0][0]
    m = np.asarray(g.masses)*g.mass_to_light
    return m/m.sum(), np.asarray(g.sigmas)/pc


def iso_q(f, s, q, R):
    """Isophotal axis ratio b/R of sum_k f_k/(2 pi s_k^2) exp(-(x^2 + y^2/q_k^2)/(2 s_k'^2)), s_k' = s_k/sqrt(q_k)."""
    surf, sm = f/(2*np.pi*s**2), s/np.sqrt(q)
    S = lambda x, y: np.sum(surf*np.exp(-0.5*(x**2+(y/q)**2)/sm**2))
    out = []
    for r in np.atleast_1d(R):
        t = S(r, 0.)
        out.append(brentq(lambda b: S(0., b)-t, 1e-3*r, 3*r)/r)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmin", type=float, default=70.); ap.add_argument("--rmax", type=float, default=1300.)
    ap.add_argument("--qlo", type=float, default=float(np.cos(np.radians(44)))+0.01)
    ap.add_argument("--tie-in", type=float, default=80.)
    a = ap.parse_args()
    S = json.loads((ROOT/SHAPE).read_text())
    Rm, qm, em = (np.array(S[k]) for k in ("R_iso", "q_iso", "q_iso_err"))
    f, s = stellar_mge()
    free = np.where((s >= a.tie_in) & (f > 1e-4))[0]                 # fitted Gaussians
    def expand(p):
        q = np.interp(np.arange(len(s)), free, p)                    # ties: inner -> first free, outer -> last free
        return q
    k = (Rm >= a.rmin) & (Rm <= a.rmax)
    res = lambda p: (iso_q(f, s, expand(p), Rm[k])-qm[k])/em[k]
    sol = least_squares(res, np.full(len(free), 0.88), bounds=(a.qlo, 1.0), x_scale=0.05)
    qk = expand(sol.x)
    J = sol.jac; act = (sol.x <= a.qlo+1e-6) | (sol.x >= 1-1e-6)
    err = np.full(len(free), np.nan)
    if (~act).any():
        Ji = J[:, ~act]; err[~act] = np.sqrt(np.diag(np.linalg.inv(Ji.T@Ji)))
    chi2 = float(np.sum(sol.fun**2))
    # comparison models
    Rg = Rm
    q_match, q_glob = iso_q(f, s, qk, Rg), iso_q(f, s, np.full(len(s), 0.891), Rg)
    P = H.QProfile.from_json(ROOT/SHAPE, qin=1.0); q_old = iso_q(f, s, H.qvals(P, s), Rg)
    c2 = lambda q: float(np.sum(((q-qm)/em)[k]**2))
    inc = np.array([44., 48., 50., 55.])
    qint = {f"{i:g}": [float(np.sqrt(max(x**2-np.cos(np.radians(i))**2, 0))/np.sin(np.radians(i))) for x in qk] for i in inc}
    out = dict(src=SHAPE, rmin=a.rmin, rmax=a.rmax, qlo=a.qlo, tie_in=a.tie_in, sigma_arcsec=s.tolist(), light_frac=f.tolist(),
               q_k=qk.tolist(), free_index=free.tolist(), q_free=sol.x.tolist(), q_free_err=err.tolist(),
               at_bound=act.tolist(), chi2=chi2, npts=int(k.sum()),
               chi2_same_range=dict(matched=c2(q_match), global_0891=c2(q_glob), qprofile_old=c2(q_old)),
               q_intrinsic=qint, R=Rg.tolist(), q_iso_matched=q_match.tolist(), q_iso_global=q_glob.tolist(),
               q_iso_qprofile_old=q_old.tolist())
    tag = f"r{a.rmax:g}"
    (ROOT/f"results/shape/q_matched_{tag}.json").write_text(json.dumps(out, indent=1))

    fig, ax = plt.subplots(2, 1, figsize=(8, 7.5), sharex=True, gridspec_kw=dict(height_ratios=(2, 1)))
    ax[0].fill_between(Rm, qm-em, qm+em, color="k", alpha=.15, label="measured q_iso (star counts) ±1σ")
    ax[0].plot(Rm, qm, "k-", lw=1)
    for q, lab, c in ((q_match, f"matched q'_k (chi2 {c2(q_match):.1f})", "C3"), (q_glob, f"global 0.891 ({c2(q_glob):.0f})", "C0"),
                      (q_old, f"old QProfile ({c2(q_old):.0f})", "C2")):
        ax[0].plot(Rg, q, "-", c=c, lw=1.5, label=lab)
        ax[1].plot(Rg, (q-qm)/em, "-", c=c, lw=1.2)
    for sk, qq, fk in zip(s, qk, f):
        if fk > 1e-4:
            ax[0].plot(sk/np.sqrt(qq), qq, "v", c="C3", ms=4+12*fk**.5)
    ax[0].axvspan(a.rmin, a.rmax, color="C1", alpha=.06, label=f"matched range {a.rmin:g}-{a.rmax:g}\"")
    ax[0].axhline(np.cos(np.radians(44)), c="r", ls=":", lw=.8, label="cos 44°")
    ax[0].set_ylabel("isophotal q'"); ax[0].legend(fontsize=7); ax[0].set_ylim(0.7, 1.1); ax[0].grid(alpha=.3)
    ax[0].set_title("JAM stellar MGE isophotes vs measured; triangles: q'_k at sigma_major (size ~ light)", fontsize=9)
    ax[1].axhline(0, c="k", lw=.5); ax[1].set_ylim(-9, 9); ax[1].set_ylabel("(model-meas)/err")
    ax[1].set_xscale("log"); ax[1].set_xlabel("R (arcsec, major axis)"); ax[1].grid(alpha=.3)
    plt.tight_layout(); p = ROOT/f"plots/q_matched_{tag}_20261115.png"; plt.savefig(p, dpi=100)
    print("sigma\"", np.round(s, 1).tolist()); print("light", np.round(f, 4).tolist())
    print("q'_k", np.round(qk, 3).tolist(), "free err", np.round(err, 3).tolist(), "at bound", act.tolist())
    print("chi2 (same range, n=%d):" % k.sum(), {kk: round(v, 1) for kk, v in out["chi2_same_range"].items()})
    print("intrinsic q:", {kk: np.round(v, 2).tolist() for kk, v in qint.items()}); print(p)


if __name__ == "__main__":
    main()
