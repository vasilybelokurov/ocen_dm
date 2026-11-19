#!/usr/bin/env python3
"""Local (restarted Nelder-Mead) MAP fits of the spherical Jeans spline families (2026-11-19; docs/PLAN_2026-11-19_limits.md).

Objective: -[lnL_data + ln N(D; 5.43, 0.05)] in the sampling coordinates (log10 for log-uniform priors, linear for
uniform), with the family's prior box as hard bounds. Parameters can be held fixed (profile likelihood). Data: the
snapshot of the production spline run (results/fits/recover_real_K1_spline_nw: corrected rotation terms, membership-
weighted Gaia nodes), or any KinematicData (mocks). Restarts continue until a round gains < 0.01.

Module use:  L = Local("NFW", inner_node=2.0); res = L.fit(x0, fixed={"M_dm_100": 1e6, "r_s": 300})
Returns dict(x, lnL (data only), lnprior_D, lnpost, chi2 per dataset).
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np
from scipy.optimize import minimize

from ocen_dm.kinematics import FitProblem
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config
from recover_tuned_mock_nested import FAMILIES, spline_config

DATA_RUN = ROOT/"results/fits/recover_real_K1_spline_nw"
D0, SD = 5.43, 0.05


def real_data():
    return problem_for(load_run_metadata(DATA_RUN)).data


def family(kind: str, inner_node: float | None = 2.0):
    cfg = family_config(problem_for(load_run_metadata(ROOT/"results/fits"/FAMILIES[kind])).family)
    return family_from_config(spline_config(cfg, inner_node))


class Local:
    def __init__(self, kind: str, inner_node: float | None = 2.0, data=None, prior_overrides: dict | None = None):
        self.fam = family(kind, inner_node)
        self.P = FitProblem(self.fam, real_data() if data is None else data)
        self.names = list(self.fam.names)
        self.log = np.array([p.prior.kind == "loguniform" for p in self.fam.parameters])
        lo, hi = [], []
        for p in self.fam.parameters:
            a, b = (prior_overrides or {}).get(p.name, (p.prior.lo, p.prior.hi))
            if p.prior.kind == "normal":
                a, b = 4.0, 7.0
            lo.append(np.log10(a) if p.prior.kind == "loguniform" else a)
            hi.append(np.log10(b) if p.prior.kind == "loguniform" else b)
        self.lo, self.hi = np.array(lo), np.array(hi)
        self.iD = self.names.index("distance")

    def z(self, x):
        x = np.asarray(x, float); return np.where(self.log, np.log10(np.maximum(x, 1e-300)), x)

    def x(self, z):
        z = np.asarray(z, float); return np.where(self.log, 10**z, z)

    def terms(self, x):
        lnl = float(np.sum(self.P.loglike_vector(np.asarray(x, float))))
        lpd = -0.5*((x[self.iD]-D0)/SD)**2
        return lnl, lpd

    def fit(self, x0, fixed: dict | None = None, max_rounds: int = 8, tol: float = 0.01, verbose=False):
        fixed = fixed or {}
        fi = np.array([self.names.index(k) for k in fixed], int)
        zfull = self.z(x0).copy()
        for k, v in fixed.items():
            i = self.names.index(k); zfull[i] = np.log10(v) if self.log[i] else v
        free = np.array([i not in fi for i in range(len(self.names))])
        lo, hi = self.lo[free], self.hi[free]

        def f(zf):
            if np.any(zf < lo) or np.any(zf > hi):
                return 1e10
            zz = zfull.copy(); zz[free] = zf
            lnl, lpd = self.terms(self.x(zz))
            return -(lnl + lpd) if np.isfinite(lnl) else 1e10

        zf = np.clip(zfull[free], lo + 1e-9, hi - 1e-9)
        best = f(zf)
        for rnd in range(max_rounds):
            r = minimize(f, zf, method="Nelder-Mead",
                         options=dict(maxiter=6000, xatol=1e-6, fatol=1e-5, adaptive=True))
            gain, zf, best = best - r.fun, r.x, min(best, r.fun)
            if verbose:
                print(f"  round {rnd} obj {best:.4f} gain {gain:.4f} nfev {r.nfev}", flush=True)
            if gain < tol:
                break
        zz = zfull.copy(); zz[free] = zf
        x = self.x(zz)
        lnl, lpd = self.terms(x)
        chi2 = {k: float(v) for k, (v, _) in self.P.chi2(x).items()}
        return dict(x=x, lnL=lnl, lnprior_D=lpd, lnpost=lnl+lpd, chi2=chi2, chi2_total=sum(chi2.values()),
                    theta=dict(zip(self.names, map(float, x))), rounds=rnd+1)
