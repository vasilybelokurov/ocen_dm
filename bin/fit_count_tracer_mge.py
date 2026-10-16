#!/usr/bin/env python3
"""Fit a multi-Gaussian (MGE) tracer profile to the HST + Gaia star counts, for the Jeans models.

Data: ocen_counts_hst_f625w19_c (2-250") and ocen_counts_gaia_g17_c (300-2400"), each with its per-bin detection
completeness, and the amplitude tie ocen_counts_tie_c.json (complete-sample ratio N(F625W<19)/N(G<17) with a
Gaussian prior), all from bin/build_count_completeness.py. Likelihood: kinematics.counts.tied_compare (joint Poisson
deviance of both products + the ratio prior; a common amplitude, the ratio and the Gaia uniform field profiled).

Model: Sigma(R) = sum_k f_k exp(-R^2 / 2 s_k^2) / (2 pi s_k^2), with N_GAUSS widths s_k fixed on a log grid and the
fractions f_k >= 0 (softmax of free logits). Fitted by L-BFGS-B from several starts; the best is kept.
The fit is in arcsec only, so it does not depend on the distance.

Output: data/processed/kinematics/tracer_mge_counts_c.json (MGEFit fields + fit record),
plots/tracer_mge_counts_c.png (counts with the model, residuals, and the composite HST/Trager MGE of the Jeans
reference model for comparison).
Usage: PYTHONPATH=src python bin/fit_count_tracer_mge.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from ocen_dm.kinematics.counts import CountProfile, CountTie, tied_compare
from ocen_dm.light_model import MGEFit
from ocen_dm.paths import processed_dir

KIN = processed_dir()/"kinematics"
N_GAUSS = 12
SIGMAS = np.geomspace(1.0, 1500., N_GAUSS)
REFERENCE_JEANS = ROOT/"results/launches/rung2_turnover_20260921/rung2_K1_turnover.json"


def load():
    num = CountProfile.from_dict(json.loads((KIN/"ocen_counts_hst_f625w19_c.json").read_text()))
    den = CountProfile.from_dict(json.loads((KIN/"ocen_counts_gaia_g17_c.json").read_text()))
    tie = CountTie.from_dict(json.loads((KIN/"ocen_counts_tie_c.json").read_text()))
    return num, den, tie


def fractions(z):
    w = np.exp(z-z.max())
    return w/w.sum()


def sigma_at(fr, R):
    R = np.asarray(R, float)[..., None]
    return np.sum(fr/(2*np.pi*SIGMAS**2)*np.exp(-0.5*(R/SIGMAS)**2), axis=-1)


def objective(z, num, den, tie):
    fr = fractions(z)
    try:
        tn, td, tt = tied_compare(num, den, sigma_at(fr, num.r_nodes), sigma_at(fr, den.r_nodes), tie)
    except ValueError:
        return 1e12
    return tn["deviance"]+td["deviance"]+tt["chi2"]


def main():
    num, den, tie = load()
    rng = np.random.default_rng(2026)
    starts = [np.zeros(N_GAUSS)]+[rng.normal(0, 2, N_GAUSS) for _ in range(7)]
    ref = json.loads(REFERENCE_JEANS.read_text())["model"]["mge_fit"]
    ref_fit = MGEFit(np.array(ref["sigmas_arcsec"]), np.array(ref["fractions"]), ref["rms_mag"], ref["max_abs_resid_mag"], ref["n_grid"])
    z0 = np.log(np.maximum([ref_fit.surface_intensity(s)[0]*s**2 for s in SIGMAS], 1e-12))   # start near the composite shape
    starts.append(z0-z0.max())
    results = []
    for z in starts:
        r = minimize(objective, z, args=(num, den, tie), method="L-BFGS-B", bounds=[(-25, 5)]*N_GAUSS,
                     options=dict(maxiter=4000, maxfun=200000))
        results.append(r)
        print(f"start: objective {r.fun:.3f} ({r.nit} it, {r.message})", flush=True)
    best = min(results, key=lambda r: r.fun)
    fr = fractions(best.x)
    tn, td, tt = tied_compare(num, den, sigma_at(fr, num.r_nodes), sigma_at(fr, den.r_nodes), tie)
    keep = fr > 1e-8
    mge = MGEFit(SIGMAS[keep], fr[keep]/fr[keep].sum(), 0.0, 0.0, 64)
    record = dict(sigmas_arcsec=mge.sigmas_arcsec.tolist(), fractions=mge.fractions.tolist(), rms_mag=0.0,
                  max_abs_resid_mag=0.0, n_grid=64,
                  fit=dict(created_utc=datetime.now(timezone.utc).isoformat(), n_gauss=N_GAUSS, sigma_grid_arcsec=SIGMAS.tolist(),
                           objective=float(best.fun), objectives_all_starts=[float(r.fun) for r in results],
                           deviance_hst=tn["deviance"], deviance_gaia=td["deviance"], n_hst=num.n, n_gaia=den.n,
                           tie_ratio=tt["ratio"], tie_chi2=tt["chi2"], tie_prior=[tie.ratio, tie.sigma],
                           amplitude_hst=tn["amplitude"], amplitude_gaia=td["amplitude"],
                           gaia_field_per_arcmin2=td["field_density_per_arcmin2"],
                           inputs=["ocen_counts_hst_f625w19_c", "ocen_counts_gaia_g17_c", "ocen_counts_tie_c"],
                           note="rms_mag/max_abs_resid_mag are not used (fit is Poisson on counts); see deviances"))
    (KIN/"tracer_mge_counts_c.json").write_text(json.dumps(record, indent=1)+"\n")

    # figure: data (completeness-corrected, field-subtracted) vs model, and the composite MGE scaled to the counts
    fig, ax = plt.subplots(2, 1, figsize=(6.0, 6.0), sharex=True, gridspec_kw=dict(hspace=0.05, height_ratios=[2, 1]))
    Rg = np.geomspace(1, 3000, 400)
    a_h = tn["amplitude"]
    def corrected(p, t, field):
        c = np.ones(p.n) if p.completeness is None else np.asarray(p.completeness)
        d = (p.counts/p.area_arcsec2-field)/c*3600.
        e = np.sqrt(p.counts)/p.area_arcsec2/c*3600.
        return np.sqrt(p.r_lower*p.r_upper), d, e
    for p, t, fld, scale, col, lab in ((num, tn, 0., 1., "C0", "HST F625W < 19"),
                                       (den, td, td["field"], tt["ratio"], "C3", r"Gaia $G<17$ $\times$ ratio")):
        R, d, e = corrected(p, t, fld)
        ax[0].errorbar(R, d*scale, e*scale, fmt="o", ms=3.5, color=col, label=lab+" (completeness-corrected)")
        model = a_h*p.bin_average(sigma_at(fr, p.r_nodes))*3600.
        ax[1].errorbar(R, d*scale/model-1, e*scale/model, fmt="o", ms=3.5, color=col)
        ref_model = a_h*p.bin_average(ref_fit.surface_intensity(p.r_nodes.ravel()).reshape(p.r_nodes.shape))*3600.
        ax[1].plot(R, ref_model/model-1, "--", color=col, lw=0.9)
    ref_scale = np.sum(num.counts)/np.sum(num.area_arcsec2*num.bin_average(ref_fit.surface_intensity(num.r_nodes.ravel()).reshape(num.r_nodes.shape))*np.asarray(num.completeness))
    ax[0].plot(Rg, a_h*sigma_at(fr, Rg)*3600, "k-", lw=1.1, label="MGE fitted to the counts")
    ax[0].plot(Rg, ref_scale*ref_fit.surface_intensity(Rg)*3600, "k--", lw=0.9, label="composite HST/Trager MGE (Jeans reference), scaled to HST")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_ylabel(r"$\Sigma$ [stars arcmin$^{-2}$, $F625W<19$]")
    ax[0].legend(fontsize=6.8, loc="lower left")
    ax[1].axhline(0, color="k", lw=0.5); ax[1].set_ylim(-0.4, 0.4); ax[1].set_ylabel("data / MGE − 1")
    ax[1].set_xlabel("R [arcsec]"); ax[1].text(0.02, 0.06, "dashed: composite MGE / count MGE − 1 (same amplitude)", transform=ax[1].transAxes, fontsize=6.8)
    fig.savefig(ROOT/"plots/tracer_mge_counts_c.png", dpi=150, bbox_inches="tight")
    print(json.dumps({k: v for k, v in record["fit"].items() if k not in ("sigma_grid_arcsec",)}, indent=1))
    print("fractions:", np.round(fr, 5).tolist())


if __name__ == "__main__":
    main()
