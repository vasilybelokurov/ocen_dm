#!/usr/bin/env python3
"""Data / model / residual plots and the per-bin halo gain (NFW - K1) for a set of real-data fits (2026-11-15;
docs/JAM_FIT_PLAN.md step 3 and 5).

Sets:
  base : JAM K1/NFW with one global q' (all bins; R <= 1100"), spherical Jeans K1/NFW (UltraNest ML, results/fits/
         recover_real_{K1,NFW}_base) -- identical data (one rotation field, all stars per bin).
  qpvar: q'(R) variants on R <= 1100" (remnants round, outer round, i = 48).
  qp   : JAM K1/NFW with the measured q'(R) profile (all bins; R <= 1100"), with the base JAM fits for comparison.
Per-bin halo gain: dlnL_i = lnL_i(NFW) - lnL_i(K1) (split-normal per bin, as the fit), summed per dataset and for
R <= 1100" / R > 1100". Radii, values and errors from tune_jam_to_data.Setup (checked against each fit's data_used).

Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/report_jam_halo_gain.py --set base
Output: plots/jam_<set>_fits_20261115.png, plots/jam_<set>_halo_gain_20261115.png, results/rotation/halo_gain_<set>.json
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

import jam_halo_models as H
from ocen_dm.kinematics import FitProblem
from ocen_dm.kinematics.likelihood import ProfileLikelihood
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config
from tune_jam_to_data import Setup

R = ROOT/"results/rotation"
SHORT = {"hst_pm_radial_ours": "HST-R", "hst_pm_tangential_ours": "HST-T", "muse_los_dispersion": "MUSE",
         "gaia_edr3_ours_radial": "Gaia-R", "gaia_edr3_ours_tangential": "Gaia-T"}
SETS = {
    "base": {"JAM all": ("jamfit_real_K1_lsq_base.json", "jamfit_real_NFW_lsq_base.json"),
             "JAM R<=1100": ("jamfit_real_K1_rmax1100_lsq_base.json", "jamfit_real_NFW_rmax1100_lsq_base.json")},
    "qp": {"JAM q'(R) all": ("jamfit_real_K1_lsq_qprofhold_remstars_qp.json", "jamfit_real_NFW_lsq_qprofhold_remstars_qp.json"),
           "JAM q'(R) R<=1100": ("jamfit_real_K1_rmax1100_lsq_qprofhold_remstars_qp.json",
                                 "jamfit_real_NFW_rmax1100_lsq_qprofhold_remstars_qp.json"),
           "JAM global q' all": ("jamfit_real_K1_lsq_base.json", "jamfit_real_NFW_lsq_base.json"),
           "JAM global q' R<=1100": ("jamfit_real_K1_rmax1100_lsq_base.json", "jamfit_real_NFW_rmax1100_lsq_base.json")},
    "qpvar": {"q'(R) i44": ("jamfit_real_K1_rmax1100_lsq_qprofhold_remstars_qp.json",
                            "jamfit_real_NFW_rmax1100_lsq_qprofhold_remstars_qp.json"),
              "q'(R) remnants round": ("jamfit_real_K1_rmax1100_lsq_qprofhold_remround_qp.json",
                                       "jamfit_real_NFW_rmax1100_lsq_qprofhold_remround_qp.json"),
              "q'(R) outer round": ("jamfit_real_K1_rmax1100_lsq_qprofround_remstars_qp.json",
                                    "jamfit_real_NFW_rmax1100_lsq_qprofround_remstars_qp.json"),
              "q'(R) i48": ("jamfit_real_K1_rmax1100_lsq_qprofhold_remstars_i48_qp.json",
                            "jamfit_real_NFW_rmax1100_lsq_qprofhold_remstars_i48_qp.json")},
    "qm": {"global q' turnover i44 (base)": ("jamfit_real_K1_lsq_base.json", None),
           "global q' spline i44": ("jamfit_real_K1_lsq_bspline_qm.json", None),
           "matched turnover i44": ("jamfit_real_K1_lsq_qmatch_remstars_qm.json", None),
           "matched spline i44": ("jamfit_real_K1_lsq_qmatch_remstars_bspline_qm.json", None),
           "matched spline i50": ("jamfit_real_K1_lsq_qmatch_remstars_i50_bspline_qm.json", None),
           "matched spline i55": ("jamfit_real_K1_lsq_qmatch_remstars_i55_bspline_qm.json", None)},
    "basepair": {"no halo (K1), global q' i44": ("jamfit_real_K1_lsq_base.json", None),
                 "NFW halo, global q' i44": ("jamfit_real_NFW_lsq_base.json", None)},
    "qmnfw": {"K1 matched spline i44": ("jamfit_real_K1_lsq_qmatch_remstars_bspline_qm.json", None),
              "NFW matched spline i44": ("jamfit_real_NFW_lsq_qmatch_remstars_bspline_qm.json", None),
              "K1 matched spline i50": ("jamfit_real_K1_lsq_qmatch_remstars_i50_bspline_qm.json", None),
              "NFW matched spline i50": ("jamfit_real_NFW_lsq_qmatch_remstars_i50_bspline_qm.json", None)},
}
SPH = {"K1": ("jeans_K1_turnover_counts_c", "recover_real_K1_base"), "NFW": ("jeans_K2_nfw_turnover_counts_c", "recover_real_NFW_base")}


def spherical(data):
    out = {}
    for k, (run, d) in SPH.items():
        fam = family_from_config(family_config(problem_for(load_run_metadata(ROOT/"results/fits"/run)).family))
        pred = FitProblem(fam, data).predict(np.load(ROOT/"results/fits"/d/"ml_x.npy"))
        out[k] = dict(pred={n: np.asarray(v) for n, v in pred.items()},
                      summary=json.loads((ROOT/"results/fits"/d/"summary.json").read_text()))
    return out


def lnl_bins(pred, p):
    return np.asarray(ProfileLikelihood._split_normal_lnlike(np.asarray(pred), p))


def k1_only(name, profs):
    """No-halo (K1) fits only: data/model/residual per dataset, one colour per fit; chi2 table."""
    fits = {lab: json.loads((R/fk).read_text()) for lab, (fk, _) in SETS[name].items() if (R/fk).exists()}
    fig, ax = plt.subplots(2, len(profs), figsize=(4.2*len(profs), 7.5), sharex="col", gridspec_kw=dict(height_ratios=(2.2, 1)))
    for j, p in enumerate(profs):
        ax[0, j].errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt="ko", ms=3, lw=.8)
        for i, (lab, d) in enumerate(fits.items()):
            m = np.asarray(d["pred"][p.name]); res = (m-p.value)/np.where(m > p.value, p.err_hi, p.err_lo)
            ls = "--" if "turnover" in lab and "base" not in lab else "-"
            ax[0, j].plot(p.r, m, ls, c=f"C{i}", lw=1.3, label=f"{lab}: chi2 {np.sum(res**2):.0f}")
            ax[1, j].plot(p.r, res, ls, c=f"C{i}", lw=1)
        ax[0, j].set_title(p.name, fontsize=9); ax[0, j].legend(fontsize=5.5); ax[0, j].set_xscale("log")
        ax[1, j].axhline(0, c="k", lw=.5); ax[1, j].set_ylim(-6, 6); ax[1, j].set_xlabel("R (arcsec)")
    ax[0, 0].set_ylabel("dispersion (mas/yr or km/s)"); ax[1, 0].set_ylabel("(model-data)/err")
    fig.suptitle(f"set '{name}': JAM fits, all bins" + ("" if name in ("basepair", "qmnfw") else " (K1, no DM)"), fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, .96)); f1 = ROOT/f"plots/jam_{name}_fits_20261116.png"; plt.savefig(f1, dpi=90); plt.close()
    for lab, d in fits.items():
        print(f"{lab:32s} lnpost {-d['neg_lnpost']:7.1f} chi2 {d['chi2_total']:6.1f} | "
              + " ".join(f"{SHORT[k]} {v['chi2']:.0f}" for k, v in d["stats"].items()))
    print(f1)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--set", choices=sorted(SETS), default="base"); a = ap.parse_args()
    S = Setup(run="jeans_K1_turnover_counts_c", models=H.jam_models)
    profs = S.data.profiles
    if all(fn is None for _, fn in SETS[a.set].values()):
        return k1_only(a.set, profs)
    pairs = {}
    for lab, (fk, fn) in SETS[a.set].items():
        pk, pn = (R/fk), (R/fn)
        if not (pk.exists() and pn.exists()):
            print("missing", lab); continue
        dk, dn = json.loads(pk.read_text()), json.loads(pn.read_text())
        for d in (dk, dn):                                   # same data as the fit?
            for p in profs:
                du = d["data_used"][p.name]
                assert np.allclose(du["value"], p.value) and np.allclose(du["err_lo"], p.err_lo), (lab, p.name)
        pairs[lab] = (dk, dn)
    sph = spherical(S.data) if a.set == "base" else {}
    if sph:
        pairs["spherical (all)"] = ({"pred": sph["K1"]["pred"], "args": {"rmax": None}, "sph": sph["K1"]["summary"]},
                                    {"pred": sph["NFW"]["pred"], "args": {"rmax": None}, "sph": sph["NFW"]["summary"]})

    # 1. data / model / residuals (all-bin fits, K1 solid/dotted, NFW)
    fig, ax = plt.subplots(2, len(profs), figsize=(4.2*len(profs), 7.5), sharex="col", gridspec_kw=dict(height_ratios=(2.2, 1)))
    sty = {0: "-", 1: "--", 2: ":", 3: "-."}
    for j, p in enumerate(profs):
        a0, a1 = ax[0, j], ax[1, j]
        a0.errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt="ko", ms=3, lw=.8)
        for i, (lab, (dk, dn)) in enumerate(pairs.items()):
            for d, c, fam in ((dk, "C0", "K1"), (dn, "C3", "NFW")):
                m = np.asarray(d["pred"][p.name]); res = (m-p.value)/np.where(m > p.value, p.err_hi, p.err_lo)
                rm = d["args"].get("rmax")
                k = np.ones(len(m), bool) if rm is None else np.asarray(p.r) <= rm
                a0.plot(np.asarray(p.r)[k], m[k], sty[i % 4], c=c, lw=1.2, label=f"{lab} {fam}: chi2 {np.sum(res[k]**2):.1f}")
                a1.plot(np.asarray(p.r)[k], res[k], sty[i % 4], c=c, lw=1)
        a0.set_title(p.name, fontsize=9); a0.legend(fontsize=5.5); a0.set_xscale("log")
        a1.axhline(0, c="k", lw=.5); a1.set_ylim(-4, 4); a1.set_xlabel("R (arcsec)")
    ax[0, 0].set_ylabel("dispersion (mas/yr or km/s)"); ax[1, 0].set_ylabel("(model-data)/err")
    fig.suptitle(f"set '{a.set}': blue K1 (no DM), red NFW; line style = fit (see legends)", fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, .96)); f1 = ROOT/f"plots/jam_{a.set}_fits_20261115.png"; plt.savefig(f1, dpi=90); plt.close()

    # 2. per-bin halo gain
    fig, ax = plt.subplots(1, len(profs), figsize=(4.2*len(profs), 3.8), sharey=True)
    table = {}
    for i, (lab, (dk, dn)) in enumerate(pairs.items()):
        rm = dk["args"].get("rmax"); row = {}
        tot_in = tot_out = 0.
        for j, p in enumerate(profs):
            g = lnl_bins(dn["pred"][p.name], p)-lnl_bins(dk["pred"][p.name], p)
            r = np.asarray(p.r); k = np.ones(len(r), bool) if rm is None else r <= rm
            ax[j].plot(r[k], g[k], sty[i % 4]+"o", ms=3, lw=1, label=lab, c=f"C{i}")
            gi, go = float(g[k & (r <= 1100)].sum()), float(g[k & (r > 1100)].sum())
            row[p.name] = dict(total=round(gi+go, 3), R_le_1100=round(gi, 3), R_gt_1100=round(go, 3))
            tot_in += gi; tot_out += go
        row["sum"] = dict(total=round(tot_in+tot_out, 3), R_le_1100=round(tot_in, 3), R_gt_1100=round(tot_out, 3))
        table[lab] = row
    for j, p in enumerate(profs):
        ax[j].axhline(0, c="k", lw=.5); ax[j].axvline(1100, c="gray", ls=":", lw=.8); ax[j].set_xscale("log")
        ax[j].set_title(p.name, fontsize=9); ax[j].set_xlabel("R (arcsec)")
    ax[0].set_ylabel("lnL(NFW) - lnL(K1) per bin"); ax[0].legend(fontsize=7)
    plt.tight_layout(); f2 = ROOT/f"plots/jam_{a.set}_halo_gain_20261115.png"; plt.savefig(f2, dpi=90); plt.close()
    (R/f"halo_gain_{a.set}.json").write_text(json.dumps(table, indent=1))
    names = [p.name for p in profs]
    print("halo gain dlnL (NFW - K1): total [R<=1100 | R>1100]")
    for lab, row in table.items():
        print(f"{lab:22s} " + "  ".join(f"{SHORT.get(n, n)} {row[n]['total']:+.2f}" for n in names)
              + f"  SUM {row['sum']['total']:+.2f} [{row['sum']['R_le_1100']:+.2f} | {row['sum']['R_gt_1100']:+.2f}]")
    print(f1); print(f2)


if __name__ == "__main__":
    main()
