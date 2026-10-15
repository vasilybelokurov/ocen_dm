#!/usr/bin/env python3
"""Model-level validation of the central point mass in the regularized DF, at the base best-fit parameters.

Checks (nothing is refitted):
  1. M_bh = 0 reproduces the saved base fit exactly (chi2_kin, count deviances).
  2. For each M_bh: equilibrium convergence and closure; central isotropy; the inner stellar
     density slope against the isotropic Kepler prediction rho ~ r^-3/2 for a DF finite at J -> 0;
     projected moments at the HST/MUSE data nodes against AGAMA's direct GalaxyModel moments;
     sensitivity of the bin predictions to velocity_nodes (32 -> 64) and r_min (1e-4 -> 1e-5 pc).
Output: results/plot_data/validate_df_point_mass.json.

Usage: python bin/validate_df_point_mass.py results/df/dftwo_counts_area_20261013 observed_no_halo_start3
"""
from __future__ import annotations

import json
import os
import sys
from dataclasses import replace
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np

from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.counts import CountProfile
from ocen_dm.kinematics.df_fit import DFJointProblem, build_df_model, model_config_from_dict
from ocen_dm.kinematics.positive_df import agama_pc
from ocen_dm.kinematics.run_io import read_data_snapshot


def main():
    batch, job = Path(sys.argv[1]), sys.argv[2]
    s = json.loads((batch/"fits"/job/"summary.json").read_text())
    report = json.loads((batch/"report.json").read_text())
    saved = [j for j in report["jobs"] if j["job"] == job][0]
    saved = saved.get("gates") or saved["status"]["gates"]
    base = refined_config(model_config_from_dict(s["best"]["config"]))
    od = batch/"observed"
    rec = json.loads((od/"problem.json").read_text())
    counts = [CountProfile.from_dict(json.loads((od/f"{n}.json").read_text())) for n in rec["counts"]]
    prob = DFJointProblem(read_data_snapshot(od, rec["data_snapshot"]), None, counts)
    pc = base.distance_kpc*1e3/206264.806
    G = agama_pc().G
    out = {}

    def evaluate(cfg):
        m = build_df_model(cfg)
        return m, prob.evaluate(cfg, m)

    m0, ev0 = evaluate(base)
    out["mbh0_reproduction"] = dict(
        chi2_kin=ev0["chi2_kinematic"], saved_chi2_kin=saved["chi2_kinematic"],
        dev_hst=ev0["counts"]["ocen_counts_hst_f625w19"]["deviance"], saved_dev_hst=saved["deviance_counts"]["ocen_counts_hst_f625w19"],
        dev_gaia=ev0["counts"]["ocen_counts_gaia_g17"]["deviance"], saved_dev_gaia=saved["deviance_counts"]["ocen_counts_gaia_g17"])
    print("M_bh = 0:", out["mbh0_reproduction"], flush=True)

    for mbh in (8200., 4.0e4):
        cfg = replace(base, matter=replace(base.matter, M_bh=mbh))
        m, ev = evaluate(cfg)
        sigma0 = np.sqrt(m.intrinsic_moments(np.array([0.5]))["radial_pressure"][0]/m.intrinsic_moments(np.array([0.5]))["rho"][0])
        r_infl = G*mbh/sigma0**2
        r_lo = 1.5*m.moment_table["r"][0]
        r = np.geomspace(max(1e-3*r_infl, r_lo), r_infl, 7)
        mom = m.intrinsic_moments(r)
        slope = np.gradient(np.log(mom["rho"]), np.log(r))
        rec_m = dict(r_infl_pc=r_infl, diagnostics={k: v for k, v in m.diagnostics.items() if isinstance(v, (int, float))},
                     chi2_kin=ev["chi2_kinematic"], dev_hst=ev["counts"]["ocen_counts_hst_f625w19"]["deviance"],
                     dev_gaia=ev["counts"]["ocen_counts_gaia_g17"]["deviance"],
                     inner_r_pc=r.tolist(), inner_slope=slope.tolist(), inner_beta=mom["beta"].tolist())
        # projected moments vs AGAMA direct, at HST/MUSE nodes inside 2 pc and a few outside
        R = np.unique(np.concatenate([p.r*pc for p in prob.data.profiles[:3]]))
        R = R[R < 10.]
        pm = m.projected_moments(R)
        dm = m.direct_projected_moments(R)
        rel = {k: float(np.max(np.abs(np.sqrt(pm[k]/pm["Sigma"])/np.sqrt(dm[k]/dm["Sigma"])-1))) for k in ("los", "pmr", "pmt")}
        rel["Sigma"] = float(np.max(np.abs(pm["Sigma"]/dm["Sigma"]-1)))
        rec_m["projected_vs_direct_max_rel"] = rel
        # numerical sensitivity of the bin predictions
        for label, num in (("velocity_nodes_64", replace(cfg.numerics, velocity_nodes=2*cfg.numerics.velocity_nodes)),
                           ("r_min_1e-5", replace(cfg.numerics, r_min=cfg.numerics.r_min/10))):
            m2, ev2 = evaluate(replace(cfg, numerics=num))
            shift = max(float(np.max(np.abs(ev2["predictions"][p.name]-ev["predictions"][p.name])/(.5*(p.err_lo+p.err_hi))))
                        for p in prob.data.profiles)
            rec_m[f"max_shift_sigma_{label}"] = shift
            rec_m[f"chi2_kin_{label}"] = ev2["chi2_kinematic"]
        out[f"mbh_{mbh:g}"] = rec_m
        print(f"M_bh = {mbh:g}:", json.dumps(rec_m)[:1500], flush=True)
    path = ROOT/"results/plot_data/validate_df_point_mass.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    print(path)


if __name__ == "__main__":
    main()
