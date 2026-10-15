#!/usr/bin/env python3
"""Does AGAMA's spherical action finder bias point-mass DF models? Compare with exact quadrature actions.

For each case (stellar DF + remnants + central point mass), two comparisons:
  frozen   same potential and DF (built self-consistently with AGAMA actions); only the action
           calculation in the moment table changes. Reports the shift of every kinematic bin in units
           of its error, D^2 = sum (dsigma/err)^2 + sum (dmu)^2/mu over count bins, and the change of
           chi2_kin and of the count deviances.
  rebuilt  the whole self-consistent model rebuilt with numerics.exact_actions = True.
Cases: base best fit + M_bh 8,200 and 4e4 Msun (nothing refitted), and one or more fit starts given as
BATCH::INDEX (the AGAMA build of a start may fail; the exact rebuild is then reported alone).
Output: results/plot_data/compare_df_actions[_starts].json (--starts-only skips the base cases;
--only LABEL runs one case, e.g. base+BH8200 or dfbh8200_20261015::2, into its own file).

Usage: python bin/compare_df_actions.py results/df/dftwo_counts_area_20261013 observed_no_halo_start3 \
           [results/df/dfbh8200_20261015::2 ...]
"""
from __future__ import annotations

import json
import os
import sys
import time
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
from ocen_dm.kinematics.positive_df import DFConvergenceError
from ocen_dm.kinematics.regularized_df import ExactSphericalActions
from ocen_dm.kinematics.run_io import read_data_snapshot


def summary(ev):
    return dict(chi2_kin=ev["chi2_kinematic"], deviance_counts=ev["deviance_counts"],
                objective=ev["objective"], deviances={k: v["deviance"] for k, v in ev["counts"].items()})


def shifts(prob, ev_a, ev_b):
    kin = {p.name: (ev_b["predictions"][p.name]-ev_a["predictions"][p.name])/(.5*(p.err_lo+p.err_hi))
           for p in prob.data.profiles}
    d2_kin = float(sum(np.sum(x**2) for x in kin.values()))
    d2_counts = float(sum(np.sum((ev_b["counts"][k]["mu"]-ev_a["counts"][k]["mu"])**2/ev_a["counts"][k]["mu"])
                          for k in ev_a["counts"]))
    return dict(max_abs_shift_sigma=float(max(np.max(np.abs(x)) for x in kin.values())),
                D2_kin=d2_kin, D2_counts=d2_counts,
                delta_chi2_kin=ev_b["chi2_kinematic"]-ev_a["chi2_kinematic"],
                delta_deviance_counts=ev_b["deviance_counts"]-ev_a["deviance_counts"],
                delta_objective=ev_b["objective"]-ev_a["objective"])


def run_case(prob, cfg, label):
    rec = dict(label=label, M_bh=cfg.matter.M_bh)
    t = time.time()
    try:
        m = build_df_model(cfg)
        ev_a = prob.evaluate(cfg, m)
        rec["agama"] = dict(build_s=time.time()-t, diagnostics={k: v for k, v in m.diagnostics.items()
                                                               if isinstance(v, (int, float))}, **summary(ev_a))
    except DFConvergenceError as e:
        m, ev_a = None, None
        rec["agama"] = dict(build_s=time.time()-t, failed=repr(e)[:300])
    print(label, "AGAMA:", json.dumps(rec["agama"])[:600], flush=True)
    if m is not None:
        # frozen potential and DF: swap only the action finder of the moment table
        table_a = m.moment_table
        m.af = ExactSphericalActions(m.potential)
        for key in ("moment_table", "_splines"):
            m.__dict__.pop(key, None)
        t = time.time()
        table_x = m.moment_table
        rec["frozen_moment_table_s"] = time.time()-t
        ev_x = prob.evaluate(cfg, m)
        inner = table_a["r"] < 100.
        rec["frozen"] = dict(**shifts(prob, ev_a, ev_x), **summary(ev_x),
                             max_rel_moment_change_r_lt_100pc={
                                 k: float(np.max(np.abs(table_x[k][inner]/table_a[k][inner]-1)))
                                 for k in ("rho", "radial_pressure", "tangential_pressure")})
        print(label, "frozen:", json.dumps(rec["frozen"])[:700], flush=True)
    cfg_x = replace(cfg, numerics=replace(cfg.numerics, exact_actions=True))
    t = time.time()
    try:
        mx = build_df_model(cfg_x)
        ev_r = prob.evaluate(cfg_x, mx)
        rec["rebuilt"] = dict(build_s=time.time()-t, diagnostics={k: v for k, v in mx.diagnostics.items()
                                                                 if isinstance(v, (int, float))}, **summary(ev_r))
        if ev_a is not None:
            rec["rebuilt"].update(shifts(prob, ev_a, ev_r))
    except DFConvergenceError as e:
        rec["rebuilt"] = dict(build_s=time.time()-t, failed=repr(e)[:300])
    print(label, "rebuilt:", json.dumps(rec["rebuilt"])[:700], flush=True)
    return rec


def main():
    batch, job = Path(sys.argv[1]), sys.argv[2]
    s = json.loads((batch/"fits"/job/"summary.json").read_text())
    base = refined_config(model_config_from_dict(s["best"]["config"]))
    od = batch/"observed"
    rec = json.loads((od/"problem.json").read_text())
    counts = [CountProfile.from_dict(json.loads((od/f"{n}.json").read_text())) for n in rec["counts"]]
    prob = DFJointProblem(read_data_snapshot(od, rec["data_snapshot"]), None, counts)
    only = sys.argv[sys.argv.index("--only")+1] if "--only" in sys.argv else None
    starts = [a for a in sys.argv[3:] if "::" in a and a != only]
    cases = [] if "--starts-only" in sys.argv else [
        (replace(base, matter=replace(base.matter, M_bh=mbh)), f"base+BH{mbh:g}") for mbh in (8200., 4e4)]
    for spec in starts:
        b, i = spec.split("::")
        start = json.loads((Path(b)/"batch.json").read_text())["jobs"][int(i)]
        cases.append((refined_config(model_config_from_dict(start["config"])), f"{Path(b).name}::{i}"))
    out = [run_case(prob, cfg, label) for cfg, label in cases if only in (None, label)]
    tag = ("_"+only.replace("::", "_").replace("+", "_")) if only else ("_starts" if "--starts-only" in sys.argv else "")
    path = ROOT/f"results/plot_data/compare_df_actions{tag}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    print(path)


if __name__ == "__main__":
    main()
