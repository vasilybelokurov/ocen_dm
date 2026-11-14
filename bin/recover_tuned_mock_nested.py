#!/usr/bin/env python3
"""Step 5 with the production sampler (2026-11-13): UltraNest fit of one spherical Jeans family to a noiseless tuned
flattened JAM mock.

The local optimiser in recover_tuned_mock.py gave start-dependent likelihoods on the flattened mocks, so the fits
are redone with the same nested sampler and settings as the production runs (fit.run_nested: 400 live points,
dlogz 0.5, slice sampler with the production walk length). Mock data: the saved tuned prediction (noiseless) with
the corrected real errors and streaming terms, built exactly as in recover_tuned_mock.py.

Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/recover_tuned_mock_nested.py \
           --truth q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1 --family K1
       --truth real fits the real data with the corrected HST rotation terms (same data as the JAM tuning).
Output: results/fits/recover_<truth>_<family>/ (summary.json, posterior.ecsv, ...)
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np
import yaml

from ocen_dm.kinematics import FitProblem, KinematicData
from ocen_dm.kinematics.fit import run_nested
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config
from tune_jam_to_data import Setup

FAMILIES = {"K1": "jeans_K1_turnover_counts_c", "NFW": "jeans_K2_nfw_turnover_counts_c"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="jam_tuned tag, or 'real' for the corrected real data")
    ap.add_argument("--family", choices=sorted(FAMILIES), required=True)
    args = ap.parse_args()
    S = Setup(per_bin=10)                     # corrected real errors and streaming terms; values replaced below
    if args.truth == "real":                  # the real data with the corrected HST rotation terms (no mock)
        mdata = S.data
        prov = {"kind": "real", "hst_rotation": "results/rotation/rotation_field_cylindrical_20261113.json (MAP), "
                "ring <m^2>; HST-T errors with the field-chain spread (tune_jam_to_data.Setup._corrected_data)",
                "script": "bin/recover_tuned_mock_nested.py --truth real"}
    else:
        T = json.loads((ROOT/f"results/rotation/jam_tuned_{args.truth}.json").read_text())
        mdata = KinematicData(tuple(replace(p, value=np.array(T["pred"][p.name])) for p in S.data.profiles))
        prov = {"kind": "mock", "truth": f"results/rotation/jam_tuned_{args.truth}.json",
                "noise": "none", "script": "bin/recover_tuned_mock_nested.py"}
    prod = ROOT/"results/fits"/FAMILIES[args.family]
    fam = family_from_config(family_config(problem_for(load_run_metadata(prod)).family))
    smp = yaml.safe_load((prod/"run.yaml").read_text())["sampler"]
    out = ROOT/f"results/fits/recover_{args.truth}_{args.family}"
    s = run_nested(FitProblem(fam, mdata), out, n_live=smp["n_live"], dlogz=smp["dlogz"], seed=smp["seed"],
                   step_sampler=smp["step_sampler"], slice_nsteps=smp["slice_nsteps"],
                   data_provenance=prov)
    print(f"{args.truth} {args.family}: lnL_max {s['lnL_max']:.2f}  logZ {s['logz']:.2f}+-{s['logzerr']:.2f}  "
          f"chi2 {s['chi2_ml_total']}  ({s['elapsed_s']/60:.0f} min)", flush=True)


if __name__ == "__main__":
    main()
