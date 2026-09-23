#!/usr/bin/env python3
"""Measure evaluation noise in compact-DF finite-difference Jacobians.

At a saved fit's best point, compares residual vectors and forward-difference
columns under different stellar-potential seeds and iteration tolerances:

  chained  - each probe seeded by the previous probe (the v2 behaviour)
  base     - every probe seeded by the base point's converged potential
  cold     - every evaluation from the default Plummer seed
  *_tight  - the same with iteration_tolerance divided by 10

Noise is reported as max |r_a - r_b| in adopted errors for repeated base-point
evaluations, and as the spread of Jacobian columns between modes.

Usage: python bin/probe_compact_jacobian_noise.py results/df/compact_recovery_20260923_v2 cored_dm_free_halo_start0
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
import sys
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np

from ocen_dm.kinematics.compact_recovery import residual_vector
from ocen_dm.kinematics.df_fit import (DFJointProblem, FitCoordinate, PhotometricData, config_at,
                                       model_config_from_dict)
from ocen_dm.kinematics.regularized_df import RegularizedDFModel
from ocen_dm.kinematics.run_io import read_data_snapshot

STEP = .002
PROBED = ("M_star", "stellar.J0", "matter.rho20")


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batch", type=Path)
    parser.add_argument("job")
    args = parser.parse_args()
    out = args.batch.resolve()
    job = next(j for j in read(out/"batch.json")["jobs"] if j["id"] == args.job)
    mock = out/"mocks"/job["case"]
    problem = DFJointProblem(read_data_snapshot(mock, read(mock/"mock.json")["data_snapshot"]),
                             PhotometricData.from_dict(read(mock/"photometry.json")))
    coords = [FitCoordinate(**c) for c in job["coordinates"]]
    low, high = np.array([c.bounds for c in coords]).T
    best = read(out/"fits"/args.job/"best.json")
    x0 = np.array(best["x"])
    base = model_config_from_dict(best["config"])
    probed = [i for i, c in enumerate(coords) if c.path in PROBED]

    def evaluate(x, seed, tight):
        config = config_at(base, coords, low+x*(high-low))
        if tight:
            config = replace(config, numerics=replace(config.numerics,
                             iteration_tolerance=config.numerics.iteration_tolerance/10))
        t = time.monotonic()
        model = RegularizedDFModel(config, initial_stellar_potential=seed)
        r = residual_vector(problem, problem.evaluate(config, model))
        return r, model.stellar_potential, time.monotonic()-t

    result = dict(batch=str(out), job=args.job, step=STEP, probed=[coords[i].path for i in probed], modes={})
    for tight in (False, True):
        r_cold, pot_base, sec = evaluate(x0, None, tight)
        r_base, _, _ = evaluate(x0, pot_base, tight)
        columns = {"chained": [], "base": [], "cold": []}
        warm = pot_base
        for i in probed:
            x = x0.copy()
            x[i] += STEP if x[i]+STEP <= 1 else -STEP
            s = x[i]-x0[i]
            r, warm, _ = evaluate(x, warm, tight)
            columns["chained"].append(((r-r_cold)/s).tolist())
            r, _, _ = evaluate(x, pot_base, tight)
            columns["base"].append(((r-r_base)/s).tolist())
            r, _, _ = evaluate(x, None, tight)
            columns["cold"].append(((r-r_cold)/s).tolist())
        name = "tight" if tight else "nominal"
        J = {k: np.array(v) for k, v in columns.items()}
        scale = np.linalg.norm(J["cold"], axis=1)
        result["modes"][name] = dict(
            seconds_per_cold_eval=sec,
            objective_cold=float(r_cold@r_cold), objective_base_reseeded=float(r_base@r_base),
            base_repeat_max_shift_sigma=float(np.max(abs(r_cold-r_base))),
            column_norm_cold=scale.tolist(),
            chained_vs_cold_rel=(np.linalg.norm(J["chained"]-J["cold"], axis=1)/scale).tolist(),
            base_vs_cold_rel=(np.linalg.norm(J["base"]-J["cold"], axis=1)/scale).tolist(),
            columns=columns)
        print(name, {k: v for k, v in result["modes"][name].items() if k != "columns"}, flush=True)
    dest = ROOT/"results/maintenance"/f"jacobian_noise_{out.name}_{args.job}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result))
    print(dest)


if __name__ == "__main__":
    main()
