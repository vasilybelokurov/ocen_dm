#!/usr/bin/env python3
"""Run the independent fixed-potential DF flexibility challenge in a fresh directory."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np
from scipy.optimize import least_squares
from scipy.special import softmax

from ocen_dm.kinematics.df_capacity import (
    CapacityObservations, FixedPotentialProjector, capacity_bounds, PARAMETER_NAMES,
)
from ocen_dm.kinematics.df_mock import challenge_mock
from ocen_dm.kinematics.df_fit import PhotometricData
from ocen_dm.kinematics.positive_df import agama_pc
from ocen_dm.kinematics.run_io import sha256, read_data_snapshot, write_data_snapshot


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")


def starting_vector(count, start, previous, rng):
    """First start nests the previous winner; other starts explore distinct shapes."""
    if previous is not None and len(previous) == 7*(count-1)-1 and start == 0:
        old = count-1
        shapes = np.asarray(previous[:6*old]).reshape(old, 6)
        weights = softmax(np.r_[0., previous[6*old:]])
        split = int(np.argmax(weights))
        shapes = np.vstack((shapes, shapes[split]))
        weights[split] /= 2
        weights = np.r_[weights, weights[split]]
        # Small asymmetric perturbations release the duplicated component.
        shapes[split, 0] -= .04
        shapes[-1, 0] += .04
        return np.r_[shapes.ravel(), np.log(weights[1:]/weights[0])]
    shapes = np.tile([np.log(100.), 0., 8., np.log(2.), 1.5, 1.1], (count, 1))
    shapes[:, 0] += np.linspace(-.4, .8, count)+.35*(start-1)
    shapes += rng.normal(size=shapes.shape)*np.array([.25, .2, .5, .15, .12, .15])
    return np.r_[shapes.ravel(), np.zeros(count-1)]


class Objective:
    def __init__(self, kernel, observations, count, history, label):
        self.kernel, self.observations, self.count = kernel, observations, count
        self.history, self.label = history, label
        self.x = None
        self.calls = 0
        self.best = None
        self.started = time.monotonic()

    def evaluate(self, x):
        if self.x is None or not np.array_equal(self.x, x):
            self.calls += 1
            projected, jac = self.kernel.project(x, self.count)
            self.residual, self.deriv = self.observations.residual(projected, jac)
            self.x = x.copy()
            score = float(self.residual @ self.residual)
            row = dict(evaluation=self.calls, score=score, x=x.tolist(),
                       elapsed_seconds=time.monotonic()-self.started)
            with self.history.open("a") as stream:
                stream.write(json.dumps(row, allow_nan=False)+"\n")
            if self.best is None or score < self.best["score"]:
                self.best = row
            if self.calls % 50 == 0:
                print(f"{self.label}: {self.calls} evaluations; best Q={self.best['score']:.6g}", flush=True)
        return self.residual, self.deriv


def metrics(residual, nkin):
    residual = np.asarray(residual)
    return dict(score=float(residual @ residual),
                score_kinematic=float(residual[:nkin] @ residual[:nkin]),
                score_photometric=float(residual[nkin:] @ residual[nkin:]),
                rms_sigma=float(np.sqrt(np.mean(residual**2))),
                max_sigma=float(np.max(abs(residual))),
                max_kinematic_sigma=float(np.max(abs(residual[:nkin]))),
                max_photometric_sigma=float(np.max(abs(residual[nkin:]))))


def validate(kernel, observations, truth, x, count):
    settings = dict(kernel.settings)
    for k in ("radial_nodes", "velocity_nodes", "projection_nodes"):
        settings[k] *= 2
    refined = FixedPotentialProjector(truth.agama_potential(800), observations.radii, **settings)
    base = kernel.project(x, count, False)
    fine = refined.project(x, count, False)
    shift = observations.residual(fine)-observations.residual(base)
    # Double the outer integration domain separately from resolution.
    wide_settings = dict(settings, r_max=2*settings["r_max"])
    wide = FixedPotentialProjector(refined.potential, observations.radii, **wide_settings)
    wide_shift = observations.residual(wide.project(x, count, False))-observations.residual(fine)
    # Native AGAMA adaptive integration does not use our cached quadrature.
    r = np.geomspace(.045, 68., 14)
    independent = FixedPotentialProjector(refined.potential, r, **settings)
    df = refined.native_df(x, count)
    density, tensor = agama_pc().GalaxyModel(refined.potential, df).moments(np.column_stack((r, r*0)))
    native = np.column_stack((density, density*tensor[:, 2], density*tensor[:, 0], density*tensor[:, 1]))
    own = independent.project(x, count, False)
    relative = np.max(abs(own/native-1), axis=0)
    rforce = np.geomspace(.01, truth.r_t*2, 100)
    xyz = np.column_stack((rforce, rforce*0, rforce*0))
    force_error = np.max(abs(-kernel.potential.force(xyz)[:, 0]*rforce*rforce/truth.G/truth.enclosed_mass(rforce)-1))
    # Independent mock line-of-sight quadrature refinement.
    mock_fine = truth.projected_moments(observations.radii, order=320)
    mock_array = np.column_stack([mock_fine[k] for k in ("Sigma", "los", "pmr", "pmt")])
    mock_shift = observations.residual(mock_array)
    answer = dict(refinement_max_shift_sigma=float(np.max(abs(shift))),
                  radial_domain_max_shift_sigma=float(np.max(abs(wide_shift))),
                  native_projection_max_relative=relative.tolist(),
                  mock_projection_max_shift_sigma=float(np.max(abs(mock_shift))),
                  force_max_relative=float(force_error), refined_settings=settings,
                  refined_metrics=metrics(observations.residual(fine), observations.nkin))
    answer["passed"] = bool(max(abs(shift)) < .05 and max(abs(wide_shift)) < .05
                             and max(relative) < .005 and max(abs(mock_shift)) < .001
                             and force_error < 1e-5)
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["single_isotropic", "mixed_no_dm", "mixed_with_dm"], required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--max-evaluations", type=int, default=600)
    parser.add_argument("--seed", type=int, default=220926)
    parser.add_argument("--wide-bounds", action="store_true", help="explicit boundary sensitivity control")
    parser.add_argument("--components", type=int, nargs="+", default=[1,2,3])
    parser.add_argument("--initial-run", type=Path, help="reuse the corresponding saved best fit as start zero")
    parser.add_argument("--data-run", type=Path, default=ROOT/"results/df/pilot_no_dm_20260922")
    args = parser.parse_args()
    if args.starts < 2 or args.max_evaluations < 10:
        parser.error("need at least two starts and ten evaluations")
    if args.components != sorted(set(args.components)) or any(n not in (1,2,3) for n in args.components):
        parser.error("components must be an increasing subset of 1 2 3")
    args.out.mkdir(parents=True, exist_ok=False)
    began = time.monotonic()
    manifest = dict(status="running", case=args.case, started_utc=datetime.now(timezone.utc).isoformat(),
                    seed=args.seed, starts=args.starts, max_evaluations=args.max_evaluations,
                    stage="fixed known potential; noise-free capacity test; no mass inference",
                    pass_criterion="RMS discrepancy <0.1 error and max discrepancy <0.3 error; numerical validation passes",
                    parameter_names=list(PARAMETER_NAMES), wide_bounds=args.wide_bounds,
                    bounds={str(n):[v.tolist() for v in capacity_bounds(n, args.wide_bounds)] for n in args.components})
    write_json(args.out/"run.json", manifest)
    try:
        code = args.out/"code"
        code.mkdir()
        sources = [Path(__file__), ROOT/"src/ocen_dm/kinematics/df_capacity.py",
                   ROOT/"src/ocen_dm/kinematics/df_mock.py", ROOT/"src/ocen_dm/kinematics/positive_df.py",
                   ROOT/"src/ocen_dm/kinematics/likelihood.py", ROOT/"src/ocen_dm/kinematics/df_fit.py",
                   ROOT/"src/ocen_dm/kinematics/run_io.py", ROOT/"tests/test_df_capacity.py"]
        manifest["source_sha256"] = {str(p.relative_to(ROOT)):sha256(p) for p in sources}
        for p in sources:
            shutil.copy2(p, code/p.name)
        agama = agama_pc()
        manifest["agama"] = dict(version=agama.__version__, module=agama.__file__,
             binary_sha256={p.name:sha256(p) for p in Path(agama.__file__).parent.glob("*.so")})
        metadata = json.loads((args.data_run/"summary.json").read_text())
        data = read_data_snapshot(args.data_run, metadata["data_snapshot"])
        manifest["observational_template"] = write_data_snapshot(data, args.out)
        photo = PhotometricData.from_dict(json.loads((args.data_run/"photometry.json").read_text()))
        write_json(args.out/"photometry_template.json", photo.to_dict())
        truth = challenge_mock(args.case)
        write_json(args.out/"truth.json", truth.to_dict())
        observations = CapacityObservations(data, photo, truth)
        write_json(args.out/"mock_observations.json", observations.to_dict())
        kernel = FixedPotentialProjector(truth.agama_potential(), observations.radii)
        manifest["quadrature"] = kernel.settings
        manifest["anchor_radius_pc"] = kernel.anchor_radius
        write_json(args.out/"run.json", manifest)
        initial_fits = {}
        if args.initial_run:
            seed_file = args.initial_run/"summary.json"
            seed_summary = json.loads(seed_file.read_text())
            if seed_summary["case"] != args.case:
                raise ValueError("initial fit belongs to a different mock case")
            initial_fits = {e["components"]:e["x"] for e in seed_summary["ladder"]}
            manifest["initial_run"] = dict(path=str(args.initial_run.resolve()), summary_sha256=sha256(seed_file))
            shutil.copy2(seed_file, args.out/"initial_summary.json")
            write_json(args.out/"run.json", manifest)
        rng = np.random.default_rng(args.seed)
        ladder, previous = [], None
        for count in args.components:
            fits = []
            for start in range(args.starts):
                label = f"{args.case}/{count}df/start{start}"
                objective = Objective(kernel, observations, count, args.out/f"{count}df_start{start}.jsonl", label)
                initial = starting_vector(count, start, previous, rng)
                if start == 0 and count in initial_fits:
                    initial = np.array(initial_fits[count])
                lower, upper = capacity_bounds(count, args.wide_bounds)
                initial = np.clip(initial, lower+1e-7, upper-1e-7)
                initial_score = float(np.sum(objective.evaluate(initial)[0]**2))
                print(f"{label}: initial Q={initial_score:.6g}", flush=True)
                result = least_squares(lambda x: objective.evaluate(x)[0], initial,
                         jac=lambda x: objective.evaluate(x)[1], bounds=(lower, upper),
                         max_nfev=args.max_evaluations, x_scale="jac", ftol=1e-8, xtol=1e-8, gtol=1e-7)
                x = np.array(objective.best["x"])
                residual = observations.residual(kernel.project(x, count, False))
                boundary = np.flatnonzero(np.minimum((x-lower)/(upper-lower), (upper-x)/(upper-lower)) < 1e-3)
                row = dict(start=start, x=x.tolist(), initial_score=initial_score, metrics=metrics(residual, observations.nkin),
                           success=bool(result.success), status=int(result.status), message=str(result.message),
                           optimality=float(result.optimality), evaluations=objective.calls,
                           near_bound_indices=boundary.tolist(), elapsed_seconds=time.monotonic()-objective.started)
                fits.append(row)
                write_json(args.out/f"{count}df_start{start}.json", row)
                print(f"{label}: finished Q={row['metrics']['score']:.6g}; {result.message}", flush=True)
            best = min(fits, key=lambda f:f["metrics"]["score"])
            previous = np.array(best["x"])
            validation = validate(kernel, observations, truth, previous, count)
            projected = kernel.project(previous, count, False)
            residual = observations.residual(projected)
            row = dict(components=count, best_start=best["start"], fits=fits, x=previous.tolist(),
                       metrics=best["metrics"], validation=validation,
                       projected=projected.tolist(), predictions=observations.reduce(projected).tolist(),
                       residual_sigma=residual.tolist())
            fine = validation["refined_metrics"]
            row["capacity_demonstrated"] = bool(validation["passed"] and fine["rms_sigma"] < .1 and fine["max_sigma"] < .3)
            ladder.append(row)
            write_json(args.out/"summary.json", dict(case=args.case, status="running", ladder=ladder))
            print(f"{args.case}/{count}df: validation={validation['passed']}; capacity={row['capacity_demonstrated']}", flush=True)
        summary = dict(case=args.case, status="completed", ladder=ladder,
                       elapsed_seconds=time.monotonic()-began,
                       interpretation="Positive light DF capacity in fixed true gravity; not self-consistent fitted mass models or posterior constraints")
        write_json(args.out/"summary.json", summary)
        manifest.update(status="completed", elapsed_seconds=summary["elapsed_seconds"])
        write_json(args.out/"run.json", manifest)
    except BaseException as exc:
        manifest.update(status="failed", error=repr(exc), elapsed_seconds=time.monotonic()-began)
        write_json(args.out/"run.json", manifest)
        raise


if __name__ == "__main__":
    main()
