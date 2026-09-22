#!/usr/bin/env python3
"""Evaluate or locally optimize positive AGAMA stellar DFs on the saved data.

No posterior sampler is launched. Fresh result directories are mandatory.
Use --photometry-snapshot for exact replay, or explicitly choose --sigma-mag.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from ocen_dm.kinematics.df_fit import (
    DFJointProblem, FitCoordinate, PhotometricData, config_at, model_config_from_dict,
)
from ocen_dm.kinematics.likelihood import ARCSEC_PER_RAD
from ocen_dm.kinematics.positive_df import agama_pc
from ocen_dm.kinematics.run_io import code_state, read_data_snapshot, sha256, write_data_snapshot
from ocen_dm.light_model import load_tracer_profile


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")


def plot_result(problem, result, outpath):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    rows = []
    for ax, p in zip(axes.flat, problem.data.profiles):
        prediction = result["predictions"][p.name]
        ax.errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi], fmt=".", color="k", label="Data")
        ax.plot(p.r, prediction, "o-", ms=3, color="C0", label="Positive DF pilot")
        ax.set(xscale="log", xlabel="Projected radius [arcsec]",
               ylabel="Dispersion [km/s]" if p.kind == "los" else "Dispersion [mas/yr]",
               title=p.name.replace("_", " "))
        ax.text(.04, .06, f"χ² = {result['terms'][p.name]['chi2']:.1f} ({p.n} bins)",
                transform=ax.transAxes, fontsize=9)
        rows.append(dict(name=p.name, r_arcsec=p.r.tolist(), data=p.value.tolist(),
                         err_lo=p.err_lo.tolist(), err_hi=p.err_hi.tolist(),
                         prediction=prediction.tolist()))
    axes.flat[0].legend(fontsize=8)
    ax = axes.flat[-1]
    photo = problem.photometry
    ax.errorbar(photo.r_arcsec, photo.mu, yerr=photo.sigma_mag, fmt=".", color="k")
    ax.plot(photo.r_arcsec, result["photometry"]["prediction"], color="C0")
    ax.set(xscale="log", xlabel="Projected radius [arcsec]", ylabel="Surface brightness [mag]",
           title="Photometric shape (adopted errors)")
    ax.invert_yaxis()
    fig.suptitle("AGAMA positive stellar DF — diagnostic pilot, not a posterior fit")
    fig.savefig(outpath, dpi=170)
    plt.close(fig)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--data-run", type=Path, default=ROOT/"results/fits/rung2_K1_turnover")
    photo_group = parser.add_mutually_exclusive_group(required=True)
    photo_group.add_argument("--sigma-mag", type=float)
    photo_group.add_argument("--photometry-snapshot", type=Path)
    parser.add_argument("--max-evaluations", type=int, default=0,
                        help="0 evaluates the supplied model; positive values run bounded Powell")
    args = parser.parse_args()
    if args.max_evaluations < 0:
        parser.error("--max-evaluations must be non-negative")
    specification = json.loads(args.config.read_text())
    config = model_config_from_dict(specification["model"])
    coordinates = [FitCoordinate(**p) for p in specification.get("fit_coordinates", [])]
    if args.max_evaluations and not coordinates:
        parser.error("optimization requires explicit fit_coordinates")
    plotpath = ROOT/"plots"/f"df_{args.out.name}.png"
    tablepath = ROOT/"results/plot_data"/f"df_{args.out.name}.json"
    if plotpath.exists() or tablepath.exists():
        raise FileExistsError("choose a fresh run name; matching plot/table already exists")
    metadata = json.loads((args.data_run/"summary.json").read_text())
    data = read_data_snapshot(args.data_run, metadata["data_snapshot"])
    # This compact report layout is specifically for the adopted five profiles.
    if len(data.profiles) != 5:
        raise ValueError("pilot report expects the adopted five-profile dataset")
    photometry = (PhotometricData.from_dict(json.loads(args.photometry_snapshot.read_text()))
                  if args.photometry_snapshot else PhotometricData.from_relative_weights(
                      load_tracer_profile("composite"), args.sigma_mag))
    problem = DFJointProblem(data, photometry)
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    manifest = dict(status="running", started_utc=datetime.now(timezone.utc).isoformat(),
                    specification=specification, data_source=str(args.data_run.resolve()),
                    max_evaluations=args.max_evaluations, code=code_state(),
                    interpretation="local diagnostic objective; no evidence or posterior",
                    streaming="published moment-level subtraction; odd rotating DF not fitted")
    try:
        manifest["data_snapshot"] = write_data_snapshot(data, args.out)
        agama = agama_pc()
        manifest["agama"] = dict(version=getattr(agama, "__version__", None),
                                  module=agama.__file__,
                                  binary_sha256={p.name: sha256(p) for p in
                                      Path(agama.__file__).parent.glob("*.so")})
        photo_file = args.out/"photometry.json"
        write_json(photo_file, photometry.to_dict())
        manifest["photometry"] = dict(file=photo_file.name, sha256=sha256(photo_file))
        code_dir = args.out/"code"
        code_dir.mkdir()
        sources = [Path(__file__).resolve(),
                   ROOT/"src/ocen_dm/kinematics/positive_df.py",
                   ROOT/"src/ocen_dm/kinematics/regularized_df.py",
                   ROOT/"src/ocen_dm/kinematics/lowered_isothermal.py",
                   ROOT/"src/ocen_dm/kinematics/df_fit.py",
                   ROOT/"src/ocen_dm/kinematics/likelihood.py"]
        for source in sources:
            shutil.copy2(source, code_dir/source.name)
        manifest["entrypoint_sha256"] = sha256(Path(__file__))
        write_json(args.out/"run.json", manifest)
        best = problem.evaluate(config)
        initial = best["objective"]
        initial_kinematic = best["chi2_kinematic"]
        calls = 0
        history = args.out/"evaluations.jsonl"

        def objective(x):
            nonlocal best, calls
            calls += 1
            row = dict(evaluation=calls, coordinates=np.asarray(x).tolist())
            try:
                candidate = config_at(config, coordinates, x)
                result = problem.evaluate(candidate)
                score = result["objective"]
                row.update(objective=score, chi2_kinematic=result["chi2_kinematic"],
                           chi2_photometric=result["photometry"]["chi2"])
                if score < best["objective"]:
                    best = result
                    write_json(args.out/"best_model.json", candidate.to_dict())
            except (ValueError, FloatingPointError) as exc:
                score = 1e100
                row["rejected"] = str(exc)
            with history.open("a") as stream:
                stream.write(json.dumps(row, allow_nan=False)+"\n")
            if calls % 10 == 0:
                print(f"{calls} evaluations; best joint objective {best['objective']:.2f}; "
                      f"kinematic χ² {best['chi2_kinematic']:.2f}", flush=True)
            return score

        optimizer = None
        print(f"Initial: kinematic χ² {best['chi2_kinematic']:.2f}, "
              f"photometric χ² {best['photometry']['chi2']:.2f}", flush=True)
        if args.max_evaluations:
            result = minimize(objective, [p.get(config) for p in coordinates], method="Powell",
                              bounds=[p.bounds for p in coordinates],
                              options=dict(maxfev=args.max_evaluations, xtol=.002, ftol=.001))
            optimizer = dict(success=bool(result.success), message=str(result.message),
                             evaluations=int(result.nfev))
        model = best["model"]
        write_json(args.out/"best_model.json", model.config.to_dict())
        # Independent native AGAMA projection plus a complete refined equilibrium.
        R = np.geomspace(.02, 80., 16)
        fast, direct = model.projected_moments(R), model.direct_projected_moments(R)
        projection_error = {k: float(np.max(np.abs(fast[k]/direct[k]-1))) for k in fast}
        n = model.config.numerics
        fine = replace(n, potential_nodes=2*n.potential_nodes, velocity_nodes=2*n.velocity_nodes,
                       moment_nodes=2*n.moment_nodes, projection_nodes=2*n.projection_nodes,
                       iteration_tolerance=n.iteration_tolerance/2)
        refined = problem.evaluate(replace(model.config, numerics=fine))
        prediction_shift = {p.name: float(np.max(np.abs(
            refined["predictions"][p.name]-best["predictions"][p.name])/np.minimum(p.err_lo, p.err_hi)))
            for p in data.profiles}
        validation = dict(direct_projection_relative_error=projection_error,
                          refinement_max_shift_in_data_sigma=prediction_shift,
                          refined_chi2_kinematic=refined["chi2_kinematic"],
                          refined_objective=refined["objective"],
                          refined_diagnostics=refined["model"].diagnostics)
        validation["passed"] = bool(max(projection_error.values()) < .005 and
                                    max(prediction_shift.values()) < .1)
        plotpath.parent.mkdir(parents=True, exist_ok=True)
        tablepath.parent.mkdir(parents=True, exist_ok=True)
        rows = plot_result(problem, best, plotpath)
        r = np.geomspace(.001, 300., 180)
        intrinsic = model.intrinsic_moments(r)
        write_json(tablepath, dict(run=str(args.out.resolve()), kinematics=rows,
                                  photometry=photometry.to_dict(),
                                  photometric_prediction=best["photometry"]["prediction"].tolist(),
                                  r_pc=r.tolist(), intrinsic={k: v.tolist() for k, v in intrinsic.items()}))
        summary = dict(status="pilot_finished" if optimizer else "evaluated", optimizer=optimizer,
                       initial_objective=initial, initial_chi2_kinematic=initial_kinematic,
                       objective=best["objective"], chi2_kinematic=best["chi2_kinematic"],
                       chi2_photometric=best["photometry"]["chi2"],
                       loglike_kinematic=best["loglike_kinematic"],
                       loglike_profiled=best["loglike_profiled"], photometric_zero_point=best["photometry"]["offset"],
                       terms=best["terms"], n_kinematic=data.n_points,
                       n_photometric=len(photometry.mu), diagnostics=model.diagnostics,
                       validation=validation, model=model.config.to_dict(),
                       data_snapshot=manifest["data_snapshot"], photometric_errors=photometry.error_model,
                       elapsed_seconds=time.monotonic()-start, plot=str(plotpath), plot_data=str(tablepath))
        write_json(args.out/"summary.json", summary)
        manifest.update(status=summary["status"], elapsed_seconds=summary["elapsed_seconds"])
        write_json(args.out/"run.json", manifest)
        print(json.dumps({k: summary[k] for k in ("status", "chi2_kinematic", "chi2_photometric", "validation")}, indent=2))
    except BaseException as exc:
        manifest.update(status="failed", error=repr(exc), elapsed_seconds=time.monotonic()-start)
        write_json(args.out/"run.json", manifest)
        raise


if __name__ == "__main__":
    main()
