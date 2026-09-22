#!/usr/bin/env python3
"""Prepare, run/resume, and report bounded noiseless compact-DF recovery.

Workers execute frozen source. Each evaluated improvement is saved atomically.
Resume restarts the optimizer from the best saved point, not its lost memory.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np
from scipy.optimize import least_squares

from ocen_dm.kinematics.compact_recovery import (
    mass_profiles, mock_problem, recovery_gates, recovery_metrics,
    refined_config, residual_vector,
)
from ocen_dm.kinematics.df_fit import (
    DFJointProblem, FitCoordinate, PhotometricData, build_df_model, config_at,
    model_config_from_dict,
)
from ocen_dm.kinematics.regularized_df import RegularizedDFModel
from ocen_dm.kinematics.run_io import code_state, read_data_snapshot, sha256, write_data_snapshot


def now():
    return datetime.now(timezone.utc).isoformat()


def dump(path, row):
    temp = path.with_suffix(path.suffix+".tmp")
    with temp.open("w") as stream:
        json.dump(row, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temp.replace(path)


def read(path):
    return json.loads(path.read_text())


def load_problem(directory):
    record = read(directory/"mock.json")
    return DFJointProblem(read_data_snapshot(directory, record["data_snapshot"]),
                          PhotometricData.from_dict(read(directory/"photometry.json")))


def prepare(out, project):
    out.mkdir(parents=True, exist_ok=False)
    frozen = out/"code"
    (frozen/"bin").mkdir(parents=True)
    shutil.copytree(project/"src", frozen/"src", ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".DS_Store"))
    shutil.copy2(Path(__file__), frozen/"bin"/Path(__file__).name)
    shutil.copytree(project/"configs/df", frozen/"configs/df")
    shutil.copytree(project/"tests", frozen/"tests",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".DS_Store"))
    source = project/"results/df/pilot_no_dm_20260922"
    inputs = out/"inputs"
    inputs.mkdir()
    for name in ("summary.json", "data_snapshot.json", "photometry.json"):
        shutil.copy2(source/name, inputs/name)
    spec = read(frozen/"configs/df/regularized_no_dm.json")
    base = model_config_from_dict(spec["model"])
    template = read_data_snapshot(inputs, read(inputs/"summary.json")["data_snapshot"])
    photo = PhotometricData.from_dict(read(inputs/"photometry.json"))
    manifest = dict(schema_version=1, status="preparing", created_utc=now(),
                    runtime=code_state(),
                    project_root=str(project), source_commit=subprocess.check_output(
                        ["git", "rev-parse", "HEAD"], cwd=project, text=True).strip(),
                    scope="Matched regularized stellar family, noiseless, cored halo, common mass/light DF; diagnostic objective only",
                    fixed=dict(distance_kpc=base.distance_kpc, gamma=0, halo_taper_pc=base.matter.r_t),
                    gates=dict(rms_sigma=.1, max_sigma=.3, M_star_fraction=.05,
                               total_mass_profile_fraction=.1, rho20_fraction_if_DM=.1,
                               rho20_absolute_if_no_DM=.1, numerical_shift_sigma=.1),
                    source_sha256={str(p.relative_to(frozen)): sha256(p)
                                   for p in frozen.rglob("*") if p.is_file()},
                    input_sha256={p.name: sha256(p) for p in inputs.iterdir()}, jobs=[])
    from ocen_dm.kinematics.positive_df import agama_pc
    ag = agama_pc()
    manifest["agama"] = dict(version=getattr(ag, "__version__", None),
                             binary_sha256={p.name: sha256(p) for p in
                                            Path(ag.__file__).parent.glob("*.so")})
    dump(out/"batch.json", manifest)
    # All five stellar and both remnant coordinates; halo allows rho20=0 exactly.
    coords = spec["fit_coordinates"]+[
        dict(path="matter.rho20", lower=0., upper=10.),
        dict(path="matter.r_s", lower=5., upper=150., log=True)]
    starts = [(2.4e6, 260., 1.3, .25, 180., 2e5, 2., .7, 20.),
              (3.6e6, 340., 1.8, .7, 450., 4.5e5, 5., 4., 60.)]
    for case, rho in (("no_dm", 0.), ("cored_dm", 2.)):
        d = out/"mocks"/case
        d.mkdir(parents=True)
        c = replace(base, matter=replace(base.matter, rho20=rho))
        t = time.monotonic()
        print(f"Building refined {case} truth", flush=True)
        truth = build_df_model(refined_config(c))
        problem = mock_problem(template, photo, truth)
        snapshot = write_data_snapshot(problem.data, d)
        dump(d/"photometry.json", problem.photometry.to_dict())
        nominal = problem.evaluate(c)
        nominal_r = residual_vector(problem, nominal)
        # Error is an independently refined truth vs a cold nominal equilibrium.
        resolution_ok = bool(np.max(abs(nominal_r)) < .1)
        record = dict(config=truth.config.to_dict(), nominal_config=c.to_dict(),
                      diagnostics=truth.diagnostics, data_snapshot=snapshot,
                      profiles=mass_profiles(truth, np.geomspace(.5, 80., 64)),
                      nominal_residual_sigma=nominal_r.tolist(),
                      nominal_resolution_passed=resolution_ok,
                      seconds=time.monotonic()-t)
        dump(d/"mock.json", record)
        if not resolution_ok:
            manifest.update(status="truth_resolution_failed", failed_case=case)
            dump(out/"batch.json", manifest)
            raise ValueError(f"{case} numerical shift exceeds 0.1 adopted errors")
        print(f"{case} truth saved; max numerical shift {max(abs(nominal_r)):.5g} errors", flush=True)
        for branch in ("free_halo", "no_halo"):
            if branch == "no_halo" and case != "no_dm":
                continue
            chosen = coords if branch == "free_halo" else coords[:-2]
            cc = [FitCoordinate(**row) for row in chosen]
            for index, values in enumerate(starts):
                start = config_at(base, cc, [p.encode(v) for p, v in zip(cc, values)])
                manifest["jobs"].append(dict(id=f"{case}_{branch}_start{index}", case=case,
                                              branch=branch, start=index, config=start.to_dict(),
                                              coordinates=chosen))
        dump(out/"batch.json", manifest)
    manifest["jobs"].sort(key=lambda j: (j["branch"] == "no_halo", j["start"], j["case"]))
    manifest.update(status="prepared", prepared_utc=now(),
                    mock_sha256={str(p.relative_to(out)): sha256(p)
                                 for p in (out/"mocks").rglob("*.json")})
    dump(out/"batch.json", manifest)


class EvaluationBudget(Exception):
    pass


def fit(out, job_id, max_calls):
    manifest = read(out/"batch.json")
    job = next(j for j in manifest["jobs"] if j["id"] == job_id)
    d = out/"fits"/job_id
    d.mkdir(parents=True, exist_ok=True)
    if (d/"summary.json").exists():
        print(f"{job_id}: already finished", flush=True)
        return
    problem = load_problem(out/"mocks"/job["case"])
    truth = read(out/"mocks"/job["case"]/"mock.json")
    coords = [FitCoordinate(**c) for c in job["coordinates"]]
    base = model_config_from_dict(job["config"])
    low, high = np.array([c.bounds for c in coords]).T
    x0 = (np.array([c.get(base) for c in coords])-low)/(high-low)
    checkpoint = d/"best.json"
    best = read(checkpoint) if checkpoint.exists() else None
    if best:
        x0 = np.array(best["x"])
    attempt = now()
    calls, warm, started = 0, None, time.monotonic()
    cache = {}
    status = dict(status="running", pid=os.getpid(), attempt_started_utc=attempt,
                  max_calls_this_attempt=max_calls,
                  resume="new optimizer from saved best; no restored Hessian" if best else None)
    dump(d/"status.json", status)

    def fun(x):
        nonlocal calls, warm, best
        key = tuple(x)
        if key in cache:
            return cache[key].copy()
        if calls >= max_calls:
            raise EvaluationBudget()
        calls += 1
        row = dict(attempt=attempt, call=calls, x=x.tolist())
        config = config_at(base, coords, low+x*(high-low))
        try:
            model = RegularizedDFModel(config, initial_stellar_potential=warm)
            ev = problem.evaluate(config, model)
            residual = residual_vector(problem, ev)
            score = float(residual @ residual)
            warm = model.stellar_potential
            row.update(score=score, iterations=model.diagnostics["iterations"],
                       rho20=config.matter.rho20, M_star=config.M_star)
            if best is None or score < best["score"]:
                best = dict(score=score, x=x.tolist(), config=config.to_dict(),
                            residual_sigma=residual.tolist(), diagnostics=model.diagnostics,
                            profiles=mass_profiles(model, truth["profiles"]["r_pc"]),
                            attempt=attempt, call=calls, saved_utc=now())
                dump(checkpoint, best)
        except (ValueError, FloatingPointError) as exc:
            row["rejected"] = repr(exc)
            residual = np.full(problem.data.n_points+len(problem.photometry.mu), 1e6)
        row["seconds"] = time.monotonic()-started
        with (d/"evaluations.jsonl").open("a") as stream:
            stream.write(json.dumps(row, allow_nan=False)+"\n")
            stream.flush()
            os.fsync(stream.fileno())
        dump(d/"status.json", dict(status, calls=calls, seconds=row["seconds"],
                                   best_score=None if best is None else best["score"]))
        if calls == 1 or calls % 10 == 0:
            print(f"{job_id}: {calls} calls; best {None if best is None else best['score']}", flush=True)
        cache[key] = residual.copy()
        return residual

    def jacobian(x):
        # An absolute step in normalized coordinates stays resolved at rho20=0.
        base_residual = fun(x)
        columns = []
        for j in range(len(x)):
            trial = x.copy()
            step = .002 if x[j] + .002 <= 1 else -.002
            trial[j] += step
            columns.append((fun(trial)-base_residual)/step)
        return np.column_stack(columns)

    try:
        try:
            opt = least_squares(fun, x0, jac=jacobian,
                                bounds=(np.zeros(len(x0)), np.ones(len(x0))),
                                x_scale="jac", ftol=1e-5, xtol=2e-4,
                                gtol=2e-3, max_nfev=max_calls)
            optimizer = dict(success=bool(opt.success), message=opt.message,
                             nfev=opt.nfev, optimality=float(opt.optimality))
        except EvaluationBudget:
            optimizer = dict(success=False, message="actual equilibrium evaluation budget reached")
        if best is None:
            raise ValueError("no valid equilibrium was evaluated")
        dump(d/"status.json", dict(status, status="validating", calls=calls))
        c = model_config_from_dict(best["config"])
        cold = problem.evaluate(c)
        fine = problem.evaluate(refined_config(c))
        cr, fr = residual_vector(problem, cold), residual_vector(problem, fine)
        R = np.geomspace(.05, 80., 12)
        fast = cold["model"].projected_moments(R)
        direct = cold["model"].direct_projected_moments(R)
        projection = {k: float(np.max(abs(fast[k]/direct[k]-1))) for k in fast}
        cold_shift = float(np.max(abs(cr-np.array(best["residual_sigma"]))))
        fine_shift = float(np.max(abs(fr-cr)))
        numerical = bool(cold_shift < .1 and fine_shift < .1 and max(projection.values()) < .005)
        profiles = mass_profiles(fine["model"], truth["profiles"]["r_pc"])
        metrics = recovery_metrics(truth["profiles"], profiles)
        gates = recovery_gates(optimizer["success"], numerical, fr, metrics)
        summary = dict(job=job, optimizer=optimizer, calls=calls, seconds=time.monotonic()-started,
                       best=best, cold_score=float(cr@cr), refined_score=float(fr@fr),
                       refined_residual_sigma=fr.tolist(), refined_profiles=profiles,
                       validation=dict(cold_shift_sigma=cold_shift, refinement_shift_sigma=fine_shift,
                                       direct_projection_error=projection, passed=numerical),
                       recovery=metrics, gates=gates, finished_utc=now())
        dump(d/"summary.json", summary)
        dump(d/"status.json", dict(status, status="finished", calls=calls, gates=gates))
        print(f"{job_id}: {json.dumps(gates)}", flush=True)
    except BaseException as exc:
        dump(d/"status.json", dict(status, status="interrupted" if isinstance(exc, KeyboardInterrupt)
                                   else "failed", error=repr(exc), calls=calls))
        raise


def run(out, workers, max_calls):
    from concurrent.futures import ThreadPoolExecutor, as_completed
    manifest = read(out/"batch.json")
    if manifest["status"] not in ("prepared", "running", "finished", "interrupted"):
        raise ValueError("batch preparation has not passed")
    for relative, expected in manifest["source_sha256"].items():
        if sha256(out/"code"/relative) != expected:
            raise ValueError(f"frozen source changed: {relative}")
    for relative, expected in manifest["mock_sha256"].items():
        if sha256(out/relative) != expected:
            raise ValueError(f"saved mock changed: {relative}")
    # Prevent duplicate controllers; OS releases this lock when the process dies.
    import fcntl
    with (out/"controller.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest.update(status="running", controller_pid=os.getpid(), workers=workers,
                        max_calls_per_attempt=max_calls, resumed_utc=now())
        dump(out/"batch.json", manifest)

        def worker(job):
            cmd = [sys.executable, str(out/"code/bin/run_compact_df_recovery.py"), "fit",
                   "--out", str(out), "--job", job["id"], "--max-calls", str(max_calls)]
            with (out/(job["id"]+".log")).open("a") as stream:
                process = subprocess.Popen(cmd, stdout=stream, stderr=subprocess.STDOUT)
                print(f"Started {job['id']} pid={process.pid}", flush=True)
                return job["id"], process.wait()
        completed = {}
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(worker, j) for j in manifest["jobs"]]
            for future in as_completed(futures):
                name, rc = future.result()
                completed[name] = rc
                manifest["returncodes"] = completed
                dump(out/"batch.json", manifest)
                report(out)
                print(f"Finished {name}, exit={rc}", flush=True)
        manifest.update(status="finished", finished_utc=now())
        dump(out/"batch.json", manifest)
        report(out)


def report(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    manifest = read(out/"batch.json")
    project = Path(manifest["project_root"])
    rows = []
    for job in manifest["jobs"]:
        d = out/"fits"/job["id"]
        status = read(d/"status.json") if (d/"status.json").exists() else dict(status="queued")
        row = dict(job=job["id"], case=job["case"], branch=job["branch"], start=job["start"], status=status)
        if (d/"summary.json").exists():
            s = read(d/"summary.json")
            row.update(gates=s["gates"], recovery=s["recovery"], profiles=s["refined_profiles"],
                       score=s["refined_score"])
        elif (d/"best.json").exists():
            b = read(d/"best.json")
            row.update(profiles=b["profiles"], score=b["score"])
        rows.append(row)
    cases = list(dict.fromkeys(j["case"] for j in manifest["jobs"]))
    fig, axes = plt.subplots(1, len(cases), figsize=(11, 4), squeeze=False, constrained_layout=True)
    truths = {}
    for ax, case in zip(axes[0], cases):
        truth = read(out/"mocks"/case/"mock.json")["profiles"]
        truths[case] = truth
        ax.axhline(1, color="black", lw=1, label="Injected total mass")
        for row in rows:
            if row["case"] != case or "profiles" not in row:
                continue
            p = row["profiles"]
            label = row["job"].removeprefix(case+"_")+f"; rho20={p['rho20']:.3g}"
            if "gates" not in row:
                label += " (interim)"
            ax.plot(p["r_pc"], np.array(p["total"])/truth["total"], label=label)
        ax.set(xscale="log", xlabel="Radius [pc]", ylabel="Recovered / injected enclosed mass",
               title=case.replace("_", " "))
        ax.legend(fontsize=7)
    fig.suptitle("Noiseless compact DF recovery; interim values are not mass measurements")
    plot = project/"plots"/(out.name+".png")
    fig.savefig(plot, dpi=160)
    plt.close(fig)
    (project/"results/plot_data").mkdir(parents=True, exist_ok=True)
    dump(project/"results/plot_data"/(out.name+".json"),
         dict(batch=str(out), updated_utc=now(), batch_status=manifest["status"],
              truths=truths, jobs=rows, source_commit=manifest["source_commit"], plot_sha256=sha256(plot)))
    dump(out/"report.json", dict(updated_utc=now(), jobs=rows))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("prepare", "run", "fit", "report"))
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--project", type=Path, default=ROOT)
    p.add_argument("--job")
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--max-calls", type=int, default=180)
    p.add_argument("--detach", action="store_true")
    args = p.parse_args()
    if not 1 <= args.workers <= 2 or args.max_calls < 1:
        p.error("use one or two workers and a positive evaluation budget")
    out = args.out.resolve()
    if args.command == "prepare":
        prepare(out, args.project.resolve())
    elif args.command == "fit":
        if args.job is None:
            p.error("fit requires --job")
        import fcntl
        fitdir = out/"fits"/args.job
        fitdir.mkdir(parents=True, exist_ok=True)
        with (fitdir/"worker.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            fit(out, args.job, args.max_calls)
    elif args.command == "report":
        report(out)
    elif args.detach:
        cmd = [sys.executable, str(out/"code/bin/run_compact_df_recovery.py"), "run",
               "--out", str(out), "--workers", str(args.workers), "--max-calls", str(args.max_calls)]
        with (out/"controller.log").open("a") as log:
            process = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        print(f"Controller PID {process.pid}; log {out/'controller.log'}")
    else:
        run(out, args.workers, args.max_calls)


if __name__ == "__main__":
    main()
