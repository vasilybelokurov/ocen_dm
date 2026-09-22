"""Launch preserved-source turnover or single-transition fits from rung-1 inputs.

Run with the project Python and PYTHONPATH=src. `launch` starts background
workers; `worker BATCH LABEL` is used only by the saved launcher. Existing
fit directories or batch directories are never reused.
Use ``launch BATCH --profile simple`` for the matched-prior single transition;
the default ``turnover`` preserves the original rung-2 launch recipe.
"""
from dataclasses import replace
from datetime import datetime, timezone
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import yaml

from ocen_dm.kinematics.fit import DarkMatterModel, FitProblem, Prior, run_nested
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config, data_fingerprint, read_data_snapshot, write_data_snapshot, sha256

ROOT = Path(os.environ.get("OCEN_DM_ROOT", Path(__file__).resolve().parents[1]))
SOURCES = {"rung2_K1_turnover": "rung1_K1_slice",
           "rung2_K2_nfw_turnover": "rung1_K2_nfw",
           "rung2_K2_cored_turnover": "rung1_K2_cored_slice"}


def companion_family(template, profile):
    """Keep the saved mass/distance model and select an explicit anisotropy prior.

    The simple case nests the diagnostic constant-beta family at beta_inf=beta_0.
    Custom priors are stored in the resolved model snapshot used by the worker
    and reporting code; no global constructor defaults are changed.
    """
    if profile not in ("simple", "turnover"):
        raise ValueError("profile must be simple or turnover")
    if (template.backend != "jeans" or not template.constant_beta
            or template.nuisance_instruments or "beta_0" not in template.names):
        raise ValueError("parent must be a constant-beta Jeans fit without instrument scales")
    kw = dict(mge_fit=template.mge_fit, tracer=template.tracer,
              distance_prior=template.distance_prior, distance_kpc=template.fixed_distance_kpc,
              fixed=template.fixed, instruments=(), backend="jeans",
              anisotropy_profile="monotonic" if profile == "simple" else "turnover",
              beta0_max=.5 if profile == "simple" else -.5)
    if isinstance(template, DarkMatterModel):
        kw.update(gamma=template.gamma, r_t=template.r_t)
    family = type(template)(**kw)
    original = {p.name: p for p in template.parameters}
    if profile == "simple" and original["beta_0"].prior != Prior("uniform", -1., .5):
        raise ValueError("simple comparison requires the agreed rung-1 beta prior [-1, 0.5]")
    # Preserve even non-default mass/distance priors from the parent snapshot.
    family.parameters = tuple(
        original[p.name] if p.name in original and p.name != "beta_0" else
        replace(p, prior=Prior("uniform", -1., .5)) if profile == "simple"
        and p.name in ("beta_0", "beta_inf") else p
        for p in family.parameters)
    family.names = tuple(p.name for p in family.parameters)
    return family


def check_constant_limit(old, problem, summary):
    """Verify the parent's best fit remains exactly representable at three scales."""
    theta = {name: summary["parameters"][name]["ml"] for name in old.family.names}
    reference = old.loglike_vector(np.array([theta[n] for n in old.family.names]))
    checks = {}
    for radius in (.5, 10., 100.):
        trial = {**theta, "beta_inf": theta["beta_0"], "r_beta": radius}
        loglike = problem.loglike_vector(np.array([trial[n] for n in problem.family.names]))
        if not np.isclose(loglike, reference, rtol=0., atol=1e-8):
            raise ValueError(f"constant-beta limit failed at {radius} pc: {loglike} vs {reference}")
        checks[str(radius)] = dict(loglike=float(loglike), difference=float(loglike-reference))
    return checks


def now():
    return datetime.now(timezone.utc).isoformat()


def worker(batch, label):
    config_path = batch/f"{label}.json"
    launch = json.loads((batch/"launch.json").read_text())
    assert sha256(config_path) == launch["job_config_sha256"][label]
    for name, expected in launch["source_sha256"].items():
        assert sha256(batch/"code"/name) == expected, name
    config = json.loads(config_path.read_text())
    family = family_from_config(config["model"])
    data = read_data_snapshot(batch, config["data_snapshot"])
    assert data_fingerprint(data) == config["data_fingerprint"]
    summary = run_nested(FitProblem(family, data), Path(config["output"]),
                         n_live=400, dlogz=.5, seed=42, step_sampler=True, verbose=True,
                         data_provenance={"kind": "real", "parent_rung1": config["parent_rung1"],
                                          "launch_config_sha256": sha256(config_path)},
                         dataset_options=config["dataset_options"])
    print(json.dumps({k: summary[k] for k in ("logz", "logzerr", "chi2_ml_total", "n_calls", "elapsed_s")}), flush=True)


def launch(batch, profile="turnover"):
    if profile not in ("simple", "turnover"):
        raise ValueError("profile must be simple or turnover")
    sources = {label.replace("_turnover", "_simple_beta") if profile == "simple" else label: source
               for label, source in SOURCES.items()}
    if batch.exists():
        raise FileExistsError(f"Batch already exists: {batch}")
    for label in sources:
        if (ROOT/"results/fits"/label).exists():
            raise FileExistsError(f"Fit already exists: {label}")
    jobs = {}
    shared_data = None
    for label, source in sources.items():
        # The report loader verifies each completed rung-1 maximum likelihood
        # against its saved observations and full resolved model.
        metadata = load_run_metadata(ROOT/"results/fits"/source)
        old = problem_for(metadata)
        family = companion_family(old.family, profile)
        problem = FitProblem(family, old.data)
        constant_limit = check_constant_limit(old, problem, metadata["summary"]) if profile == "simple" else None
        # Exercise the same snapshot reconstruction the background worker uses.
        restored = family_from_config(family_config(family))
        assert family_config(restored) == family_config(family)
        rng = np.random.default_rng(123)
        cubes = rng.uniform(.001, .999, size=(256, len(family.names)))
        cubes = np.vstack((cubes, np.full(len(family.names), .01), np.full(len(family.names), .99)))
        loglike = np.array([problem.loglike_unit(u) for u in cubes])
        assert np.all(np.isfinite(loglike)) and np.min(loglike) > -1e50
        fingerprint = data_fingerprint(old.data)
        assert old.data.n_points == 89
        if shared_data is not None:
            assert fingerprint == data_fingerprint(shared_data)
        shared_data = old.data
        previous = yaml.safe_load((ROOT/"results/fits"/source/"run.yaml").read_text())
        options = {**previous["dataset_options"], "constant_beta": False, "isotropic": False,
                   "anisotropy_profile": family.anisotropy_profile, "no_scales": True}
        jobs[label] = dict(model=family_config(family), parent_rung1=source,
                           comparison_profile=profile,
                           parent_manifest_sha256=sha256(ROOT/"results/fits"/source/"run.yaml"),
                           data_fingerprint=fingerprint, dataset_options=options,
                           output=str(ROOT/"results/fits"/label),
                           preflight=dict(prior_points=len(cubes), all_finite=True,
                                          constant_limit=constant_limit,
                                          loglike_min=float(loglike.min()), loglike_max=float(loglike.max())))
        print(label, len(family.names), "parameters;", len(cubes), "finite prior checks; original likelihood replayed", flush=True)
    batch.mkdir(parents=True)
    snapshot = write_data_snapshot(shared_data, batch)
    code = batch/"code"
    for path in list((ROOT/"src").rglob("*.py"))+[ROOT/"pyproject.toml", Path(__file__).resolve()]:
        relative = path.relative_to(ROOT.resolve()) if path.is_relative_to(ROOT.resolve()) else path.relative_to(ROOT)
        destination = code/relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    source_hashes = {str(p.relative_to(code)): sha256(p) for p in code.rglob("*") if p.is_file()}
    for label, config in jobs.items():
        config["data_snapshot"] = snapshot
        (batch/f"{label}.json").write_text(json.dumps(config, indent=2)+"\n")
    env = {"PYTHONPATH": str(code/"src"), "OCEN_DM_ROOT": str(ROOT),
           "PYTHONDONTWRITEBYTECODE": "1", "PYTHONUNBUFFERED": "1",
           "MPLCONFIGDIR": "/private/tmp/ocen-nbody-mpl",
           **{name: "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")}}
    record = dict(created_utc=now(), status="prepared", comparison_profile=profile,
                  source_sha256=source_hashes,
                  environment=env, jobs={},
                  job_config_sha256={label: sha256(batch/f"{label}.json") for label in jobs})
    path = batch/"launch.json"
    path.write_text(json.dumps(record, indent=2)+"\n")
    for label in jobs:
        command = [sys.executable, "-u", str(code/"bin/run_rung2_fits.py"), "worker", str(batch), label]
        with (batch/f"{label}.log").open("x") as log:
            process = subprocess.Popen(command, cwd=ROOT, env={**os.environ, **env},
                                       stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
        caffeine = subprocess.Popen(["/usr/bin/caffeinate", "-i", "-w", str(process.pid)],
                                     stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                     stderr=subprocess.DEVNULL, start_new_session=True)
        record["jobs"][label] = dict(pid=process.pid, caffeinate_pid=caffeine.pid, command=command)
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(record, indent=2)+"\n")
        temporary.replace(path)
        print(label, "PID", process.pid, flush=True)
    record.update(status="launched", launched_utc=now())
    # Workers read the immutable source/config hashes already saved above.
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(record, indent=2)+"\n")
    temporary.replace(path)
    print("Launch record:", path, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["launch", "worker"])
    parser.add_argument("batch", type=Path)
    parser.add_argument("label", nargs="?")
    parser.add_argument("--profile", choices=["turnover", "simple"], default="turnover")
    args = parser.parse_args()
    if args.action == "launch":
        launch(args.batch.resolve(), args.profile)
    else:
        worker(args.batch.resolve(), args.label)
