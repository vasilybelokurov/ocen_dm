#!/usr/bin/env python3
"""Launch a Jeans reference model with the star-count tracer, as a twin of a J-Trager run.

``--source`` picks the J-Trager parent: rung2_K1_turnover (no DM; twin jeans_K1_turnover_counts_c, J-counts),
rung2_K2_nfw_turnover or rung2_K2_cored_turnover (with a dark halo; twins jeans_K2_{nfw,cored}_turnover_counts_c).

J-Trager is rung2_K1_turnover (no dark matter, free black hole, turnover anisotropy), whose tracer and stellar-mass
shape is the MGE of the composite HST-count/Trager light profile. J-counts copies that job exactly -- the same frozen
89-bin data snapshot (results/launches/rung2_turnover_20260921), parameters, priors, sampler settings (400 live
points, dlogz 0.5, seed 42, slice steps, one thread) -- and replaces only the MGE by the fit to the
completeness-corrected HST + Gaia counts with the tied amplitude ratio (bin/fit_count_tracer_mge.py,
data/processed/kinematics/tracer_mge_counts_c.json).

The worker is bin/run_rung2_fits.py `worker`, run from a frozen copy of the code in the batch directory.
Usage: PYTHONPATH=src python bin/run_jeans_counts_tracer.py BATCH_DIR [--source LABEL] [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import numpy as np

from ocen_dm.kinematics.fit import FitProblem
from ocen_dm.kinematics.run_io import data_fingerprint, family_config, family_from_config, read_data_snapshot, sha256

SOURCE_BATCH = ROOT/"results/launches/rung2_turnover_20260921"
SOURCES = {"rung2_K1_turnover": "jeans_K1_turnover_counts_c",
           "rung2_K2_nfw_turnover": "jeans_K2_nfw_turnover_counts_c",
           "rung2_K2_cored_turnover": "jeans_K2_cored_turnover_counts_c"}
TRACER = ROOT/"data/processed/kinematics/tracer_mge_counts_c.json"


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("batch", type=Path)
    p.add_argument("--source", choices=sorted(SOURCES), default="rung2_K1_turnover")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    SOURCE_LABEL, LABEL = args.source, SOURCES[args.source]
    batch = args.batch.resolve()
    if batch.exists() or (ROOT/"results/fits"/LABEL).exists():
        raise FileExistsError("batch or fit directory already exists")
    source = json.loads((SOURCE_BATCH/f"{SOURCE_LABEL}.json").read_text())
    data = read_data_snapshot(SOURCE_BATCH, source["data_snapshot"])
    assert data_fingerprint(data) == source["data_fingerprint"] and data.n_points == 89
    tracer = json.loads(TRACER.read_text())
    model = dict(source["model"], tracer="counts_c",
                 mge_fit={k: tracer[k] for k in ("sigmas_arcsec", "fractions", "rms_mag", "max_abs_resid_mag", "n_grid")})
    family = family_from_config(model)
    assert family_config(family_from_config(family_config(family))) == family_config(family)
    reference = family_from_config(source["model"])
    assert [(q.name, q.prior) for q in family.parameters] == [(q.name, q.prior) for q in reference.parameters]
    problem = FitProblem(family, data)
    rng = np.random.default_rng(123)
    cubes = rng.uniform(.001, .999, size=(256, len(family.names)))
    cubes = np.vstack((cubes, np.full(len(family.names), .01), np.full(len(family.names), .99)))
    loglike = np.array([problem.loglike_unit(u) for u in cubes])
    assert np.all(np.isfinite(loglike)) and loglike.min() > -1e50
    # the J-Trager maximum-likelihood point evaluated with the count tracer (diagnostic only)
    ml = json.loads((ROOT/"results/fits"/SOURCE_LABEL/"summary.json").read_text())["parameters"]
    theta = np.array([ml[n]["ml"] for n in family.names])
    ll_ref, ll_new = float(FitProblem(reference, data).loglike_vector(theta)), float(problem.loglike_vector(theta))
    print(f"{LABEL}: {len(family.names)} parameters, {len(cubes)} finite prior checks; at the J-Trager ML point "
          f"lnL {ll_ref:.2f} (Trager tracer) -> {ll_new:.2f} (count tracer)", flush=True)
    if args.dry_run:
        return
    batch.mkdir(parents=True)
    shutil.copy2(SOURCE_BATCH/source["data_snapshot"]["file"], batch/source["data_snapshot"]["file"])
    shutil.copy2(TRACER, batch/TRACER.name)
    code = batch/"code"
    for path in list((ROOT/"src").rglob("*.py"))+[ROOT/"pyproject.toml", ROOT/"bin/run_rung2_fits.py", Path(__file__).resolve()]:
        dest = code/path.resolve().relative_to(ROOT.resolve())
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, dest)
    config = dict(model=family_config(family), parent_rung1=source["parent_rung1"], twin_of=SOURCE_LABEL,
                  tracer_file=TRACER.name, tracer_sha256=sha256(TRACER), data_snapshot=source["data_snapshot"],
                  data_fingerprint=source["data_fingerprint"], dataset_options=dict(source["dataset_options"], tracer="counts_c"),
                  output=str(ROOT/"results/fits"/LABEL),
                  preflight=dict(prior_points=len(cubes), all_finite=True, loglike_min=float(loglike.min()),
                                 loglike_max=float(loglike.max()), loglike_at_jtrager_ml=dict(trager=ll_ref, counts=ll_new)))
    (batch/f"{LABEL}.json").write_text(json.dumps(config, indent=2)+"\n")
    env = {"PYTHONPATH": str(code/"src"), "OCEN_DM_ROOT": str(ROOT), "PYTHONDONTWRITEBYTECODE": "1", "PYTHONUNBUFFERED": "1",
           **{k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                               "NUMEXPR_NUM_THREADS")}}
    record = dict(created_utc=now(), status="prepared", source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                  text=True).strip(), source_sha256={str(q.relative_to(code)): sha256(q) for q in code.rglob("*") if q.is_file()},
                  environment=env, job_config_sha256={LABEL: sha256(batch/f"{LABEL}.json")}, jobs={})
    command = [sys.executable, "-u", str(code/"bin/run_rung2_fits.py"), "worker", str(batch), LABEL]
    with (batch/f"{LABEL}.log").open("x") as log:
        proc = subprocess.Popen(command, cwd=ROOT, env={**os.environ, **env}, stdin=subprocess.DEVNULL, stdout=log,
                                stderr=subprocess.STDOUT, start_new_session=True)
    caff = subprocess.Popen(["/usr/bin/caffeinate", "-i", "-w", str(proc.pid)], stdin=subprocess.DEVNULL,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    record["jobs"][LABEL] = dict(pid=proc.pid, caffeinate_pid=caff.pid, command=command)
    record.update(status="launched", launched_utc=now())
    (batch/"launch.json").write_text(json.dumps(record, indent=2)+"\n")
    print(LABEL, "PID", proc.pid, flush=True)


if __name__ == "__main__":
    main()
