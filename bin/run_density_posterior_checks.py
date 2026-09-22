#!/usr/bin/env python3
"""Launch independent free-density and rho20=2 conditional posterior checks.

Sources, observations, resolved priors and sampler settings are frozen per
batch. Original runs and existing output directories are never reused.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback

import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.density_conditioning import DensityConditionedModel
from ocen_dm.kinematics.fit import FitProblem, run_nested
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import (data_fingerprint, family_config, family_from_config,
                                      read_data_snapshot, sha256, write_data_snapshot)


ROOT=Path(os.environ.get("OCEN_DM_ROOT",Path(__file__).resolve().parents[1])).resolve()
PARENT="rung2_K2_cored_turnover"


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path,payload):
    temporary=path.with_name(f"{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload,indent=2)+"\n")
    temporary.replace(path)


def labels_for(batch):
    return {kind:f"rung2_K2_cored_turnover_{kind}_{batch.name}" for kind in ("freecheck","rho2")}


def derived_report(out,family,summary):
    post=Table.read(out/"posterior.ecsv")
    if isinstance(family,DensityConditionedModel):
        rho=np.full(len(post),family.density)
        post["M_dm_100"]=[family.complete(dict(zip(post.colnames,row))) ["M_dm_100"] for row in post]
    else:
        rho=np.array([float(family.halo({"M_dm_100":m,"r_s":rs}).density(20.).item())
                      for m,rs in zip(post["M_dm_100"],post["r_s"])])
    post["rho_dm_20"]=rho
    post["rho_dm_20"].unit="solMass / pc3"
    post.meta.update(density_definition="local point density at 20 pc",weights="equal posterior weights")
    post.write(out/"density_posterior.ecsv",overwrite=False)
    report=dict(created_utc=now(),rho_quantiles=dict(zip(("p16","p50","p84","p95","p99"),
                    map(float,np.percentile(rho,[16,50,84,95,99])))),
                tail_probabilities={str(t):float(np.mean(rho>t)) for t in (1.03,2.,3.,5.)},
                logz=summary["logz"],logzerr=summary["logzerr"],lnL_max=summary["lnL_max"],
                conditional_density=isinstance(family,DensityConditionedModel))
    if isinstance(family,DensityConditionedModel):
        report["parent_prior_pdf_at_fixed_density"]=family.prior_density_at_constraint
    atomic_json(out/"density_diagnostics.json",report)


def compare_if_complete(batch):
    labels=labels_for(batch)
    paths={k:ROOT/"results/fits"/label/"density_diagnostics.json" for k,label in labels.items()}
    if not all(p.exists() for p in paths.values()):
        return
    runs={k:json.loads(p.read_text()) for k,p in paths.items()}
    free,conditional=runs["freecheck"],runs["rho2"]
    logratio=conditional["logz"]-free["logz"]
    result=dict(created_utc=now(),runs=runs,
        conditional_minus_free_logz=logratio,
        independent_logz_error=float(np.hypot(free["logzerr"],conditional["logzerr"])),
        core_posterior_pdf_at_rho2=float(conditional["parent_prior_pdf_at_fixed_density"]*np.exp(logratio)),
        interpretation="p(rho=2 | data, core) = p_prior(rho=2 | core) * Z_conditional / Z_free; "
                       "this is a density per (Msun pc^-3), not a probability atom or DM/no-DM model odds",
        input_sha256={str(p.relative_to(ROOT)):sha256(p) for p in paths.values()})
    atomic_json(batch/"density_comparison.json",result)


def worker(batch,label):
    status=batch/f"{label}.status.json"
    try:
        manifest=json.loads((batch/"launch.json").read_text())
        config_path=batch/f"{label}.json"
        assert sha256(config_path)==manifest["job_config_sha256"][label]
        for name,digest in manifest["source_sha256"].items():
            assert sha256(batch/"code"/name)==digest,name
        config=json.loads(config_path.read_text())
        for name,digest in config["parent_sha256"].items():
            assert sha256(ROOT/name)==digest,name
        family=family_from_config(config["model"])
        data=read_data_snapshot(batch,config["data_snapshot"])
        assert data_fingerprint(data)==config["data_fingerprint"]
        atomic_json(status,dict(status="running",started_utc=now(),pid=os.getpid()))
        out=Path(config["output"])
        summary=run_nested(FitProblem(family,data),out,**config["sampler"],verbose=True,
            data_provenance={"kind":"real","parent_run":PARENT,
                             "launch_config_sha256":sha256(config_path),
                             "density_condition":config["density_condition"]},
            dataset_options=config["dataset_options"])
        derived_report(out,family,summary)
        atomic_json(status,dict(status="complete",finished_utc=now(),pid=os.getpid(),
                    logz=summary["logz"],logzerr=summary["logzerr"],
                    lnL_max=summary["lnL_max"],elapsed_s=summary["elapsed_s"]))
        compare_if_complete(batch)
    except Exception:
        atomic_json(status,dict(status="failed",failed_utc=now(),pid=os.getpid(),
                                traceback=traceback.format_exc()))
        raise


def launch(batch):
    labels=labels_for(batch)
    if batch.exists():
        raise FileExistsError(f"Batch already exists: {batch}")
    for label in labels.values():
        if (ROOT/"results/fits"/label).exists():
            raise FileExistsError(f"Fit already exists: {label}")
    parent_path=ROOT/"results/fits"/PARENT
    original=load_run_metadata(parent_path)
    old=problem_for(original)  # exact saved ML replay before making any output
    assert old.data.n_points==89
    witness=json.loads((ROOT/"results/diagnostics/dm_density_profile_20260922/core.json").read_text())
    states={"freecheck":witness["free"],
            "rho2":next(r for r in witness["profile"] if np.isclose(r["rho"],2.))}
    families={"freecheck":family_from_config(family_config(old.family)),
              "rho2":DensityConditionedModel(family_from_config(family_config(old.family)),2.)}
    jobs={}
    for i,(kind,family) in enumerate(families.items()):
        model=family_config(family)
        assert family_config(family_from_config(model))==model
        if kind=="freecheck":
            assert model==original["summary"]["model"]
        problem=FitProblem(family,old.data)
        rng=np.random.default_rng(987+i)
        cubes=rng.uniform(.001,.999,size=(128,len(family.names)))
        likelihoods=np.array([problem.loglike_unit(u) for u in cubes])
        assert np.isfinite(likelihoods).all() and likelihoods.min()>-1e50
        theta=states[kind]["parameters"]
        x=np.array([theta[n] for n in family.names])
        expected=old.loglike_vector(np.array([theta[n] for n in old.family.names]))
        actual=problem.loglike_vector(x)
        np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-7)
        if kind=="rho2":
            for u in cubes[:32]:
                t=family.complete(family.to_dict(family.transform(u)))
                np.testing.assert_allclose(old.family.halo(t).density(20.),2.,rtol=1e-12)
        label=labels[kind]
        jobs[label]=dict(kind=kind,model=model,output=str(ROOT/"results/fits"/label),
            data_fingerprint=data_fingerprint(old.data),dataset_options=original["summary"]["dataset_options"],
            density_condition=None if kind=="freecheck" else dict(radius_pc=20.,density=2.,
                units="Msun pc^-3",mass_support_envelope=family.mass_support_envelope,
                parent_prior_pdf=family.prior_density_at_constraint,
                prior="exact conditional of archived loguniform halo-mass prior"),
            sampler=dict(n_live=800,dlogz=.25,seed=1729+i,step_sampler=True,
                         slice_nsteps=4*len(family.names),resume="overwrite",n_profile_samples=300),
            preflight=dict(prior_likelihoods=len(cubes),all_finite=True,witness_loglike=actual,
                           witness_replay_difference=actual-expected),
            parent_sha256={str((parent_path/n).relative_to(ROOT)):sha256(parent_path/n)
                           for n in ("summary.json","run.yaml","data_snapshot.json","posterior.ecsv","ml_x.npy")})
        print(kind,len(family.names),"parameters; prior and best-fit replay checks passed",flush=True)
    batch.mkdir(parents=True,exist_ok=False)
    snapshot=write_data_snapshot(old.data,batch)
    code=batch/"code"
    for path in list((ROOT/"src").rglob("*.py"))+[ROOT/"pyproject.toml",Path(__file__).resolve()]:
        target=code/path.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target)
    for label,config in jobs.items():
        config["data_snapshot"]=snapshot
        (batch/f"{label}.json").write_text(json.dumps(config,indent=2)+"\n")
    env={"PYTHONPATH":str(code/"src"),"OCEN_DM_ROOT":str(ROOT),"PYTHONDONTWRITEBYTECODE":"1",
         "PYTHONUNBUFFERED":"1","MPLCONFIGDIR":"/private/tmp/ocen-nbody-mpl",
         **{name:"1" for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
                                  "VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS")}}
    manifest=dict(created_utc=now(),status="prepared",parent=PARENT,environment=env,jobs={},
        source_sha256={str(p.relative_to(code)):sha256(p) for p in code.rglob("*") if p.is_file()},
        job_config_sha256={label:sha256(batch/f"{label}.json") for label in jobs})
    path=batch/"launch.json"
    atomic_json(path,manifest)
    for label in jobs:
        command=[sys.executable,"-u",str(code/"bin/run_density_posterior_checks.py"),"worker",str(batch),label]
        with (batch/f"{label}.log").open("x") as log:
            proc=subprocess.Popen(command,cwd=ROOT,env={**os.environ,**env},stdin=subprocess.DEVNULL,
                stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        caffeine=subprocess.Popen(["/usr/bin/caffeinate","-i","-w",str(proc.pid)],stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True)
        manifest["jobs"][label]=dict(pid=proc.pid,caffeinate_pid=caffeine.pid,command=command)
        atomic_json(path,manifest)
        print(label,"PID",proc.pid,flush=True)
    manifest.update(status="launched",launched_utc=now())
    atomic_json(path,manifest)
    print("Launch record:",path,flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action",choices=["launch","worker","compare"])
    parser.add_argument("batch",type=Path)
    parser.add_argument("label",nargs="?")
    args=parser.parse_args()
    batch=args.batch.resolve()
    if args.action=="launch":
        launch(batch)
    elif args.action=="worker":
        worker(batch,args.label)
    else:
        compare_if_complete(batch)
