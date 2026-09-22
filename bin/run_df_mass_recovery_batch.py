#!/usr/bin/env python3
"""Preserved, bounded mass-recovery batch with a noiseless gate before noisy fits.

Two worker processes by default. Every job uses a frozen copy of the scientific
source. A noisy job is queued only after a numerically validated noiseless fit
passes the previously defined capacity threshold for its injected cluster.
"""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def dump(path,row):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
    temp.replace(path)


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--noise-realizations',type=int,default=8)
    p.add_argument('--max-evaluations',type=int,default=60)
    args=p.parse_args()
    if not 1<=args.workers<=3 or args.noise_realizations<0 or args.max_evaluations<1:
        p.error('invalid worker count/noise count/evaluation budget')
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    code=out/'code';(code/'bin').mkdir(parents=True)
    shutil.copytree(ROOT/'src',code/'src',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for file in ['run_df_mass_recovery.py','run_df_mass_recovery_batch.py']:
        shutil.copy2(ROOT/'bin'/file,code/'bin'/file)
    # Preserve the tests used to validate this implementation as well.
    (code/'tests').mkdir()
    for file in ['test_df_mass_recovery.py','test_positive_df.py','test_df_capacity.py']:
        shutil.copy2(ROOT/'tests'/file,code/'tests'/file)
    inputs=out/'inputs';(inputs/'data_run').mkdir(parents=True)
    for name in ['summary.json','data_snapshot.json','photometry.json']:
        shutil.copy2(ROOT/'results/df/pilot_no_dm_20260922'/name,inputs/'data_run'/name)
    for case in ['mixed_no_dm','mixed_with_dm']:
        (inputs/case).mkdir()
        shutil.copy2(ROOT/'results/df'/f'capacity_20260922_{case}'/'summary.json',inputs/case/'summary.json')
    manifest=dict(status='running',controller_pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),
                  workers=args.workers,noise_realizations_per_case=args.noise_realizations,
                  seed_base=23092026,max_evaluations_per_stage=args.max_evaluations,
                  source_sha256={str(f.relative_to(code)):checksum(f) for f in code.rglob('*') if f.is_file()},
                  input_sha256={str(f.relative_to(inputs)):checksum(f) for f in inputs.rglob('*') if f.is_file()},
                  gate='best noiseless fit: optimizer terminated, numerical validation passed, RMS<0.1 and max<0.3 adopted errors',
                  uncertainty_scope='pilot noise-realization recovery scatter; no posterior interval coverage claim')
    dump(out/'batch.json',manifest)
    env=os.environ.copy()
    for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        env[name]='1'
    env['PYTHONDONTWRITEBYTECODE']='1';env['MPLCONFIGDIR']='/private/tmp/ocen-df-mpl'
    cases=['mixed_no_dm','mixed_with_dm']
    jobs=[]
    def job(case,start,seed=None,initial=None,shared=False):
        tag=f'{case}_'+('noiseless' if seed is None else f'noise{seed}')+f'_start{start}'+('_shared_ml' if shared else '')
        directory=out/tag
        cmd=[sys.executable,str(code/'bin/run_df_mass_recovery.py'),'--case',case,'--out',str(directory),
             '--data-run',str(inputs/'data_run'),
             '--capacity-run',str(inputs/case),
             '--start',str(start),'--max-evaluations',str(args.max_evaluations)]
        if seed is not None:cmd+=['--noise-seed',str(seed)]
        if initial is not None:cmd+=['--initial-model',str(initial)]
        if shared:cmd+=['--shared-ml']
        return dict(id=tag,case=case,start=start,noise_seed=seed,shared_ml=shared,
                    directory=str(directory),command=cmd,status='queued')
    # Interleave cases so both have early results.
    for start in range(3):
        jobs.extend(job(case,start) for case in cases)
    active={}
    def run_queue():
        while any(j['status']=='queued' for j in jobs) or active:
            for j in jobs:
                if len(active)>=args.workers:break
                if j['status']!='queued':continue
                stream=(out/f"{j['id']}.log").open('w')
                process=subprocess.Popen(j['command'],cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
                stream.close()
                j.update(status='running',pid=process.pid,started_utc=datetime.now(timezone.utc).isoformat())
                active[j['id']]=process
                print(f"Started {j['id']} pid={process.pid}",flush=True)
            for j in jobs:
                if j['id'] not in active:continue
                rc=active[j['id']].poll()
                if rc is not None:
                    j.update(status='completed' if rc==0 else 'failed',exit_code=rc,finished_utc=datetime.now(timezone.utc).isoformat())
                    active.pop(j['id'])
                    print(f"Finished {j['id']} exit={rc}",flush=True)
            dump(out/'jobs.json',jobs)
            if active:time.sleep(3)
    try:
        run_queue()
        gates={}
        for case in cases:
            candidates=[]
            for j in jobs:
                if j['case']!=case or j['status']!='completed':continue
                s=json.loads((Path(j['directory'])/'summary.json').read_text())
                last=s['stages'][-1]
                passed=(last['stage']=='full' and last['optimizer']['success'] and s['validation']['passed']
                        and s['validation']['refined_rms_residual_sigma']<.1
                        and s['validation']['refined_max_residual_sigma']<.3)
                if passed:candidates.append((s['best']['score'],j))
            gates[case]=dict(passed=bool(candidates))
            if not candidates:continue
            winner=min(candidates,key=lambda item:item[0])[1]
            initial=Path(winner['directory'])/'best_model.json'
            gates[case].update(best_run=winner['directory'],initial_model=str(initial))
            # Shared-M/L control uses the same starting shapes and mass scales.
            jobs.append(job(case,0,initial=initial,shared=True))
            for realization in range(args.noise_realizations):
                seed=manifest['seed_base']+1000*cases.index(case)+realization
                jobs.extend(job(case,start,seed,initial) for start in range(2))
        manifest['noiseless_gate']=gates
        dump(out/'batch.json',manifest)
        run_queue()
        final_status=('completed_with_failed_jobs' if any(j['status']=='failed' for j in jobs) else 'completed')
        manifest.update(status=final_status if all(g['passed'] for g in gates.values()) else 'noiseless_gate_not_passed',
                        completed_jobs=sum(j['status']=='completed' for j in jobs),failed_jobs=sum(j['status']=='failed' for j in jobs),
                        finished_utc=datetime.now(timezone.utc).isoformat())
        dump(out/'batch.json',manifest)
    except BaseException as exc:
        # Workers keep their own provenance; do not silently terminate other fits.
        manifest.update(status='controller_failed',error=repr(exc))
        dump(out/'batch.json',manifest)
        raise


if __name__=='__main__':main()
