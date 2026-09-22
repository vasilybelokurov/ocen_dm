#!/usr/bin/env python3
"""Preserved, staged lifetime experiment queue. A failed numerical gate stops it."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import yaml
from astropy.table import Table

from ocen_dm.dynamics import TIME_UNIT_MYR
from ocen_dm.dynamics.config import Experiment, Orbit
from ocen_dm.dynamics.experiment import prepare, run
from ocen_dm.dynamics.orbits import end_time


def stamp():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def specifications(root):
    base = Experiment.load(root/'configs/dynamics/remnant_cusp.yaml')
    passage = Experiment.load(root/'results/dynamics/remnant_passage_20260921T114939Z/remnant_cusp_passage/config.yaml')
    jobs = []
    bridge_duration = end_time(replace(base, orbit=Orbit(duration_myr=500.),
                                       integrator=replace(base.integrator, max_step_myr=.01)))*TIME_UNIT_MYR
    def add(stage, label, config, engine='live'):
        jobs.append(dict(stage=stage, label=label, config=config.to_dict(), engine=engine, status='pending'))
    for shape, gamma in [('cusp', 1.), ('core', 0.)]:
        c = replace(base, components=(replace(base.components[0], gamma=gamma),))
        for control, count, step in [('base', 100000, .02), ('dt_half', 100000, .01), ('n2', 200000, .01)]:
            # Identical full-passage endpoints, independent of rounded timestep.
            for host in ('isolation', 'class1'):
                orbit = replace(passage.orbit, kind=host, duration_myr=end_time(passage)*TIME_UNIT_MYR,
                                initial_xv=passage.orbit.initial_xv if host == 'class1' else None)
                case = replace(c, orbit=orbit, components=(replace(c.components[0], n_particles=count),),
                               integrator=replace(c.integrator, max_step_myr=step, output_step_myr=.5))
                add('passage', f'{shape}_{control}_{host}', case)
        # Particle mass and timestep controls of a responsive nucleus, plus seed spread.
        for control, count, seed, step in [('base', 100000, 42, .01), ('n2', 200000, 42, .01),
                                          ('seed43', 100000, 43, .01), ('dt_half', 100000, 42, .005)]:
            for host in ('isolation', 'class1'):
                case = replace(c, seed=seed,
                               components=(replace(c.components[0], n_particles=count),),
                               nucleus=replace(c.nucleus, live=True, beta=-.1, n_particles=count),
                               orbit=Orbit(kind=host, duration_myr=bridge_duration, sample_step_myr=.02,
                                           end_at_present=host == 'class1'),
                               integrator=replace(c.integrator, max_step_myr=step, levels=5,
                                                  output_step_myr=2., accuracy=1e-11))
                add('bridge', f'{shape}_{control}_{host}_live_nucleus', case)
        # Conditional fixed-potential long controls. A live lifetime grid is
        # intentionally left for review of the measured bridge and heating.
        for history, duration in [('static', 10000.), ('bar', 8*TIME_UNIT_MYR)]:
            for host in ('isolation', 'class1' if history == 'static' else 'class3'):
                case = replace(c, orbit=Orbit(kind=host, duration_myr=duration, sample_step_myr=.02,
                                              end_at_present=host != 'isolation'),
                               integrator=replace(c.integrator, max_step_myr=.01, output_step_myr=5.,
                                                  accuracy=1e-12))
                add('lifetime_frozen', f'{shape}_{history}_{host}_frozen', case, 'frozen')
    # Within a stage alternate cusp/core so concurrent slots cover both.
    return sorted(jobs, key=lambda j: (['passage', 'bridge', 'lifetime_frozen'].index(j['stage']),
                                      j['label'].split('_', 1)[1], j['label']))


def load_table(batch, label):
    path = batch/'runs'/label/'evolution_live'
    record = yaml.safe_load((path/'run.yaml').read_text())
    data = path/'diagnostics.ecsv'
    if record['status'] != 'complete' or hashlib.sha256(data.read_bytes()).hexdigest() != record['diagnostics_sha256']:
        raise ValueError(f'Incomplete or changed diagnostics: {label}')
    return Table.read(data)


def gate(batch, stage):
    """Conservative numerical screening, not observational survival thresholds."""
    checks = []
    def check(name, value, limit):
        checks.append(dict(name=name, value=float(value), limit=limit,
                           passed=bool(np.isfinite(value) and abs(value) <= limit)))
    for shape in ('cusp', 'core'):
        suffix = '_live_nucleus' if stage == 'bridge' else ''
        comparisons = ['base', 'dt_half', 'n2'] + (['seed43'] if stage == 'bridge' else [])
        changes = {}
        for control in comparisons:
            iso = load_table(batch, f'{shape}_{control}_isolation{suffix}')
            tide = load_table(batch, f'{shape}_{control}_class1{suffix}')
            # Broad late phase window reduces dependence on a single particle crossing.
            key = 'dm_20pc_rho_msun_pc3'
            tail = np.asarray(tide['elapsed_myr']) >= .8*float(tide['elapsed_myr'][-1])
            iso_rho = np.interp(tide['elapsed_myr'], iso['elapsed_myr'], iso[key])
            changes[control] = float(np.median(np.asarray(tide[key])[tail]/iso_rho[tail]-1))
            check(f'{shape}/{control}: minimum DM20 shell count shortfall',
                  max(0., 500-float(np.min(tide['dm_20pc_shell_neff']))), 0.)
            check(f'{shape}/{control}: isolation DM20 drift',
                  np.median(iso[key][np.asarray(iso['elapsed_myr']) >= .8*iso['elapsed_myr'][-1]])/iso[key][0]-1, .10)
            if stage == 'bridge':
                for key, tolerance in [('nucleus_half_mass_radius_pc', .05), ('nucleus_bound_mass_msun', .02),
                                       ('nucleus_sigma_x_kms', .05), ('nucleus_sigma_y_kms', .05),
                                       ('nucleus_sigma_z_kms', .05)]:
                    check(f'{shape}/{control}: isolation {key} drift',
                          np.median(iso[key][-max(3, len(iso)//5):])/iso[key][0]-1, tolerance)
        for control in comparisons[1:]:
            check(f'{shape}: DM20 tidal contrast {control} minus base',
                  changes[control]-changes['base'], .02 if control == 'dt_half' else .05)
        if stage == 'bridge':
            base = load_table(batch, f'{shape}_base_class1{suffix}')
            for control in ['n2', 'dt_half']:
                other = load_table(batch, f'{shape}_{control}_class1{suffix}')
                for key in ['nucleus_half_mass_radius_pc', 'nucleus_bound_mass_msun']:
                    check(f'{shape}: tidal {key} {control}/base-1',
                          np.median(other[key][-20:])/np.median(base[key][-20:])-1, .05)
    return dict(stage=stage, checked_utc=stamp(), passed=all(c['passed'] for c in checks), checks=checks,
                interpretation='Numerical screening only; passing does not establish observational compatibility.')


def worker(batch, label):
    record = json.loads((batch/'batch.json').read_text())
    job = next(j for j in record['jobs'] if j['label'] == label)
    from ocen_dm.kinematics.run_io import read_data_snapshot
    data = read_data_snapshot(batch, record['stellar_observations'])
    c = Experiment.from_dict(job['config'])
    prepared = prepare(c, batch/'runs'/label, stellar_data=data if c.nucleus.live else None)
    run(prepared, job['engine'])


def execute(batch, workers):
    manifest = batch/'batch.json'
    record = json.loads(manifest.read_text())
    if record['status'] != 'ready':
        raise ValueError('Queue already started; preserve it and inspect before a new launch')
    for name, digest in record['source_sha256'].items():
        if hashlib.sha256((batch/'code'/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Archived queue source changed: {name}')
    record.update(status='running', pid=os.getpid(), started_utc=stamp())
    save(manifest, record)
    def child(job):
        env = os.environ.copy()
        env['PYTHONPATH'] = str(batch/'code/src')
        env['OCEN_DM_ROOT'] = record['project_root']
        with (batch/'logs'/f"{job['label']}.log").open('xb') as log:
            return subprocess.run([sys.executable, batch/'code/bin/run_lifetime_batch.py', 'worker',
                                   str(batch), '--label', job['label']], env=env, stdout=log,
                                  stderr=subprocess.STDOUT).returncode
    try:
        for stage in ('passage', 'bridge', 'lifetime_frozen'):
            record['active_stage'] = stage
            with ThreadPoolExecutor(max_workers=workers) as pool:
                pending = {pool.submit(child, j): j for j in record['jobs'] if j['stage'] == stage}
                for job in pending.values():
                    job['status'] = 'submitted'
                save(manifest, record)
                for future in as_completed(pending):
                    job = pending[future]
                    code = future.result()
                    job.update(status='complete' if code == 0 else 'failed', returncode=code, finished_utc=stamp())
                    save(manifest, record)
                    print(stamp(), job['label'], job['status'], flush=True)
            if any(j['status'] == 'failed' for j in record['jobs'] if j['stage'] == stage):
                raise RuntimeError(f'{stage} contains failed runs; later stages not launched')
            if stage != 'lifetime_frozen':
                assessment = gate(batch, stage)
                save(batch/f'{stage}_gate.json', assessment)
                if not assessment['passed']:
                    raise RuntimeError(f'{stage} numerical gate failed; later stages not launched')
        record.update(status='complete', completed_utc=stamp())
    except BaseException as exc:
        record.update(status='stopped', error=f'{type(exc).__name__}: {exc}', stopped_utc=stamp())
        raise
    finally:
        save(manifest, record)


def initialise(batch, root):
    batch.mkdir(parents=True, exist_ok=False)
    for name in ('runs', 'logs', 'code/bin'):
        (batch/name).mkdir(parents=True)
    shutil.copytree(root/'src', batch/'code/src', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copyfile(root/'pyproject.toml', batch/'code/pyproject.toml')
    for filename in ['run_lifetime_batch.py', 'nemo_env.sh', 'build_nbody_keys.sh', 'build_nbody_kernel.py', 'ocen_keep_keys.cc']:
        shutil.copyfile(root/'bin'/filename, batch/'code/bin'/filename)
    shutil.copytree(root/'bin/nbody_compat', batch/'code/bin/nbody_compat')
    source = {str(p.relative_to(batch/'code')): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted((batch/'code').rglob('*')) if p.is_file()}
    fit = root/'results/fits/rung1_K2_cored_slice'
    observations = yaml.safe_load((fit/'run.yaml').read_text())['data_snapshot']
    shutil.copyfile(fit/observations['file'], batch/observations['file'])
    save(batch/'batch.json', dict(schema_version=1, status='ready', created_utc=stamp(),
                                project_root=str(root), source_sha256=source,
                                stellar_observations=observations, jobs=specifications(root)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['init', 'execute', 'worker'])
    parser.add_argument('batch', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--label')
    args = parser.parse_args()
    batch = args.batch.resolve()
    if args.workers < 1 or args.workers > 8:
        parser.error('workers must be between 1 and 8')
    if args.action == 'init':
        initialise(batch, Path(__file__).resolve().parents[1])
    elif args.action == 'worker':
        worker(batch, args.label)
    else:
        execute(batch, args.workers)
