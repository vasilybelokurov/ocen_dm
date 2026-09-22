#!/usr/bin/env python3
"""Report measured lifetime-batch results, orbital forcing and local heating."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yaml
from astropy.table import Table

from ocen_dm.dynamics import TIME_UNIT_MYR
from ocen_dm.dynamics.config import Experiment
from ocen_dm.dynamics.experiment import load_prepared
from ocen_dm.dynamics.lifetime import heating_time_gyr
from ocen_dm.dynamics.models import agama_module, xyz
from ocen_dm.dynamics.orbits import build_track, Track


def report(batch, output):
    state = json.loads((batch/'batch.json').read_text())
    pause_path = batch/'pause.json'
    pause = json.loads(pause_path.read_text()) if pause_path.exists() else {}
    paused = pause.get('status') == 'paused'
    queue_status = 'paused by user (SIGSTOP; see pause.json)' if paused else state['status']
    output.mkdir(parents=True, exist_ok=True)
    assets = Path(__file__).resolve().parents[1]/'plots'
    assets.mkdir(exist_ok=True)
    data_output = Path(__file__).resolve().parents[1]/'results/plot_data'
    data_output.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat()
    records, completed, hashes = [], {}, {}
    for job in state['jobs']:
        prepared = batch/'runs'/job['label']
        evolution = prepared/('evolution_'+job['engine'])
        manifest = evolution/'run.yaml'
        status = job['status']
        if manifest.exists():
            record = yaml.safe_load(manifest.read_text())
            status = record['status']
            if status == 'complete':
                path = evolution/'diagnostics.ecsv'
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                if digest != record['diagnostics_sha256']:
                    raise ValueError(f'Changed diagnostics: {path}')
                hashes[str(path)] = digest
                completed[job['label']] = Table.read(path)
        elif (prepared/'run.yaml').exists():
            status = yaml.safe_load((prepared/'run.yaml').read_text())['status']
        elif status == 'submitted':
            status = 'queued'
        if paused and status == 'running':
            status = 'paused'
        records.append((job['stage'], job['label'], status))

    # Reference tracks only: do not confuse long orbital histories with evolved particles.
    history_records = {}
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, constrained_layout=True)
    for history, host, color in [('static', 'class1', '#2465a5'), ('bar', 'class3', '#c85925')]:
        label = f'cusp_{history}_{host}_frozen'
        c = Experiment.from_dict(next(j['config'] for j in state['jobs'] if j['label'] == label))
        path = batch/'reference_histories'/history
        if not path.exists():
            path.mkdir(parents=True)
            track, record = build_track(c, path)
            (path/'history.json').write_text(json.dumps(record, indent=2)+'\n')
        else:
            track, record = Track.load(path), json.loads((path/'history.json').read_text())
        time = (track.time-track.time[-1])*TIME_UNIT_MYR/1000
        radius = np.linalg.norm(track.xv[:, :3], axis=1)
        peri = np.where((radius[1:-1] < radius[:-2]) & (radius[1:-1] < radius[2:]))[0]+1
        z = track.xv[:, 2]
        crossings = np.where(z[1:]*z[:-1] < 0)[0]
        record.update(pericentres=len(peri), disc_crossings=len(crossings),
                      minimum_radius_kpc=float(radius.min()), maximum_radius_kpc=float(radius.max()))
        history_records[history] = record
        axes[0].plot(time[peri], radius[peri], '.', color=color,
                     label=f'{history}: {len(peri)} pericentres')
        axes[1].step(time[crossings], np.arange(1, len(crossings)+1), where='post', color=color,
                     label=f'{history}: {len(crossings)} crossings')
        hashes[str(path/'orbit.npz')] = hashlib.sha256((path/'orbit.npz').read_bytes()).hexdigest()
    axes[0].set(ylabel='Pericentre radius [kpc]', title='Reference orbital histories ending at the present epoch')
    axes[1].set(xlabel='Time relative to present [Gyr]', ylabel='Cumulative disc crossings')
    for ax in axes:
        ax.legend(); ax.grid(alpha=.2)
    fig.savefig(assets/'lifetime_orbital_forcing.png', dpi=180)
    plt.close(fig)

    # Local physical relaxation estimate in the actual initial DM equilibrium.
    # The physical mean stellar mass and effective mass are explicit assumptions.
    heating_rows = []
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    agama = agama_module()
    for shape, color in [('cusp', '#2465a5'), ('core', '#c85925')]:
        prepared = batch/'runs'/f'{shape}_base_isolation'
        if not (prepared/'equilibrium.json').exists():
            continue
        c, _, _, _ = load_prepared(prepared)
        potential = agama.Potential(str(prepared/'satellite_initial.ini'))
        dm = agama.Potential(str(prepared/'dm_potential.ini'))
        stars = agama.Potential(str(prepared/'nucleus.ini'))
        df = agama.DistributionFunction(type='QuasiSpherical', density=dm, potential=potential, beta0=0.)
        radius = np.unique(np.r_[np.geomspace(.3, 70., 150), c.radii_pc])
        _, moments = agama.GalaxyModel(potential, df).moments(xyz(radius/1000), dens=True, vel=False, vel2=True)
        vrms = np.sqrt(moments.sum(axis=1))
        density = stars.density(xyz(radius/1000))/1e9
        coulomb = np.log(.4*np.array([stars.enclosedMass(r/1000) for r in radius])/.5)
        times = [heating_time_gyr(vrms, density, mass, coulomb) for mass in [.3, 1., 3.]]
        ax.fill_between(radius, times[2], times[0], alpha=.15, color=color)
        ax.loglog(radius, times[1], color=color, label=shape)
        for i, r in enumerate(radius):
            heating_rows.append(dict(shape=shape, radius_pc=float(r), vrms_dm_kms=float(vrms[i]),
                 stellar_density_msun_pc3=float(density[i]), coulomb_log=float(coulomb[i]),
                 heat_gyr_meff03=float(times[0][i]), heat_gyr_meff1=float(times[1][i]), heat_gyr_meff3=float(times[2][i])))
    ax.axhspan(7.82234, 10, color='grey', alpha=.15, label='7.82–10 Gyr exposure')
    ax.axvline(20., color='grey', ls=':', lw=1)
    ax.set(xlabel='Radius [pc]', ylabel='Local DM heating time [Gyr]',
           title='Initial equilibrium; physical stellar masses, not simulation particles')
    ax.text(.03, .97, r'Line: effective stellar mass = 1 M$_\odot$'+'\n'+r'Band: 0.3–3 M$_\odot$; mean mass = 0.5 M$_\odot$',
            va='top', transform=ax.transAxes, fontsize=10)
    ax.legend(loc='lower right'); ax.grid(alpha=.2)
    if heating_rows:
        fig.savefig(assets/'lifetime_physical_heating.png', dpi=180)
        Table(rows=heating_rows).write(data_output/'lifetime_physical_heating.ecsv', overwrite=True)
    plt.close(fig)

    if completed:
        fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
        for label, table in completed.items():
            t = table['elapsed_myr']/1000
            # Only base cases in summary plots; all controls remain in tables/gates.
            if '_base_' not in label and not label.endswith('_frozen'):
                continue
            line, = axes[0, 0].plot(t, table['dm_20pc_rho_msun_pc3'], label=label)
            color = line.get_color()
            axes[0, 1].plot(t, table['dm_20pc_rho_msun_pc3']/table['dm_20pc_rho_msun_pc3'][0], color=color)
            if 'nucleus_bound_fraction' in table.colnames:
                axes[1, 0].plot(t, table['nucleus_bound_mass_msun']/1e6, color=color)
                axes[1, 1].plot(t, table['nucleus_half_mass_radius_pc'], color=color)
        for ax, y in zip(axes.ravel(), [r'DM shell density at 20 pc [M$_\odot$/pc$^3$]',
                'DM shell density / initial', r'Bound stellar mass [10$^6$ M$_\odot$]', 'Stellar 3D half-mass radius [pc]']):
            ax.set(xlabel='Elapsed time [Gyr]', ylabel=y); ax.grid(alpha=.2)
        axes[0, 0].legend(fontsize=6)
        fig.suptitle('Completed integrations: DM response and stellar evolution')
        fig.savefig(assets/'lifetime_measured_evolution.png', dpi=180)
        plt.close(fig)

    lines = ['# Lifetime evolution: implementation and running batch', '', f'Updated {now}.', '',
        f'Batch: `{batch.name}`. Queue state: **{queue_status}**; active stage: **{state.get("active_stage", "not started")}**.', '',
        'The new calculations follow accumulated forcing without resetting the particle distribution between passages. '
        'The live nucleus can expand and lose stars; diagnostics follow its measured centre. '
        'These runs are controlled experiments. No multi-Gyr particle result or observationally accepted survival history is available yet.' if not any('frozen' in k for k in completed) else
        'Only completed, checksum-verified integrations enter the measured plots. Frozen-potential histories remain conditional controls.', '',
        '| Stage | Cases | Purpose |', '|---|---:|---|',
        '| Full passage, 88.102 Myr | 12 | Cusp/core × tide/isolation × original timestep, half timestep, twice the particles |',
        '| Live nucleus, 500.004 Myr | 16 | Cusp/core × tide/isolation × baseline, twice the particles, second seed, half timestep |',
        '| Frozen potential, 7.822/10 Gyr | 8 | Cusp/core × bar/static × tide/matched isolation; conditional on both earlier numerical checks |', '',
        'Four workers run concurrently. Baseline live-nucleus cases contain 100,000 DM particles and 100,000 stellar particles; '
        'the particle-number controls double both. Stellar mass is 3.55 million solar masses and the initial Plummer scale is 5.4 pc. '
        'The stellar DF has β = −0.1 in both halos: the isotropic Plummer nucleus fails the signed-DF check in the cusp. '
        'This simulation choice is separate from the fitted rung-1 anisotropy.', '',
        'Passing a numerical gate only permits the next computation. It does not establish stellar survival or agreement with the observations. '
        'The queue checks matched density contrasts, isolation drift, shell counts and stellar mass/size/dispersion stability. '
        'A failed integration or failed gate stops advancement and preserves its products.', '',
        'Validation: 54 focused tests pass across the dynamics, kernel and lifetime test modules, including compiled N-body '
        'checks in isolation and an external host, nucleus gravity counted once, stellar-centre translation/boost invariance, '
        'DF admissibility, physical-time units, sparse observational bins, and rejection of drift or changed gate inputs. '
        'Long-term convergence is what the running batch will measure; these tests establish the implementation checks.', '',
        '## Present-epoch orbital forcing', '',
        f'![Orbital forcing]({os.path.relpath(assets/"lifetime_orbital_forcing.png", output)})', '',
        '| History | Duration [Gyr] | Pericentres | Disc crossings | Smallest pericentre [kpc] | Endpoint error [pc; km/s] |',
        '|---|---:|---:|---:|---:|---|']
    for name, r in history_records.items():
        lines.append(f'| {name} | {r["effective_duration_myr"]/1000:.5f} | {r["pericentres"]} | {r["disc_crossings"]} | '
            f'{r["minimum_radius_kpc"]:.3f} | {r["present_endpoint_position_error_pc"]:.4g}; {r["present_endpoint_velocity_error_kms"]:.4g} |')
    lines += ['', 'These are integrated reference orbits, not particle-density histories. The static host repeats today’s potential; '
              'the bar uses its stored absolute clock and Ω_final = 24 km/s/kpc. Neither includes a host wake, progenitor mass loss or dynamical friction.', '',
              '## Physical stellar heating', '',
              f'![Heating estimate]({os.path.relpath(assets/"lifetime_physical_heating.png", output)})', '',
              'The local estimate uses equation 11 of [Bertone & Fairbairn (2008)](https://doi.org/10.1103/PhysRevD.77.043515): '
              'T_heat = 0.814 v_rms³ / (G² m_eff ρ_star ln Λ). The initial DM equilibrium supplies v_rms; '
              'the Plummer nucleus supplies ρ_star. We use ln Λ = ln[0.4 M_star(<r)/0.5 M_sun] and bracket '
              'm_eff = <m²>/<m> between 0.3 and 3 M_sun. These are illustrative mass-spectrum assumptions, not measurements.', '',
              '| Halo | Radius [pc] | Heating time, m_eff=1 [Gyr] | Range, m_eff=3 to 0.3 [Gyr] |', '|---|---:|---:|---|']
    for r in heating_rows:
        if r['radius_pc'] in [3., 10., 20.]:
            lines.append(f'| {r["shape"]} | {r["radius_pc"]:g} | {r["heat_gyr_meff1"]:.1f} | {r["heat_gyr_meff3"]:.1f}–{r["heat_gyr_meff03"]:.1f} |')
    lines += ['', 'This formula assumes comparable stellar and DM velocity scales. It is a local timescale, not a predicted density loss. '
              'Particles observed at 20 pc can visit the denser centre, so orbital averaging and a calibrated collisional treatment remain necessary '
              'before claiming that physical relaxation is negligible. Simulation-particle scattering is a numerical effect; '
              'the particle-number and isolation controls test it separately.', '',
              '## Stellar-survival constraint', '',
              'Each live-nucleus evolution records bound stellar mass, 3D and three projected half-mass radii, velocity dispersions, '
              'central density and displacement from the reference orbit. Initial and final projected kinematic profiles are also binned '
              'against the archived HST, MUSE and Gaia observations; bins with fewer than 30 effective particles remain unmeasured.', '',
              'The three Cartesian views are diagnostic projections with complete annular selection and constant stellar M/L. '
              'They do not replace the measured selection function, surface-brightness fit or a rotation model. '
              'The current Plummer initial conditions are not fitted to those observations. We therefore record response and residuals '
              'without assigning an arbitrary “survived” label. Observationally selected long-lived histories still require a calibrated '
              'stellar initial-condition family, physical stellar/remnant evolution and collisional evolution where needed.', '',
              '## Run state', '', '| Stage | Case | State |', '|---|---|---|']
    lines += [f'| {stage} | `{label}` | {status} |' for stage, label, status in records]
    if completed:
        lines += ['', f'![Measured evolution]({os.path.relpath(assets/"lifetime_measured_evolution.png", output)})']
    lines += ['', 'Refresh this report with:', '', '```sh',
              f'PYTHONPATH=src python bin/report_lifetime_batch.py {batch.relative_to(Path.cwd()) if batch.is_relative_to(Path.cwd()) else batch}',
              '```', '', 'Raw configurations, source hashes, logs and gate decisions live in the batch directory. '
              'The controller executes its archived source copy, so later edits do not alter queued physics.', '']
    (output/'LIFETIME_BATCH_REPORT.md').write_text('\n'.join(lines))
    (data_output/'lifetime_provenance.json').write_text(json.dumps(dict(updated_utc=now, batch=str(batch),
        input_sha256=hashes, reference_histories=history_records,
        report_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')
    print(output/'LIFETIME_BATCH_REPORT.md')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batch', type=Path)
    parser.add_argument('--output', type=Path, default=Path('docs'))
    args = parser.parse_args()
    report(args.batch.resolve(), args.output.resolve())
