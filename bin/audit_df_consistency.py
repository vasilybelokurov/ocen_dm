#!/usr/bin/env python3
"""Audit archived Jeans fits without changing runs or assigning a unique DF.

Outputs are necessary-condition screens, never certificates of DF positivity.
Run again after active fits finish to include their completed posteriors.
"""
from datetime import datetime, timezone
import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
import yaml
from astropy.table import Table

from ocen_dm.kinematics.df_consistency import necessary_profiles, summarize_profiles
from ocen_dm.kinematics.fit import FitProblem
from ocen_dm.kinematics.report import load_run_metadata, _family_for, data_for
from ocen_dm.kinematics.run_io import family_from_config, read_data_snapshot, sha256

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/plot_data'
ARCSEC = 206264.806


def aperture(data):
    representatives = np.concatenate([p.r for p in data.profiles])
    nodes = np.concatenate([p.r_nodes.ravel() if p.r_nodes is not None else p.r
                            for p in data.profiles])
    outer = max(float(np.max(p.r_upper if p.has_edges else p.r)) for p in data.profiles)
    inner = min(float(np.min(p.r_lower if p.has_edges else p.r)) for p in data.profiles)
    return dict(representative_min_arcsec=float(representatives.min()),
                representative_max_arcsec=float(representatives.max()),
                node_min_arcsec=float(nodes[nodes > 0].min()),
                node_max_arcsec=float(nodes.max()), inner_edge_arcsec=inner,
                outer_edge_arcsec=outer)


def resolve(run):
    s = run['summary']
    if run['dir'].name.startswith('preset_') and 'model' not in s:
        from ocen_dm.kinematics.presets import build_preset
        family = build_preset(run['dir'].name.removeprefix('preset_'))[1]
    else:
        options = dict(run.get('run', {}).get('dataset_options', {}))
        # This legacy run predates saved resolved models; its name and the absence
        # of all s_instrument parameters specify the no-scale variant.
        if run['dir'].name == 'K1_noscale':
            options['no_scales'] = True
        family = _family_for(s, options=options)
    if set(family.names) != set(s['parameters']):
        raise ValueError('saved parameter names do not match reconstructed family')
    data, _ = data_for(s, run_dir=run['dir'])
    return family, data


def profile_replay(family, directory, table):
    """Verify the physical reconstruction against *every* archived profile draw.

    A legacy likelihood may have changed while density/force/anisotropy remain
    reproducible. Do not bypass its failed likelihood check for a fit comparison.
    """
    path = directory/'profiles.npz'
    if not path.exists():
        return {'status': 'no_saved_profiles'}
    errors = []
    with np.load(path) as saved:
        for j, index in enumerate(saved['sample_index']):
            theta = {n: float(table[n][int(index)]) for n in family.names}
            model, _, _ = family.build(theta)
            actual = model.mass.radial_profile(saved['r'])
            for key in ('M_stars', 'M_remnants', 'M_imbh', 'M_total', 'M_dm', 'rho_total'):
                expected = saved[key][j]
                err = float(np.max(np.abs(actual[key]-expected)/np.maximum(np.abs(expected), 1e-30)))
                errors.append(err)
                np.testing.assert_allclose(actual[key], expected, rtol=2e-7, atol=1e-10,
                                           err_msg=f'{directory.name} {key} profile {j}')
        count = len(saved['sample_index'])
    return {'status': 'verified', 'draws': count, 'max_relative_error': max(errors),
            'sha256': sha256(path)}


def audit_point(family, theta, bounds, ngrid):
    model, distance, _ = family.build(theta)
    full = family.complete(theta)
    beta0 = float(model.anisotropy.beta_0)
    if float(full['M_bh']) <= 0:
        raise ValueError('this central test requires a positive, unsoftened point mass')
    rmin = bounds['node_min_arcsec']*distance*1000/ARCSEC
    rmax = bounds['outer_edge_arcsec']*distance*1000/ARCSEC
    r = np.geomspace(rmin, rmax, ngrid)
    profiles = necessary_profiles(model, r)
    detailed = summarize_profiles(profiles, beta0)
    fine = summarize_profiles(necessary_profiles(model, np.geomspace(rmin, rmax, 2*ngrid)), beta0)
    agree = all(detailed[k]['fails_necessary_condition'] == fine[k]['fails_necessary_condition']
                for k in ('B1','B2','R2_over_R'))
    if not agree:
        raise ValueError('finite-radius failure classification changes on grid refinement')
    inner = summarize_profiles(necessary_profiles(model, np.geomspace(1e-5, rmax, 2*ngrid)), beta0)
    point = dict(beta0=beta0, M_bh=float(full['M_bh']), gamma0=0.,
                 central_margin=-beta0-.5, central_condition_fails=bool(beta0 > -.5),
                 distance_kpc=distance, radial_scale_interval_pc=[rmin, rmax],
                 spatial_density_nonnegative=bool(np.all(profiles['rho_total'] >= 0)
                                                    and np.all(profiles['rho_stars'] >= 0)),
                 resolved_scale_screen=detailed, inner_to_edge_screen=inner,
                 refinement_classification_agrees=agree,
                 df_nonnegative_certified=False)
    return point, profiles


def audit_posterior(family, table, bounds, ngrid):
    """Audit every unique draw; preserve the equal-weight posterior multiplicities."""
    samples = np.column_stack([np.asarray(table[n], float) for n in family.names])
    if not np.all(np.isfinite(samples)):
        raise ValueError('nonfinite posterior parameters')
    unique, multiplicities = np.unique(samples, axis=0, return_counts=True)
    flags = np.zeros((len(unique), 4), bool)
    minima = np.zeros((len(unique), 3))
    central_beta = np.zeros(len(unique))
    grid = np.geomspace(bounds['node_min_arcsec'], bounds['outer_edge_arcsec'], ngrid)
    positive = True
    for i, sample in enumerate(unique):
        theta = dict(zip(family.names, sample))
        model, d, _ = family.build(theta)
        full = family.complete(theta)
        if full['M_bh'] <= 0:
            raise ValueError('positive BH required for central screen')
        central_beta[i] = model.anisotropy.beta_0
        p = necessary_profiles(model, grid*d*1000/ARCSEC)
        positive &= bool(np.all(p['rho_total'] >= 0) and np.all(p['rho_stars'] >= 0))
        minima[i] = [p[k].min() for k in ('B1','B2','R2_over_R')]
        flags[i] = [central_beta[i] > -.5,
                    central_beta[i] <= .5 and minima[i,0] < -1e-8,
                    central_beta[i] <= -.5 and minima[i,1] < -1e-8,
                    minima[i,2] < -1e-8]
    fraction = lambda b: float(multiplicities[b].sum()/multiplicities.sum())
    summary = dict(draws=len(samples), unique_draws=len(unique),
                   fraction_failing_central_condition=fraction(flags[:,0]),
                   fraction_failing_separable_B1=fraction(flags[:,1]),
                   fraction_failing_separable_B2=fraction(flags[:,2]),
                   fraction_failing_separable_R2=fraction(flags[:,3]),
                   fraction_failing_any_separable_test=fraction(flags[:,1:].any(axis=1)),
                   fraction_failing_either_screen=fraction(flags.any(axis=1)),
                   beta0_range=[float(central_beta.min()), float(central_beta.max())],
                   all_spatial_densities_nonnegative=positive,
                   df_nonnegative_certified=False)
    return summary, dict(theta=unique, multiplicity=multiplicities, failure_flags=flags,
                         minima=minima, parameter_names=np.array(family.names),
                         flag_names=np.array(['central','separable_B1','separable_B2','separable_R2']))


def scan_run(path, ngrid, posterior_grid):
    directory = path.parent
    label = str(directory.relative_to(ROOT/'results'))
    raw = path.read_bytes()
    metadata = yaml.safe_load(raw)
    entry = dict(label=label, run_sha256=hashlib.sha256(raw).hexdigest(),
                 status=metadata.get('status', 'complete_legacy'),
                 excluded_from_science='_invalid_' in label or '/diagnostics/' in str(directory),
                 posterior_available=(directory/'posterior.ecsv').exists())
    if not (directory/'summary.json').exists():
        family = family_from_config(metadata['model'])
        data = read_data_snapshot(directory, metadata['data_snapshot'])
        entry['aperture'] = aperture(data)
        entry['model_definition'] = 'saved_resolved_model'
        checkpoint = directory/'ultranest/extra/sampling-stuck-it%d.npz'
        if checkpoint.exists():
            with np.load(checkpoint) as z:
                theta = dict(zip(family.names, z['sample_v'][0]))
                # Verify physical vs unit-cube parameter conventions before use.
                np.testing.assert_allclose(family.transform(z['sample_u'][0]), z['sample_v'][0])
                u = z['u']
                if np.any((u < 0) | (u > 1)):
                    raise ValueError('checkpoint live points outside unit cube')
                values = np.array([family.transform(row) for row in u])
            entry['checkpoint_point'], _ = audit_point(family, theta, entry['aperture'], ngrid)
            entry['checkpoint_scope'] = 'stuck proposal and 400 saved live points; not a posterior or best fit'
            entry['checkpoint_sha256'] = sha256(checkpoint)
            live = Table(values, names=family.names)
            entry['checkpoint_live_screen'], _ = audit_posterior(family, live, entry['aperture'], posterior_grid)
        else:
            entry['pending'] = 'No completed posterior or saved physical checkpoint; audit after completion.'
        return entry, None
    run = load_run_metadata(directory)
    family, data = resolve(run)
    bounds = aperture(data)
    theta = {n: run['summary']['parameters'][n]['ml'] for n in family.names}
    table = Table.read(directory/'posterior.ecsv')
    entry.update(aperture=bounds,
                 aperture_source='saved_snapshot' if 'data_snapshot' in run['summary'] else 'current_legacy_data_loader',
                 model_definition='saved_resolved_model' if 'model' in run['summary'] else 'legacy_reconstruction_checked_against_archived_mass_profiles',
                 summary_sha256=sha256(directory/'summary.json'),
                 posterior_sha256=sha256(directory/'posterior.ecsv'))
    entry['physical_profile_replay'] = profile_replay(family, directory, table)
    ll = FitProblem(family, data).loglike_vector(np.array([theta[n] for n in family.names]))
    entry['likelihood_replay_delta'] = float(ll-run['summary']['lnL_max'])
    entry['likelihood_replays'] = bool(np.isclose(ll,run['summary']['lnL_max'],rtol=1e-9,atol=1e-7))
    if 'model' not in run['summary'] and entry['physical_profile_replay']['status'] != 'verified':
        raise ValueError('legacy physical reconstruction unverified')
    entry['best_saved_sample'], curves = audit_point(family, theta, bounds, ngrid)
    entry['posterior'], arrays = audit_posterior(family, table, bounds, posterior_grid)
    cache = ROOT/'results/diagnostics/df_consistency'
    cache.mkdir(parents=True, exist_ok=True)
    filename = label.replace('/','__')+'.npz'
    np.savez_compressed(cache/filename, **arrays)
    entry['posterior_screen_file'] = str((cache/filename).relative_to(ROOT))
    return entry, curves


def scan_optimizations(ngrid):
    entries = []
    base = ROOT/'results/diagnostics/dm_density_profile_20260922'
    restarts = json.loads((base/'high_density_multistart_check.json').read_text())
    for case in ('core','cusp','core_lower_stellar_floor'):
        path = base/(case+'.json')
        state = json.loads(path.read_text())
        family = family_from_config(state['model'])
        parent = load_run_metadata(ROOT/'results/fits'/state['label'])
        data, _ = data_for(parent['summary'], run_dir=parent['dir'])
        bounds = aperture(data)
        records = [('free',state['free'])]+[(f'rho{r["rho"]:g}',r) for r in state['profile']]
        records.append(('rho5_restart',restarts[case]['check']))
        for kind, record in records:
            point, _ = audit_point(family,record['parameters'],bounds,ngrid)
            source = base/'high_density_multistart_check.json' if kind == 'rho5_restart' else path
            entries.append(dict(label=f'{case}/{kind}',source=str(source.relative_to(ROOT)),
                                source_sha256=sha256(source), model_source=str(path.relative_to(ROOT)),
                                model_source_sha256=sha256(path), rho20=record['rho'], point=point))
    return entries


def make_plot(entries, curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
    labels = ['rung0_K1','rung1_K1_slice','rung2_K1_simple_beta','rung2_K1_turnover']
    names = ['Rung 0: isotropic','Rung 1: constant beta','Rung 2: simple beta(r)','Rung 2: turnover beta(r)']
    for label,name in zip(labels,names):
        p=curves['fits/'+label]
        axes[0].plot(p['r_pc'],p['B1'],label=name)
    for label,name in [('rung2_K1_turnover','No DM'),('rung2_K2_nfw_turnover','Cusped DM'),('rung2_K2_cored_turnover','Cored DM')]:
        p=curves['fits/'+label]
        axes[1].plot(p['r_pc'],p['B2'],label=name)
    finished=[e for e in entries if e['label'].startswith('fits/rung') and 'best_saved_sample' in e]
    for i,e in enumerate(finished):
        p=e['best_saved_sample']; label=e['label'].removeprefix('fits/')
        axes[2].plot(p['central_margin'],i,'o',color='#c54438' if p['central_condition_fails'] else '#297f56')
    axes[2].set_yticks(range(len(finished)),[e['label'].removeprefix('fits/').replace('rung','R').replace('_',' ') for e in finished],fontsize=8)
    axes[2].invert_yaxis();axes[2].axvline(0,color='k',ls=':',lw=1)
    axes[2].set_xlabel(r'Central margin $-\beta_0-1/2$');axes[2].set_title('Central black-hole condition')
    axes[0].set_title('No-DM best fits: first derivative')
    axes[0].set_ylabel(r'$B_1=\gamma-2\beta$ (sign of $P^\prime$)')
    axes[1].set_title('Turnover best fits: second derivative')
    axes[1].set_ylabel(r'$B_2$ (sign of $P^{\prime\prime}$)')
    for ax in axes[:2]:
        ax.axhline(0,color='k',ls=':',lw=1);ax.set_xscale('log');ax.set_xlabel('Radius [pc]')
        ax.set_xlim(.08,10);ax.set_ylim(-.2,1.2 if ax is axes[0] else .3)
        ax.legend(fontsize=8,loc='upper left');ax.grid(alpha=.2)
    fig.suptitle('Stellar DF necessary conditions: finite-radius panels assume separable augmented density',fontsize=12)
    fig.savefig(ROOT/'plots/df_consistency_audit.png',dpi=180)
    plt.close(fig)


def write_table(result):
    lines=['# DF audit: run inventory','',f'Generated {result["created_utc"]}.',
           '', 'A failed central condition rules out the unmodified central extrapolation. '
           'Finite-radius failures below refer **only to separable augmented-density DFs**. '
           'No entry is a certificate of a non-negative DF. Posterior percentages use all saved equal-weight draws.', '',
           '| Run | State | Central ML condition | Separable ML screen | Posterior central failures | Posterior separable failures |',
           '|---|---|---|---|---:|---:|']
    for e in result['runs']:
        if 'best_saved_sample' in e:
            p=e['best_saved_sample'];post=e['posterior']
            lines.append(f'| {e["label"]} | {e["status"]} | {"FAIL" if p["central_condition_fails"] else "passes necessary condition"} | '
                         f'{"FAIL" if p["resolved_scale_screen"]["separable_df_ruled_out"] else "no failure found"} | '
                         f'{100*post["fraction_failing_central_condition"]:.3f}% | {100*post["fraction_failing_any_separable_test"]:.3f}% |')
        else:
            lines.append(f'| {e["label"]} | {e["status"]} | no completed fit | no completed fit | — | — |')
    lines += ['', 'The invalid-prefix-solver archive and diagnostic smoke test are inventoried but excluded from science conclusions. '
              'Interrupted-run checkpoints are screened separately in the JSON; live points are not posterior samples. '
              'Active runs remain pending until they produce completed posteriors.', '',
              'Legacy physical reconstructions are checked against every archived mass-profile draw. '
              'Several legacy likelihoods do not reproduce under current code/data; this audit does not validate their likelihoods or evidence comparisons. '
              'Legacy finite-radius bounds use current loaders; modern runs use frozen data snapshots.', '',
              '## Density-profile optimizations and repeat starts', '',
              '| Optimization | rho_DM(20 pc) | Central condition | Separable screen on observed radial scales |',
              '|---|---:|---|---|']
    for e in result['optimizations']:
        p=e['point']
        lines.append(f'| {e["label"]} | {e["rho20"]:.4g} | {"FAIL" if p["central_condition_fails"] else "passes necessary condition"} | '
                     f'{"FAIL" if p["resolved_scale_screen"]["separable_df_ruled_out"] else "no failure found"} |')
    (OUT/'df_consistency_inventory.md').write_text('\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grid',type=int,default=2048)
    parser.add_argument('--posterior-grid',type=int,default=512)
    args=parser.parse_args()
    if min(args.grid,args.posterior_grid)<64:
        parser.error('grids must have at least 64 points')
    warnings.filterwarnings('ignore',message='Legacy run has no data snapshot')
    OUT.mkdir(parents=True,exist_ok=True)
    paths=sorted((ROOT/'results/fits').rglob('run.yaml'))
    paths += sorted((ROOT/'results/diagnostics').glob('density_conditioning_smoke_*/run.yaml'))
    entries, curves = [], {}
    for path in paths:
        if 'code' in path.relative_to(ROOT).parts:
            continue
        print('Auditing',path.parent.name,flush=True)
        entry,p=scan_run(path,args.grid,args.posterior_grid)
        entries.append(entry)
        if p is not None:
            curves[entry['label']]=p
        (OUT/'df_consistency_audit.partial.json').write_text(json.dumps(entries,indent=2,allow_nan=False)+'\n')
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),
                scope='stellar tracer in fitted spherical potential; necessary conditions, not a unique DF inversion',
                assumptions=dict(stellar_density='positive MGE; central logarithmic slope zero',
                    potential='saved unsoftened point mass plus positive extended components',
                    finite_radius='separable augmented density P(Psi) R(r^2)',
                    spatial_scope='spherical radii anchored to projected radial scales; central limit treated separately',
                    dm_remnant_df='unspecified by these Jeans fits; not certified'),
                grid=args.grid,posterior_grid=args.posterior_grid,tolerance=1e-8,
                source_sha256={str(p.relative_to(ROOT)):sha256(p) for p in
                               [Path(__file__).resolve(),ROOT/'src/ocen_dm/kinematics/df_consistency.py']},
                runs=entries,optimizations=scan_optimizations(args.grid))
    (OUT/'df_consistency_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_table(result)
    make_plot(entries,curves)
    print('Wrote audit, inventory and plots/df_consistency_audit.png',flush=True)


if __name__=='__main__':
    main()
