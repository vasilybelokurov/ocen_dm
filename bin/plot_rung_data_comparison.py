#!/usr/bin/env python3
"""Observed kinematics and the two lowest-chi-square ML models in each rung."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullFormatter, StrMethodFormatter
import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import data_fingerprint


def main():
    root = Path(__file__).resolve().parents[1]
    prefix = root/'plots/rung0_rung1_data_models'
    prefix.parent.mkdir(exist_ok=True)
    data_prefix = root/'results/plot_data/rung0_rung1_data_models'
    data_prefix.parent.mkdir(parents=True, exist_ok=True)
    groups = [['rung0_K1', 'rung0_K2_nfw', 'rung0_K2_cored'],
              ['rung1_K1_slice', 'rung1_K2_nfw', 'rung1_K2_cored_slice']]
    selected, provenance, hashes = [], [], {}
    fingerprint = None
    profiles = None
    for rung, group in enumerate(groups):
        runs = [load_run_metadata(root/'results/fits'/name) for name in group]
        runs.sort(key=lambda run: run['summary']['chi2_ml_total'])
        for run in runs[:2]:
            problem = problem_for(run)  # Enforces replay of the saved maximum likelihood.
            identity = data_fingerprint(problem.data)
            if fingerprint is not None and identity != fingerprint:
                raise ValueError('Selected runs do not use identical observations')
            fingerprint = identity
            profiles = problem.data.profiles
            summary = run['summary']
            x = np.array([summary['parameters'][name]['ml'] for name in problem.family.names])
            prediction = problem.predict(x)
            chi2 = problem.chi2(x)
            total = sum(value[0] for value in chi2.values())
            if not np.isclose(total, summary['chi2_ml_total'], rtol=0, atol=.011):
                raise ValueError('Replayed chi-square differs from the stored total')
            halo = 'NFW' if 'nfw' in run['label'] else 'Core'
            selected.append(dict(label=run['label'], rung=rung, halo=halo,
                color='#0072B2' if halo == 'NFW' else '#D55E00',
                linestyle=(0, (5, 3)) if rung == 0 else '-', prediction=prediction, chi2=total))
            provenance.append(dict(label=run['label'], rung=rung, halo=halo, chi2=total,
                lnL_saved=summary['lnL_max'], data_fingerprint=identity,
                parameters_ml=dict(zip(problem.family.names, x.tolist())),
                dataset_chi2={name: value[0] for name, value in chi2.items()},
                data_replay='saved snapshot' if 'data_snapshot' in summary else 'legacy inputs; verified likelihood replay'))
            for name in ('summary.json', 'run.yaml', 'data_snapshot.json'):
                path = run['dir']/name
                if path.exists():
                    hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    selected.sort(key=lambda row: (row['rung'], row['halo'] != 'NFW'))
    names = ['hst_pm_radial_ours', 'hst_pm_tangential_ours', 'muse_los_dispersion',
             'gaia_edr3_ours_radial', 'gaia_edr3_ours_tangential']
    titles = ['HST: radial proper motion', 'HST: tangential proper motion', 'MUSE: line of sight',
              'Gaia: radial proper motion', 'Gaia: tangential proper motion']
    by_name = {p.name: p for p in profiles}
    plt.rcParams.update({'font.size': 10.5, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
    fig = plt.figure(figsize=(14, 9.2))
    grid = fig.add_gridspec(2, 3, left=.066, right=.985, top=.89, bottom=.075,
                          hspace=.39, wspace=.29)
    exported = []
    for index, (name, title) in enumerate(zip(names, titles)):
        profile = by_name[name]
        sub = grid[index//3, index%3].subgridspec(2, 1, height_ratios=[3, 1.15], hspace=.04)
        ax = fig.add_subplot(sub[0])
        residual_ax = fig.add_subplot(sub[1], sharex=ax)
        lo = np.maximum(profile.r_lower, profile.r.min()*.75)
        ax.errorbar(profile.r, profile.value, yerr=[profile.err_lo, profile.err_hi],
                    xerr=[profile.r-lo, profile.r_upper-profile.r], fmt='o', ms=3.8,
                    mfc='white', mec='#222222', mew=.85, ecolor='#989898',
                    elinewidth=.8, capsize=0, zorder=5)
        max_residual = 0.
        for model in selected:
            prediction = model['prediction'][name]
            # These are the likelihood's measured-star/bin averages with the
            # streaming correction; joining them does not add a smooth prediction.
            style = dict(color=model['color'], linestyle=model['linestyle'], lw=1.9,
                         marker='o' if model['halo'] == 'NFW' else 's', ms=2.5,
                         markerfacecolor='white', markeredgewidth=.7)
            ax.plot(profile.r, prediction, **style)
            error = np.where(prediction > profile.value, profile.err_hi, profile.err_lo)
            residual = (profile.value-prediction)/error
            residual_ax.plot(profile.r, residual, **style)
            max_residual = max(max_residual, float(np.max(np.abs(residual))))
        for i in range(profile.n):
            exported.append(dict(dataset=name, radius_arcsec=float(profile.r[i]),
                bin_lower_arcsec=float(profile.r_lower[i]), bin_upper_arcsec=float(profile.r_upper[i]),
                value=float(profile.value[i]), err_lo=float(profile.err_lo[i]), err_hi=float(profile.err_hi[i]),
                **{model['label']: float(model['prediction'][name][i]) for model in selected}))
        ax.set_title(f'{title}  ({profile.n} bins)', fontsize=11.5, pad=9)
        ax.set_ylabel(r'$\sigma_{\rm LOS}$ [km s$^{-1}$]' if profile.kind == 'los' else
                      (r'$\sigma_{\mu,R}$ [mas yr$^{-1}$]' if profile.kind == 'pmr' else
                       r'$\sigma_{\mu,T}$ [mas yr$^{-1}$]'))
        ax.set_xscale('log')
        ax.tick_params(labelbottom=False)
        residual_ax.axhspan(-1, 1, color='#999999', alpha=.14, lw=0, zorder=0)
        residual_ax.axhline(0, color='#777777', lw=.7, zorder=0)
        limit = max(3, np.ceil(max_residual*1.12))
        residual_ax.set_ylim(-limit, limit)
        residual_ax.set_yticks([-int(limit//2), 0, int(limit//2)])
        residual_ax.set_ylabel(r'$\Delta/\epsilon$')
        residual_ax.set_xlabel('Projected radius [arcsec]')
        for a in (ax, residual_ax):
            a.grid(axis='y', alpha=.15)
            a.tick_params(which='both', direction='out', labelsize=9.5)
        ax.set_xlim(max(.1, profile.r.min()*.76), profile.r_upper.max()*1.05)
        if name.startswith('gaia_'):
            residual_ax.xaxis.set_major_locator(FixedLocator([400, 600, 1000, 2000]))
            residual_ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
            residual_ax.xaxis.set_minor_formatter(NullFormatter())

    legend_ax = fig.add_subplot(grid[1, 2])
    legend_ax.set_axis_off()
    handles = [Line2D([], [], color=m['color'], ls=m['linestyle'], lw=2.2,
                      marker='o' if m['halo'] == 'NFW' else 's', ms=4, mfc='white',
                      label=f"Rung {m['rung']}  {m['halo']}     $\\chi^2 = {m['chi2']:.2f}$") for m in selected]
    handles.append(Line2D([], [], color='#222222', ls='', marker='o', mfc='white', ms=5,
                          label='Observed data with 1-sigma errors'))
    legend_ax.legend(handles=handles, loc='upper left', frameon=False, fontsize=11,
                     handlelength=3.5, labelspacing=1.2, borderaxespad=0)
    legend_ax.text(0, .33, 'Dashed: rung 0, isotropic\nSolid: rung 1, constant anisotropy\n\n'
        'Points on model lines are the predictions\nevaluated in the likelihood; lines join bins.\n\n'
        r'Residuals: $\Delta/\epsilon =$ (data - model)/error'+'\nShading marks the +/-1-sigma interval.',
        transform=legend_ax.transAxes, va='top', fontsize=10.3, linespacing=1.5)
    fig.suptitle('Data and the two best-fitting models from each rung', fontsize=18, y=.973)
    fig.text(.5, .929, 'Maximum-likelihood samples selected by lowest total chi-square; same 89 measurements in all four fits',
             ha='center', fontsize=11, color='#444444')
    fig.savefig(prefix.with_suffix('.png'), dpi=180)
    plt.close(fig)
    table = Table(rows=exported)
    table.meta.update(data_fingerprint=fingerprint, prediction='exact likelihood bin prediction including streaming correction',
                      parameter_choice='maximum likelihood; two lowest chi-square families in each rung',
                      units='radius: arcsec; HST/Gaia value and predictions: mas/yr; MUSE: km/s')
    table.write(data_prefix.with_suffix('.ecsv'), overwrite=True)
    record = dict(selection='two smallest chi2_ml_total among K1, NFW and cored in each rung',
                  models=provenance, inputs_sha256=hashes,
                  generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    record['outputs_sha256'] = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in (prefix.with_suffix('.png'), data_prefix.with_suffix('.ecsv'))}
    data_prefix.with_suffix('.json').write_text(json.dumps(record, indent=2)+'\n')
    print(prefix.with_suffix('.png'))


if __name__ == '__main__':
    main()
