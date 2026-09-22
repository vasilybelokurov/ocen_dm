#!/usr/bin/env python3
"""Report a preserved free-potential DF recovery batch, including active progress."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
import hashlib
import os
from pathlib import Path
import sys
import time

os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ocen-df-mpl')
os.environ.setdefault('OMP_NUM_THREADS','1')
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ocen_dm.kinematics.df_mock import LoweredIsothermalMock,MockPopulation
from ocen_dm.kinematics.run_io import sha256


def report(batch,prefix):
    meta=json.loads((batch/'batch.json').read_text())
    jobs=json.loads((batch/'jobs.json').read_text())
    cases=['mixed_no_dm','mixed_with_dm']
    records=[]
    truths={}
    for job in jobs:
        directory=Path(job['directory'])
        if not (directory/'best.json').exists():continue
        best_bytes=(directory/'best.json').read_bytes()
        best=json.loads(best_bytes)
        truth_row=json.loads((directory/'truth.json').read_text())
        if job['case'] not in truths:
            m=LoweredIsothermalMock([MockPopulation(**p) for p in truth_row['populations']],W0=truth_row['W0'],
                                    radius_pc=truth_row['radius_pc'],velocity_kms=truth_row['velocity_kms'])
            r=np.array(best['r_pc'])
            dm=[i for i,p in enumerate(m.populations) if p.name=='extended_dark_component']
            truths[job['case']]=dict(r_pc=r.tolist(),total_mass=m.enclosed_mass(r).tolist(),
                                    total_mass_20=float(m.enclosed_mass(np.array([20.]))[0]),
                                    stellar_mass=float(sum(v for v,p in zip(m.masses,m.populations) if p.light_per_mass>0)),
                                    halo_mass=0. if not dm else float(m.masses[dm[0]]),
                                    rho_dm_20=0. if not dm else float(m.intrinsic_components(np.array([20.]))[0,dm[0],0]),
                                    physical_specification=truth_row)
        stage='initial'
        history=[];offset=0
        for name in ['weights','scales','full']:
            file=directory/f'{name}_evaluations.jsonl'
            if not file.exists():continue
            stage=name
            lines=file.read_text().splitlines();last=0
            for i,line in enumerate(lines):
                try:row=json.loads(line)
                except json.JSONDecodeError:
                    if i==len(lines)-1:break
                    raise
                last=max(last,row['elapsed_seconds'])
                if not row['probe'] and 'score' in row:
                    history.append([offset+row['elapsed_seconds'],row['score']])
            offset+=last
        summary=json.loads((directory/'summary.json').read_text()) if (directory/'summary.json').exists() else None
        qualified=bool(summary and summary['status']=='completed' and summary['validation']['passed'] and
                       summary['stages'][-1]['stage']=='full' and summary['stages'][-1]['optimizer']['success'])
        records.append(dict(job=job,best=best,stage=stage,qualified=qualified,history=history,
                            summary=summary,source_sha256={'best.json':hashlib.sha256(best_bytes).hexdigest()}))
    if not records:return False
    now=datetime.now(timezone.utc).isoformat()
    plotdir=ROOT/'plots';datadir=ROOT/'results/plot_data'
    plotdir.mkdir(exist_ok=True);datadir.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(2,2,figsize=(12,8),height_ratios=[2,1],constrained_layout=True)
    for col,case in enumerate(cases):
        if case not in truths:continue
        t=truths[case];r=np.array(t['r_pc']);mass=np.array(t['total_mass'])
        axes[0,col].plot(r,mass,color='black',lw=2,label='Injected true mass')
        for record in records:
            j,b=record['job'],record['best']
            if j['case']!=case or j['noise_seed'] is not None:continue
            label=f"Start {j['start']}: {record['stage']}"+(' (shared M/L)' if j['shared_ml'] else '')
            color=f"C{j['start']}"
            ls=':' if j['shared_ml'] else '-' if record['qualified'] else '--'
            m=np.array(b['total_mass_profile'])
            axes[0,col].plot(r,m,color=color,ls=ls,label=label)
            axes[1,col].plot(r,m/mass-1,color=color,ls=ls)
        axes[0,col].set(xscale='log',yscale='log',ylabel='Enclosed total mass [M☉]',title=case.replace('_',' '))
        axes[1,col].set(xscale='log',xlabel='Radius [pc]',ylabel='Recovered / true − 1')
        axes[1,col].axhline(0,color='black',lw=.7)
        axes[1,col].axvline(20,color='grey',ls=':',lw=.7)
        axes[0,col].legend(fontsize=8)
    fig.suptitle(f"Free-potential DF recovery — {meta['status'].replace('_',' ')}\nDashed curves are interim fits; stellar mass/light weights are independent")
    fig.savefig(plotdir/f'{prefix}_mass_profiles.png',dpi=170);plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    for col,case in enumerate(cases):
        for record in records:
            j=record['job']
            if j['case']!=case or j['noise_seed'] is not None:continue
            h=np.array(record['history'])
            if not len(h):continue
            axes[col].plot(h[:,0]/60,np.minimum.accumulate(h[:,1]),color=f"C{j['start']}",
                           ls=':' if j['shared_ml'] else '-',label=f"Start {j['start']}"+(' shared M/L' if j['shared_ml'] else ''))
        axes[col].set(yscale='log',xlabel='Elapsed fitting time [min]',ylabel='Best joint approximation score Q',title=case.replace('_',' '))
        axes[col].legend(fontsize=8)
    fig.suptitle('Noiseless recovery optimization — lower Q alone does not establish unbiased component masses')
    fig.savefig(plotdir/f'{prefix}_progress.png',dpi=170);plt.close(fig)

    noisy_best={}
    for record in records:
        j=record['job']
        if j['noise_seed'] is None or not record['qualified']:continue
        key=(j['case'],j['noise_seed'])
        if key not in noisy_best or record['best']['score']<noisy_best[key]['best']['score']:
            noisy_best[key]=record
    noise_statistics={}
    if noisy_best:
        fig,axes=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
        for row,case in enumerate(cases):
            selected=sorted([r for (c,seed),r in noisy_best.items() if c==case],key=lambda r:r['job']['noise_seed'])
            if case not in truths:continue
            t=truths[case]
            arrays=[np.array([r['best']['total_mass_20']/t['total_mass_20']-1 for r in selected]),
                    np.array([r['best']['config']['mass']['M_star']/t['stellar_mass']-1 for r in selected]),
                    np.array([r['best']['config']['dark']['halo_mass'] for r in selected]),
                    np.array([r['best']['rho_dm_20'] for r in selected])]
            targets=[0.,0.,t['halo_mass'],t['rho_dm_20']]
            labels=['M(<20 pc) / true − 1','Stellar mass / true − 1','Halo mass [M☉]','Halo density at 20 pc [M☉ pc⁻³]']
            noise_statistics[case]=dict(realizations=len(selected),seeds=[r['job']['noise_seed'] for r in selected],
                columns={label:dict(values=a.tolist(),injected=target,median=None if not len(a) else float(np.median(a)),
                         standard_deviation=None if len(a)<2 else float(np.std(a,ddof=1)))
                         for label,a,target in zip(labels,arrays,targets)})
            for ax,a,target,label in zip(axes[row],arrays,targets,labels):
                ax.axhline(target,color='black',ls=':',label='Injected truth')
                if len(a):ax.plot(np.arange(1,len(a)+1),a,'o',color=f'C{row}')
                else:ax.text(.5,.5,'Awaiting qualified noisy fits',transform=ax.transAxes,ha='center',fontsize=8)
                ax.set(xlabel='Noise realization',ylabel=label,title=case.replace('_',' '))
                ax.legend(fontsize=8)
        fig.suptitle('Pilot noise recovery — best qualified result per realization; scatter is not a posterior interval')
        fig.savefig(plotdir/f'{prefix}_noise_recovery.png',dpi=170);plt.close(fig)

    lines=['# Free-potential DF recovery: current status','',f'Snapshot: {now}. Batch status: **{meta["status"]}**.','',
           'The table reports the latest saved solution in each noiseless fit. Interim and unconverged solutions are not mass measurements.','',
           '| Mock | Start | Stage | Job | Q | Stellar mass | Halo mass | DM density at 20 pc | Qualified final fit |',
           '|---|---:|---|---|---:|---:|---:|---:|---|']
    for record in records:
        j,b=record['job'],record['best']
        if j['noise_seed'] is not None:continue
        lines.append(f"| {j['case']}"+(' shared M/L' if j['shared_ml'] else '')+
                     f" | {j['start']} | {record['stage']} | {j['status']} | {b['score']:.4g} | {b['config']['mass']['M_star']:.0f} | {b['config']['dark']['halo_mass']:.0f} | {b['rho_dm_20']:.5g} | {record['qualified']} |")
    lines+=['','Masses are in solar masses; densities are in solar masses per cubic parsec. A qualified final fit has finished the full shape stage, met an optimizer termination criterion and passed numerical validation.','',
            f"Jobs: {sum(j['status']=='running' for j in jobs)} running, {sum(j['status']=='queued' for j in jobs)} queued, {sum(j['status']=='completed' for j in jobs)} completed, {sum(j['status']=='failed' for j in jobs)} failed.",'',
            'Noisy fits and shared-M/L controls are added only after the noiseless accuracy gate passes. The planned eight noise realizations per mock assess recovery scatter; they are not a posterior-coverage certification.','',
            f'![Current mass profiles](../plots/{prefix}_mass_profiles.png)','',f'![Optimization progress](../plots/{prefix}_progress.png)','',
            '[Methods and scope](DF_MASS_RECOVERY.md).','',f'Batch: `{batch}`.']
    if noisy_best:
        lines+=['','## Noisy recovery pilot','',
                'Only completed, numerically validated full-stage fits with optimizer termination enter this summary. The lower-Q qualified start is selected per realization. This small ensemble measures recovery scatter, not posterior coverage.','']
        for case,statistics in noise_statistics.items():
            lines.append(f"- {case}: {statistics['realizations']} qualified noise realizations.")
        lines+=['',f'![Noisy recovery](../plots/{prefix}_noise_recovery.png)']
    statusfile=ROOT/'docs/DF_MASS_RECOVERY_STATUS.md'
    statusfile.write_text('\n'.join(lines)+'\n')
    payload=dict(snapshot_utc=now,batch=str(batch),batch_manifest=meta,jobs=jobs,truths=truths,records=records,
                 noise_recovery_statistics=noise_statistics,
                 generator_sha256=sha256(Path(__file__)))
    (datadir/f'{prefix}.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print(f"Reported {len(records)} saved fits; batch {meta['status']}",flush=True)
    return meta['status']=='running'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',type=Path,required=True)
    parser.add_argument('--prefix',default='df_mass_recovery_20260922')
    parser.add_argument('--watch',action='store_true')
    args=parser.parse_args()
    if Path(args.prefix).name!=args.prefix:parser.error('prefix must be a filename stem')
    while True:
        active=report(args.batch,args.prefix)
        if not args.watch or not active:break
        time.sleep(60)


if __name__=='__main__':main()
