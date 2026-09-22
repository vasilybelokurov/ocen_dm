#!/usr/bin/env python3
"""Fit independent mass/light DFs to one preserved mock realization.

Fresh output directories are required. Every stage and rejected trial is saved.
This is optimization and a mock recovery experiment, not a posterior sampler.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys
import time

for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name,"1")
os.environ.setdefault("MPLCONFIGDIR","/private/tmp/ocen-df-mpl")
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))

import numpy as np
from scipy.optimize import least_squares
from scipy.special import softmax

from ocen_dm.kinematics.df_capacity import CapacityObservations,FixedPotentialProjector
from ocen_dm.kinematics.df_fit import PhotometricData
from ocen_dm.kinematics.df_mass_recovery import (
    DensityTemplate,MassLightConfig,MassLightDFModel,TemplateMassConfig,RecoveryProblem,gaussian_mock_noise,
)
from ocen_dm.kinematics.df_mock import challenge_mock
from ocen_dm.kinematics.positive_df import ActionComponent,DFModelConfig,DFNumerics,DFConvergenceError,agama_pc
from ocen_dm.kinematics.run_io import read_data_snapshot,write_data_snapshot,sha256


def write_json(path,row):
    # Atomic replacement makes progress files readable while optimization runs.
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(row,indent=2,allow_nan=False)+"\n")
    temporary.replace(path)


def initial_config(capacity_run,truth,start):
    summary=json.loads((capacity_run/"summary.json").read_text())
    if summary['case'] != ('mixed_with_dm' if len(truth.populations)==4 else 'mixed_no_dm'):
        raise ValueError("capacity seed belongs to the wrong mock")
    entry=next(e for e in summary['ladder'] if e['components']==3)
    x=np.array(entry['x'])
    kernel=FixedPotentialProjector(truth.agama_potential(),np.array([1.,20.]))
    weights=softmax(np.r_[0.,x[18:]])
    light=np.array([w*kernel.native_df(p,1).totalMass() for w,p in zip(weights,x[:18].reshape(3,6))])
    light/=sum(light)
    mass_weights=softmax(np.log(light)+np.array([.25,-.1,-.1])*(start-1))
    components=tuple(ActionComponent(fraction=w,J0=np.exp(p[0]),slope_in=p[1],slope_out=p[2],
                      steepness=np.exp(p[3]),h_r=p[4],g_r=p[5])
                      for w,p in zip(mass_weights,x[:18].reshape(3,6)))
    m=DFModelConfig(M_star=[3e6,2.4e6,3.6e6][start%3],M_bh=0.,M_rem=0.,M_dm_100=0.,components=components,
                   numerics=DFNumerics(potential_nodes=160,velocity_nodes=32,moment_nodes=180,
                    projection_nodes=80,iteration_tolerance=1e-5,max_iterations=80))
    dark=TemplateMassConfig(remnant_mass=[5e4,1e5,2e4][start%3],halo_mass=[3e5,8e5,5e4][start%3])
    return MassLightConfig(m,tuple(light),dark)


def summarize_model(result,templates):
    model=result['model'];config=result['config']
    radii=np.geomspace(.05,100.,100)
    stellar=-model.stellar_potential.force(np.column_stack((radii,radii*0,radii*0)))[:,0]*radii**2/agama_pc().G
    halo_rho20=float(templates['halo'].density(np.array([20.]),config.dark.halo_mass,config.dark.halo_scale)[0])
    return dict(score=result['score'],config=config.to_dict(),diagnostics=model.diagnostics,
                r_pc=radii.tolist(),total_mass_profile=model.mass_profile(radii).tolist(),
                stellar_mass_profile=stellar.tolist(),rho_dm_20=halo_rho20,
                total_mass_20=float(model.mass_profile(np.array([20.]))[0]),
                projected=result['projected'].tolist(),prediction=result['prediction'].tolist(),
                residual_sigma=result['residual'].tolist(),
                max_residual_sigma=float(np.max(abs(result['residual']))),
                rms_residual_sigma=float(np.sqrt(np.mean(result['residual']**2))))


class ScaledObjective:
    """Explicit finite steps resolve the accuracy of the equilibrium solve.

Accepted points use cold, reproducible equilibria. Derivative probes share the
same converged starting potential, with the same strict convergence/closure
guards. No invalid DF or moment is clipped into a valid likelihood.
    """
    def __init__(self,problem,out,stage,step=5e-4):
        self.problem,self.out,self.stage,self.step=problem,out,stage,step
        self.span=problem.upper-problem.lower
        self.cache_y=None
        self.cache_result=None
        self.calls=0
        self.accepted_evaluations=0
        self.rejections=0
        self.best=None
        self.started=time.monotonic()

    def evaluate(self,y,probe=False,seed=None):
        if not probe and self.cache_y is not None and np.array_equal(y,self.cache_y):
            return self.cache_result
        self.calls+=1
        x=self.problem.lower+self.span*y
        row=dict(call=self.calls,probe=probe,coordinates=x.tolist(),elapsed_seconds=time.monotonic()-self.started)
        try:
            self.problem.last_model=seed
            result=self.problem.evaluate(x,warm=probe and seed is not None)
            row.update(score=result['score'],iterations=result['model'].diagnostics['iterations'])
            if not probe:
                self.accepted_evaluations+=1
                if self.best is None or result['score']<self.best['score']:
                    self.best=result
                    write_json(self.out/'best_model.json',result['config'].to_dict())
                    write_json(self.out/'best.json',summarize_model(result,self.problem.templates))
            answer=result
        except (DFConvergenceError,FloatingPointError,ValueError) as exc:
            self.rejections+=1
            row['rejected']=str(exc)
            answer=None
        with (self.out/f'{self.stage}_evaluations.jsonl').open('a') as stream:
            stream.write(json.dumps(row,allow_nan=False)+'\n')
        if not probe:
            self.cache_y=np.array(y).copy()
            self.cache_result=answer
        if self.calls%25==0:
            print(f'{self.stage}: {self.calls} model calls, {self.rejections} rejected; best Q='+
                  ('none' if self.best is None else f"{self.best['score']:.5g}"),flush=True)
        return answer

    def fun(self,y):
        result=self.evaluate(y)
        return np.full(len(self.problem.target),1e6) if result is None else result['residual']

    def jac(self,y):
        base=self.evaluate(y)
        if base is None:
            raise DFConvergenceError('optimizer requested derivatives at invalid equilibrium')
        output=np.empty((len(self.problem.target),len(y)))
        for i in range(len(y)):
            direction=1 if y[i]+self.step<=1 else -1
            trial=y.copy();trial[i]+=direction*self.step
            changed=self.evaluate(trial,True,base['model'])
            if changed is None:
                direction=-direction
                trial=y.copy();trial[i]+=direction*self.step
                if trial[i]<0 or trial[i]>1:
                    raise DFConvergenceError(f'no valid derivative probe for {self.problem.names[i]}')
                changed=self.evaluate(trial,True,base['model'])
            if changed is None:
                raise DFConvergenceError(f'both derivative probes failed for {self.problem.names[i]}')
            output[:,i]=(changed['residual']-base['residual'])/(direction*self.step)
        return output


def validate(problem,result):
    config=result['config']
    n=config.mass.numerics
    refined_cfg=replace(config,mass=replace(config.mass,numerics=replace(n,
        potential_nodes=2*n.potential_nodes,velocity_nodes=2*n.velocity_nodes,
        moment_nodes=2*n.moment_nodes,projection_nodes=2*n.projection_nodes,
        iteration_tolerance=n.iteration_tolerance/2)))
    refined_problem=RecoveryProblem(problem.observations,refined_cfg,shapes='full',
         templates=problem.templates,noise=problem.noise,shared_ml=problem.shared_ml)
    refined=refined_problem.evaluate(refined_problem.encode())
    shift=refined['residual']-result['residual']
    r=np.geomspace(.05,68.,14)
    fast,direct=refined['model'].projected_moments(r),refined['model'].direct_projected_moments(r)
    relative={k:float(np.max(abs(fast[k]/direct[k]-1))) for k in fast}
    cold=problem.evaluate(problem.encode(config))
    replay=float(np.max(abs(cold['residual']-result['residual'])))
    return dict(passed=bool(max(abs(shift))<.05 and max(relative.values())<.005 and replay<1e-8),
                refinement_max_shift_sigma=float(np.max(abs(shift))),native_projection_max_relative=relative,
                cold_replay_max_shift_sigma=replay,refined_score=refined['score'],
                refined_rms_residual_sigma=float(np.sqrt(np.mean(refined['residual']**2))),
                refined_max_residual_sigma=float(np.max(abs(refined['residual']))),
                refined_diagnostics=refined['model'].diagnostics)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=['mixed_no_dm','mixed_with_dm'],required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--data-run',type=Path,required=True)
    p.add_argument('--capacity-run',type=Path,required=True)
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--noise-seed',type=int)
    p.add_argument('--initial-model',type=Path)
    p.add_argument('--shared-ml',action='store_true')
    p.add_argument('--stages',nargs='+',choices=['weights','scales','full'],default=['weights','scales','full'])
    p.add_argument('--max-evaluations',type=int,default=60)
    p.add_argument('--evaluate-only',action='store_true')
    args=p.parse_args()
    if args.max_evaluations<1 or args.start<0: p.error('invalid fitting budget/start')
    args.out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    manifest=dict(status='running',case=args.case,started_utc=datetime.now(timezone.utc).isoformat(),
        start=args.start,noise_seed=args.noise_seed,shared_ml=args.shared_ml,stages=args.stages,
        max_evaluations_per_stage=args.max_evaluations,
        scope='Self-consistent positive stellar mass/light DFs; matched dark density families with free masses and radial scales; optimization, no posterior',
        dark_positivity='spatial density certified; dark-component DFs not certified in fitted potentials',
        initialization='action shapes from known-potential capacity test; different stellar/dark mass and M/L starts')
    write_json(args.out/'run.json',manifest)
    try:
        sources=list((ROOT/'src/ocen_dm/kinematics').glob('df_*.py'))+[ROOT/'src/ocen_dm/kinematics/positive_df.py',Path(__file__)]
        code=args.out/'code';code.mkdir()
        manifest['source_sha256']={str(f.relative_to(ROOT)):sha256(f) for f in sources}
        for f in sources: shutil.copy2(f,code/f.name)
        agama=agama_pc()
        manifest['agama']=dict(version=agama.__version__,binary_sha256={f.name:sha256(f) for f in Path(agama.__file__).parent.glob('*.so')})
        truth=challenge_mock(args.case)
        halo_truth=truth if args.case=='mixed_with_dm' else challenge_mock('mixed_with_dm')
        templates=dict(remnant=DensityTemplate.from_mock(truth,'dark_remnants'),halo=DensityTemplate.from_mock(halo_truth,'extended_dark_component'))
        write_json(args.out/'templates.json',{k:v.to_dict() for k,v in templates.items()})
        write_json(args.out/'truth.json',truth.to_dict())
        meta=json.loads((args.data_run/'summary.json').read_text())
        data=read_data_snapshot(args.data_run,meta['data_snapshot'])
        manifest['data_snapshot']=write_data_snapshot(data,args.out)
        photo=PhotometricData.from_dict(json.loads((args.data_run/'photometry.json').read_text()))
        write_json(args.out/'photometry.json',photo.to_dict())
        obs=CapacityObservations(data,photo,truth)
        noise=np.zeros(len(obs.truth)) if args.noise_seed is None else gaussian_mock_noise(obs,args.noise_seed)
        write_json(args.out/'mock_observations.json',dict(**obs.to_dict(),noise_values=noise.tolist(),target_values=(obs.truth+noise).tolist()))
        if args.initial_model:
            cfg=MassLightConfig.from_dict(json.loads(args.initial_model.read_text()))
            manifest['initial_model']=dict(path=str(args.initial_model.resolve()),sha256=sha256(args.initial_model))
            if args.start:
                fractions=softmax(np.log([c.fraction for c in cfg.mass.components])+
                                  .08*np.cos(np.arange(len(cfg.mass.components))+args.start))
                cfg=replace(cfg,mass=replace(cfg.mass,M_star=np.clip(cfg.mass.M_star*.96,5.01e5,6.99e6),
                    components=tuple(replace(c,fraction=float(w)) for c,w in zip(cfg.mass.components,fractions))),
                    dark=replace(cfg.dark,halo_mass=min(2.99e6,cfg.dark.halo_mass*1.3+3e4),
                                           remnant_mass=min(7.99e5,cfg.dark.remnant_mass*1.15)))
        else:
            cfg=initial_config(args.capacity_run,truth,args.start)
        write_json(args.out/'initial_model.json',cfg.to_dict())
        manifest['capacity_seed']=dict(path=str(args.capacity_run.resolve()),summary_sha256=sha256(args.capacity_run/'summary.json'))
        write_json(args.out/'run.json',manifest)
        stages=[]
        for stage in args.stages:
            problem=RecoveryProblem(obs,cfg,templates=templates,shapes=stage,noise=noise,shared_ml=args.shared_ml)
            objective=ScaledObjective(problem,args.out,stage)
            y=(problem.encode()-problem.lower)/(problem.upper-problem.lower)
            y=np.clip(y,1e-8,1-1e-8)
            initial=objective.evaluate(y)
            if initial is None: raise DFConvergenceError('initial equilibrium fails accuracy guards')
            print(f'{args.case}/{stage}: initial Q={initial["score"]:.5g}, {len(y)} free parameters',flush=True)
            if args.evaluate_only:
                optimizer=None
            else:
                optimized=least_squares(objective.fun,y,jac=objective.jac,bounds=(np.zeros(len(y)),np.ones(len(y))),
                    max_nfev=args.max_evaluations,x_scale='jac',ftol=2e-5,xtol=2e-5,gtol=1e-5)
                optimizer=dict(success=bool(optimized.success),message=str(optimized.message),nfev=int(optimized.nfev),
                               njev=int(optimized.njev),optimality=float(optimized.optimality))
            best=objective.best
            cfg=best['config']
            record=dict(stage=stage,initial_score=initial['score'],best=summarize_model(best,templates),
                         optimizer=optimizer,model_calls=objective.calls,rejected_trials=objective.rejections,
                         parameter_names=problem.names,bounds=[problem.lower.tolist(),problem.upper.tolist()])
            stages.append(record)
            write_json(args.out/f'{stage}_result.json',record)
            write_json(args.out/'summary.json',dict(status='running',case=args.case,stages=stages))
            print(f'{args.case}/{stage}: Q={best["score"]:.5g}; optimizer={optimizer}',flush=True)
            if args.evaluate_only: break
        validation=validate(problem,best)
        record=dict(status='completed',case=args.case,noise_seed=args.noise_seed,shared_ml=args.shared_ml,stages=stages,
                    best=summarize_model(best,templates),validation=validation,elapsed_seconds=time.monotonic()-started)
        write_json(args.out/'summary.json',record)
        manifest.update(status='completed',validation_passed=validation['passed'],elapsed_seconds=record['elapsed_seconds'])
        write_json(args.out/'run.json',manifest)
        print(json.dumps(dict(score=best['score'],validation=validation),indent=2),flush=True)
    except BaseException as exc:
        manifest.update(status='failed',error=repr(exc),elapsed_seconds=time.monotonic()-started)
        write_json(args.out/'run.json',manifest)
        raise


if __name__=='__main__':
    main()
