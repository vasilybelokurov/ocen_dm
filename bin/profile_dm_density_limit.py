#!/usr/bin/env python3
"""Refit nuisance parameters at fixed DM density to explain the posterior limit.

The objective is kinematic chi-square plus the archived Gaussian distance
constraint. This profile is not a marginal posterior or a calibrated confidence
interval. The relaxed-stellar-floor case is a labelled prior-support diagnostic.
"""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from astropy.table import Table
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import data_fingerprint


ROOT = Path(__file__).resolve().parents[1]
TARGETS = [.1,.3,.6,1.03,1.5,2.,3.,5.]
CASES = [("cusp","rung2_K2_nfw_turnover",1e6),
         ("core","rung2_K2_cored_turnover",1e6),
         ("core_lower_stellar_floor","rung2_K2_cored_turnover",1e4)]
OUTPUT = ROOT/"results/diagnostics/dm_density_profile_20260922"
NODES,WEIGHTS = leggauss(64)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mass100_for_density(rho,rs,gamma,rt):
    r = 50*(NODES+1)
    x = r/rs
    integral = 4*np.pi*50*np.sum(WEIGHTS*r*r*x**(-gamma)*(1+x)**(gamma-3)*(1+(r/rt)**2)**-2)
    x20 = 20/rs
    shape20 = x20**(-gamma)*(1+x20)**(gamma-3)*(1+(20/rt)**2)**-2
    return float(rho*integral/shape20)


class ProfileProblem:
    def __init__(self,problem,stellar_floor):
        self.problem = problem
        self.family = problem.family
        self.names = self.family.names
        self.logarithmic = np.array([p.prior.kind=="loguniform" for p in self.family.parameters])
        self.lower,self.upper = [],[]
        self.priors = {p.name:p.prior for p in self.family.parameters}
        for p in self.family.parameters:
            q = p.prior
            if q.kind=="normal":
                lo,hi = q.mu-8*q.sigma,q.mu+8*q.sigma
            else:
                lo,hi = q.lo,q.hi
            if p.name=="M_star":
                lo = stellar_floor
            if p.name in ("beta_inf","beta_mid"):
                hi = min(hi,1-1e-8)
            if q.kind=="loguniform":
                lo,hi = np.log10(lo),np.log10(hi)
            self.lower.append(lo);self.upper.append(hi)
        self.lower,self.upper = np.array(self.lower),np.array(self.upper)
        self.dm_index = self.names.index("M_dm_100")
        self.keep = np.arange(len(self.names))!=self.dm_index
        self.calls = 0

    def encode(self,x,rho):
        v = np.array(x,float)
        v[self.logarithmic] = np.log10(v[self.logarithmic])
        v = np.clip(v,self.lower+1e-9,self.upper-1e-9)
        return v if rho is None else v[self.keep]

    def decode(self,v,rho):
        encoded = np.empty(len(self.names))
        if rho is None:
            encoded[:] = v
        else:
            encoded[self.keep] = v;encoded[self.dm_index] = 5
        x = encoded.copy();x[self.logarithmic] = 10**x[self.logarithmic]
        if rho is not None:
            x[self.dm_index] = mass100_for_density(rho,x[self.names.index("r_s")],self.family.gamma,self.family.r_t)
        return x

    def residual(self,v,rho):
        self.calls += 1
        x = self.decode(v,rho)
        theta = self.family.to_dict(x)
        m = self.priors["M_dm_100"]
        if not m.lo <= theta["M_dm_100"] <= m.hi:
            raise ValueError("Fixed density violates archived DM-mass support")
        prediction = self.problem.predict(x)
        residual = [(prediction[p.name]-p.value)/np.where(prediction[p.name]>p.value,p.err_hi,p.err_lo)
                    for p in self.problem.data.profiles]
        d = self.priors["distance"]
        residual.append(np.array([(theta["distance"]-d.mu)/d.sigma]))
        return np.concatenate(residual)

    def solve(self,starts,rho):
        fits = []
        bounds = (self.lower,self.upper) if rho is None else (self.lower[self.keep],self.upper[self.keep])
        for i,x in enumerate(starts):
            before = self.calls;start = time.monotonic()
            fit = least_squares(lambda v:self.residual(v,rho),self.encode(x,rho),bounds=bounds,
                                x_scale="jac",diff_step=2e-5,ftol=1e-9,xtol=1e-9,gtol=1e-6,max_nfev=300)
            xx = self.decode(fit.x,rho)
            residual = self.residual(fit.x,rho)
            fits.append(dict(x=xx,objective=float(residual@residual),success=bool(fit.success),
                             message=fit.message,calls=self.calls-before,nfev=fit.nfev,
                             optimality=float(fit.optimality),seconds=time.monotonic()-start))
        best = min(fits,key=lambda f:f["objective"])
        # Polish with a smaller differentiation step to detect discretisation sensitivity.
        fit = least_squares(lambda v:self.residual(v,rho),self.encode(best["x"],rho),bounds=bounds,
                            x_scale="jac",diff_step=2e-6,ftol=1e-10,xtol=1e-10,gtol=1e-7,max_nfev=200)
        residual = self.residual(fit.x,rho)
        if residual@residual < best["objective"]:
            xx = self.decode(fit.x,rho)
        else:
            xx = best["x"]
        theta = self.family.to_dict(xx)
        prediction = self.problem.predict(xx)
        chi = self.problem.chi2(xx)
        encoded = self.encode(xx,None)
        fraction = (encoded-self.lower)/(self.upper-self.lower)
        boundaries = {n:dict(side="lower" if f<.001 else "upper",fraction=float(f),value=theta[n])
                      for n,f in zip(self.names,fraction) if f<.001 or f>.999}
        density = float(np.asarray(self.family.halo(theta).density(20)).reshape(-1)[0])
        if rho is not None:
            np.testing.assert_allclose(density,rho,rtol=1e-6)
        prior = self.priors["distance"]
        dchi = ((theta["distance"]-prior.mu)/prior.sigma)**2
        return dict(rho=density,parameters=theta,chi2_kinematics=float(sum(c[0] for c in chi.values())),
                    chi2_distance=float(dchi),objective=float(sum(c[0] for c in chi.values())+dchi),
                    chi2_by_dataset={k:float(v[0]) for k,v in chi.items()},
                    prediction={k:np.asarray(v).tolist() for k,v in prediction.items()},active_boundaries=boundaries,
                    start_diagnostics=[{k:v for k,v in f.items() if k!="x"} for f in fits],
                    polish=dict(success=bool(fit.success),objective=float(residual@residual),nfev=fit.nfev,
                                optimality=float(fit.optimality)),calls=self.calls)


def run_case(case):
    name,label,floor = case
    path = OUTPUT/f"{name}.json"
    if path.exists():
        raise FileExistsError(f"Preserve existing diagnostic: {path}")
    run = load_run_metadata(ROOT/"results/fits"/label)
    problem = problem_for(run);profile = ProfileProblem(problem,floor)
    ml = np.array([run["summary"]["parameters"][n]["ml"] for n in problem.family.names])
    posterior = Table.read(run["dir"]/"posterior.ecsv")
    samples = np.column_stack([posterior[n] for n in problem.family.names])
    # Use the existing likelihood cache for informed, separated starting points.
    caches = sorted((ROOT/"results/diagnostics/dm_density_constraints").glob(label+"_*.npz"))
    cache = caches[-1]
    meta = json.loads(cache.with_suffix(".json").read_text())
    assert meta["cache_sha256"]==digest(cache)
    for f,h in meta["input_sha256"].items():
        assert digest(ROOT/f)==h
    with np.load(cache) as a:
        points,rhos,ll = a["theta"],a["rho"],a["loglike"].sum(axis=1)
    starts = [ml,points[np.argmax(ll)],samples[len(samples)//2]]
    print(name,"optimising free density",flush=True)
    free = profile.solve(starts,None)
    previous = np.array([free["parameters"][n] for n in problem.family.names])
    records = []
    bins = []
    for p in problem.data.profiles:
        bins.append(dict(dataset=p.name,radius_arcsec=p.r.tolist(),value=p.value.tolist(),
                         err_lo=p.err_lo.tolist(),err_hi=p.err_hi.tolist()))
    state = dict(case=name,label=label,stellar_floor=floor,model=run["summary"]["model"],
        data_fingerprint=data_fingerprint(problem.data),bins=bins,
        objective="kinematic chi-square plus ((distance-5.43)/0.05)^2; other priors impose support only",
        input_sha256=meta["input_sha256"],source_sha256=meta["source_sha256"],
        script_sha256=digest(Path(__file__)),created_utc=datetime.now(timezone.utc).isoformat(),
        status="running",free=free,profile=records)
    OUTPUT.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(state,indent=2)+"\n")
    print(name,'free',round(free['rho'],3),round(free['objective'],3),free['active_boundaries'],flush=True)
    for rho in TARGETS:
        distance = abs(np.log(rhos/rho))
        near = np.where(distance<=np.quantile(distance,.15))[0]
        best_near = points[near[np.argmax(ll[near])]]
        result = profile.solve([previous,ml,best_near],rho)
        previous = np.array([result["parameters"][n] for n in problem.family.names])
        records.append(result)
        path.write_text(json.dumps(state,indent=2)+"\n")
        print(name,'rho',rho,'objective',round(result['objective'],3),'chi2',round(result['chi2_kinematics'],3),
              'boundaries',result['active_boundaries'],flush=True)
    # A fixed-density optimum can reveal a better basin for the free model.
    best = min(records,key=lambda r:r["objective"])
    check = profile.solve([np.array([best["parameters"][n] for n in problem.family.names])],None)
    if check["objective"] < free["objective"]:
        state["free"] = check
    state["status"] = "complete"
    state["finished_utc"] = datetime.now(timezone.utc).isoformat()
    path.write_text(json.dumps(state,indent=2)+"\n")
    return str(path)


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=2) as pool:
        print(list(pool.map(run_case,CASES)),flush=True)
