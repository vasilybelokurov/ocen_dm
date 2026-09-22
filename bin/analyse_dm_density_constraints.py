#!/usr/bin/env python3
"""Trace data contributions to the 20-pc density posterior using saved fits.

Likelihood tempering is a sensitivity calculation, not a replacement for
independent fits with datasets removed. Cache all unique posterior predictions
and count duplicate samples when computing the original posterior weights.
"""
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import numpy as np
from astropy.table import Table
from scipy.special import logsumexp

from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import data_fingerprint
from plot_dm_density_posterior import density_at_radius, model_probabilities, weighted_quantiles


ROOT = Path(__file__).resolve().parents[1]
GROUPS = {"simple":"simple_beta", "flexible":"turnover"}
FAMILIES = ("K1","K2_nfw","K2_cored")
PRIORS = np.array([.5,.25,.25])
DATA_NAMES = {"hst_pm_radial_ours":"HST radial", "hst_pm_tangential_ours":"HST tangential",
              "muse_los_dispersion":"MUSE line of sight", "gaia_edr3_ours_radial":"Gaia radial",
              "gaia_edr3_ours_tangential":"Gaia tangential"}
CACHE = ROOT/"results/diagnostics/dm_density_constraints"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate(label):
    directory = ROOT/"results/fits"/label
    inputs = {str((directory/f).relative_to(ROOT)):digest(directory/f)
              for f in ("summary.json","run.yaml","posterior.ecsv","data_snapshot.json")}
    source = {str(p.relative_to(ROOT)):digest(p) for p in (ROOT/"src/ocen_dm").rglob("*.py")}
    # Include the conversion implementation as well as all likelihood source.
    source["bin/plot_dm_density_posterior.py"] = digest(ROOT/"bin/plot_dm_density_posterior.py")
    key = hashlib.sha256(json.dumps(dict(inputs=inputs,source=source,cache_version=1),sort_keys=True).encode()).hexdigest()[:16]
    path = CACHE/f"{label}_{key}.npz"
    meta_path = path.with_suffix(".json")
    if path.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        if digest(path) != meta["cache_sha256"]:
            raise ValueError(f"Changed cache: {path}")
        print(f"{label}: verified cached predictions",flush=True)
        return str(path)
    run = load_run_metadata(directory)
    if run["run"]["status"] != "complete":
        raise ValueError(f"Incomplete run: {label}")
    problem = problem_for(run)
    table = Table.read(directory/"posterior.ecsv")
    theta, counts = np.unique(np.column_stack([table[n] for n in problem.family.names]),axis=0,return_counts=True)
    loglike = np.empty((len(theta),problem.data.n_points))
    prediction = np.empty_like(loglike)
    bins = []
    for p in problem.data.profiles:
        for i in range(p.n):
            bins.append(dict(dataset=p.name, index=i, radius_arcsec=float(p.r[i]),
                             radius_pc_reference=float(p.r[i]*5360/206264.806),
                             value=float(p.value[i]),err_lo=float(p.err_lo[i]),err_hi=float(p.err_hi[i])))
    print(f"{label}: evaluating {len(theta)} unique draws ({len(table)} posterior samples)",flush=True)
    for i,x in enumerate(theta):
        pred = problem.predict(x)
        prediction[i] = np.concatenate([pred[p.name] for p in problem.data.profiles])
        loglike[i] = np.concatenate([problem.likelihood._split_normal_lnlike(pred[p.name],p) for p in problem.data.profiles])
        if i in (0,len(theta)//2,len(theta)-1):
            np.testing.assert_allclose(loglike[i].sum(),problem.loglike_vector(x),atol=1e-8,rtol=0)
        if (i+1)%1000 == 0:
            print(f"{label}: {i+1}/{len(theta)}",flush=True)
    if not np.all(np.isfinite(loglike)):
        raise ValueError(f"Non-finite posterior likelihood: {label}")
    s = run["summary"]
    if "M_dm_100" in problem.family.names:
        rho = density_at_radius(theta[:,problem.family.names.index("M_dm_100")],
                                theta[:,problem.family.names.index("r_s")],s["model"]["gamma"],s["model"]["r_t"])
    else:
        rho = np.zeros(len(theta))
    CACHE.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,theta=theta,counts=counts,loglike=loglike,prediction=prediction,rho=rho)
    meta = dict(label=label,logz=s["logz"],logzerr=s["logzerr"],bins=bins,parameter_names=problem.family.names,
                data_fingerprint=data_fingerprint(problem.data),input_sha256=inputs,source_sha256=source,
                cache_sha256=digest(path),created_utc=datetime.now(timezone.utc).isoformat())
    meta_path.write_text(json.dumps(meta,indent=2)+"\n")
    print(f"{label}: finished",flush=True)
    return str(path)


def tempered_weights(base_weights, removed_loglike, fraction):
    """Change L_group into L_group**(1-fraction); return weights and ln(Znew/Z)."""
    if not 0 <= fraction <= 1:
        raise ValueError("Removal fraction must lie in [0,1]")
    logw = np.log(base_weights/base_weights.sum())-fraction*np.asarray(removed_loglike)
    correction = logsumexp(logw)
    return np.exp(logw-correction),float(correction)


def distribution(rows,mask,fraction,split=None):
    corrections, weights, diagnostics, densities = [], [], [], []
    for row in rows:
        keep = np.ones(len(row["rho"]),bool) if split is None else row["split"]==split
        weight,correction = tempered_weights(row["counts"][keep].astype(float),row["loglike"][keep][:,mask].sum(axis=1),fraction)
        corrections.append(correction)
        weights.append(weight)
        densities.append(row["rho"][keep])
        diagnostics.append(dict(effective_unique_draws=float(1/np.sum(weight**2)),largest_unique_weight=float(weight.max())))
    model_p = model_probabilities(np.array([r["meta"]["logz"] for r in rows])+corrections,PRIORS)
    rho = np.concatenate(densities)
    weight = np.concatenate([w*p for w,p in zip(weights,model_p)])
    positive = rho>0
    q = weighted_quantiles(rho,weight)
    qc = weighted_quantiles(rho[positive],weight[positive])
    return dict(model_probabilities=model_p.tolist(),probability_dm=float(1-model_p[0]),
                rho_quantiles=dict(zip(("p16","p50","p84","p95"),map(float,q))),
                rho_conditional_dm_quantiles=dict(zip(("p16","p50","p84","p95"),map(float,qc))),
                probability_rho_gt_1=float(weight[rho>1].sum()),
                effective_unique_draws=float(1/np.sum(weight**2)),
                within_model_diagnostics=diagnostics)


def local_tail_sensitivity(rows):
    """Derivative of P(rho>1) when weakening each bin, including model weights."""
    model_p = model_probabilities([r["meta"]["logz"] for r in rows],PRIORS)
    rho = np.concatenate([r["rho"] for r in rows])
    w = np.concatenate([p*r["counts"]/r["counts"].sum() for p,r in zip(model_p,rows)])
    ll = np.concatenate([r["loglike"] for r in rows])
    def derivative(keep):
        ww = w[keep]/w[keep].sum()
        indicator = (rho[keep]>1).astype(float)
        p = ww@indicator
        return -((ww*(indicator-p))@ll[keep])
    return derivative(np.ones(len(rho),bool)),derivative(rho>0)


def posterior_branches(row):
    """Expose nuisance-parameter tradeoffs within the flexible cored family."""
    result = {}
    names = row["meta"]["parameter_names"]
    for branch,keep in [("low",row["rho"]<.1),("high",row["rho"]>1.)]:
        w = row["counts"][keep].astype(float)
        w /= w.sum()
        prediction = row["prediction"][keep]
        result[branch] = dict(posterior_fraction=float(row["counts"][keep].sum()/row["counts"].sum()),
            parameter_medians={n:float(weighted_quantiles(row["theta"][keep,i],w,[.5])[0]) for i,n in enumerate(names)},
            prediction_mean=np.average(prediction,axis=0,weights=w).tolist(),
            prediction_p16_p84=[weighted_quantiles(prediction[:,i],w,[.16,.84]).tolist() for i in range(prediction.shape[1])])
    return result


def figures(results,point_rows,branches,bins):
    assets = ROOT/"plots"
    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
    names = list(DATA_NAMES.values())
    dataset_colors = ["#4477AA","#66AADD","#777777","#D55E00","#EEAA33"]
    fig,axes = plt.subplots(2,2,figsize=(13,9),gridspec_kw={"height_ratios":[1,1.2]})
    fig.subplots_adjust(left=.12,right=.98,top=.87,bottom=.10,wspace=.34,hspace=.44)
    fig.suptitle("Which measurements limit the DM density at 20 pc?",fontsize=19,y=.98)
    fig.text(.5,.936,"Sensitivity of the existing posterior; no-DM, cusp and core weights are updated together",ha="center",fontsize=11)
    for col,group in enumerate(GROUPS):
        result = results[group]
        baseline = result["baseline"]["rho_quantiles"]["p95"]
        ax = axes[0,col]
        for i,(name,color) in enumerate(zip(names,dataset_colors)):
            r = result["experiments"][name]["0.2"]
            if not r["overlap_screen_passed"]:
                raise ValueError("Do not plot a failed overlap screen as a measured upper limit")
            value = r["rho_quantiles"]["p95"]
            ax.plot([baseline,value],[i,i],color=color,lw=3)
            ax.plot(value,i,"o",color=color,ms=7)
            side = 1 if value>=baseline else -1
            ax.text(value+side*.008,i,f"{value:.3f}",va="center",ha="left" if side>0 else "right",fontsize=10)
        ax.axvline(baseline,color="black",ls=":",lw=1.4)
        ax.set(yticks=range(len(names)),yticklabels=names,ylim=(4.6,-.7),xlim=(.86,1.36),
               xlabel=r"95% upper limit on $\rho_{\rm DM}(20\,\mathrm{pc})$ [$M_\odot\,\mathrm{pc}^{-3}$]",
               title=f"{group.capitalize()} anisotropy: weaken one dataset by 20%")
        ax.text(.98,.98,f"All data: {baseline:.3f} (dotted)",transform=ax.transAxes,ha="right",va="top",fontsize=10)
        ax.grid(axis="x",alpha=.15)
        ax = axes[1,col]
        for (key,name),color,marker in zip(DATA_NAMES.items(),dataset_colors,["o","s","+","o","^"]):
            points = [r for r in point_rows if r["group"]==group and r["dataset"]==key]
            x = np.array([r["radius_pc_reference"] for r in points])
            y = np.array([r["tail_probability_derivative"] for r in points])*10
            ax.plot(x,y,marker=marker,color=color,lw=.8,ms=4,label=name)
        ax.axhline(0,color="black",lw=.8)
        ax.axvline(20,color="grey",ls=":",lw=.8)
        ax.set(xscale="log",xlim=(.07,70),xlabel="Projected radius [pc; distance = 5.36 kpc]",
               ylabel="Change in P(density > 1) [percentage points]",
               title="Individual bins: response to 10% weakening (linear)")
        ax.grid(alpha=.15)
        ax.legend(fontsize=8.5,frameon=True,facecolor="white",edgecolor="none",loc="lower left")
    lo = min(a.get_ylim()[0] for a in axes[1]);hi=max(a.get_ylim()[1] for a in axes[1])
    for ax in axes[1]:
        ax.set_ylim(lo,hi)
    fig.text(.5,.035,"Positive lower-panel values identify bins that suppress high DM density; negative values identify bins that favour it.",ha="center",fontsize=10)
    paths = [assets/"dm_density_data_constraints.png"]
    fig.savefig(paths[-1],dpi=200);plt.close(fig)

    fig,axes = plt.subplots(2,2,figsize=(11.5,7.5),sharex="col",gridspec_kw={"height_ratios":[2,1]})
    fig.subplots_adjust(left=.09,right=.98,top=.81,bottom=.12,wspace=.19,hspace=.10)
    fig.suptitle("How high-density solutions change the Gaia fit",fontsize=18,y=.975)
    fig.text(.5,.935,"Flexible anisotropy, cored halo: posterior subsets with all other parameters free",ha="center",fontsize=11)
    colors = {"low":"#0072B2","high":"#D55E00"}
    legends = [Line2D([],[],color="black",marker="o",mfc="white",ls="",label="Gaia data")]
    legends += [Line2D([],[],color=colors[b],lw=2,label=label) for b,label in [("low",r"$\rho_{20}<0.1$"),("high",r"$\rho_{20}>1$ [$M_\odot\,\mathrm{pc}^{-3}$]")]]
    fig.legend(handles=legends,loc="upper center",bbox_to_anchor=(.5,.916),ncol=3,frameon=False,fontsize=11)
    for col,key in enumerate(["gaia_edr3_ours_radial","gaia_edr3_ours_tangential"]):
        ids = np.array([i for i,b in enumerate(bins) if b["dataset"]==key])
        x = np.array([bins[i]["radius_pc_reference"] for i in ids])
        data = np.array([bins[i]["value"] for i in ids])
        lo = np.array([bins[i]["err_lo"] for i in ids]);hi=np.array([bins[i]["err_hi"] for i in ids])
        ax = axes[0,col];res = axes[1,col]
        ax.errorbar(x,data,yerr=[lo,hi],fmt="o",mfc="white",mec="black",color="black",ms=4,lw=.9,zorder=5)
        for branch in ["low","high"]:
            mean = np.array(branches[branch]["prediction_mean"])[ids]
            q = np.array(branches[branch]["prediction_p16_p84"])[ids]
            ax.plot(x,mean,color=colors[branch],lw=1.8)
            ax.fill_between(x,q[:,0],q[:,1],color=colors[branch],alpha=.14,lw=0)
            residual = (data-mean)/np.where(mean>data,hi,lo)
            res.plot(x,residual,"o-",color=colors[branch],lw=1.4,ms=4)
        ax.set(title=DATA_NAMES[key],ylabel="Dispersion [mas/yr]",ylim=(.17,.58))
        res.axhspan(-1,1,color="grey",alpha=.12,lw=0);res.axhline(0,color="grey",lw=.8)
        res.set(xlabel="Projected radius [pc; distance = 5.36 kpc]",ylabel="(Data − model) / error",ylim=(-2.3,3.8))
        for a in [ax,res]:
            a.set(xscale="log",xlim=(8,63));a.grid(alpha=.12)
        res.set_xticks([10,20,30,50],labels=["10","20","30","50"])
        res.xaxis.set_minor_formatter(NullFormatter())
    fig.text(.5,.029,"Lines: conditional posterior mean predictions. Bands: 16th–84th percentiles within each density subset.",ha="center",fontsize=10)
    paths.append(assets/"dm_density_gaia_tradeoff.png")
    fig.savefig(paths[-1],dpi=200);plt.close(fig)
    return paths


def analyse(paths):
    rng = np.random.default_rng(57)
    loaded = {}
    for path in paths:
        path = Path(path)
        with np.load(path) as a:
            row = {k:a[k] for k in a.files}
        row["meta"] = json.loads(path.with_suffix(".json").read_text())
        row["cache"] = str(path.relative_to(ROOT))
        row["split"] = rng.integers(0,2,len(row["rho"]))
        loaded[row["meta"]["label"]] = row
    assert len({r["meta"]["data_fingerprint"] for r in loaded.values()})==1
    first = next(iter(loaded.values()))
    bins = first["meta"]["bins"]
    assert all(r["meta"]["bins"]==bins for r in loaded.values())
    dataset = np.array([b["dataset"] for b in bins])
    radius = np.array([b["radius_pc_reference"] for b in bins])
    masks = {name:dataset==key for key,name in DATA_NAMES.items()}
    masks.update({"All HST":np.char.startswith(dataset,"hst"),"All Gaia":np.char.startswith(dataset,"gaia"),
                  "Gaia below 20 pc":np.char.startswith(dataset,"gaia")&(radius<20),
                  "Gaia above 20 pc":np.char.startswith(dataset,"gaia")&(radius>=20)})
    results,point_rows = {},[]
    for group,suffix in GROUPS.items():
        rows = [loaded[f"rung2_{f}_{suffix}"] for f in FAMILIES]
        base = distribution(rows,np.zeros(len(bins),bool),0.)
        results[group] = dict(baseline=base,experiments={})
        full,conditional = local_tail_sensitivity(rows)
        for i,b in enumerate(bins):
            point_rows.append(dict(group=group,**b,tail_probability_derivative=float(full[i]),
                                   conditional_tail_probability_derivative=float(conditional[i])))
        for name,mask in masks.items():
            out = {}
            for fraction in [.1,.2,1.]:
                result = distribution(rows,mask,fraction)
                halves = [distribution(rows,mask,fraction,split=i) for i in [0,1]]
                q95 = result["rho_quantiles"]["p95"]
                result["split_half_relative_q95_difference"] = abs(halves[0]["rho_quantiles"]["p95"]-halves[1]["rho_quantiles"]["p95"])/q95
                result["split_half_probability_dm_difference"] = abs(halves[0]["probability_dm"]-halves[1]["probability_dm"])
                result["overlap_screen_passed"] = bool(all(d["effective_unique_draws"]>=200 and d["largest_unique_weight"]<=.02 for d in result["within_model_diagnostics"])
                    and result["split_half_relative_q95_difference"]<.15 and result["split_half_probability_dm_difference"]<.1)
                out[str(fraction)] = result
                if fraction in [.2,1.]:
                    print(group,name,fraction,'rho95',round(q95,3),'conditional95',round(result['rho_conditional_dm_quantiles']['p95'],3),
                          'PDM',round(result['probability_dm'],3),'minESS',round(min(d['effective_unique_draws'] for d in result['within_model_diagnostics'])),
                          'screen',result['overlap_screen_passed'],flush=True)
            results[group]["experiments"][name] = out
    branches = posterior_branches(loaded["rung2_K2_cored_turnover"])
    plot_paths = figures(results,point_rows,branches,bins)
    output = ROOT/"results/plot_data/dm_density_data_constraints.json"
    output.write_text(json.dumps(dict(created_utc=datetime.now(timezone.utc).isoformat(),
        method="Posterior importance reweighting after likelihood tempering; no independent refits",
        model_priors=PRIORS.tolist(),tail_threshold_msun_pc3=1.,radius_reference_distance_kpc=5.36,
        cache_sha256={r["cache"]:digest(ROOT/r["cache"]) for r in loaded.values()},
        script_sha256=digest(Path(__file__)),results=results,flexible_core_density_subsets=branches,
        figure_sha256={str(p.relative_to(ROOT)):digest(p) for p in plot_paths}),indent=2)+"\n")
    Table(rows=point_rows).write(output.with_suffix(".ecsv"),overwrite=True)
    return results,point_rows


def main():
    labels = [f"rung2_{f}_{suffix}" for suffix in GROUPS.values() for f in FAMILIES]
    with ProcessPoolExecutor(max_workers=2) as pool:
        paths = list(pool.map(evaluate,labels))
    analyse(paths)


if __name__ == "__main__":
    main()
