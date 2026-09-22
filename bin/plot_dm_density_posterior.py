#!/usr/bin/env python3
"""Posterior local DM density at 20 pc, including evidence-weighted no-DM models.

The ECSV posteriors are already equally weighted posterior draws. Never apply
the likelihood to them a second time. Model probabilities use Z times the
explicit model prior, with a separate probability atom at zero for no DM.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp

from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import data_fingerprint


ROOT = Path(__file__).resolve().parents[1]
GROUPS = {"simple": "simple_beta", "flexible": "turnover"}
FAMILIES = ("K1", "K2_nfw", "K2_cored")
NAMES = ("No DM", "Cusp", "Core")
COLORS = ("#777777", "#0072B2", "#D55E00")
MODEL_PRIORS = np.array([.5, .25, .25])
RADIUS_PC = 20.


def model_probabilities(logz, priors):
    logz, priors = np.asarray(logz, float), np.asarray(priors, float)
    if logz.shape != priors.shape or logz.ndim != 1:
        raise ValueError("One evidence and prior are required per model")
    if not np.all(np.isfinite(logz)) or not np.all(np.isfinite(priors)):
        raise ValueError("Non-finite evidence or prior")
    if np.any(priors <= 0) or not np.isclose(priors.sum(), 1.):
        raise ValueError("Model priors must be positive and sum to one")
    logw = logz + np.log(priors)
    return np.exp(logw-logsumexp(logw))


def density_at_radius(m100, rs, gamma, rt, radius=RADIUS_PC, order=96):
    """Normalise the fitted truncated gNFW density by M(<100 pc)."""
    m100, rs = np.broadcast_arrays(np.atleast_1d(m100), np.atleast_1d(rs))
    if np.any(m100 <= 0) or np.any(rs <= 0) or rt <= 0 or radius <= 0:
        raise ValueError("Masses and length scales must be positive")
    nodes, weights = leggauss(order)
    r = 50*(nodes+1)
    x = r[None, :]/rs[:, None]
    shape = x**(-gamma)*(1+x)**(gamma-3)*(1+(r[None, :]/rt)**2)**-2
    unit_m100 = 4*np.pi*50*((r[None, :]**2*shape)@weights)
    x0 = radius/rs
    return m100/unit_m100*x0**(-gamma)*(1+x0)**(gamma-3)*(1+(radius/rt)**2)**-2


def mixture_samples(values, probabilities):
    """Give each posterior draw its model's probability divided by draw count."""
    probabilities = np.asarray(probabilities, float)
    if len(values) != len(probabilities) or np.any(probabilities < 0) or not np.isclose(probabilities.sum(), 1):
        raise ValueError("Invalid mixture probabilities")
    if any(np.asarray(v).ndim != 1 or len(v) == 0 or not np.all(np.isfinite(v)) or np.any(v < 0) for v in values):
        raise ValueError("Density samples must be finite, non-negative arrays")
    sample = np.concatenate(values)
    weight = np.concatenate([np.full(len(v), p/len(v)) for v, p in zip(values, probabilities)])
    return sample, weight


def weighted_quantiles(values, weights, levels=(.16, .5, .84, .95)):
    """Inverse empirical CDF: retains a discrete atom at zero exactly."""
    order = np.argsort(values)
    cdf = np.cumsum(np.asarray(weights)[order])
    cdf /= cdf[-1]
    idx = np.minimum(np.searchsorted(cdf, levels, side="left"), len(order)-1)
    return np.asarray(values)[order[idx]]


def summarise(values, probabilities):
    sample, weight = mixture_samples(values, probabilities)
    q = weighted_quantiles(sample, weight)
    positive = sample > 0
    conditional = weighted_quantiles(sample[positive], weight[positive])
    return dict(probability_zero=float(weight[~positive].sum()),
                probability_dm=float(weight[positive].sum()),
                quantiles=dict(zip(("p16", "p50", "p84", "p95"), map(float, q))),
                conditional_on_dm_quantiles=dict(zip(("p16", "p50", "p84", "p95"), map(float, conditional))))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    figure = ROOT/"plots/dm_density_20pc_posterior.png"
    prefix = ROOT/"results/plot_data/dm_density_20pc_posterior"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    figure.parent.mkdir(exist_ok=True)
    runs, exports, hashes = {}, [], {}
    identity = None
    for group, suffix in GROUPS.items():
        runs[group] = []
        for family, name in zip(FAMILIES, NAMES):
            label = f"rung2_{family}_{suffix}"
            run = load_run_metadata(ROOT/"results/fits"/label)
            if run["run"]["status"] != "complete":
                raise ValueError(f"Incomplete run: {label}")
            problem = problem_for(run)  # Replays the saved maximum likelihood.
            fingerprint = data_fingerprint(problem.data)
            if identity is not None and fingerprint != identity:
                raise ValueError("Evidence comparisons require identical observations")
            identity = fingerprint
            s = run["summary"]
            posterior = Table.read(run["dir"]/"posterior.ecsv")
            if family == "K1":
                rho = np.zeros(1)  # One atom, independent of posterior sample count.
                error = 0.
            else:
                config = s["model"]
                rho = density_at_radius(posterior["M_dm_100"], posterior["r_s"], config["gamma"], config["r_t"])
                # An independent route through the production mass component.
                errors = []
                for i in np.linspace(0, len(posterior)-1, 9, dtype=int):
                    theta = {n: float(posterior[n][i]) for n in problem.family.names}
                    expected = float(np.asarray(problem.family.halo(theta).density(RADIUS_PC)).reshape(-1)[0])
                    errors.append(abs(rho[i]/expected-1))
                error = max(errors)
                if error > 1e-5:
                    raise ValueError(f"Density disagrees with production halo: {label}")
            record = dict(label=label, name=name, logz=s["logz"], logzerr=s["logzerr"],
                          lnL_max=s["lnL_max"], posterior_samples=len(posterior),
                          resolved_model=s["model"], local_density_model_check_max_relative_error=error,
                          conditional_quantiles=dict(zip(("p16", "p50", "p84", "p95"), map(float, np.quantile(rho, [.16,.5,.84,.95])))))
            runs[group].append(dict(record=record, rho=rho))
            exports.append(Table(dict(group=np.repeat(group, len(rho)), model=np.repeat(label, len(rho)),
                                      sample_index=np.arange(len(rho)), rho_dm_20pc=rho)))
            for filename in ("summary.json", "posterior.ecsv", "run.yaml", "data_snapshot.json"):
                path = run["dir"]/filename
                hashes[str(path.relative_to(ROOT))] = digest(path)

    rng = np.random.default_rng(42)
    summaries = {}
    for group, rows in runs.items():
        logz = np.array([r["record"]["logz"] for r in rows])
        sigma = np.array([r["record"]["logzerr"] for r in rows])
        probability = model_probabilities(logz, MODEL_PRIORS)
        values = [r["rho"] for r in rows]
        result = summarise(values, probability)
        # Numerical evidence uncertainty only; the posterior shapes are held fixed.
        trials = rng.normal(logz, sigma, size=(20000, 3))+np.log(MODEL_PRIORS)
        trials = np.exp(trials-logsumexp(trials, axis=1)[:,None])
        result.update(model_probabilities=probability.tolist(),
                      model_probability_numerical_p16_p84=np.quantile(trials,[.16,.84],axis=0).tolist(),
                      bayes_factor_dm_vs_no_dm=float(np.exp(logsumexp(logz[1:])-np.log(2)-logz[0])),
                      equal_three_model_prior_probabilities=model_probabilities(logz,np.ones(3)/3).tolist())
        rows[0]["probability"] = probability
        summaries[group] = result

    allrows = [r for rows in runs.values() for r in rows]
    joint = model_probabilities([r["record"]["logz"] for r in allrows], np.tile(MODEL_PRIORS/2, 2))
    summaries["joint"] = summarise([r["rho"] for r in allrows], joint)
    summaries["joint"].update(model_probabilities=joint.tolist(),
                              anisotropy_probabilities={g:float(joint[3*i:3*i+3].sum()) for i,g in enumerate(GROUPS)})

    # Point densities for the simulation ICs, not finite-shell estimates.
    simulation = {}
    for shape in ("cusp", "core"):
        path = ROOT/"results/dynamics/lifetime_20260921T162000Z/runs"/f"{shape}_base_class1_live_nucleus/equilibrium.json"
        eq = json.loads(path.read_text())["components"]["dm"]
        gamma = 1 if shape == "cusp" else 0
        rho = eq["density_norm_msun_kpc3"]/1e9*2**(gamma-3)*np.exp(-(.2)**2)
        simulation[shape] = dict(initial_rho_dm_20pc=float(rho), initial_mass_inside_20pc=1e5)
        hashes[str(path.relative_to(ROOT))] = digest(path)
    for result, rows, probability in [(summaries[g],runs[g],runs[g][0]["probability"]) for g in GROUPS]+[(summaries["joint"],allrows,joint)]:
        sample, weight = mixture_samples([r["rho"] for r in rows],probability)
        result["probability_above_simulation_initial_density"] = {s:float(weight[sample>=v["initial_rho_dm_20pc"]].sum()) for s,v in simulation.items()}

    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
    fig, axes = plt.subplots(2,2,figsize=(11.5,8),gridspec_kw={"height_ratios":[1,2]})
    fig.subplots_adjust(left=.09,right=.98,bottom=.12,top=.87,hspace=.40,wspace=.24)
    fig.suptitle("Local dark-matter density at 20 pc",fontsize=19,y=.975)
    fig.text(.5,.932,"Prior: 50% no DM, 25% cusp, 25% core, within each anisotropy family",ha="center")
    positive = np.concatenate([r["rho"] for r in allrows if r["record"]["name"] != "No DM"])
    bins = np.linspace(np.floor(np.log10(positive.min())),np.ceil(np.log10(positive.max())),64)
    for col,(group,rows) in enumerate(runs.items()):
        prob = rows[0]["probability"]
        result = summaries[group]
        ax = axes[0,col]
        ax.bar(NAMES,100*prob,color=COLORS,width=.6)
        ax.set(ylim=(0,75),ylabel="Posterior model probability [%]",title=f"Rung 2: {group} anisotropy")
        low,high = np.array(result["model_probability_numerical_p16_p84"])*100
        ax.errorbar(np.arange(3),100*prob,yerr=[100*prob-low,high-100*prob],fmt="none",color="black",capsize=3)
        for i,p in enumerate(prob):
            ax.text(i,max(100*p,high[i])+2,f"{100*p:.1f}%",ha="center",fontsize=11)
        ax.grid(axis="y",alpha=.15)
        ax = axes[1,col]
        total = np.zeros(len(bins)-1)
        for i in [1,2]:
            counts,_ = np.histogram(np.log10(rows[i]["rho"]),bins=bins)
            density = counts/len(rows[i]["rho"])*prob[i]/np.diff(bins)
            total += density
            ax.stairs(density,10**bins,color=COLORS[i],lw=1.6,label=f"{NAMES[i]} contribution")
        ax.stairs(total,10**bins,color="#242424",lw=2,label="Combined positive density")
        np.testing.assert_allclose(np.sum(total*np.diff(bins)),1-prob[0],atol=1e-12)
        median = result["quantiles"]["p50"]
        upper = result["quantiles"]["p95"]
        ax.axvline(median,color="#7B3294",ls="--",lw=1.6)
        ax.axvline(upper,color="#242424",ls=":",lw=1.2)
        ax.text(.025,.96,f"P(density = 0) = {prob[0]:.1%} (bar above)\n"
                f"Median including zero (dashed) = {median:.3g}\n"
                f"95% upper limit (dotted) = {upper:.3g} solar masses/pc³",
                transform=ax.transAxes,va="top",fontsize=10,
                bbox=dict(facecolor="white",edgecolor="none",alpha=.95,pad=2))
        ax.set(xscale="log",xlim=(10**bins[0],10**bins[-1]),
               xlabel=r"Local $\rho_{\rm DM}(20\,\mathrm{pc})$ [$M_\odot\,\mathrm{pc}^{-3}$]",
               ylabel="Posterior probability density per dex")
        ax.set_ylim(0,max(total)*1.55)
        ax.legend(loc="upper left",bbox_to_anchor=(0,.75),fontsize=9,
                  frameon=True,facecolor="white",edgecolor="none",framealpha=.95)
        ax.grid(alpha=.15)
    fig.text(.5,.044,"Bar errors propagate numerical uncertainty in ln Z; they do not include model uncertainty.",ha="center",fontsize=9)
    fig.text(.5,.021,"The zero-density probability is discrete. Positive curves integrate to P(DM), not to one.",ha="center",fontsize=10)
    fig.savefig(figure,dpi=200)
    plt.close(fig)

    samples = vstack(exports,metadata_conflicts="silent")
    samples["rho_dm_20pc"].unit = "solMass / pc3"
    samples.meta["meaning"] = "Equally weighted within-model posterior samples; each no-DM model is represented by one zero atom. Model weights are in the JSON."
    samples.write(prefix.with_suffix(".ecsv"),overwrite=True)
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(),radius_pc=RADIUS_PC,
        density_definition="Local point density, not M(<20 pc)/(4 pi 20^3/3) or a finite-shell average",
        model_priors=dict(zip(NAMES,map(float,MODEL_PRIORS))),anisotropy_priors={"simple":.5,"flexible":.5},
        numerical_probability_intervals="Independent Gaussian lnZ errors; fixed posterior density shapes; 20000 draws, seed 42. Not model-systematic uncertainty.",
        data_fingerprint=identity,models={g:[r["record"] for r in rows] for g,rows in runs.items()},
        summaries=summaries,simulation_initial_conditions=simulation,input_sha256=hashes,
        script_sha256=digest(Path(__file__)),figure_sha256=digest(figure),samples_sha256=digest(prefix.with_suffix(".ecsv")))
    prefix.with_suffix(".json").write_text(json.dumps(provenance,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(summaries=summaries,simulation_initial_conditions=simulation),indent=2))


if __name__ == "__main__":
    main()
