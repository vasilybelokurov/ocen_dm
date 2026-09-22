#!/usr/bin/env python3
"""Validate and illustrate fixed-density refits without changing archived fits."""
from datetime import datetime, timezone
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter
import numpy as np
from astropy.table import Table

from profile_dm_density_limit import ROOT, OUTPUT, CASES, digest
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import data_fingerprint


COLORS = {"cusp":"#2476b5", "core":"#d9791e", "core_lower_stellar_floor":"#9053a6"}
LABELS = {"cusp":"Cusp: original bounds", "core":"Core: original bounds",
          "core_lower_stellar_floor":r"Core: stellar floor $10^4\,M_\odot$"}
DS_LABELS = {"hst_pm_radial_ours":"HST radial", "hst_pm_tangential_ours":"HST tangential",
             "muse_los_dispersion":"MUSE LOS", "gaia_edr3_ours_radial":"Gaia radial",
             "gaia_edr3_ours_tangential":"Gaia tangential", "distance":"Distance constraint"}
RHO_LABEL = r"Local DM density at 20 pc [$M_\odot\,\mathrm{pc}^{-3}$]"
PLOTS = ROOT/"plots"
PRODUCTS = ROOT/"results/plot_data"


def residual(record, b):
    pred = np.asarray(record["prediction"][b["dataset"]])
    y = np.asarray(b["value"])
    return (pred-y)/np.where(pred>y, b["err_hi"], b["err_lo"])


def fixed(state, rho):
    return next(r for r in state["profile"] if np.isclose(r["rho"], rho))


def validate(states):
    """Replay every solution, check original inputs, and record per-bin evidence."""
    rows, checks, inputs = [], {}, {}
    for name, state in states.items():
        assert state["status"] == "complete"
        for path, expected in (state["input_sha256"] | state["source_sha256"]).items():
            assert digest(ROOT/path) == expected, path
        assert digest(ROOT/"bin/profile_dm_density_limit.py") == state["script_sha256"]
        inputs.update(state["input_sha256"])
        problem = problem_for(load_run_metadata(ROOT/"results/fits"/state["label"]))
        assert data_fingerprint(problem.data) == state["data_fingerprint"]
        errors = []
        for kind, records in [("free", [state["free"]]), ("fixed", state["profile"])]:
            for r in records:
                x = np.array([r["parameters"][n] for n in problem.family.names])
                prediction = problem.predict(x)
                chi = problem.chi2(x)
                lnlike = problem.loglike_vector(x)
                lognorm = sum(np.log(np.sqrt(2/np.pi)/(p.err_lo+p.err_hi)).sum()
                              for p in problem.data.profiles)
                np.testing.assert_allclose(-2*(lnlike-lognorm), r["chi2_kinematics"], atol=1e-8)
                distance = r["parameters"]["distance"]
                np.testing.assert_allclose(r["objective"], r["chi2_kinematics"]+((distance-5.43)/.05)**2)
                assert 5.03 < distance < 5.83  # numerical normal-prior bounds never active
                density = float(np.asarray(problem.family.halo(r["parameters"]).density(20.)).item())
                np.testing.assert_allclose(density, r["rho"], rtol=1e-8)
                for b in state["bins"]:
                    ds = b["dataset"]
                    np.testing.assert_allclose(prediction[ds], r["prediction"][ds], atol=1e-10)
                    errors.append(abs(chi[ds][0]-r["chi2_by_dataset"][ds]))
                    res = residual(r, b)
                    base = residual(state["free"], b)
                    for i in range(len(res)):
                        rows.append((name, kind, r["rho"], ds, i+1,
                                     b["radius_arcsec"][i], b["radius_arcsec"][i]*5360/206264.806,
                                     b["value"][i], b["err_lo"][i], b["err_hi"][i],
                                     float(prediction[ds][i]), float(res[i]), float(res[i]**2-base[i]**2)))
        spread = [np.ptp([d["objective"] for d in r["start_diagnostics"]]) for r in state["profile"]]
        checks[name] = dict(maximum_replayed_chi2_difference=max(errors),
                            maximum_three_start_objective_spread=float(max(spread)),
                            all_profile_polishes_successful=all(r["polish"]["success"] for r in state["profile"]))
    table = Table(rows=rows, names=["case", "fit_type", "rho20", "dataset", "bin",
                  "radius_arcsec", "radius_pc_reference", "observed", "err_lo", "err_hi",
                  "prediction", "residual_sigma", "delta_chi2_from_free"])
    table.meta.update(radius_reference_distance_kpc=5.36, density_units="Msun pc^-3",
                      residual_convention="(model-data)/error on the model side",
                      reference="Each case's free-density optimum; component changes are signed")
    table.write(PRODUCTS/"dm_density_profile_residuals.ecsv", overwrite=True)
    return checks, inputs


def profile_figure(states, limit):
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2))
    ax, zoom, mass, bars = axes.flat
    for name, s in states.items():
        records = sorted(s["profile"]+[s["free"]], key=lambda r:r["rho"])
        rho = np.array([r["rho"] for r in records])
        cost = np.array([r["objective"]-s["free"]["objective"] for r in records])
        style = "--" if "lower" in name else "-"
        for a in (ax, zoom):
            a.plot(rho, cost, style, color=COLORS[name], marker="o", ms=4, lw=2, label=LABELS[name])
        mass.plot(rho, [r["parameters"]["M_star"]/1e6 for r in records], style,
                  color=COLORS[name], marker="o", ms=4, lw=2)
    ax.set(xlim=(0,5.2), ylim=(-.5,33), title="(a) High densities eventually worsen the best fit",
           ylabel=r"$\Delta Q$ from each halo's free-density optimum")
    zoom.set(xlim=(0,2.2), ylim=(-.08,2.8), title="(b) Little fit penalty near the quoted upper limit",
             ylabel=r"$\Delta Q$ (enlarged view)")
    zoom.annotate(r"At $\rho=2$: $\Delta Q=1.15$ (core)", xy=(2,1.146), xytext=(.48,2.17),
                  arrowprops=dict(arrowstyle="->", color=".35"), fontsize=10)
    mass.set(xlim=(0,5.2), ylim=(-.08,2.65), title="(c) Extra DM is offset by fewer light-tracing stars",
             ylabel=r"Fitted stellar mass [$10^6\,M_\odot$]")
    mass.axhline(1., color=".5", lw=1, ls=":")
    mass.text(3.05,1.09,"Original stellar-mass floor", fontsize=10, color=".35")
    for a in (ax, zoom, mass):
        a.axvline(limit, color=".2", ls=":", lw=1.5)
        a.set_xlabel(RHO_LABEL)
        a.grid(alpha=.18)
    ax.text(limit+.12,30.7,"Bayesian 95% upper limit = 1.026", fontsize=10)
    names = list(states["core"]["free"]["chi2_by_dataset"])+["distance"]
    for offset, name in [(-.19,"core"), (.19,"core_lower_stellar_floor")]:
        s = states[name]; r = fixed(s,5.)
        values = [r["chi2_by_dataset"][ds]-s["free"]["chi2_by_dataset"][ds] for ds in names[:-1]]
        values += [r["chi2_distance"]-s["free"]["chi2_distance"]]
        bars.barh(np.arange(6)+offset, values, height=.35, color=COLORS[name], alpha=.85)
    bars.set_yticks(np.arange(6), [DS_LABELS[n] for n in names])
    bars.invert_yaxis(); bars.axvline(0,color=".3",lw=1)
    bars.set(xlim=(-7,14), xlabel=r"Change in $\chi^2$ relative to free-density fit",
             title=r"(d) Where the core fit worsens at $\rho=5$")
    bars.grid(axis="x",alpha=.18)
    fig.suptitle("Testing densities above the posterior upper limit", fontsize=17, y=.985)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5,.955), ncol=3, frameon=False, fontsize=10)
    fig.text(.5,.017, r"All remaining parameters refitted.  $Q=\chi^2_{\rm kinematics}+[(D-5.43)/0.05]^2$."
             "  These curves are not posterior probabilities.", ha="center",fontsize=10)
    fig.subplots_adjust(top=.885,bottom=.09,left=.075,right=.98,hspace=.36,wspace=.36)
    path=PLOTS/"dm_density_profile_limit.png"; fig.savefig(path,dpi=200);plt.close(fig)
    return path


def data_figure(states):
    s = states["core"]
    choices = [(s["free"],"#222222", "-", r"Free density: $\rho=0.48$"),
               (fixed(s,2.),"#2476b5", "--", r"$\rho=2$: original bounds"),
               (fixed(s,5.),"#ce453b", "-", r"$\rho=5$: original bounds"),
               (fixed(states["core_lower_stellar_floor"],5.),"#9053a6", ":",
                r"$\rho=5$: lower stellar floor")]
    datasets = ["gaia_edr3_ours_radial", "gaia_edr3_ours_tangential", "hst_pm_radial_ours"]
    fig, axes = plt.subplots(2,3,figsize=(14.5,7.7),sharex="col",gridspec_kw={"height_ratios":[1.5,1]})
    for j,ds in enumerate(datasets):
        b = next(b for b in s["bins"] if b["dataset"]==ds)
        radius = np.asarray(b["radius_arcsec"])*5360/206264.806
        ax, res = axes[:,j]
        ax.errorbar(radius,b["value"],yerr=[b["err_lo"],b["err_hi"]],fmt="o",color="k",ms=4,
                    elinewidth=1,capsize=2,label="Data",zorder=5)
        for r,color,style,label in choices:
            ax.plot(radius,r["prediction"][ds],style,color=color,lw=2,marker=".",label=label)
            res.plot(radius,residual(r,b),style,color=color,lw=2,marker=".")
        res.axhspan(-1,1,color=".93"); res.axhline(0,color=".5",lw=.8)
        res.set_ylim(-4.3,4.5)
        ax.set_title(DS_LABELS[ds]);ax.set_ylabel(r"Dispersion [mas yr$^{-1}$]")
        res.set_ylabel("(Model − data) / error");res.set_xlabel("Projected radius [pc]")
        for a in (ax,res):
            a.set_xscale("log");a.grid(alpha=.18);a.set_xticks([10,20,30,50] if j<2 else [.1,.3,1,3,8])
            a.xaxis.set_major_formatter(ScalarFormatter());a.xaxis.set_minor_formatter(NullFormatter())
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,.948),ncol=3,frameon=False,fontsize=10)
    fig.suptitle("Which measurements resist very high DM density?  Flexible cored halo",fontsize=16,y=.99)
    fig.text(.5,.02,"Predictions use each fit's distance and the archived bin averaging / rotation correction. "
             "Display radii use D = 5.36 kpc.\nDensity units: solar masses per cubic parsec. "
             "Lines connect bin predictions; they are not continuous model profiles.",ha="center",fontsize=10)
    fig.subplots_adjust(top=.82,bottom=.12,left=.065,right=.99,hspace=.09,wspace=.29)
    path=PLOTS/"dm_density_profile_data.png";fig.savefig(path,dpi=200);plt.close(fig)
    return path


def main():
    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
    states = {name:json.loads((OUTPUT/f"{name}.json").read_text()) for name,_,_ in CASES}
    checks, inputs = validate(states)
    posterior = ROOT/"results/plot_data/dm_density_20pc_posterior.json"
    limit = json.loads(posterior.read_text())["summaries"]["flexible"]["quantiles"]["p95"]
    figures = [profile_figure(states,limit),data_figure(states)]
    supplement = OUTPUT/"high_density_multistart_check.json"
    multistart = json.loads(supplement.read_text())
    for name, check in multistart.items():
        assert check["input_sha256"] == digest(OUTPUT/f"{name}.json")
        # A materially better solution must be incorporated before reporting this profile.
        assert abs(check["check"]["objective"]-fixed(states[name],5.)["objective"]) < .01
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        density_units="Msun pc^-3",model_averaged_bayesian_upper95=limit,
        objective="Q = chi2_kinematics + ((D-5.43)/0.05)^2; other priors are bounds only",
        interpretation="Multi-start profile objective, not a marginal posterior or calibrated confidence interval",
        validation=checks,high_density_multistart_check=multistart,
        inputs_sha256=inputs,
        profile_sha256={str((OUTPUT/f"{n}.json").relative_to(ROOT)):digest(OUTPUT/f"{n}.json") for n in states},
        source_sha256=states["core"]["source_sha256"],
        posterior_sha256=digest(posterior),script_sha256=digest(Path(__file__)),
        figure_sha256={str(p.relative_to(ROOT)):digest(p) for p in figures},
        table_sha256=digest(PRODUCTS/"dm_density_profile_residuals.ecsv"),cases={})
    for name,s in states.items():
        summary["cases"][name]=dict(stellar_floor=s["stellar_floor"],free=s["free"],profile=s["profile"])
    (PRODUCTS/"dm_density_profile_limit.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(checks,indent=2))
    print("\n".join(str(p) for p in figures))


if __name__=="__main__":
    main()
