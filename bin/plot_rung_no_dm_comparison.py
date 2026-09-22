#!/usr/bin/env python3
"""Data, no-DM and lowest-chi-square DM fits across the completed ladder.

Figures belong in plots/; numerical exports and provenance in results/plot_data/.
Model curves join the exact predictions evaluated in the likelihood.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, NullFormatter, StrMethodFormatter
import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import data_fingerprint, family_config


ROOT = Path(__file__).resolve().parents[1]
STAGES = [
    ("Rung 0", "Isotropic", "rung0_K1", ["rung0_K2_nfw", "rung0_K2_cored"]),
    ("Rung 1", "Constant beta", "rung1_K1_slice", ["rung1_K2_nfw", "rung1_K2_cored_slice"]),
    ("Rung 2: simple", "One smooth transition", "rung2_K1_simple_beta",
     ["rung2_K2_nfw_simple_beta", "rung2_K2_cored_simple_beta"]),
    ("Rung 2: flexible", "Two smooth transitions", "rung2_K1_turnover",
     ["rung2_K2_nfw_turnover", "rung2_K2_cored_turnover"]),
]
DATASETS = [
    ("hst_pm_radial_ours", "HST: radial proper motion", r"$\sigma_{\mu,R}$ [mas yr$^{-1}$]"),
    ("hst_pm_tangential_ours", "HST: tangential proper motion", r"$\sigma_{\mu,T}$ [mas yr$^{-1}$]"),
    ("muse_los_dispersion", "MUSE: line of sight", r"$\sigma_{\rm LOS}$ [km s$^{-1}$]"),
    ("gaia_edr3_ours_radial", "Gaia: radial proper motion", r"$\sigma_{\mu,R}$ [mas yr$^{-1}$]"),
    ("gaia_edr3_ours_tangential", "Gaia: tangential proper motion", r"$\sigma_{\mu,T}$ [mas yr$^{-1}$]"),
]
STYLES = {
    "no_dm": dict(color="#0072B2", linestyle=(0, (5, 2.5)), lw=1.9,
                  marker="o", ms=2.7, markerfacecolor="white", markeredgewidth=.7),
    "dm": dict(color="#D55E00", linestyle="-", lw=1.9,
               marker="s", ms=2.7, markerfacecolor="white", markeredgewidth=.7),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_models():
    models, provenance, selections, hashes = [], [], [], {}
    fingerprint, data = None, None
    for stage, _, no_dm, halos in STAGES:
        candidates = [load_run_metadata(ROOT/"results/fits"/label) for label in halos]
        chosen = min(candidates, key=lambda row: row["summary"]["chi2_ml_total"])
        selections.append(dict(stage=stage, selected_dm=chosen["label"],
                               candidate_dm_chi2={row["label"]: row["summary"]["chi2_ml_total"]
                                                 for row in candidates}))
        pair = {}
        for category, run in [("no_dm", load_run_metadata(ROOT/"results/fits"/no_dm)), ("dm", chosen)]:
            problem = problem_for(run)  # Checks the saved best likelihood before plotting.
            identity = data_fingerprint(problem.data)
            if fingerprint is not None and identity != fingerprint:
                raise ValueError("Selected runs use different observations")
            fingerprint, data = identity, problem.data
            summary = run["summary"]
            x = np.array([summary["parameters"][name]["ml"] for name in problem.family.names])
            prediction, chi2 = problem.predict(x), problem.chi2(x)
            total = sum(value[0] for value in chi2.values())
            if not np.isclose(total, summary["chi2_ml_total"], rtol=0, atol=.011):
                raise ValueError(f"Chi-square replay failed: {run['label']}")
            residuals = {}
            for p in data.profiles:
                error = np.where(prediction[p.name] > p.value, p.err_hi, p.err_lo)
                residuals[p.name] = (p.value-prediction[p.name])/error
                np.testing.assert_allclose(np.sum(residuals[p.name]**2), chi2[p.name][0], rtol=1e-12)
            pair[category] = dict(label=run["label"], prediction=prediction, residuals=residuals,
                                  chi2=total, dataset_chi2={key: value[0] for key, value in chi2.items()})
            provenance.append(dict(stage=stage, category=category, label=run["label"], chi2=total,
                                   lnL_saved=summary["lnL_max"], data_fingerprint=identity,
                                   resolved_model=family_config(problem.family),
                                   parameters_ml=problem.family.to_dict(x),
                                   dataset_chi2=pair[category]["dataset_chi2"],
                                   replay="saved snapshot" if "data_snapshot" in summary else
                                          "legacy inputs; verified likelihood replay"))
        models.append(pair)
        for label in [no_dm, *halos]:
            for name in ("summary.json", "run.yaml", "data_snapshot.json"):
                path = ROOT/"results/fits"/label/name
                if path.exists():
                    hashes[str(path.relative_to(ROOT))] = digest(path)
    assert data.n_points == 89 and len(provenance) == 8
    return models, provenance, selections, hashes, data, fingerprint


def main():
    figure_path = ROOT/"plots/rung0_rung1_rung2_no_dm_vs_dm.png"
    data_prefix = ROOT/"results/plot_data/rung0_rung1_rung2_no_dm_vs_dm"
    figure_path.parent.mkdir(exist_ok=True)
    data_prefix.parent.mkdir(parents=True, exist_ok=True)
    models, provenance, selections, hashes, data, fingerprint = load_models()
    profiles = {p.name: p for p in data.profiles}
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(18.5, 17.2))
    grid = fig.add_gridspec(5, 4, left=.065, right=.987, bottom=.075, top=.857,
                          hspace=.34, wspace=.14)
    fig.suptitle("Data and best-fit no-DM / DM models across the ladder", fontsize=21, y=.987)
    fig.text(.5, .963, "Same 89 measurements in every column; best saved maximum-likelihood samples",
             ha="center", fontsize=12, color="#444444")
    handles = [Line2D([], [], color="#242424", linestyle="", marker="o", mfc="white", ms=5,
                      label="Data with 1-sigma errors")]
    handles += [Line2D([], [], **{k: v for k, v in STYLES[key].items() if k not in ("ms",)},
                       label=label) for key, label in [("no_dm", "No DM"), ("dm", "DM: cored halo")]]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .950), ncol=3,
               frameon=False, fontsize=12, handlelength=3.6, columnspacing=3)
    assert all("cored" in row["selected_dm"] for row in selections), "Update the halo legend for this selection"
    for col, (stage, description, _, _) in enumerate(STAGES):
        box = grid[0, col].get_position(fig)
        center = (box.x0+box.x1)/2
        fig.text(center, .913, stage, ha="center", fontsize=15, fontweight="bold")
        fig.text(center, .896, description, ha="center", fontsize=10.5, color="#444444")
        fig.text(center, .881, rf"$\chi^2$ (no DM / DM): {models[col]['no_dm']['chi2']:.2f} / {models[col]['dm']['chi2']:.2f}",
                 ha="center", fontsize=10.5)

    exports = []
    for row, (name, title, unit) in enumerate(DATASETS):
        p = profiles[name]
        values = [p.value-p.err_lo, p.value+p.err_hi]
        values += [pair[category]["prediction"][name] for pair in models for category in STYLES]
        ymin, ymax = min(np.min(v) for v in values), max(np.max(v) for v in values)
        padding = .12*(ymax-ymin)
        maximum = max(np.max(np.abs(pair[category]["residuals"][name]))
                      for pair in models for category in STYLES)
        limit = max(3., np.ceil(maximum*1.1))
        lower = np.maximum(p.r_lower, p.r.min()*.75)
        for col, pair in enumerate(models):
            sub = grid[row, col].subgridspec(2, 1, height_ratios=[3, 1.05], hspace=.055)
            ax = fig.add_subplot(sub[0])
            residual_ax = fig.add_subplot(sub[1], sharex=ax)
            ax.errorbar(p.r, p.value, yerr=[p.err_lo, p.err_hi],
                        xerr=[p.r-lower, p.r_upper-p.r], fmt="o", ms=3.8,
                        mfc="white", mec="#242424", mew=.8, ecolor="#969696", elinewidth=.75,
                        capsize=0, zorder=5)
            for category, style in STYLES.items():
                ax.plot(p.r, pair[category]["prediction"][name], **style)
                residual_ax.plot(p.r, pair[category]["residuals"][name], **style)
            ax.set_xscale("log")
            ax.set_ylim(max(0, ymin-padding), ymax+padding)
            ax.set_xlim(max(.1, p.r.min()*.75), p.r_upper.max()*1.05)
            ax.tick_params(labelbottom=False)
            residual_ax.axhspan(-1, 1, color="#b5b5b5", alpha=.2, zorder=0, lw=0)
            residual_ax.axhline(0, color="#737373", lw=.65, zorder=0)
            residual_ax.set_ylim(-limit, limit)
            half = max(1, int(limit//2))
            residual_ax.set_yticks([-half, 0, half])
            for a in (ax, residual_ax):
                a.grid(axis="y", alpha=.13)
                a.tick_params(which="both", direction="out", labelsize=9.2)
            if col == 0:
                ax.set_ylabel(unit, fontsize=11)
                ax.set_title(f"{title} ({p.n} bins)", fontsize=10.5, loc="left", pad=7)
                residual_ax.set_ylabel(r"$\Delta/\epsilon$", fontsize=11)
            else:
                ax.tick_params(labelleft=False)
                residual_ax.tick_params(labelleft=False)
            if row == 4:
                residual_ax.set_xlabel("Projected radius [arcsec]", fontsize=10.5)
            if name.startswith("gaia_"):
                residual_ax.xaxis.set_major_locator(FixedLocator([400, 600, 1000, 2000]))
            else:
                ticks = [v for v in [5, 10, 30, 100, 300] if ax.get_xlim()[0] <= v <= ax.get_xlim()[1]]
                residual_ax.xaxis.set_major_locator(FixedLocator(ticks))
            residual_ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
            residual_ax.xaxis.set_minor_formatter(NullFormatter())
        for i in range(p.n):
            record = dict(dataset=name, radius_arcsec=float(p.r[i]), bin_lower_arcsec=float(p.r_lower[i]),
                          bin_upper_arcsec=float(p.r_upper[i]), value=float(p.value[i]),
                          err_lo=float(p.err_lo[i]), err_hi=float(p.err_hi[i]))
            for pair in models:
                for model in pair.values():
                    record[model["label"]] = float(model["prediction"][name][i])
                    record[model["label"]+"_residual"] = float(model["residuals"][name][i])
            exports.append(record)

    fig.text(.065, .027, "Lines join exact likelihood-bin predictions, including streaming corrections. "
             "Residuals = (data - model) / asymmetric 1-sigma error; grey bands mark +/-1.", fontsize=10)
    fig.text(.065, .014, "Axis ranges are shared across each row. The cored halo has the lowest DM chi-square in every column. "
             "Flexible rung 2 uses a different central-anisotropy prior.", fontsize=10, color="#444444")
    fig.savefig(figure_path, dpi=220)
    plt.close(fig)
    table = Table(rows=exports)
    table.meta.update(data_fingerprint=fingerprint,
                      prediction="Exact likelihood-bin predictions including streaming correction; lines join bins",
                      selection="K1 no-DM and the lowest-chi-square of NFW/cored within each stage",
                      units="radius: arcsec; HST/Gaia: mas/yr; MUSE: km/s; residuals: dimensionless")
    table_path = data_prefix.with_suffix(".ecsv")
    table.write(table_path, overwrite=True)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), models=provenance,
                  selections=selections, inputs_sha256=hashes, data_fingerprint=fingerprint,
                  generator_sha256=digest(Path(__file__)),
                  outputs_sha256={str(p.relative_to(ROOT)): digest(p) for p in [figure_path, table_path]})
    data_prefix.with_suffix(".json").write_text(json.dumps(record, indent=2)+"\n")
    print(figure_path)
    for stage, pair in zip(STAGES, models):
        print(stage[0], {category: (model["label"], round(model["chi2"], 2)) for category, model in pair.items()})


if __name__ == "__main__":
    main()
