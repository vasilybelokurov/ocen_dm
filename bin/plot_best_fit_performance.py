#!/usr/bin/env python3
"""Performance report of selected compact-DF fits to the observed Omega Cen data.

Figure 1 (<name>_data.png): each dataset with the models and normalized
residuals; the rung-2 Jeans turnover fit (no DM) is shown dashed for the
kinematics (it was fitted to the same 89 bins; its photometry is a separate MGE).
Figure 2 (<name>_diagnostics.png): residual histogram vs N(0,1); chi2/n per
dataset with the chi2 survival probability; Wald-Wolfowitz runs test on the
residual signs per dataset (z < -2 flags coherent structure); intrinsic beta(r);
enclosed-mass components. Also writes a JSON table of all numbers.

Residuals use the adopted (raw published) errors; bins are treated as
independent, so p-values are indicative, not calibrated.

Usage: python bin/plot_best_fit_performance.py --name best_fit_performance_20260923 \
    "results/df/dftwo_20260923::observed_no_halo_start0::two transitions, no halo" \
    "results/df/dftwo_halo_20260923::observed_free_halo_start0::two transitions + halo"
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scipy import stats

from ocen_dm.kinematics.compact_recovery import mass_profiles, refined_config, residual_vector
from ocen_dm.kinematics.counts import CountProfile
from ocen_dm.kinematics.df_fit import DFJointProblem, PhotometricData, build_df_model, model_config_from_dict
from ocen_dm.kinematics.run_io import read_data_snapshot

KINDS = [("phot", r"$\mu$ [mag]", "surface brightness"), ("los", r"$\sigma_{\rm los}$ [km/s]", "line of sight"),
         ("pmr", r"$\sigma_{\rm pm,R}$ [mas/yr]", "PM radial"), ("pmt", r"$\sigma_{\rm pm,T}$ [mas/yr]", "PM tangential")]
JEANS = "rung2_K1_turnover"


def read(path):
    return json.loads(Path(path).read_text())


def runs_z(residual):
    """Wald-Wolfowitz runs test on signs; strongly negative z = too few runs = coherent residuals."""
    s = np.sign(residual)
    s = s[s != 0]
    n1, n2 = np.sum(s > 0), np.sum(s < 0)
    if min(n1, n2) == 0:
        return float("-inf")
    runs = 1+np.sum(s[1:] != s[:-1])
    mu = 2*n1*n2/(n1+n2)+1
    var = 2*n1*n2*(2*n1*n2-n1-n2)/((n1+n2)**2*(n1+n2-1))
    return float((runs-mu)/np.sqrt(var))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--name", required=True)
    parser.add_argument("specs", nargs="+", help="BATCH::JOB::LABEL")
    args = parser.parse_args()
    entries = [tuple(e.split("::")) for e in args.specs]
    d = Path(entries[0][0])/"observed"
    record0 = read(d/"problem.json")
    photometry = (PhotometricData.from_dict(read(d/"photometry.json"))
                  if (d/"photometry.json").exists() and record0.get("photometry", True) else None)
    counts = [CountProfile.from_dict(read(d/f"{n}.json")) for n in record0.get("counts", [])]
    problem = DFJointProblem(read_data_snapshot(d, record0["data_snapshot"]), photometry, counts)
    pc = 5.43e3*np.pi/(180*3600)
    jeans = Table.read(ROOT/"results/plot_data/rung0_rung1_rung2_no_dm_vs_dm.ecsv")
    names = [p.name for p in problem.data.profiles]+(["photometry"] if problem.photometry is not None else [])+[c.name for c in problem.counts]
    radii = np.geomspace(.5, 80., 64)
    fits = []
    for batch, job, label in entries:
        s = read(Path(batch)/"fits"/job/"summary.json")
        config = refined_config(model_config_from_dict(s["best"]["config"]))
        model = build_df_model(config)
        ev = problem.evaluate(config, model)
        res = -residual_vector(problem, ev)  # residual_vector is (model-data)/err; plot (data-model)/err
        split, i = {}, 0
        for p in problem.data.profiles:
            split[p.name] = res[i:i+p.n]; i += p.n
        if problem.photometry is not None:
            # Radius order (the photometric splice is interleaved) so the runs test sees the profile.
            npho = len(problem.photometry.mu)
            split["photometry"] = res[i:i+npho][np.argsort(problem.photometry.r_arcsec)]; i += npho
        for c in problem.counts:
            split[c.name] = -res[i:i+c.n]; i += c.n   # deviance residuals are already (data-model)-signed
        fits.append(dict(label=label, config=s["best"]["config"], ev=ev, res=res, split=split,
                         beta=model.intrinsic_moments(np.geomspace(.05, 80., 200))["beta"],
                         mass=mass_profiles(model, radii)))
        print(label, "chi2_kin %.1f density term %.1f" % (ev["chi2_kinematic"], (ev["photometry"] or {}).get("chi2", 0.)+ev.get("deviance_counts", 0.)), flush=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # ---------------- figure 1: data, models, residuals
    fig, axes = plt.subplots(2, 4, figsize=(16, 6.8), sharex="col", constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[3, 1.5]))
    for col, (kind, ylabel, title) in enumerate(KINDS):
        ax, rx = axes[0, col], axes[1, col]
        if kind == "phot" and problem.photometry is not None:
            ph = problem.photometry
            ax.errorbar(ph.r_arcsec*pc, ph.mu, yerr=ph.sigma_mag, fmt=".", color="0.35", ms=3, lw=.6, label="observed")
            order = np.argsort(ph.r_arcsec)
            for i, f in enumerate(fits):
                ax.plot(ph.r_arcsec[order]*pc, f["ev"]["photometry"]["prediction"][order], color=colors[i], lw=1.2,
                        label=f["label"])
                rx.plot(np.sort(ph.r_arcsec)*pc, f["split"]["photometry"], ".", color=colors[i], ms=3)
            ax.invert_yaxis()
        elif kind == "phot":
            for c, mk in zip(problem.counts, ("o", "s")):
                ax.errorbar(c.r_median*pc, c.counts/c.area_arcsec2*3600., np.sqrt(np.maximum(c.counts, 1))/c.area_arcsec2*3600.,
                            fmt=mk, color="0.35", ms=3, mfc="white", lw=.6, label=c.name.split("_")[2]+" counts")
                for i, f in enumerate(fits):
                    out = f["ev"]["counts"][c.name]
                    ax.plot(c.r_median*pc, out["mu"]/c.area_arcsec2*3600., color=colors[i], lw=1.2,
                            label=f["label"] if c is problem.counts[0] else None)
                    rx.plot(c.r_median*pc, f["split"][c.name], ".-", color=colors[i], ms=3, lw=.6)
            ax.set_yscale("log")
        profiles = [p for p in problem.data.profiles if p.kind == kind]
        for j, p in enumerate(profiles):
            ax.errorbar(p.r*pc, p.value, yerr=[p.err_lo, p.err_hi], fmt=".", color="0.35", ms=3, lw=.6,
                        label="observed" if j == 0 else None)
            for i, f in enumerate(fits):
                ax.plot(p.r*pc, f["ev"]["predictions"][p.name], color=colors[i], lw=1.2,
                        label=f["label"] if j == 0 else None)
                rx.plot(p.r*pc, f["split"][p.name], ".-", color=colors[i], ms=3, lw=.6)
            rows = jeans[jeans["dataset"] == p.name]
            ax.plot(np.asarray(rows["radius_arcsec"])*pc, rows[JEANS], color="0.2", ls="--", lw=1,
                    label="Jeans turnover + BH (no DM)" if j == 0 else None)
            rx.plot(np.asarray(rows["radius_arcsec"])*pc, rows[JEANS+"_residual"], color="0.2", ls="--", lw=.8)
        for y in (-2, 0, 2):
            rx.axhline(y, color="black" if y == 0 else "0.7", lw=.7)
        ax.set(xscale="log", ylabel=ylabel if not (kind == "phot" and problem.photometry is None) else r"stars per arcmin$^2$", title=title)
        rx.set(xscale="log", xlabel="projected radius [pc]",
               ylabel="deviance residual" if (kind == "phot" and problem.photometry is None) else "(data-model)/err")
        ax.legend(fontsize=6)
    fig.suptitle("Best-fit compact DF models vs observed Omega Cen data (published errors)")
    plot1 = ROOT/"plots"/(args.name+"_data.png")
    fig.savefig(plot1, dpi=160)

    # ---------------- figure 2: diagnostics
    table = {}
    fig, axes = plt.subplots(1, 5, figsize=(20, 4.2), constrained_layout=True)
    ax = axes[0]
    bins = np.linspace(-5, 5, 26)
    for i, f in enumerate(fits):
        kin = f["res"][:problem.data.n_points]
        ax.hist(kin, bins=bins, density=True, histtype="step", lw=2, color=colors[i], label=f["label"])
    x = np.linspace(-5, 5, 200)
    ax.plot(x, stats.norm.pdf(x), color="black", ls=":", label="N(0,1)")
    ax.set(xlabel="kinematic residual / err", ylabel="density", title="(a) kinematic residuals (89 bins)")
    ax.legend(fontsize=6)
    width = .8/len(fits)
    for panel, key, title in ((axes[1], "chi2n", r"(b) $\chi^2/n$ per dataset"),
                              (axes[2], "runs", "(c) runs test z (coherence)")):
        for i, f in enumerate(fits):
            vals = []
            for k, n in enumerate(names):
                r = f["split"][n]
                c2 = float(r @ r)
                entry = table.setdefault(f["label"], {}).setdefault(n, {})
                entry.update(n=int(r.size), chi2=c2, chi2_per_n=c2/r.size, p_value=float(stats.chi2.sf(c2, r.size)),
                             runs_z=runs_z(r), max_abs=float(np.max(abs(r))))
                vals.append(entry["chi2_per_n"] if key == "chi2n" else entry["runs_z"])
            panel.bar(np.arange(len(names))+i*width, vals, width, color=colors[i], label=f["label"])
        panel.set_xticks(np.arange(len(names))+width*(len(fits)-1)/2)
        panel.set_xticklabels([n.replace("_ours", "").replace("_dispersion", "").replace("gaia_edr3", "gaia")
                               .replace("hst_pm_", "hst ") for n in names], rotation=35, ha="right", fontsize=7)
        panel.set_title(title)
    axes[1].axhline(1, color="black", lw=.7); axes[1].axhline(2, color="0.6", ls="--", lw=.7)
    axes[2].axhline(0, color="black", lw=.7); axes[2].axhline(-2, color="0.6", ls="--", lw=.7)
    axes[1].legend(fontsize=6)
    ax = axes[3]
    rr = np.geomspace(.05, 80., 200)
    for i, f in enumerate(fits):
        ax.plot(rr, f["beta"], color=colors[i], lw=2, label=f["label"])
    ax.axvspan(63, 80, color="0.9", label="beyond data")
    ax.axhline(0, color="black", lw=.7)
    ax.set(xscale="log", xlabel="r [pc]", ylabel=r"$\beta(r)$", title="(d) intrinsic anisotropy")
    ax.legend(fontsize=6)
    ax = axes[4]
    for i, f in enumerate(fits):
        m = f["mass"]
        ax.plot(m["r_pc"], m["total"], color=colors[i], lw=2, label=f["label"]+" total")
        ax.plot(m["r_pc"], m["stars"], color=colors[i], ls="--", lw=1)
        ax.plot(m["r_pc"], m["remnants"], color=colors[i], ls="-.", lw=1)
        if max(m["halo"]) > 0:
            ax.plot(m["r_pc"], m["halo"], color=colors[i], ls=":", lw=1.5)
    ax.set(xscale="log", yscale="log", ylim=(1e3, 6e6), xlabel="r [pc]", ylabel=r"$M(<r)$ [$M_\odot$]",
           title="(e) enclosed mass: total, stars --, remnants -., halo :")
    ax.legend(fontsize=6)
    plot2 = ROOT/"plots"/(args.name+"_diagnostics.png")
    fig.savefig(plot2, dpi=160)

    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), fits=[e for e in args.specs], table=table,
                  totals={f["label"]: dict(chi2_kin=f["ev"]["chi2_kinematic"],
                                           chi2_phot=(f["ev"]["photometry"] or {}).get("chi2"),
                                           deviance_counts=f["ev"].get("deviance_counts"), config=f["config"]) for f in fits},
                  plots={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (plot1, plot2)})
    (ROOT/"results/plot_data"/(args.name+".json")).write_text(json.dumps(record))
    for label, per in table.items():
        print(label)
        for n, e in per.items():
            print(f"   {n:28s} chi2/n {e['chi2_per_n']:5.2f}  p {e['p_value']:.2g}  runs z {e['runs_z']:+.1f}  max|r| {e['max_abs']:.1f}")
    print(plot1, plot2)


if __name__ == "__main__":
    main()
