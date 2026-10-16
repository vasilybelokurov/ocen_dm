#!/usr/bin/env python3
"""Is the central deficit of the HST F625W<19 counts crowding incompleteness, or Poisson noise?

A (data only). oMEGACat (Haeberle et al. 2025) radial number-density profiles in F625W slices, on a finer
radial grid in the core. Crowding removes faint stars first, so incompleteness of the bright slices would
show as a central DROP of faint/bright ratios; a flat ratio means the slices lose stars in the same proportion
(or none). Same centre and pixel scale as bin/build_count_profiles.py.
B (models). Observed counts in the inner HST bins against each fitted model's expected counts (mu, with its
profiled amplitude), with Poisson tail probabilities P(N <= n_obs | mu).

Outputs: plots/hst_counts_completeness_20261015.png, results/plot_data/hst_counts_completeness_20261015.json.
Usage: python bin/check_hst_count_completeness.py BATCH::JOB::LABEL [...]
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.stats import poisson

from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.counts import CountProfile
from ocen_dm.kinematics.df_fit import DFJointProblem, build_df_model, model_config_from_dict
from ocen_dm.kinematics.run_io import read_data_snapshot

CATALOGUE = ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits"
CENTRE, SCALE = 15000., 0.04          # pixel, arcsec/px (as in build_count_profiles.py)
SLICES = [(None, 17.), (17., 18.), (18., 19.), (19., 20.), (20., 21.)]
NAME = "ocen_counts_hst_f625w19"


def data_profiles():
    with fits.open(CATALOGUE, memmap=True) as h:
        t = h[1].data
        x, y, mag = (np.asarray(t[k], float) for k in ("x", "y", "f625w"))
    r = np.hypot(x-CENTRE, y-CENTRE)*SCALE
    edges = np.geomspace(1., 120., 17)          # arcsec; inside the full-coverage region (30") and beyond
    area = np.pi*np.diff(edges**2)
    out = {}
    for lo, hi in SLICES:
        sel = np.isfinite(mag) & (mag < hi) & ((mag >= lo) if lo is not None else True)
        n, _ = np.histogram(r[sel], edges)
        out[f"{'' if lo is None else lo}-{hi}"] = n
    return edges, area, out


def model_counts(specs):
    rows = {}
    for spec in specs:
        batch, job, label = spec.split("::")
        batch = Path(batch)
        s = json.loads((batch/"fits"/job/"summary.json").read_text())
        cfg = refined_config(model_config_from_dict(s["best"]["config"]))
        od = batch/"observed"
        rec = json.loads((od/"problem.json").read_text())
        counts = [CountProfile.from_dict(json.loads((od/f"{n}.json").read_text())) for n in rec["counts"]]
        prob = DFJointProblem(read_data_snapshot(od, rec["data_snapshot"]), None, counts)
        ev = prob.evaluate(cfg, build_df_model(cfg))
        c = [p for p in counts if p.name == NAME][0]
        rows[label] = dict(mu=np.asarray(ev["counts"][NAME]["mu"], float), obs=np.asarray(c.counts, float),
                           r_lower=np.asarray(c.r_lower), r_upper=np.asarray(c.r_upper))
        print(label, "built", flush=True)
    return rows


def main():
    edges, area, prof = data_profiles()
    rmid = np.sqrt(edges[:-1]*edges[1:])
    rows = model_counts(sys.argv[1:])
    fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
    outer = (rmid > 30) & (rmid < 120)
    for k, n in prof.items():
        dens = n/area
        norm = np.sum(n[outer])/np.sum(area[outer])
        err = np.sqrt(np.maximum(n, 1))/area/norm
        ax[0].errorbar(rmid, dens/norm, err, marker="o", ms=3, capsize=2, label=f"F625W {k}")
    ax[0].set(xscale="log", yscale="log", xlabel="R [arcsec]", ylabel="density / mean density at 30-120\"",
              title="oMEGACat number density by magnitude slice")
    ax[0].legend(fontsize=8)
    ref = prof["18.0-19.0"]
    for k in ("-17.0", "17.0-18.0", "19.0-20.0", "20.0-21.0"):
        n = prof[k]
        ratio = n/np.maximum(ref, 1)
        err = ratio*np.sqrt(1/np.maximum(n, 1)+1/np.maximum(ref, 1))
        ax[1].errorbar(rmid, ratio/np.mean(ratio[outer]), err/np.mean(ratio[outer]), marker="o", ms=3, capsize=2,
                       label=f"N({k}) / N(18-19)")
    ax[1].axhline(1, color="k", lw=.8)
    ax[1].set(xscale="log", xlabel="R [arcsec]", ylabel="ratio / its mean at 30-120\"", ylim=(0, 2),
              title="faint/bright ratios (crowding would drop them inwards)")
    ax[1].legend(fontsize=8)
    table = {}
    for i, (label, row) in enumerate(rows.items()):
        rc = np.sqrt(row["r_lower"]*row["r_upper"])
        p_low = poisson.cdf(row["obs"], row["mu"])
        table[label] = dict(r_arcsec=rc[:8].tolist(), obs=row["obs"][:8].tolist(), mu=row["mu"][:8].tolist(),
                            p_le_obs=p_low[:8].tolist())
        ax[2].plot(rc, row["mu"]/(row["r_upper"]**2-row["r_lower"]**2)/np.pi, color=f"C{i}", label=label)
    r0 = list(rows.values())[0]
    rc = np.sqrt(r0["r_lower"]*r0["r_upper"])
    a = np.pi*(r0["r_upper"]**2-r0["r_lower"]**2)
    ax[2].errorbar(rc, r0["obs"]/a, np.sqrt(np.maximum(r0["obs"], 1))/a, fmt="ko", ms=3, capsize=2,
                   label="HST F625W<19 (Poisson)")
    ax[2].set(xscale="log", yscale="log", xlabel="R [arcsec]", ylabel="stars per arcsec$^2$",
              title="inner HST counts vs models (expected mu)")
    ax[2].legend(fontsize=8)
    png = ROOT/"plots/hst_counts_completeness_20261015.png"
    fig.savefig(png, dpi=140)
    out = dict(slices={k: v.tolist() for k, v in prof.items()}, edges_arcsec=edges.tolist(), models=table)
    (ROOT/"results/plot_data/hst_counts_completeness_20261015.json").write_text(json.dumps(out, indent=1))
    for label, t in table.items():
        print(label)
        for r, o, m, p in zip(t["r_arcsec"], t["obs"], t["mu"], t["p_le_obs"]):
            print(f"   R {r:6.2f}\"  obs {o:5.0f}  mu {m:7.1f}  P(N<=obs) {p:.3g}")
    print(png)


if __name__ == "__main__":
    main()
