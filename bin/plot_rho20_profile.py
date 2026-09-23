#!/usr/bin/env python3
"""Profile likelihood in rho20 from the fixed-rho20 scan batches.

Reads results/df/rho20free_<TAG> (free halo), results/df/rho20scan_<RHO>_<TAG>
(fixed rho20) and the no-halo reference batch, takes the best converged job of
each, and plots (a) the refined objective minus its minimum against rho20, with
the kinematic, count and prior components; (b) the profiled r_s and the halo
mass inside 20 and 63 pc. The objective is chi2_kin + Poisson deviances +
prior, so Delta values are on a chi2 scale only if the likelihood is
calibrated; the plot marks Delta = 1 and 3.84 for reference, not as intervals.

Usage: python bin/plot_rho20_profile.py [--tag 20260923]
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import glob
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text())


def best_job(batch):
    rows = []
    for f in glob.glob(f"{batch}/fits/*/summary.json"):
        s = read(f)
        if "best" in s:
            rows.append((s["refined_score"], f, s))
    if not rows:
        return None
    rows.sort(key=lambda r: r[0])
    return rows[0]


def halo_mass(profiles, r):
    return float(np.interp(r, profiles["r_pc"], profiles["halo"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", default="20260923")
    parser.add_argument("--reference", default="results/df/dftwo_counts_20260923", help="no-halo batch (rho20 = 0)")
    args = parser.parse_args()
    points = []
    ref = best_job(ROOT/args.reference)
    if ref:
        points.append(dict(rho20=0., label="no halo", **_extract(ref)))
    free = best_job(ROOT/f"results/df/rho20free_{args.tag}")
    if free:
        points.append(dict(rho20=free[2]["best"]["config"]["matter"]["rho20"], label="free halo", **_extract(free)))
    for b in sorted(glob.glob(str(ROOT/f"results/df/rho20scan_*_{args.tag}"))):
        rho = float(Path(b).name.split("_")[1])
        bj = best_job(b)
        if bj:
            points.append(dict(rho20=rho, label=f"fixed {rho:g}", **_extract(bj)))
    points.sort(key=lambda p: p["rho20"])
    if not points:
        raise SystemExit("no finished scan batches")
    obj = np.array([p["objective"] for p in points]); rho = np.array([p["rho20"] for p in points])
    base = obj.min()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    ax = axes[0]
    for key, lab in (("objective", "total"), ("chi2_kin", "kinematics"), ("deviance_counts", "star counts"), ("prior", "distance prior")):
        v = np.array([p[key] for p in points]); v0 = v[np.argmin(obj)]
        ax.plot(rho, v-v0, marker="o", label=f"{lab}: value - value at minimum")
    for y, txt in ((1., r"$\Delta=1$"), (3.84, r"$\Delta=3.84$")):
        ax.axhline(y, color="0.6", ls=":", lw=.8); ax.text(rho.max(), y, txt, fontsize=7, va="bottom", ha="right")
    ax.set(xlabel=r"$\rho_{20}$ [$M_\odot$ pc$^{-3}$]", ylabel=r"$\Delta$ objective", title="(a) profile in rho20 (all else refitted)")
    ax.legend(fontsize=7)
    ax = axes[1]
    ax.plot(rho, [p["r_s"] for p in points], marker="o")
    ax.set(xlabel=r"$\rho_{20}$", ylabel=r"profiled $r_s$ [pc]", yscale="log", title="(b) halo scale radius (bounds 5-500 pc)")
    ax.axhline(5, color="0.6", ls=":", lw=.8); ax.axhline(500, color="0.6", ls=":", lw=.8)
    ax = axes[2]
    ax.plot(rho, [p["M_halo_20"] for p in points], marker="o", label=r"$M_{\rm DM}(<20\,{\rm pc})$")
    ax.plot(rho, [p["M_halo_63"] for p in points], marker="s", label=r"$M_{\rm DM}(<63\,{\rm pc})$")
    ax.plot(rho, [p["M_star"] for p in points], marker="^", label=r"$M_\star$")
    ax.plot(rho, [p["M_rem"] for p in points], marker="v", label=r"$M_{\rm rem}$")
    ax.set(xlabel=r"$\rho_{20}$", ylabel=r"mass [$M_\odot$]", yscale="log", title="(c) profiled masses")
    ax.legend(fontsize=7)
    fig.suptitle("rho20 profile likelihood: two-transition DF, Poisson counts, distance prior (published errors; Delta lines are reference only)")
    plot = ROOT/"plots"/f"rho20_profile_{args.tag}.png"
    fig.savefig(plot, dpi=150)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), points=points, plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/f"rho20_profile_{args.tag}.json").write_text(json.dumps(record))
    for p in points:
        print(f"rho20 {p['rho20']:5.2f} ({p['label']:10s}) obj {p['objective']:7.1f} (+{p['objective']-base:5.2f})  chi2_kin {p['chi2_kin']:6.1f}  dev {p['deviance_counts']:6.1f}  prior {p['prior']:.2f}  r_s {p['r_s']:6.1f}  M_DM(<20) {p['M_halo_20']:.2e}  M_DM(<63) {p['M_halo_63']:.2e}  M* {p['M_star']:.3e}  D {p['D']:.3f}  {'BOUND ' + ','.join(p['pinned']) if p['pinned'] else ''}  ({p['message'][:20]})")
    print(plot)


def _extract(row):
    score, f, s = row
    g = s["gates"]; c = s["best"]["config"]; job = s["job"]
    pinned = [cc["path"].split(".")[-1] for cc, v in zip(job["coordinates"], s["best"]["x"]) if v < 1e-3 or v > 1-1e-3]
    dev = sum(g.get("deviance_counts", {}).values())
    prior = ((c["distance_kpc"]-5.43)/0.05)**2
    return dict(objective=s["refined_score"], chi2_kin=g["chi2_kinematic"], deviance_counts=dev, prior=prior,
                r_s=c["matter"]["r_s"], M_star=c["M_star"], M_rem=c["matter"]["M_rem"], D=c["distance_kpc"],
                M_halo_20=halo_mass(s["refined_profiles"], 20.), M_halo_63=halo_mass(s["refined_profiles"], 63.),
                pinned=pinned, message=s["optimizer"]["message"], summary=f)


if __name__ == "__main__":
    main()
