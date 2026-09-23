#!/usr/bin/env python3
"""rho20 and the mass decomposition as functions of the assumed stellar mass per counted star.

Reads results/df/mpsscan_<M>_<TAG> (HST count amplitude pinned at M Msun per F625W < 19
star, everything else refitted with a free cored halo) and the free-amplitude reference
batch, takes the best converged job of each, and plots (a) the objective and its
components relative to the overall minimum, (b) the recovered rho20 and profiled r_s,
(c) the profiled masses (M_star, M_rem, M_DM inside 20 and 63 pc). The photometric
mass-function range for M/N19 (JOURNAL 2026-09-23) is shaded. As for the rho20 profile,
Delta values are on a chi2 scale only if the likelihood is calibrated.

Usage: python bin/plot_mass_per_star_profile.py [--tag 20260923]
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
ARCSEC_PER_RAD = 206264.806
COUNTS = "ocen_counts_hst_f625w19"
# photometric mass function, Msun per F625W<19 star: zones 30-175 arcsec (Kroupa / fitted slope,
# IMF -2.3 above the turnoff, both metallicities) and the full spread over zones and variants
MF_CORE = (8.3, 9.8)
MF_FULL = (7.7, 11.3)


def read(path):
    return json.loads(Path(path).read_text())


def best_job(batch):
    rows = [(s["refined_score"], f, s) for f in glob.glob(f"{batch}/fits/*/summary.json")
            for s in [read(f)] if "best" in s]
    return min(rows, key=lambda r: r[0]) if rows else None


def extract(row, label):
    score, f, s = row
    g = s["gates"]; c = s["best"]["config"]; job = s["job"]
    pinned = [cc["path"].split(".")[-1] for cc, v in zip(job["coordinates"], s["best"]["x"]) if v < 1e-3 or v > 1-1e-3]
    hst = s["refined_counts"][COUNTS]
    pc_per_arcsec = c["distance_kpc"]*1000/ARCSEC_PER_RAD
    mps = hst.get("mass_per_star", pc_per_arcsec**2/hst["amplitude"])
    prof = s["refined_profiles"]
    return dict(label=label, mass_per_star=mps, amplitude_fixed=hst.get("amplitude_fixed", False),
                objective=s["refined_score"], chi2_kin=g["chi2_kinematic"],
                deviance_counts=sum(g.get("deviance_counts", {}).values()),
                deviance_hst=g.get("deviance_counts", {}).get(COUNTS), prior=((c["distance_kpc"]-5.43)/0.05)**2,
                rho20=c["matter"]["rho20"], r_s=c["matter"]["r_s"], M_star=c["M_star"], M_rem=c["matter"]["M_rem"],
                a_rem=c["matter"]["a_rem"], D=c["distance_kpc"],
                M_halo_20=float(np.interp(20., prof["r_pc"], prof["halo"])), M_halo_63=float(np.interp(63., prof["r_pc"], prof["halo"])),
                M_tot_20=float(np.interp(20., prof["r_pc"], prof["total"])), pinned=pinned,
                message=s["optimizer"]["message"], summary=f, n_jobs=len(glob.glob(f"{os.path.dirname(os.path.dirname(f))}/*/summary.json")))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", default="20260923")
    parser.add_argument("--reference", default="results/df/rho20free_20260923", help="free-amplitude, free-halo batch")
    args = parser.parse_args()
    points = []
    ref = best_job(ROOT/args.reference)
    if ref:
        points.append(extract(ref, "free amplitude"))
    for b in sorted(glob.glob(str(ROOT/f"results/df/mpsscan_*_{args.tag}"))):
        bj = best_job(b)
        if bj:
            points.append(extract(bj, f"pinned {Path(b).name.split('_')[1]}"))
    if not points:
        raise SystemExit("no finished batches")
    points.sort(key=lambda p: p["mass_per_star"])
    x = np.array([p["mass_per_star"] for p in points]); obj = np.array([p["objective"] for p in points])
    free = [p for p in points if not p["amplitude_fixed"]]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for ax in axes:
        ax.axvspan(*MF_FULL, color="0.93", zorder=0); ax.axvspan(*MF_CORE, color="0.85", zorder=0)
        for p in free:
            ax.axvline(p["mass_per_star"], color="tab:red", ls="--", lw=1)
        ax.set_xlabel(r"assumed $M_\star$ per F625W$<$19 star [$M_\odot$]")
    ax = axes[0]
    for key, lab in (("objective", "total"), ("chi2_kin", "kinematics"), ("deviance_counts", "star counts"), ("prior", "distance prior")):
        v = np.array([p[key] for p in points]); ax.plot(x, v-v[np.argmin(obj)], marker="o", label=f"{lab}: value - value at best")
    for y, txt in ((1., r"$\Delta=1$"), (3.84, r"$\Delta=3.84$")):
        ax.axhline(y, color="0.6", ls=":", lw=.8); ax.text(x.max(), y, txt, fontsize=7, va="bottom", ha="right")
    ax.set(ylabel=r"$\Delta$ objective", title="(a) objective vs assumed mass per star (grey: photometric MF range)")
    ax.legend(fontsize=7)
    ax = axes[1]
    ax.plot(x, [p["rho20"] for p in points], marker="o", color="tab:blue", label=r"$\rho_{20}$ [$M_\odot$ pc$^{-3}$]")
    ax.set(ylabel=r"recovered $\rho_{20}$", title="(b) halo density at 20 pc and profiled $r_s$")
    ax2 = ax.twinx()
    ax2.plot(x, [p["r_s"] for p in points], marker="s", color="tab:green", ls="--", label=r"$r_s$ [pc]")
    ax2.set(yscale="log", ylabel=r"profiled $r_s$ [pc] (bounds 5-500)")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1+h2, l1+l2, fontsize=7)
    ax = axes[2]
    for key, mk, lab in (("M_star", "^", r"$M_\star$"), ("M_rem", "v", r"$M_{\rm rem}$"), ("M_halo_20", "o", r"$M_{\rm DM}(<20\,{\rm pc})$"),
                         ("M_halo_63", "s", r"$M_{\rm DM}(<63\,{\rm pc})$"), ("M_tot_20", "x", r"$M_{\rm tot}(<20\,{\rm pc})$")):
        ax.plot(x, [max(p[key], 1.) for p in points], marker=mk, label=lab)
    ax.set(yscale="log", ylabel=r"mass [$M_\odot$]", title="(c) profiled masses", ylim=(1e3, 1e7))
    ax.legend(fontsize=7)
    fig.suptitle("Dark matter vs the assumed stellar mass per counted star (HST count amplitude pinned; two-transition DF, "
                 "Poisson counts, distance prior; red dashed: free amplitude)")
    plot = ROOT/"plots"/f"mass_per_star_profile_{args.tag}.png"
    fig.savefig(plot, dpi=150)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), points=points, mf_core=MF_CORE, mf_full=MF_FULL,
                  plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/f"mass_per_star_profile_{args.tag}.json").write_text(json.dumps(record))
    base = obj.min()
    for p in points:
        print(f"M/N19 {p['mass_per_star']:5.2f} ({p['label']:14s}) obj {p['objective']:7.1f} (+{p['objective']-base:5.2f})  chi2_kin {p['chi2_kin']:6.1f}  "
              f"dev {p['deviance_counts']:6.1f} (HST {p['deviance_hst']:.1f})  prior {p['prior']:.2f}  rho20 {p['rho20']:.3f}  r_s {p['r_s']:6.1f}  "
              f"M* {p['M_star']:.3e}  Mrem {p['M_rem']:.2e}@{p['a_rem']:.2f}  M_DM(<20) {p['M_halo_20']:.2e}  D {p['D']:.3f}  "
              f"{'BOUND ' + ','.join(p['pinned']) if p['pinned'] else ''} ({p['message'][:20]})")
    print(plot)


if __name__ == "__main__":
    main()
