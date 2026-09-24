#!/usr/bin/env python3
"""Who carries the mass inside ~10 pc: stars, remnants or a dark halo (data and mocks).

 (a) Observed data, rho20 profile scan: enclosed mass of each component inside
     R_IN pc (stars, remnants, halo) and their sum, against the fixed rho20; the
     Delta objective of each point is annotated. The sum is what the kinematics
     fix; along the scan the stars give way to halo AND remnants together.
 (b) The same fits as enclosed-mass profiles M(<r) per component (no halo,
     rho20 = 1, rho20 = 4): the totals coincide over the data range.
 (c) Every realistic-mock fit: error of the stellar mass inside R_IN against the
     error of the dark mass (remnants + halo) inside R_IN. The anticorrelation
     is the leading degeneracy (stars <-> dark), driven by the count amplitude.
 (d) Halo error against remnant error inside R_IN for the fits with the stellar
     mass pinned at the truth: the residual degeneracy once the stars are fixed.

Usage: python bin/plot_remnant_halo_degeneracy.py [--tag 20260923]
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import glob
import hashlib
import json
import os
from pathlib import Path
import re

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
R_IN = 10.          # pc; the remnants (Plummer a ~ 2 pc) and the compact halo (r_s at 5 pc) both live inside
COUNTS = "ocen_counts_hst_f625w19"


def read(path):
    return json.loads(Path(path).read_text())


def best_job(batch):
    rows = [(s["refined_score"], s) for f in glob.glob(f"{batch}/fits/*/summary.json") for s in [read(f)] if "best" in s]
    return min(rows, key=lambda r: r[0])[1] if rows else None


def enclosed(prof, key, r):
    return float(np.interp(r, prof["r_pc"], prof[key]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", default="20260923")
    args = parser.parse_args()
    # ---- data scan
    scan = []
    s = best_job(ROOT/"results/df/dftwo_counts_20260923")
    scan.append(dict(rho20=0., s=s, label="no halo"))
    for b in sorted(glob.glob(str(ROOT/f"results/df/rho20scan_*_{args.tag}"))):
        s = best_job(b)
        if s:
            scan.append(dict(rho20=float(Path(b).name.split("_")[1]), s=s, label=f"rho20 = {Path(b).name.split('_')[1]}"))
    scan.sort(key=lambda p: p["rho20"])
    base = min(p["s"]["refined_score"] for p in scan)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    ax = axes[0, 0]
    rho = np.array([p["rho20"] for p in scan])
    for key, mk, lab in (("stars", "^", "stars (follow the light)"), ("remnants", "v", "remnants (Plummer)"), ("halo", "o", "DM halo"), ("total", "s", "total")):
        ax.plot(rho, [enclosed(p["s"]["refined_profiles"], key, R_IN) for p in scan], marker=mk, lw=1.5, label=lab)
    for p in scan:
        ax.annotate(f"+{p['s']['refined_score']-base:.1f}", (p["rho20"], enclosed(p["s"]["refined_profiles"], "total", R_IN)),
                    textcoords="offset points", xytext=(0, 6), ha="center", fontsize=7)
    ax.set(xlabel=r"fixed $\rho_{20}$ [$M_\odot$ pc$^{-3}$]", ylabel=fr"$M(<{R_IN:g}\,{{\rm pc}})$ [$M_\odot$]", yscale="log", ylim=(1e4, 3e6),
           title=f"(a) data: mass inside {R_IN:g} pc along the rho20 scan (labels: $\\Delta$ objective)")
    ax.legend(fontsize=8, loc="center right")
    # ---- enclosed-mass profiles
    ax = axes[0, 1]
    styles = {0.: "-", 1.: "--", 4.: ":"}
    colors = dict(stars="tab:blue", remnants="tab:red", halo="tab:green", total="black")
    for p in scan:
        if p["rho20"] not in styles:
            continue
        prof = p["s"]["refined_profiles"]; r = np.array(prof["r_pc"])
        for key in ("stars", "remnants", "halo", "total"):
            m = np.array(prof[key])
            if key == "halo" and m.max() <= 0:
                continue
            ax.plot(r, np.maximum(m, 1.), ls=styles[p["rho20"]], color=colors[key], lw=1.8 if key == "total" else 1.2,
                    label=f"{key}, {p['label']}" if key in ("total", "halo", "remnants") else (key if p["rho20"] == 0 else None))
    ax.axvspan(63, 80, color="0.92"); ax.axvspan(.5, 2., color="0.96")
    ax.set(xscale="log", yscale="log", xlim=(.5, 80.), ylim=(1e3, 6e6), xlabel="r [pc]", ylabel=r"$M(<r)$ [$M_\odot$]",
           title="(b) data: enclosed mass by component\n(solid: no halo; dashed: rho20 = 1; dotted: rho20 = 4; grey: outside data)")
    ax.legend(fontsize=7, ncol=2)
    # ---- mocks
    ax = axes[1, 0]
    sets = {"cov": ("coverage test (amplitude profiled)", "0.4", "x"), "truth": ("M/N19 pinned at truth", "tab:green", "o"),
            "9.0": ("M/N19 pinned at 9.0 (-16%)", "tab:red", "s"), "12.0": ("M/N19 pinned at 12.0 (+12%)", "tab:purple", "D")}
    rows = []
    for b in sorted(glob.glob(str(ROOT/f"results/df/cov_*_{args.tag}"))+glob.glob(str(ROOT/f"results/df/mpsmock_*_{args.tag}"))):
        s = best_job(b)
        if not s or not os.path.exists(f"{b}/mocks/mock/mock.json"):
            continue
        t = read(f"{b}/mocks/mock/mock.json"); c = s["best"]["config"]; tc = t["config"]
        name = os.path.basename(b); m = re.match(r"mpsmock_([^_]+)_", name)
        key = m.group(1) if m else "cov"
        rows.append(dict(name=name, set=key, truth_rho20=tc["matter"]["rho20"], rho20=c["matter"]["rho20"],
                         mrem_err=c["matter"]["M_rem"]/tc["matter"]["M_rem"]-1, mstar_err=c["M_star"]/tc["M_star"]-1,
                         halo_in=enclosed(s["refined_profiles"], "halo", R_IN)-enclosed(t["profiles"], "halo", R_IN),
                         rem_in=enclosed(s["refined_profiles"], "remnants", R_IN)-enclosed(t["profiles"], "remnants", R_IN),
                         star_in=enclosed(s["refined_profiles"], "stars", R_IN)-enclosed(t["profiles"], "stars", R_IN)))
    for key, (lab, color, mk) in sets.items():
        rr = [r for r in rows if r["set"] == key]
        if rr:
            ax.scatter([r["star_in"]/1e5 for r in rr], [(r["halo_in"]+r["rem_in"])/1e5 for r in rr], marker=mk, color=color, s=50,
                       edgecolor="black" if mk != "x" else None, label=lab, zorder=3)
    lim = max(abs(v) for r in rows for v in (r["star_in"], r["halo_in"]+r["rem_in"]))/1e5*1.1
    ax.plot([-lim, lim], [lim, -lim], color="0.6", ls="--", lw=1, label="total conserved")
    ax.axhline(0, color="black", lw=.6); ax.axvline(0, color="black", lw=.6)
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel=fr"recovered $-$ true $M_\star(<{R_IN:g}\,{{\rm pc}})$ [$10^5\,M_\odot$]",
           ylabel=fr"recovered $-$ true $M_{{\rm rem}}+M_{{\rm DM}}(<{R_IN:g}\,{{\rm pc}})$ [$10^5\,M_\odot$]",
           title="(c) mocks: stars trade against dark mass (remnants + halo)")
    ax.legend(fontsize=7, loc="upper right")
    ax = axes[1, 1]
    for key in ("cov", "truth"):
        lab, color, mk = sets[key]
        rr = [r for r in rows if r["set"] == key]
        ax.scatter([r["halo_in"]/1e5 for r in rr], [r["rem_in"]/1e5 for r in rr], marker=mk, color=color, s=50,
                   edgecolor="black" if mk != "x" else None, label=lab, zorder=3)
    lim2 = max(abs(v) for r in rows if r["set"] in ("cov", "truth") for v in (r["halo_in"], r["rem_in"]))/1e5*1.15
    ax.plot([-lim2, lim2], [lim2, -lim2], color="0.6", ls="--", lw=1, label="remnants + halo conserved")
    ax.axhline(0, color="black", lw=.6); ax.axvline(0, color="black", lw=.6)
    ax.set(xlim=(-lim2, lim2), ylim=(-lim2, lim2), xlabel=fr"recovered $-$ true $M_{{\rm DM}}(<{R_IN:g}\,{{\rm pc}})$ [$10^5\,M_\odot$]",
           ylabel=fr"recovered $-$ true $M_{{\rm rem}}(<{R_IN:g}\,{{\rm pc}})$ [$10^5\,M_\odot$]",
           title="(d) mocks with the stellar mass pinned: halo vs remnants")
    ax.legend(fontsize=7, loc="upper right")
    fig.suptitle(f"Who carries the mass inside {R_IN:g} pc: the kinematics fix the total (~1%), stars trade against remnants + halo,\n"
                 "and with the stars pinned the remnants trade against the halo")
    plot = ROOT/"plots"/f"remnant_halo_degeneracy_{args.tag}.png"
    fig.savefig(plot, dpi=150)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), r_in_pc=R_IN,
                  scan=[dict(rho20=p["rho20"], delta=p["s"]["refined_score"]-base,
                             **{k: enclosed(p["s"]["refined_profiles"], k, R_IN) for k in ("stars", "remnants", "halo", "total")},
                             M_rem=p["s"]["best"]["config"]["matter"]["M_rem"], a_rem=p["s"]["best"]["config"]["matter"]["a_rem"],
                             M_star=p["s"]["best"]["config"]["M_star"]) for p in scan],
                  mocks=rows, plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/f"remnant_halo_degeneracy_{args.tag}.json").write_text(json.dumps(record))
    for p in record["scan"]:
        print(f"rho20 {p['rho20']:4.2f} (+{p['delta']:5.2f}): inside {R_IN:g} pc stars {p['stars']:.3e} remnants {p['remnants']:.3e} halo {p['halo']:.3e} total {p['total']:.3e}  | M_rem {p['M_rem']:.3e} @ {p['a_rem']:.2f} pc, M* {p['M_star']:.3e}")
    for key in ("cov", "truth", "9.0", "12.0", "all"):
        rr = rows if key == "all" else [r for r in rows if r["set"] == key]
        h = np.array([r["halo_in"] for r in rr]); m = np.array([r["rem_in"] for r in rr]); st = np.array([r["star_in"] for r in rr])
        if len(rr) > 2:
            print(f"mocks {key:6s} n={len(rr):2d}  corr(stars, rem+halo) {np.corrcoef(st, h+m)[0,1]:+.2f}  corr(halo, rem) {np.corrcoef(h, m)[0,1]:+.2f}  "
                  f"rms inside {R_IN:g} pc: stars {np.std(st):.1e} rem {np.std(m):.1e} halo {np.std(h):.1e} total {np.std(h+m+st):.1e} Msun")
    print(plot)


if __name__ == "__main__":
    main()
