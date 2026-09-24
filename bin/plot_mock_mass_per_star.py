#!/usr/bin/env python3
"""Mock test of the pinned stellar mass per counted star (M/N19).

Collects the coverage mocks (results/df/cov_*, count amplitude profiled) and the
pinned refits of the same mock data (results/df/mpsmock_<M>_*), and plots:
 (a) recovered rho20 per seed, unpinned vs pinned at the truth M/N19, for both
     truths (A: rho20 = 0, B: rho20 = 1);
 (b) recovered rho20 against the assumed M/N19 (truth 10.67; wrong values 9, 12);
 (c) recovered remnant mass error against the assumed M/N19;
 (d) recovered total mass inside 20 pc against the assumed M/N19.
Usage: python bin/plot_mock_mass_per_star.py [--tag 20260923]
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
COUNTS = "ocen_counts_hst_f625w19"


def read(path):
    return json.loads(Path(path).read_text())


def collect(pattern):
    rows = []
    for b in sorted(glob.glob(str(ROOT/pattern))):
        fits = [read(f) for f in glob.glob(f"{b}/fits/*/summary.json")]
        fits = [s for s in fits if "best" in s]
        if not fits or not os.path.exists(f"{b}/mocks/mock/mock.json"):
            continue
        t = read(f"{b}/mocks/mock/mock.json"); s = min(fits, key=lambda s: s["refined_score"])
        c = s["best"]["config"]; tc = t["config"]
        name = os.path.basename(b)
        m = re.match(r"mpsmock_([^_]+)_", name)
        truth_mps = t.get("truth_mass_per_star", {}).get(COUNTS, 10.67)
        pinned = t.get("mass_per_star", {}).get(COUNTS)
        tot = lambda p, r: float(np.interp(r, p["r_pc"], p["total"]))
        rows.append(dict(name=name, seed=t["seed"], truth="B" if tc["matter"]["rho20"] > 0 else "A",
                         truth_rho20=tc["matter"]["rho20"], rho20=c["matter"]["rho20"], r_s=c["matter"]["r_s"],
                         assumed=pinned if pinned else truth_mps, pinned=pinned is not None, truth_mps=truth_mps,
                         mstar_err=c["M_star"]/tc["M_star"]-1, mrem_err=c["matter"]["M_rem"]/tc["matter"]["M_rem"]-1,
                         mtot20_err=tot(s["refined_profiles"], 20.)/tot(t["profiles"], 20.)-1,
                         chi2_kin=s["gates"]["chi2_kinematic"], message=s["optimizer"]["message"]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", default="20260923")
    args = parser.parse_args()
    free = collect(f"results/df/cov_*_{args.tag}")
    pinned = collect(f"results/df/mpsmock_*_{args.tag}")
    if not pinned:
        raise SystemExit("no pinned mock fits")
    col = {"A": "tab:blue", "B": "tab:orange"}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    # (a) per seed, unpinned vs truth-pinned
    ax = axes[0, 0]
    for truth in "AB":
        seeds = sorted({r["seed"] for r in free if r["truth"] == truth})
        for i, seed in enumerate(seeds):
            f = [r for r in free if r["truth"] == truth and r["seed"] == seed]
            p = [r for r in pinned if r["truth"] == truth and r["seed"] == seed and abs(r["assumed"]-r["truth_mps"]) < 1e-6]
            if f and p:
                ax.plot([0, 1], [f[0]["rho20"], p[0]["rho20"]], marker="o", color=col[truth], alpha=.8,
                        label=f"truth {truth}: rho20 = {f[0]['truth_rho20']:g}" if i == 0 else None)
        ref = next(r["truth_rho20"] for r in free if r["truth"] == truth)
        ax.axhline(ref, color=col[truth], ls=":", lw=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["amplitude profiled\n(coverage test)", "M/N19 pinned at truth"])
    ax.set(ylabel=r"recovered $\rho_{20}$ [$M_\odot$ pc$^{-3}$]", title="(a) same mock data, free vs pinned stellar mass (dotted: truth)")
    ax.legend(fontsize=8)
    # (b)-(d) vs assumed M/N19
    for ax, key, lab, title in ((axes[0, 1], "rho20", r"recovered $\rho_{20}$", "(b) rho20 vs assumed mass per star"),
                                (axes[1, 0], "mrem_err", r"$M_{\rm rem}$ / truth $-$ 1", "(c) remnant mass error"),
                                (axes[1, 1], "mtot20_err", r"$M_{\rm tot}(<20\,{\rm pc})$ / truth $-$ 1", "(d) total mass inside 20 pc")):
        for truth in "AB":
            rows = [r for r in pinned if r["truth"] == truth]
            ax.scatter([r["assumed"] for r in rows], [r[key] for r in rows], color=col[truth], s=50, edgecolor="black", zorder=3,
                       label=f"truth {truth}")
            for seed in sorted({r["seed"] for r in rows}):
                rr = sorted([r for r in rows if r["seed"] == seed], key=lambda r: r["assumed"])
                if len(rr) > 1:
                    ax.plot([r["assumed"] for r in rr], [r[key] for r in rr], color=col[truth], alpha=.5, lw=1)
            frees = [r for r in free if r["truth"] == truth]
            ax.scatter([r["assumed"] for r in frees], [r[key] for r in frees], marker="x", color=col[truth], s=40,
                       label=f"truth {truth}, amplitude profiled" if key == "rho20" else None)
        if key == "rho20":
            for truth in "AB":
                ax.axhline(next(r["truth_rho20"] for r in pinned if r["truth"] == truth), color=col[truth], ls=":", lw=1)
        else:
            ax.axhline(0, color="black", lw=.8)
        ax.axvline(pinned[0]["truth_mps"], color="0.5", ls="--", lw=1)
        ax.set(xlabel=r"assumed $M_\star$ per F625W$<$19 star [$M_\odot$] (dashed: truth)", ylabel=lab, title=title)
        ax.legend(fontsize=8)
    fig.suptitle("Pinning the stellar mass per counted star: realistic mocks (published noise, Poisson counts), one generic start each")
    plot = ROOT/"plots"/f"mock_mass_per_star_{args.tag}.png"
    fig.savefig(plot, dpi=150)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), free=free, pinned=pinned,
                  plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/f"mock_mass_per_star_{args.tag}.json").write_text(json.dumps(record))
    for truth in "AB":
        for label, rows in (("profiled", [r for r in free if r["truth"] == truth]),
                            *[(f"pinned {m:g}", [r for r in pinned if r["truth"] == truth and abs(r["assumed"]-m) < .01])
                              for m in sorted({round(r["assumed"], 2) for r in pinned})]):
            if rows:
                v = np.array([r["rho20"] for r in rows])
                print(f"truth {truth} (rho20 {rows[0]['truth_rho20']:g}) {label:14s} n={len(v)}  rho20 {v.min():.2f}-{v.max():.2f} (median {np.median(v):.2f})  "
                      f"M* {100*np.mean([r['mstar_err'] for r in rows]):+.1f}%  Mrem {100*np.mean([r['mrem_err'] for r in rows]):+.1f}%  "
                      f"Mtot(<20) {100*np.mean([r['mtot20_err'] for r in rows]):+.1f}%")
    print(plot)


if __name__ == "__main__":
    main()
