#!/usr/bin/env python3
"""Plot the star/halo trade-off seen in a finished compact-DF recovery batch.

Uses only saved per-trial logs and refined profiles; no model is evaluated.
Top row: every optimizer trial in (M_star, rho20), coloured by objective.
These are optimizer paths, not a profile likelihood or posterior.
Bottom row: refined enclosed-mass components relative to the injected total.

Usage: python bin/plot_compact_degeneracy.py results/df/compact_recovery_20260923_v2
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text())


def trials(path):
    rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
    return [r for r in rows if "score" in r]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batch", type=Path)
    out = parser.parse_args().batch.resolve()
    manifest = read(out/"batch.json")
    cases = list(dict.fromkeys(j["case"] for j in manifest["jobs"]))
    colors = dict(zip([j["id"] for j in manifest["jobs"]], plt.rcParams["axes.prop_cycle"].by_key()["color"]))
    fig, axes = plt.subplots(2, len(cases), figsize=(11, 8), squeeze=False, constrained_layout=True)
    norm = LogNorm(1e-3, 1e4)
    record = dict(batch=str(out), created_utc=datetime.now(timezone.utc).isoformat(),
                  source_commit=manifest["source_commit"], cases={})
    for col, case in enumerate(cases):
        mock = read(out/"mocks"/case/"mock.json")
        truth, tp = mock["config"], mock["profiles"]
        ax, bx = axes[0, col], axes[1, col]
        record["cases"][case] = dict(truth=dict(M_star=truth["M_star"], rho20=truth["matter"]["rho20"]), jobs={})
        for job in (j for j in manifest["jobs"] if j["case"] == case):
            d = out/"fits"/job["id"]
            rows = trials(d/"evaluations.jsonl")
            m = np.array([r["M_star"] for r in rows])/1e6
            rho = np.array([r["rho20"] for r in rows])
            score = np.array([r["score"] for r in rows])
            name = job["id"].removeprefix(case+"_")
            ax.plot(m, rho, color=colors[job["id"]], lw=.6, alpha=.5, zorder=1)
            sc = ax.scatter(m, rho, c=np.clip(score, norm.vmin, norm.vmax),
                            norm=norm, cmap="viridis_r", s=10, zorder=2)
            s = read(d/"summary.json")
            p = s["refined_profiles"]
            ax.scatter(p["M_star"]/1e6, p["rho20"], marker="X", s=90, color=colors[job["id"]],
                       edgecolor="black", zorder=4, label=f"{name}: final, obj {s['refined_score']:.3g}")
            r = np.array(p["r_pc"])
            bx.plot(r, np.array(p["stars"])/tp["total"], color=colors[job["id"]], label=f"{name} stars")
            if job["branch"] == "free_halo":
                bx.plot(r, np.array(p["halo"])/tp["total"], color=colors[job["id"]], ls="--", label=f"{name} halo")
            record["cases"][case]["jobs"][job["id"]] = dict(
                trial_M_star=(m*1e6).tolist(), trial_rho20=rho.tolist(), trial_score=score.tolist(),
                final_M_star=p["M_star"], final_rho20=p["rho20"], final_score=s["refined_score"],
                r_pc=p["r_pc"], stars=p["stars"], halo=p["halo"], remnants=p["remnants"], total=p["total"])
        ax.scatter(truth["M_star"]/1e6, truth["matter"]["rho20"], marker="*", s=250, color="red",
                   edgecolor="black", zorder=5, label="injected")
        ax.set(xlabel=r"$M_\star$ [$10^6\,M_\odot$]", ylabel=r"$\rho_{\rm DM}(20\,{\rm pc})$ [$M_\odot\,{\rm pc}^{-3}$]",
               title=case.replace("_", " ")+": optimizer trials (not a likelihood)")
        ax.legend(fontsize=7, loc="upper right")
        r = np.array(tp["r_pc"])
        bx.plot(r, np.array(tp["stars"])/np.array(tp["total"]), color="black", lw=2, label="injected stars")
        if truth["matter"]["rho20"] > 0:
            bx.plot(r, np.array(tp["halo"])/np.array(tp["total"]), color="black", lw=2, ls="--", label="injected halo")
        bx.set(xscale="log", xlabel="Radius [pc]", ylabel="Component / injected total enclosed mass",
               title=case.replace("_", " ")+": refined components")
        bx.legend(fontsize=7, ncol=2)
    fig.colorbar(sc, ax=axes[0, :], label=r"objective $\sum$ (residual/error)$^2$")
    fig.suptitle("Noiseless compact DF recovery: stellar mass vs halo density")
    plot = ROOT/"plots"/(out.name+"_degeneracy.png")
    fig.savefig(plot, dpi=160)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/(out.name+"_degeneracy.json")).write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
