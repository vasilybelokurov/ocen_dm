#!/usr/bin/env python3
"""Photometric stellar mass per counted star against the dynamical requirement.

Left: stacked bars per radial zone of the mass that follows the light per F625W < 19
star from the oMEGACat mass function (observed main sequence 0.52-0.78 Msun, main
sequence extrapolated to 0.1 Msun with the fitted slope, evolved stars, white dwarfs
from an IMF slope -2.3 above the turnoff), for the two isochrone metallicities; markers
give the alternative extrapolations (Kroupa below 0.5 Msun; flat, an illustrative
lower bound). The dynamical value from the free fit and its Delta = 1 / 3.84 ranges
from the M/N19 scan are drawn as horizontal bands.
Right: the same photometric values placed on the real-data scan: rho20 and
M_DM(<10 pc) the dynamics return when the stellar mass per counted star is pinned.

Usage: python bin/plot_mass_budget_tension.py [--tag 20260923]
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
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WD_KEY = "Kroupa/Salpeter IMF (-2.3)"


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", default="20260923")
    args = parser.parse_args()
    mfs = [read(ROOT/f"results/mass_function/mass_function_{t}.json") for t in (args.tag, f"feh153_{args.tag}")]
    scan = read(ROOT/"results/plot_data"/f"mass_per_star_profile_{args.tag}.json")
    x = np.array([p["mass_per_star"] for p in scan["points"]]); obj = np.array([p["objective"] for p in scan["points"]])
    grid = np.linspace(x.min(), x.max(), 4001); oi = np.interp(grid, x, obj)
    free = x[np.argmin(obj)]
    d1 = grid[oi <= obj.min()+1.]; d4 = grid[oi <= obj.min()+3.84]
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8), constrained_layout=True, gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax = axes[0]
    zones = mfs[0]["zones"]; nz = len(zones); width = .38
    comp_colors = dict(obs="tab:blue", ext="tab:cyan", ev="tab:orange", wd="0.5")
    rows = []
    for k, (mf, off, hatch) in enumerate(zip(mfs, (-width/2, width/2), (None, "//"))):
        for zi, z in enumerate(mf["zones"]):
            n = z["n_f625w_lt19"]; wd = z["white_dwarfs"][WD_KEY]["n_wd"]*mf["wd_mass"]/n
            parts = [("obs", z["mass_ms_observed"]/n, "MS observed 0.52-0.78 Msun"), ("ext", z["mass_ms_extrapolated"]["fitted slope"]/n, "MS extrapolated to 0.1 Msun (fitted slope)"),
                     ("ev", z["mass_evolved"]/n, "evolved stars"), ("wd", wd, "white dwarfs (IMF -2.3, 0.55 Msun)")]
            bottom = 0.
            for key, val, lab in parts:
                ax.bar(zi+off, val, width, bottom=bottom, color=comp_colors[key], hatch=hatch, edgecolor="black", lw=.5,
                       label=lab if (zi == 0 and k == 0) else None)
                bottom += val
            kro = z["mass_per_n19"]["Kroupa (-1.3 below 0.5) | "+WD_KEY]; flat = z["mass_per_n19"]["flat (alpha=0) | "+WD_KEY]
            ax.plot(zi+off, kro, "v", color="black", ms=6, label="Kroupa extrapolation below 0.5 Msun" if (zi == 0 and k == 0) else None)
            ax.plot(zi+off, flat, "_", color="black", ms=12, mew=2, label="flat extrapolation (lower bound)" if (zi == 0 and k == 0) else None)
            ax.text(zi+off, bottom+.15, f"{bottom:.1f}", ha="center", fontsize=7)
            rows.append(dict(feh=mf["feh_adopted"], zone_arcsec=z["zone_arcsec"], alpha=z["alpha"], total_fitted=bottom, total_kroupa=kro, total_flat=flat,
                             parts={p[0]: p[1] for p in parts}))
    ax.axhspan(d4.min(), d4.max(), color="tab:red", alpha=.12, label=r"dynamics: $\Delta$ objective $\leq$ 3.84")
    ax.axhspan(d1.min(), d1.max(), color="tab:red", alpha=.25, label=r"dynamics: $\Delta \leq 1$")
    ax.axhline(free, color="tab:red", lw=1.5, label=f"dynamics, free fit: {free:.2f}")
    ax.set_xticks(range(nz)); ax.set_xticklabels([f"{z['zone_arcsec'][0]:.0f}-{z['zone_arcsec'][1]:.0f}\"\n{z['zone_arcsec'][0]*0.02633:.1f}-{z['zone_arcsec'][1]*0.02633:.1f} pc"
                                                   for z in zones], fontsize=8)
    ax.set(ylabel=r"stellar mass per F625W$<$19 star [$M_\odot$]", ylim=(0, 13.5),
           title=f"(a) photometric mass budget (solid: [Fe/H] {mfs[0]['feh_adopted']:+.2f}; hatched: {mfs[1]['feh_adopted']:+.2f}) vs the dynamical requirement")
    ax.legend(fontsize=7, loc="upper left", ncol=2)
    # right: scan
    ax = axes[1]
    rho = np.array([p["rho20"] for p in scan["points"]]); m10 = np.array([p["M_halo_10"] for p in scan["points"]])
    ax.plot(x, rho, "o-", color="tab:blue", label=r"$\rho_{20}$ when M/N19 is pinned (real data)")
    ax.set(xlabel=r"stellar mass per F625W$<$19 star [$M_\odot$]", ylabel=r"$\rho_{20}$ [$M_\odot$ pc$^{-3}$]", title="(b) what the dynamics return for each assumed value")
    ax2 = ax.twinx(); ax2.plot(x, m10/1e5, "s--", color="tab:green", label=r"$M_{\rm DM}(<10\,{\rm pc})$ [$10^5 M_\odot$]"); ax2.set_ylabel(r"$M_{\rm DM}(<10\,{\rm pc})$ [$10^5\,M_\odot$]")
    core = [r["total_fitted"] for r in rows if 30 <= r["zone_arcsec"][0] < 175]+[r["total_kroupa"] for r in rows if 30 <= r["zone_arcsec"][0] < 175]
    ax.axvspan(min(core), max(core), color="0.85", label=f"photometric, 30-175\" zones: {min(core):.1f}-{max(core):.1f}")
    ax.axvspan(min(r["total_flat"] for r in rows), max(max(r["total_fitted"], r["total_kroupa"]) for r in rows), color="0.93", zorder=0, label="all zones and variants")
    ax.axvspan(d1.min(), d1.max(), color="tab:red", alpha=.25, label=r"dynamics $\Delta \leq 1$")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1+h2, l1+l2, fontsize=7, loc="upper right")
    fig.suptitle("The stellar mass per counted star: photometry (mass function) says 8-10, the dynamics want ~10.7; the gap is what the halo fills")
    plot = ROOT/"plots"/f"mass_budget_tension_{args.tag}.png"
    fig.savefig(plot, dpi=150)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), rows=rows, dynamical=dict(free=float(free), delta1=[float(d1.min()), float(d1.max())],
                  delta384=[float(d4.min()), float(d4.max())]), plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/f"mass_budget_tension_{args.tag}.json").write_text(json.dumps(record))
    print(f"dynamics: free {free:.2f}, Delta<=1 [{d1.min():.2f}, {d1.max():.2f}], Delta<=3.84 [{d4.min():.2f}, {d4.max():.2f}]")
    for r in rows:
        print(f"[Fe/H] {r['feh']:+.2f} zone {r['zone_arcsec'][0]:4.0f}-{r['zone_arcsec'][1]:4.0f}  alpha {r['alpha']:+.2f}  fitted {r['total_fitted']:5.2f}  Kroupa {r['total_kroupa']:5.2f}  flat {r['total_flat']:5.2f}   parts obs {r['parts']['obs']:.2f} ext {r['parts']['ext']:.2f} ev {r['parts']['ev']:.2f} wd {r['parts']['wd']:.2f}")
    print(plot)


if __name__ == "__main__":
    main()
