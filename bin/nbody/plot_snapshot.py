#!/usr/bin/env python3
"""State of disruption of the N-body models at one snapshot, models side by side.

For each run the latest common snapshot (or --time) is shown in the Galactocentric X-Y
plane: all particles, stars in blue, remnants in red, DM halo in green (drawn first),
the cluster orbit over +-T_ORB Myr in grey and the cluster centre as a cross. A zoomed
inset shows the inner +-ZOOM kpc. The legend gives the bound mass fraction per species
from diagnostics.jsonl at that time.

Usage: python bin/nbody/plot_snapshot.py --name nbody_state results/nbody/A_nodm/orbit results/nbody/B_dm_phot/orbit
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"bin"/"nbody"))
from run_nbody import SPECIES, cluster_centre, mw_potential  # noqa: E402

COLORS = dict(stars="tab:blue", remnants="tab:red", halo="tab:green")
T_ORB = 300.
ZOOM = 1.0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--name", required=True)
    parser.add_argument("--time", type=float, default=None)
    args = parser.parse_args()
    runs = [Path(r) for r in args.runs]
    agama, pot = mw_potential("McMillan17")
    if args.time is None:
        args.time = min(max(float(np.load(f)["t_myr"]) for f in glob.glob(str(r/"snap_*.npz"))) for r in runs)
    fig, axes = plt.subplots(1, len(runs), figsize=(8*len(runs), 8), constrained_layout=True, squeeze=False)
    for ax, run in zip(axes[0], runs):
        ics = np.load(run.parent/"ics.npz"); species = ics["species"]; mass = ics["mass"].astype(float)
        files = sorted(glob.glob(str(run/"snap_*.npz"))); times = np.array([float(np.load(f)["t_myr"]) for f in files])
        s = np.load(files[int(np.argmin(np.abs(times-args.time)))]); t = float(s["t_myr"]); pos = s["pos"].astype(float); vel = s["vel"].astype(float)
        lum = species <= 1
        centre, vcentre = cluster_centre(pos, vel, mass, lum, np.median(pos[lum], axis=0))
        diag = [json.loads(l) for l in open(run/"diagnostics.jsonl")]
        d = min(diag, key=lambda d: abs(d["t_myr"]-t)); m0 = diag[0]["bound_mass"]
        ic = np.concatenate((centre*1e-3, vcentre))
        for sign, col in ((-1, "0.6"), (1, "0.8")):
            _, o = agama.orbit(potential=pot, ic=ic, time=sign*T_ORB*1e-3, trajsize=1000)
            ax.plot(o[:, 0], o[:, 1], color=col, lw=.8, zorder=0, label="orbit, past 300 Myr" if sign < 0 else "orbit, next 300 Myr")
        for si in (2, 0, 1):
            k = species == si
            if not k.any():
                continue
            sp = SPECIES[si]
            ax.scatter(pos[k, 0]*1e-3, pos[k, 1]*1e-3, s=.3 if si != 1 else .6, color=COLORS[sp], alpha=.25 if si == 2 else .5, rasterized=True,
                       label=f"{sp}: bound {d['bound_mass'][sp]/m0[sp]*100:.1f}% of {m0[sp]:.2e} Msun")
        ax.plot(centre[0]*1e-3, centre[1]*1e-3, "k+", ms=12, mew=1.5)
        ax.set(xlim=(-9, 9), ylim=(-9, 9), aspect="equal", xlabel="X [kpc]", ylabel="Y [kpc]",
               title=f"{run.parent.name}: t = {t:.0f} Myr ({2000-t:.0f} Myr before today), r_gal = {np.linalg.norm(centre)*1e-3:.2f} kpc")
        ax.legend(fontsize=8, markerscale=10, loc="upper left")
        ins = ax.inset_axes([0.62, 0.02, 0.36, 0.36])
        for si in (2, 0, 1):
            k = species == si
            if k.any():
                ins.scatter((pos[k, 0]-centre[0])*1e-3, (pos[k, 1]-centre[1])*1e-3, s=.3, color=COLORS[SPECIES[si]], alpha=.3, rasterized=True)
        ins.set(xlim=(-ZOOM, ZOOM), ylim=(-ZOOM, ZOOM), aspect="equal", xticks=[-ZOOM, 0, ZOOM], yticks=[-ZOOM, 0, ZOOM])
        ins.set_title(f"inner {2*ZOOM:.0f} kpc around the cluster", fontsize=8); ins.tick_params(labelsize=7)
    fig.suptitle("Tidal disruption of the no-DM (A) and DM (B) models on omega Cen's orbit (McMillan17)")
    plot = ROOT/"plots"/f"{args.name}_t{int(round(args.time)):04d}.png"
    fig.savefig(plot, dpi=140)
    print(plot)


if __name__ == "__main__":
    main()
