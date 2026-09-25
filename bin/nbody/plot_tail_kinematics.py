#!/usr/bin/env python3
"""Kinematics of the stellar tidal tails of the N-body models, side by side.

At the chosen snapshot the unbound stars of each model are placed in orbit coordinates
(time offset dt along the cluster orbit, 6D-matched as in analyse_tails.py). Rows:
 (a) integrals of motion: energy and z-angular-momentum offsets from the cluster,
     dE = E - E_cluster and dLz, in the MW potential (leading debris: dE < 0);
 (b) mean velocity of the debris relative to the local orbit velocity along dt, decomposed
     into the along-track component and the two perpendicular ones;
 (c) velocity dispersions of the same three components along dt;
 (d) distribution of dE for both models (overlaid), with its rms; the energy scale of
     tidal stripping is set by the cluster mass inside the tidal radius.
Numbers are also written to results/plot_data/<name>.json.

Usage: python bin/nbody/plot_tail_kinematics.py --name nbody_tail_kinematics --time 250 \
    results/nbody/A_nodm/orbit results/nbody/B_dm_phot/orbit
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
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
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"bin"/"nbody"))
from run_nbody import AGAMA_T_MYR, bound_mask, cluster_centre, mw_potential  # noqa: E402
from analyse_tails import T_ORB, V_SCALE  # noqa: E402

DT_MAX = 150.      # Myr: range of the along-orbit plots
MODEL_COLORS = ("tab:blue", "tab:orange")


def track(agama, pot, centre_pc, vcentre):
    ic = np.concatenate((centre_pc*1e-3, vcentre))
    parts = []
    for sign in (-1, 1):
        t, o = agama.orbit(potential=pot, ic=ic, time=sign*T_ORB/AGAMA_T_MYR, trajsize=int(T_ORB*4)+1)
        o = np.asarray(o); parts.append((np.asarray(t)*AGAMA_T_MYR, o[:, :3]*1e3, o[:, 3:]))
    t_all = np.concatenate((parts[0][0][::-1], parts[1][0][1:]))
    p_all = np.vstack((parts[0][1][::-1], parts[1][1][1:])); v_all = np.vstack((parts[0][2][::-1], parts[1][2][1:]))
    return t_all, p_all, v_all


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--name", required=True)
    parser.add_argument("--time", type=float, default=None)
    parser.add_argument("--eps", type=float, default=0.3)
    args = parser.parse_args()
    runs = [Path(r) for r in args.runs]
    agama, pot = mw_potential("McMillan17")
    if args.time is None:
        args.time = min(max(float(np.load(f)["t_myr"]) for f in glob.glob(str(r/"snap_*.npz"))) for r in runs)
    n = len(runs)
    fig, axes = plt.subplots(4, n, figsize=(7.5*n, 17), constrained_layout=True, squeeze=False)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), time_myr=args.time, models={})
    edges = np.linspace(-DT_MAX, DT_MAX, 31); mids = 0.5*(edges[1:]+edges[:-1])
    dE_all = {}
    for i, run in enumerate(runs):
        label = run.parent.name
        ics = np.load(run.parent/"ics.npz"); species = ics["species"]; mass = ics["mass"].astype(float)
        files = sorted(glob.glob(str(run/"snap_*.npz"))); times = np.array([float(np.load(f)["t_myr"]) for f in files])
        s = np.load(files[int(np.argmin(np.abs(times-args.time)))]); t = float(s["t_myr"]); pos = s["pos"].astype(float); vel = s["vel"].astype(float)
        lum = species <= 1
        centre, vcentre = cluster_centre(pos, vel, mass, lum, np.median(pos[lum], axis=0))
        bound = bound_mask(pos, vel, mass, args.eps, centre, vcentre)
        k = (species == 0) & ~bound
        p, v = pos[k], vel[k]
        # integrals of motion in the MW potential (kpc, km/s)
        E = 0.5*np.sum(v**2, axis=1)+pot.potential(p*1e-3); Ec = 0.5*np.sum(vcentre**2)+pot.potential((centre*1e-3)[None, :])[0]
        Lz = p[:, 0]*v[:, 1]-p[:, 1]*v[:, 0]; Lzc = centre[0]*vcentre[1]-centre[1]*vcentre[0]
        dE = E-Ec; dLz = (Lz-Lzc)*1e-3
        # orbit coordinates and velocities relative to the local track
        t_all, p_all, v_all = track(agama, pot, centre, vcentre)
        _, idx = cKDTree(np.hstack((p_all, V_SCALE*v_all))).query(np.hstack((p, V_SCALE*v)))
        dt = t_all[idx]
        vt = v_all[idx]; e_par = vt/np.linalg.norm(vt, axis=1)[:, None]
        r_hat = p_all[idx]/np.linalg.norm(p_all[idx], axis=1)[:, None]
        e_z = np.cross(r_hat, e_par); e_z /= np.maximum(np.linalg.norm(e_z, axis=1), 1e-12)[:, None]     # normal to the local orbital plane
        e_perp = np.cross(e_z, e_par)                                                                     # in-plane, perpendicular to the track
        dv = v-vt
        comps = dict(along=np.sum(dv*e_par, axis=1), in_plane_perp=np.sum(dv*e_perp, axis=1), normal=np.sum(dv*e_z, axis=1))
        near = np.abs(dt) < DT_MAX
        ax = axes[0, i]
        sc = ax.scatter(dLz, dE, c=dt, cmap="coolwarm", vmin=-DT_MAX, vmax=DT_MAX, s=2, rasterized=True)
        fig.colorbar(sc, ax=ax, label="time offset along the orbit [Myr]")
        ax.axhline(0, color="black", lw=.6); ax.axvline(0, color="black", lw=.6)
        ax.set(xlabel=r"$\Delta L_z$ [kpc km/s]", ylabel=r"$\Delta E$ [km$^2$ s$^{-2}$]", title=f"{label}, t = {t:.0f} Myr: {k.sum()} unbound stars ({mass[k].sum():.2e} Msun)")
        ax = axes[1, i]; bx = axes[2, i]
        out = dict(n_unbound=int(k.sum()), unbound_mass=float(mass[k].sum()), dE_rms=float(np.std(dE)), dE_percentiles=np.percentile(dE, [5, 50, 95]).tolist(),
                   dLz_rms=float(np.std(dLz)), profiles=dict(dt_myr=mids.tolist()))
        for key, col in zip(comps, ("tab:blue", "tab:green", "tab:purple")):
            ib = np.digitize(dt[near], edges)-1
            mean = np.array([comps[key][near][ib == j].mean() if np.sum(ib == j) >= 30 else np.nan for j in range(len(mids))])
            sig = np.array([comps[key][near][ib == j].std() if np.sum(ib == j) >= 30 else np.nan for j in range(len(mids))])
            ax.plot(mids, mean, marker=".", color=col, label=key); bx.plot(mids, sig, marker=".", color=col, label=key)
            out["profiles"][key+"_mean"] = mean.tolist(); out["profiles"][key+"_sigma"] = sig.tolist()
        for a in (ax, bx):
            a.axvline(0, color="black", lw=.6); a.set_xlim(-DT_MAX, DT_MAX); a.legend(fontsize=7)
        ax.axhline(0, color="black", lw=.6)
        ax.set(xlabel="time offset along the orbit [Myr] (>0 leading)", ylabel="mean velocity relative to the orbit [km/s]", title="(b) debris streaming relative to the local orbit")
        bx.set(xlabel="time offset along the orbit [Myr]", ylabel="velocity dispersion [km/s]", title="(c) velocity dispersion of the debris", yscale="log")
        far = near & (np.abs(dt) > 10.)
        out["sigma_tail_10_150"] = {key: float(comps[key][far].std()) for key in comps}
        out["mean_tail_10_150"] = {key: float(comps[key][far].mean()) for key in comps}
        dE_all[label] = dE
        record["models"][label] = out
    ax = axes[3, 0]
    lim = max(np.percentile(np.abs(v), 99) for v in dE_all.values())
    for (label, dE), col in zip(dE_all.items(), MODEL_COLORS):
        ax.hist(dE, np.linspace(-lim, lim, 80), histtype="step", color=col, lw=1.5, density=True, label=f"{label}: rms {np.std(dE):.0f} km$^2$ s$^{{-2}}$")
    ax.axvline(0, color="black", lw=.6); ax.set(xlabel=r"$\Delta E$ [km$^2$ s$^{-2}$]", ylabel="density", title="(d) energy offset of the unbound stars"); ax.legend(fontsize=8)
    for j in range(1, n):
        axes[3, j].axis("off")
    if n > 1:
        ax = axes[3, 1]; ax.axis("on")
        for (label, out), col in zip(record["models"].items(), MODEL_COLORS):
            ax.bar(np.arange(3)+(-.2 if col == MODEL_COLORS[0] else .2), [out["sigma_tail_10_150"][key] for key in ("along", "in_plane_perp", "normal")], .4, color=col, label=label)
        ax.set_xticks(range(3)); ax.set_xticklabels(["along track", "in-plane perp.", "normal"]); ax.set(ylabel="dispersion [km/s]", title="(e) tail velocity dispersions, 10 < |dt| < 150 Myr"); ax.legend(fontsize=8)
    fig.suptitle(f"Kinematics of the stellar tidal tails at t = {args.time:.0f} Myr (unbound stars only)")
    plot = ROOT/"plots"/f"{args.name}_t{int(round(args.time)):04d}.png"
    fig.savefig(plot, dpi=130)
    (ROOT/"results/plot_data"/f"{args.name}_t{int(round(args.time)):04d}.json").write_text(json.dumps(record, default=float))
    for label, out in record["models"].items():
        print(f"{label}: {out['n_unbound']} unbound stars ({out['unbound_mass']:.2e} Msun); dE rms {out['dE_rms']:.0f} km2/s2 (5-95%: {out['dE_percentiles'][0]:+.0f}..{out['dE_percentiles'][2]:+.0f}); "
              f"dLz rms {out['dLz_rms']:.1f} kpc km/s; tail dispersions (10<|dt|<150 Myr): " + ", ".join(f"{k} {v:.1f}" for k, v in out["sigma_tail_10_150"].items()) + " km/s")
    print(plot)


if __name__ == "__main__":
    main()
