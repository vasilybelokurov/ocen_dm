#!/usr/bin/env python3
"""Tidal debris of the N-body models: mass-loss history and tail structure, models side by side.

Inputs: one or more run directories (results/nbody/<model>/orbit). Two products:
 1. History (<name>_history.png): bound mass per species versus time from diagnostics.jsonl,
    with pericentre passages marked, for all runs.
 2. Tails at a chosen time (<name>_tails_tXXXX.png; nearest snapshot): unbound particles
    (energy in the cluster's own softened potential, as in run_nbody) placed in orbit
    coordinates. The cluster-centre orbit is integrated +-T_ORB Myr in the same MW
    potential; each particle gets the time offset dt of the nearest orbit point (dt > 0
    leading, i.e. ahead of the cluster) and its perpendicular distance. Panels: X-Y and
    X-Z projections with the orbit; N(dt) per species; width and velocity dispersion
    versus dt; DM-to-star ratio along the tails (if a halo is present). At t = tback the
    debris is also shown on the sky (l, b, distance, v_los), converting the simulation
    frame (Baumgardt's X toward the GC from the Sun, Y, Z; U, V, W) to astropy's
    Galactocentric frame (x = -X, v_x = -U) with the solar parameters that reproduce the
    catalogue's heliocentric values for omega Cen (checked and printed).

Usage:
  python bin/nbody/analyse_tails.py --name tails_20260924 --time 2000 \
      results/nbody/A_nodm/orbit results/nbody/B_dm_phot/orbit
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import glob
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "2")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"bin"/"nbody"))
from run_nbody import AGAMA_T_MYR, G, OCEN_TODAY, SPECIES, T_UNIT, bound_mask, cluster_centre, mw_potential  # noqa: E402

COLORS = dict(stars="tab:blue", remnants="tab:red", halo="tab:green")
T_ORB = 400.       # Myr each way for the orbit-coordinate reference track
# Solar parameters reproducing Baumgardt's heliocentric -> Galactocentric conversion (checked in main)
R_SUN_KPC, Z_SUN_KPC, V_SUN = 8.178, 0.0, (11.1, 12.24+240.0, 7.25)   # reproduces l, b, d, v_los, pm to 0.02 deg / 0.06 km/s / 0.006 mas/yr
V_SCALE = 10.       # pc per km/s: phase-space metric for matching debris to the orbit track (10 km/s ~ 100 pc)


def read_diag(run):
    return [json.loads(l) for l in open(Path(run)/"diagnostics.jsonl")]


def pericentres(diag):
    r = np.array([d["r_gal_kpc"] for d in diag]); t = np.array([d["t_myr"] for d in diag])
    k = np.where((r[1:-1] < r[:-2]) & (r[1:-1] < r[2:]))[0]+1
    return t[k], r[k]


def load_snapshot(run, t_request):
    files = sorted(glob.glob(str(Path(run)/"snap_*.npz")))
    times = np.array([float(np.load(f)["t_myr"]) for f in files])
    k = int(np.argmin(np.abs(times-t_request)))
    s = np.load(files[k])
    return float(s["t_myr"]), s["pos"].astype(float), s["vel"].astype(float)


def orbit_coordinates(agama, pot, centre_pc, vcentre, pos_pc, vel):
    """Time offset [Myr] along the cluster orbit and perpendicular distance [pc] of each particle.

    Debris is matched to the orbit track in 6D (positions in pc, velocities scaled by
    V_SCALE) so that where the rosette orbit crosses itself a particle is assigned to the
    loop whose velocity it shares, not merely the nearest point in space."""
    ic = np.concatenate((centre_pc*1e-3, vcentre))
    tracks = []
    for sign in (-1, 1):
        t, o = agama.orbit(potential=pot, ic=ic, time=sign*T_ORB/AGAMA_T_MYR, trajsize=int(T_ORB*4)+1)
        o = np.asarray(o); tracks.append((np.asarray(t)*AGAMA_T_MYR, o[:, :3]*1e3, o[:, 3:]))
    t_all = np.concatenate((tracks[0][0][::-1], tracks[1][0][1:]))
    p_all = np.vstack((tracks[0][1][::-1], tracks[1][1][1:])); v_all = np.vstack((tracks[0][2][::-1], tracks[1][2][1:]))
    tree = cKDTree(np.hstack((p_all, V_SCALE*v_all)))
    _, idx = tree.query(np.hstack((pos_pc, V_SCALE*vel)))
    d = np.linalg.norm(pos_pc-p_all[idx], axis=1)
    return t_all[idx], d, (t_all, p_all)


def to_sky(pos_pc, vel):
    """Simulation (Baumgardt) frame -> heliocentric l, b [deg], distance [kpc], v_los [km/s], pm_l*, pm_b [mas/yr]."""
    import astropy.units as u
    from astropy.coordinates import Galactocentric, Galactic, CartesianDifferential, SkyCoord
    frame = Galactocentric(galcen_distance=R_SUN_KPC*u.kpc, z_sun=Z_SUN_KPC*u.kpc,
                           galcen_v_sun=CartesianDifferential(V_SUN[0]*u.km/u.s, V_SUN[1]*u.km/u.s, V_SUN[2]*u.km/u.s))
    c = SkyCoord(x=-pos_pc[:, 0]*u.pc, y=pos_pc[:, 1]*u.pc, z=pos_pc[:, 2]*u.pc,
                 v_x=-vel[:, 0]*u.km/u.s, v_y=vel[:, 1]*u.km/u.s, v_z=vel[:, 2]*u.km/u.s, frame=frame)
    g = c.transform_to(Galactic())
    return (g.l.deg, g.b.deg, g.distance.kpc, g.radial_velocity.to_value(u.km/u.s),
            g.pm_l_cosb.to_value(u.mas/u.yr), g.pm_b.to_value(u.mas/u.yr))


def check_frame():
    """Heliocentric values of omega Cen from OCEN_TODAY with the adopted solar parameters."""
    l, b, d, vr, pml, pmb = to_sky(OCEN_TODAY[None, :3]*1e3, OCEN_TODAY[None, 3:])
    ref = dict(l=309.102, b=14.968, d=5.43, vr=232.78, pml=None, pmb=None)
    print(f"frame check omega Cen: l {l[0]:.3f} (cat 309.102), b {b[0]:.3f} (14.968), d {d[0]:.3f} (5.43), v_los {vr[0]:.1f} (232.78), pm_l* {pml[0]:.2f}, pm_b {pmb[0]:.2f}")
    return dict(l=float(l[0]), b=float(b[0]), d=float(d[0]), v_los=float(vr[0]), pm_l=float(pml[0]), pm_b=float(pmb[0]), reference=ref)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--name", required=True)
    parser.add_argument("--time", type=float, default=None, help="Myr; default: latest common snapshot")
    parser.add_argument("--eps", type=float, default=0.3)
    args = parser.parse_args()
    runs = [Path(r) for r in args.runs]
    labels = [r.parent.name for r in runs]
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), runs=[str(r) for r in runs], frame_check=check_frame(), history={}, tails={})
    # ---- 1. history
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8), constrained_layout=True)
    for run, label, ls in zip(runs, labels, ("-", "--", ":")):
        diag = read_diag(run); t = np.array([d["t_myr"] for d in diag]); m0 = diag[0]["bound_mass"]
        tp, rp = pericentres(diag)
        for sp in m0:
            axes[0].plot(t, [d["bound_mass"][sp]/m0[sp] for d in diag], ls=ls, color=COLORS[sp], label=f"{label}: {sp}")
        axes[1].plot(t, [d["r_gal_kpc"] for d in diag], ls=ls, color="0.3", label=label if ls == "-" else None)
        axes[1].plot(t, [d["r_half_stars_pc"] for d in diag], ls=ls, color="tab:blue", label=f"{label}: r_half stars [pc]")
        for tt in tp:
            axes[0].axvline(tt, color="0.85", lw=.6, zorder=0)
        record["history"][label] = dict(t_myr=t.tolist(), bound={sp: [d["bound_mass"][sp] for d in diag] for sp in m0},
                                        r_gal_kpc=[d["r_gal_kpc"] for d in diag], r_half=[d["r_half_stars_pc"] for d in diag],
                                        pericentre_t=tp.tolist(), pericentre_r=rp.tolist(), final={sp: diag[-1]["bound_mass"][sp]/m0[sp] for sp in m0})
    axes[0].set(xlabel="t [Myr]", ylabel="bound mass / initial", title="(a) bound mass per species (grey: pericentres)"); axes[0].legend(fontsize=7)
    axes[1].set(xlabel="t [Myr]", ylabel="kpc (orbit) / pc (r_half)", title="(b) Galactocentric radius and stellar half-mass radius"); axes[1].legend(fontsize=7)
    fig.savefig(ROOT/"plots"/f"{args.name}_history.png", dpi=140)
    # ---- 2. tails at a time
    agama, pot = mw_potential("McMillan17")
    t_req = args.time if args.time is not None else min(max(float(np.load(f)["t_myr"]) for f in glob.glob(str(r/"snap_*.npz"))) for r in runs)
    n = len(runs)
    fig, axes = plt.subplots(3, n, figsize=(7*n, 15), constrained_layout=True, squeeze=False)
    fig2, axes2 = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    tback = (lambda j: j.get("t_today_myr", j["tback_myr"]))(json.loads((runs[0]/"run.json").read_text()))
    for i, (run, label) in enumerate(zip(runs, labels)):
        ics = np.load(run.parent/"ics.npz"); species = ics["species"]; mass = ics["mass"].astype(float)
        t_snap, pos, vel = load_snapshot(run, t_req)
        sel_lum = species <= 1
        centre, vcentre = cluster_centre(pos, vel, mass, sel_lum, np.median(pos[sel_lum], axis=0))
        bound = bound_mask(pos, vel, mass, args.eps, centre, vcentre)
        dt, dperp, (t_track, p_track) = orbit_coordinates(agama, pot, centre, vcentre, pos, vel)
        rel = pos-centre
        out = dict(t_myr=t_snap, centre_kpc=(centre*1e-3).tolist(), r_gal_kpc=float(np.linalg.norm(centre)*1e-3), species={})
        ax = axes[0, i]
        ax.plot(p_track[:, 0]*1e-3, p_track[:, 1]*1e-3, color="0.7", lw=.8, zorder=0)
        for si, sp in enumerate(SPECIES):
            k = (species == si) & ~bound
            if not k.any():
                continue
            ax.scatter(pos[k, 0]*1e-3, pos[k, 1]*1e-3, s=.5, color=COLORS[sp], alpha=.4, label=f"unbound {sp}: {mass[k].sum():.2e} Msun", rasterized=True)
            lead = k & (dt > 0); trail = k & (dt < 0)
            near = k & (np.abs(dt) < T_ORB)
            out["species"][sp] = dict(unbound_mass=float(mass[k].sum()), bound_mass=float(mass[(species == si) & bound].sum()),
                                      leading_mass=float(mass[lead].sum()), trailing_mass=float(mass[trail].sum()),
                                      dt_percentiles_myr=np.percentile(dt[near], [5, 25, 50, 75, 95]).tolist() if near.any() else None,
                                      width_median_pc=float(np.median(dperp[near])) if near.any() else None)
        ax.plot(centre[0]*1e-3, centre[1]*1e-3, "k+", ms=10)
        ax.set(xlabel="X [kpc]", ylabel="Y [kpc]", title=f"{label}, t = {t_snap:.0f} Myr, r_gal {np.linalg.norm(centre)*1e-3:.2f} kpc", aspect="equal")
        ax.legend(fontsize=7, markerscale=8)
        ax = axes[1, i]
        edges = np.linspace(-T_ORB, T_ORB, 81)
        for si, sp in enumerate(SPECIES):
            k = (species == si) & ~bound
            if k.any():
                ax.hist(dt[k], edges, weights=mass[k], histtype="step", color=COLORS[sp], label=sp)
        ax.axvline(0, color="black", lw=.7); ax.set(xlabel="time offset along the orbit [Myr] (>0 leading)", ylabel=r"unbound mass per 10 Myr [$M_\odot$]", yscale="log", title="(b) debris along the orbit")
        ax.legend(fontsize=7)
        ax = axes[2, i]
        mids = 0.5*(edges[1:]+edges[:-1])
        for si, sp in enumerate(SPECIES):
            k = (species == si) & ~bound & (np.abs(dt) < T_ORB)
            if k.sum() < 100:
                continue
            idx = np.digitize(dt[k], edges)-1
            w = np.array([np.median(dperp[k][idx == j]) if np.sum(idx == j) >= 20 else np.nan for j in range(len(mids))])
            ax.plot(mids, w, color=COLORS[sp], label=f"{sp}: median perpendicular distance")
        ax.set(xlabel="time offset along the orbit [Myr]", ylabel="width [pc]", yscale="log", title="(c) tail width"); ax.legend(fontsize=7)
        record["tails"][label] = out
        # sky view (only meaningful at t = tback = today)
        if abs(t_snap-tback) < 1.:
            k = ~bound & (species == 0)
            wrap = lambda x: np.where(x > 180., x-360., x)          # l in (-180, 180] so the debris is not cut at l = 0
            l, b, d, vr, pml, pmb = to_sky(pos[k], vel[k]); l = wrap(l)
            ax = axes2[0, i]
            if (species == 2).any():
                kh = ~bound & (species == 2); lh, bh, *_ = to_sky(pos[kh], vel[kh])
                ax.scatter(wrap(lh), bh, s=.3, color="tab:green", alpha=.15, rasterized=True, label="unbound DM", zorder=1)
            ax.scatter(l, b, s=.5, color="tab:blue", alpha=.5, rasterized=True, label=f"{label}: unbound stars", zorder=2)
            lc, bc, *_ = to_sky(centre[None, :], vcentre[None, :]); lc = wrap(lc); ax.plot(lc, bc, "k+", ms=12, zorder=3)
            ax.scatter(wrap(np.array([309.102])), [14.968], marker="x", color="red", s=6, zorder=3, label="omega Cen observed")
            ax.set_xlim(75, -75); ax.set_ylim(-60, 60)
            ax.set(xlabel="l [deg]", ylabel="b [deg]", title=f"{label}: debris on the sky today (+: model cluster)"); ax.legend(fontsize=7, markerscale=8)
            ax = axes2[1, i]
            ax.scatter(l, vr, s=.5, color="tab:blue", alpha=.4, rasterized=True); ax.set(xlabel="l [deg] (wrapped to -180..180)", ylabel=r"$v_{\rm los}$ [km/s]", title="line-of-sight velocity of the unbound stars"); ax.set_xlim(75, -75)
            out["sky"] = dict(l_percentiles=np.percentile(l, [5, 50, 95]).tolist(), b_percentiles=np.percentile(b, [5, 50, 95]).tolist(),
                              d_percentiles=np.percentile(d, [5, 50, 95]).tolist(), vlos_percentiles=np.percentile(vr, [5, 50, 95]).tolist(),
                              centre=dict(l=float(lc[0]), b=float(bc[0])))
    fig.suptitle(f"Tidal debris at t = {t_req:.0f} Myr (unbound = positive energy in the cluster's own potential)")
    fig.savefig(ROOT/"plots"/f"{args.name}_tails_t{int(round(t_req)):04d}.png", dpi=140)
    if record["tails"] and any("sky" in v for v in record["tails"].values()):
        fig2.savefig(ROOT/"plots"/f"{args.name}_sky.png", dpi=140)
    (ROOT/"results/plot_data"/f"{args.name}.json").write_text(json.dumps(record, default=float))
    for label, h in record["history"].items():
        print(f"{label}: t_end {h['t_myr'][-1]:.0f} Myr, {len(h['pericentre_t'])} pericentres (r_min {min(h['pericentre_r']) if h['pericentre_r'] else float('nan'):.2f} kpc); bound fraction now: "
              + ", ".join(f"{sp} {v:.4f}" for sp, v in h["final"].items()))
    for label, o in record["tails"].items():
        print(f"{label} at t = {o['t_myr']:.0f} Myr (r_gal {o['r_gal_kpc']:.2f} kpc):")
        for sp, v in o["species"].items():
            print(f"   {sp}: unbound {v['unbound_mass']:.3e} (leading {v['leading_mass']:.2e}, trailing {v['trailing_mass']:.2e}); dt 5-95% {v['dt_percentiles_myr'][0]:+.0f}..{v['dt_percentiles_myr'][-1]:+.0f} Myr; median width {v['width_median_pc']:.0f} pc" if v["dt_percentiles_myr"] else f"   {sp}: unbound {v['unbound_mass']:.3e}")
    print(ROOT/"plots"/f"{args.name}_history.png")


if __name__ == "__main__":
    main()
