#!/usr/bin/env python3
"""Self-gravitating N-body evolution of an omega Cen model, isolated or on its Galactic orbit.

Leapfrog (kick-drift-kick, shared fixed step) with falcON self-gravity (pyfalcon,
Plummer softening) plus, with --orbit, the static McMillan (2017) Milky Way potential
from AGAMA. In orbit mode the cluster is placed at omega Cen's position T_BACK Myr ago,
obtained by integrating the present-day phase-space coordinates (Baumgardt catalogue)
backwards in the same potential, so that it arrives at today's position at t = T_BACK.

Units: pc, km/s, Msun; time in Myr (internal time unit pc/(km/s) = 0.97779 Myr).
Diagnostics every --diag Myr: cluster centre (shrinking spheres on stars+remnants),
bound mass per species (energy in the cluster's own softened potential, one iteration),
half-mass radius of bound stars, total energy (isolated) or Galactocentric radius;
appended to <out>/diagnostics.jsonl. Snapshots every --snap Myr as float32 npz
(<out>/snap_XXXX.npz, positions and velocities in the Galactocentric frame for orbit
mode). A restart file (float64) is written with each snapshot; --resume continues.

Usage:
  python bin/nbody/run_nbody.py --ics results/nbody/A_nodm/ics.npz --out results/nbody/A_nodm/isolated \
      --eps 0.3 --dt 0.01 --tstop 20 --snap 10 --diag 1
  python bin/nbody/run_nbody.py --ics results/nbody/A_nodm/ics.npz --out results/nbody/A_nodm/orbit \
      --orbit --tback 2000 --eps 0.3 --dt 0.01 --tstop 2000 --snap 10 --diag 1
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import glob
import json
import os
from pathlib import Path
import sys
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
import numpy as np
import pyfalcon

ROOT = Path(__file__).resolve().parents[2]
G = 4.30091727067736e-3            # pc (km/s)^2 / Msun
T_UNIT = 0.977792221680356         # Myr per (pc / (km/s))
SPECIES = ("stars", "remnants", "halo")
# omega Cen today, Galactocentric Cartesian (Baumgardt catalogue orbits_table.txt): X Y Z [kpc], U V W [km/s]
OCEN_TODAY = np.array([4.87, -4.07, 1.40, -95.98, -21.69, -86.82])


def mw_potential(name):
    import agama
    agama.setUnits(mass=1, length=1, velocity=1)      # Msun, kpc, km/s (this process only uses AGAMA for the MW)
    ini = glob.glob(os.path.dirname(agama.__file__)+f"/**/{name}.ini", recursive=True)
    if not ini:
        raise FileNotFoundError(f"AGAMA potential file {name}.ini not found")
    return agama, agama.Potential(ini[0])


def mw_acceleration(pot, pos_pc):
    """Acceleration of the MW potential at positions in pc, returned in (km/s)^2/pc."""
    return pot.force(pos_pc*1e-3)*1e-3


def cluster_centre(pos, vel, mass, sel, guess, radii=(50., 20., 10., 5.)):
    """Shrinking spheres on the selected (luminous) particles from ``guess``; if the guess has
    lost the cluster (fewer than 100 particles within the first sphere) restart from the
    median position of the selected particles, which unbound debris shifts only slightly."""
    c = np.asarray(guess, float).copy()
    if np.sum(sel & (np.sum((pos-c)**2, axis=1) < radii[0]**2)) < 100:
        c = np.median(pos[sel], axis=0)
        print("  centre finder: guess lost the cluster, restarted from the luminous median", flush=True)
    for R in radii:
        k = sel & (np.sum((pos-c)**2, axis=1) < R*R)
        if k.sum() < 100:
            break
        c = np.average(pos[k], axis=0, weights=mass[k])
    k = sel & (np.sum((pos-c)**2, axis=1) < 10.**2)
    vc = np.average(vel[k], axis=0, weights=mass[k]) if k.sum() >= 100 else vel[sel].mean(axis=0)
    return c, vc


def bound_mask(pos, vel, mass, eps, centre, vcentre, iterations=2):
    bound = np.ones(len(pos), bool)
    for _ in range(iterations):
        m_eff = np.where(bound, mass, 1e-30*mass)      # unbound particles do not source the cluster potential
        _, phi = pyfalcon.gravity(pos.astype(np.float32), m_eff.astype(np.float32), eps)
        E = 0.5*np.sum((vel-vcentre)**2, axis=1)+G*phi
        bound = E < 0
    return bound


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ics", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--eps", type=float, default=0.3, help="Plummer softening [pc]")
    parser.add_argument("--dt", type=float, default=0.01, help="time step [Myr]")
    parser.add_argument("--tstop", type=float, required=True, help="[Myr]")
    parser.add_argument("--snap", type=float, default=10., help="snapshot interval [Myr]")
    parser.add_argument("--diag", type=float, default=1., help="diagnostics interval [Myr]")
    parser.add_argument("--orbit", action="store_true", help="include the MW potential and start on omega Cen's orbit")
    parser.add_argument("--mw", default="McMillan17")
    parser.add_argument("--tback", type=float, default=2000., help="start this many Myr before today (orbit mode)")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--nthreads", type=int, default=1)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    ics = np.load(args.ics)
    species = ics["species"]
    mass = ics["mass"].astype(np.float64)
    dt = args.dt/T_UNIT                                  # internal time unit
    restart = args.out/"restart.npz"
    if args.resume and restart.exists():
        st = np.load(restart)
        pos, vel, t_myr, step = st["pos"], st["vel"], float(st["t_myr"]), int(st["step"])
        print(f"resuming at t = {t_myr:.2f} Myr (step {step})", flush=True)
    else:
        pos, vel = ics["pos"].astype(np.float64), ics["vel"].astype(np.float64)
        t_myr, step = 0., 0
    pot = None
    if args.orbit:
        agama, pot = mw_potential(args.mw)
        if not (args.resume and restart.exists()):
            traj_t, traj = agama.orbit(potential=pot, ic=OCEN_TODAY, time=-args.tback*1e-3, trajsize=2)
            start = traj[-1]                             # kpc, km/s
            pos += start[:3]*1e3; vel += start[3:]
            print(f"orbit start {args.tback:.0f} Myr ago: X,Y,Z = {np.round(start[:3], 3).tolist()} kpc, V = {np.round(start[3:], 2).tolist()} km/s", flush=True)
    meta = dict(created_utc=datetime.now(timezone.utc).isoformat(), ics=str(args.ics), eps_pc=args.eps, dt_myr=args.dt, tstop_myr=args.tstop,
                orbit=args.orbit, mw=args.mw if args.orbit else None, tback_myr=args.tback if args.orbit else None, n=int(len(pos)),
                ocen_today=OCEN_TODAY.tolist(), G=G, t_unit_myr=T_UNIT)
    (args.out/"run.json").write_text(json.dumps(meta, indent=1))
    sel_lum = species <= 1                              # stars + remnants define the cluster centre
    centre = pos[sel_lum].mean(axis=0); vcentre = vel[sel_lum].mean(axis=0)
    pos32 = pos.astype(np.float32); m32 = mass.astype(np.float32)

    def accel(pos):
        acc, phi = pyfalcon.gravity(pos.astype(np.float32), m32, args.eps)
        acc = G*acc
        if pot is not None:
            acc += mw_acceleration(pot, pos)
        return acc, G*phi

    acc, phi = accel(pos)
    n_steps = int(round(args.tstop/args.dt))
    snap_every = max(1, int(round(args.snap/args.dt))); diag_every = max(1, int(round(args.diag/args.dt)))
    diag_file = open(args.out/"diagnostics.jsonl", "a")
    t_wall = time.time(); t_last = t_wall
    E0 = None

    last_diag_t = [t_myr]

    def diagnostics(step, t_myr, pos, vel, phi, write=True):
        nonlocal centre, vcentre, E0
        # predict the centre from its velocity over the interval since the last diagnostics
        guess = centre+vcentre*(t_myr-last_diag_t[0])/T_UNIT
        last_diag_t[0] = t_myr
        centre, vcentre = cluster_centre(pos, vel, mass, sel_lum, guess)
        bound = bound_mask(pos, vel, mass, args.eps, centre, vcentre)
        row = dict(step=step, t_myr=t_myr, centre_pc=centre.tolist(), vcentre_kms=vcentre.tolist(),
                   bound_mass={SPECIES[s]: float(mass[bound & (species == s)].sum()) for s in np.unique(species)},
                   n_bound=int(bound.sum()))
        rs = np.sort(np.linalg.norm(pos[bound & (species == 0)]-centre, axis=1))
        if len(rs):
            row["r_half_stars_pc"] = float(rs[len(rs)//2]); row["r_90_stars_pc"] = float(rs[int(.9*len(rs))])
        if pot is None:
            K = 0.5*np.sum(mass*np.sum(vel**2, axis=1)); W = 0.5*np.sum(mass*phi)
            E = K+W; E0 = E if E0 is None else E0
            row.update(kinetic=float(K), potential=float(W), energy=float(E), energy_error=float((E-E0)/abs(E0)), virial=float(2*K/abs(W)))
        else:
            row["r_gal_kpc"] = float(np.linalg.norm(centre)*1e-3)
        if write:
            diag_file.write(json.dumps(row)+"\n"); diag_file.flush()
        return row

    def snapshot(step, t_myr, pos, vel):
        np.savez(args.out/f"snap_{step//snap_every:04d}.npz", t_myr=t_myr, pos=pos.astype(np.float32), vel=vel.astype(np.float32))
        np.savez(restart, pos=pos, vel=vel, t_myr=t_myr, step=step)

    if step == 0:
        row = diagnostics(0, 0., pos, vel, phi)
        snapshot(0, 0., pos, vel)
        print(f"t = 0: bound {row['bound_mass']}, r_half(stars) {row.get('r_half_stars_pc', float('nan')):.2f} pc" + (f", virial {row['virial']:.4f}" if pot is None else ""), flush=True)
    while step < n_steps:
        vel += acc*(dt/2)
        pos += vel*dt
        acc, phi = accel(pos)
        vel += acc*(dt/2)
        step += 1; t_myr = step*args.dt
        if step % diag_every == 0:
            row = diagnostics(step, t_myr, pos, vel, phi)
        if step % snap_every == 0:
            snapshot(step, t_myr, pos, vel)
            now = time.time(); rate = (now-t_last)/snap_every; t_last = now
            eta = rate*(n_steps-step)/3600
            extra = f", energy error {row['energy_error']:+.2e}, virial {row['virial']:.4f}" if pot is None else f", r_gal {row['r_gal_kpc']:.2f} kpc"
            print(f"t = {t_myr:8.2f} Myr  bound: " + ", ".join(f"{k} {v:.3e}" for k, v in row["bound_mass"].items())
                  + f", r_half(stars) {row.get('r_half_stars_pc', float('nan')):.2f} pc{extra}  [{rate:.2f} s/step, ETA {eta:.1f} h]", flush=True)
    diag_file.close()
    (args.out/"run.json").write_text(json.dumps(dict(meta, finished_utc=datetime.now(timezone.utc).isoformat(), steps=step,
                                                   wall_hours=(time.time()-t_wall)/3600), indent=1))
    print("done")


if __name__ == "__main__":
    main()
