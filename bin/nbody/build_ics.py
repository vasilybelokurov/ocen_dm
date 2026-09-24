#!/usr/bin/env python3
"""Equilibrium N-body initial conditions for a fitted omega Cen DF model.

Stars are sampled from the fitted action-based DF (AGAMA GalaxyModel of the refined
self-consistent model); the dark remnants (Plummer) and the DM halo (cored, tapered)
are sampled from isotropic Eddington-inversion DFs ('QuasiSpherical') in the same
total potential, so every component is in equilibrium in the combined potential.
All particles have the same mass m_p = M_total/N (no inter-species mass segregation
from unequal particle masses). Units: pc, km/s, Msun. The model is non-rotating
(the fitted DF is the even part in L_z) and the MW tide is not included here.

Output <out>/ics.npz: pos (N,3) pc, vel (N,3) km/s, mass (N) Msun, species (N) int8
(0 stars, 1 remnants, 2 halo); <out>/ics.json: model, masses, counts, checks
(cumulative mass of the sample vs the model per species, virial ratio with the
softened self-gravity); <out>/model_profiles.json: analytic density and enclosed
mass per component on a radial grid, for later comparison with snapshots.

Usage: python bin/nbody/build_ics.py --model "results/df/dftwo_counts_20260923::observed_no_halo_start0" \
    --out results/nbody/A_nodm --n 500000 --seed 1
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "4")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np

from ocen_dm.kinematics.compact_recovery import mass_profiles, refined_config
from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict
from ocen_dm.kinematics.positive_df import agama_pc

SPECIES = ("stars", "remnants", "halo")
ACTION_FLOOR = 1e-6      # pc km/s; AGAMA's sampler probes actions ~1e-14 where the spherical frequency
                         # ratio is undefined. J0 ~ 50, so a floor of 1e-6 changes no resolvable phase space.


class FlooredDF:
    """The fitted DF evaluated at max(J, ACTION_FLOOR) for Jr and L (sampling robustness only)."""

    def __init__(self, df):
        self.df = df

    def __call__(self, actions):
        a = np.array(actions, float, copy=True)
        a[..., 0] = np.maximum(a[..., 0], ACTION_FLOOR)
        a[..., 1] = np.maximum(a[..., 1], ACTION_FLOOR)
        return self.df(a)


def read(path):
    return json.loads(Path(path).read_text())


def enclosed(potential, R, G):
    R = np.atleast_1d(np.asarray(R, float))
    return -potential.force(np.column_stack((R, 0*R, 0*R)))[:, 0]*R**2/G


def component_density(model, r):
    xyz = np.column_stack((r, 0*r, 0*r))
    m = model.config.matter
    out = dict(stars=np.asarray(model.stellar_potential.density(xyz), float),
               remnants=3*m.M_rem/(4*np.pi*m.a_rem**3)*(1+(r/m.a_rem)**2)**-2.5,
               halo=np.asarray(m.halo_density(r), float))
    out["total"] = out["stars"]+out["remnants"]+out["halo"]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--model", required=True, help="BATCH::JOB")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--n", type=int, default=500000, help="total particle number")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--eps", type=float, default=0.3, help="softening [pc] for the virial check only")
    args = parser.parse_args()
    agama = agama_pc()
    G = agama.G
    np.random.seed(args.seed)
    if hasattr(agama, "setRandomSeed"):
        agama.setRandomSeed(args.seed)
    batch, job = args.model.split("::")
    summary = read(ROOT/batch/"fits"/job/"summary.json")
    config = refined_config(model_config_from_dict(summary["best"]["config"]))
    t0 = time.time()
    model = build_df_model(config)
    print(f"model built in {time.time()-t0:.0f} s", flush=True)
    matter = config.matter
    components = {"stars": (agama.GalaxyModel(model.potential, FlooredDF(model.df), model.af), model.stellar_potential)}
    statics = list(model.static_potentials)
    if matter.M_rem:
        rem = statics.pop(0)
        components["remnants"] = (agama.GalaxyModel(model.potential, agama.DistributionFunction(
            type="QuasiSpherical", potential=model.potential, density=rem)), rem)
    if matter.rho20:
        halo = statics.pop(0)
        components["halo"] = (agama.GalaxyModel(model.potential, agama.DistributionFunction(
            type="QuasiSpherical", potential=model.potential, density=halo)), halo)
    masses = {k: float(v[1].totalMass()) for k, v in components.items()}
    total = sum(masses.values())
    m_p = total/args.n
    counts = {k: int(round(m/m_p)) for k, m in masses.items()}
    print("component masses:", {k: f"{v:.4e}" for k, v in masses.items()}, f"total {total:.4e}; particle mass {m_p:.3f} Msun; counts {counts}", flush=True)
    pos, vel, species, checks = [], [], [], {}
    radii_check = np.array([0.3, 1., 3., 10., 30., 60., 160., 500., 1000.])
    for si, name in enumerate(SPECIES):
        if name not in components:
            continue
        gm, pot = components[name]
        xv, _ = gm.sample(counts[name])
        r = np.linalg.norm(xv[:, :3], axis=1)
        model_m = enclosed(pot, radii_check, G)
        sampled = np.array([m_p*np.sum(r < R) for R in radii_check])
        E = 0.5*np.sum(xv[:, 3:]**2, axis=1)+model.potential.potential(xv[:, :3])
        checks[name] = dict(n=counts[name], mass=masses[name], r_pc=radii_check.tolist(), sampled_over_model=(sampled/model_m).tolist(),
                            unbound_fraction=float(np.mean(E > 0)), r_max_pc=float(r.max()))
        print(f"  {name}: {counts[name]} particles; M(<r) sampled/model {np.round(sampled/model_m, 3).tolist()} at {radii_check.tolist()} pc; unbound {np.mean(E > 0):.4f}", flush=True)
        pos.append(xv[:, :3]); vel.append(xv[:, 3:]); species.append(np.full(counts[name], si, np.int8))
    pos = np.vstack(pos); vel = np.vstack(vel); species = np.concatenate(species)
    mass = np.full(len(pos), m_p)
    # centre of mass frame (the samples are already centred to O(1/sqrt N))
    com = pos.mean(axis=0); cov = vel.mean(axis=0)
    pos -= com; vel -= cov
    # virial check with the softened self-gravity actually used in the simulation
    import pyfalcon
    acc, phi = pyfalcon.gravity(pos.astype(np.float32), mass.astype(np.float32), args.eps)
    W = 0.5*G*np.sum(mass*phi); K = 0.5*np.sum(mass*np.sum(vel**2, axis=1))
    virial = 2*K/abs(W)
    print(f"virial ratio 2K/|W| = {virial:.4f} (softening {args.eps} pc); centre-of-mass shift removed {np.linalg.norm(com):.3f} pc, {np.linalg.norm(cov):.3f} km/s", flush=True)
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez(args.out/"ics.npz", pos=pos, vel=vel, mass=mass, species=species)
    r = np.geomspace(0.05, 3000., 300)
    dens = component_density(model, r); enc = mass_profiles(model, r)
    (args.out/"model_profiles.json").write_text(json.dumps(dict(r_pc=r.tolist(), density={k: v.tolist() for k, v in dens.items()},
                                                                   enclosed={k: enc[k] for k in ("stars", "remnants", "halo", "total")})))
    meta = dict(created_utc=datetime.now(timezone.utc).isoformat(), model=args.model, source_commit=subprocess.check_output(
                    ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                config=config.to_dict(), objective=summary["refined_score"], n=int(len(pos)), particle_mass=m_p, masses=masses, counts=counts,
                seed=args.seed, units="pc, km/s, Msun", G=G, checks=checks, virial_ratio=float(virial), virial_softening_pc=args.eps,
                rotating=False, action_floor=ACTION_FLOOR,
                note="non-rotating even DF; remnants and halo isotropic (Eddington); halo as fitted incl. taper at r_t")
    (args.out/"ics.json").write_text(json.dumps(meta, indent=1, default=float))
    print(f"wrote {args.out/'ics.npz'} ({len(pos)} particles) in {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
