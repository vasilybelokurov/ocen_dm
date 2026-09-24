#!/usr/bin/env python3
"""Equilibrium check of an isolated N-body run against its DF model.

Compares the first and last snapshots with the analytic model profiles
(model_profiles.json written by build_ics.py): density per species, enclosed mass,
radial and tangential velocity dispersions and beta(r) of the stars, plus the time
series of energy error, virial ratio and half-mass radius from diagnostics.jsonl.
The model density is shell-averaged from the enclosed mass (comparing with the
density at the shell midpoint biases the ratio by ~3% for these steep profiles).
Acceptance (printed): |rho_end/rho_model - 1| < 5% for the stars at 0.3-30 pc in
shells with >= 2000 particles (shot noise < 2.2%), |2K/W - 1| < 2%, |dE/E| < 1e-3,
r_half drift < 2%. Remnant and halo ratios are reported, not gated (few particles).

Usage: python bin/nbody/analyse_isolation.py results/nbody/A_nodm/isolated [--name A_nodm_isolated]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SPECIES = ("stars", "remnants", "halo")
COLORS = dict(stars="tab:blue", remnants="tab:red", halo="tab:green")


def shells(pos, vel, mass, edges, centre, vcentre):
    r = np.linalg.norm(pos-centre, axis=1)
    idx = np.digitize(r, edges)-1
    n = len(edges)-1
    vol = 4/3*np.pi*(edges[1:]**3-edges[:-1]**3)
    count = np.bincount(idx[(idx >= 0) & (idx < n)], minlength=n)
    rho = np.bincount(idx[(idx >= 0) & (idx < n)], weights=mass[(idx >= 0) & (idx < n)], minlength=n)/vol
    rhat = (pos-centre)/np.maximum(r, 1e-12)[:, None]
    dv = vel-vcentre
    vr = np.sum(dv*rhat, axis=1); vt2 = np.sum(dv**2, axis=1)-vr**2
    ok = (idx >= 0) & (idx < n)
    s_r2 = np.bincount(idx[ok], weights=vr[ok]**2, minlength=n)/np.maximum(count, 1)
    s_t2 = np.bincount(idx[ok], weights=vt2[ok], minlength=n)/np.maximum(count, 1)/2
    return count, rho, np.sqrt(s_r2), np.sqrt(s_t2), 1-s_t2/np.maximum(s_r2, 1e-30)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run", type=Path)
    parser.add_argument("--name", default=None)
    args = parser.parse_args()
    name = args.name or "_".join(args.run.parts[-2:])
    ics_dir = args.run.parent
    ics = np.load(ics_dir/"ics.npz"); species = ics["species"]; mass = ics["mass"]
    model = json.loads((ics_dir/"model_profiles.json").read_text())
    snaps = sorted(glob.glob(str(args.run/"snap_*.npz")))
    first, last = np.load(snaps[0]), np.load(snaps[-1])
    diag = [json.loads(l) for l in open(args.run/"diagnostics.jsonl")]
    edges = np.geomspace(0.1, 300., 31); rmid = np.sqrt(edges[1:]*edges[:-1])
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    report = {}
    for snap, label, ls in ((first, f"t = {float(first['t_myr']):.0f} Myr", ":"), (last, f"t = {float(last['t_myr']):.0f} Myr", "-")):
        pos, vel = snap["pos"].astype(float), snap["vel"].astype(float)
        lum = species <= 1
        c = pos[lum].mean(axis=0)
        for R in (50., 20., 10., 5.):
            k = lum & (np.sum((pos-c)**2, axis=1) < R*R); c = np.average(pos[k], axis=0, weights=mass[k])
        k = lum & (np.sum((pos-c)**2, axis=1) < 100.); vc = np.average(vel[k], axis=0, weights=mass[k])
        for si, sp in enumerate(SPECIES):
            sel = species == si
            if not sel.any():
                continue
            count, rho, sr, st, beta = shells(pos[sel], vel[sel], mass[sel], edges, c, vc)
            m_enc = np.interp(edges, model["r_pc"], model["enclosed"][sp])
            rho_model = np.diff(m_enc)/(4/3*np.pi*(edges[1:]**3-edges[:-1]**3))     # shell-averaged model density
            good = count >= 500
            axes[0, 0].plot(rmid[good], rho[good], ls=ls, color=COLORS[sp], label=f"{sp}, {label}")
            axes[0, 1].plot(rmid[good], rho[good]/rho_model[good], ls=ls, color=COLORS[sp], marker="." if ls == "-" else None)
            if sp == "stars":
                axes[1, 0].plot(rmid[good], sr[good], ls=ls, color="tab:blue", label=f"$\\sigma_r$, {label}")
                axes[1, 0].plot(rmid[good], st[good], ls=ls, color="tab:orange", label=f"$\\sigma_t/\\sqrt{{2}}$, {label}")
                axes[1, 1].plot(rmid[good], beta[good], ls=ls, color="tab:blue", label=label)
            if ls == "-":
                inside = (count >= 2000) & (rmid > 0.3) & (rmid < 30.)
                report[sp] = dict(max_abs_density_dev=float(np.max(np.abs(rho[inside]/rho_model[inside]-1))) if inside.any() else None,
                                  r_pc=rmid[inside].tolist(), rho_over_model=(rho[inside]/rho_model[inside]).tolist())
    for sp in SPECIES:
        if sp in model["density"] and max(model["density"][sp]) > 0:
            axes[0, 0].plot(model["r_pc"], model["density"][sp], color=COLORS[sp], lw=.8, alpha=.5)
    axes[0, 0].set(xscale="log", yscale="log", xlim=(.1, 300), ylim=(1e-3, 3e4), xlabel="r [pc]", ylabel=r"$\rho$ [$M_\odot$ pc$^{-3}$]", title="(a) density (thin: DF model)")
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].axhline(1, color="black", lw=.7); [axes[0, 1].axhline(y, color="0.7", ls=":", lw=.7) for y in (.95, 1.05)]
    axes[0, 1].set(xscale="log", xlim=(.1, 300), ylim=(.7, 1.3), xlabel="r [pc]", ylabel="N-body / model", title="(b) density ratio to the shell-averaged model (shells with >= 500 particles)")
    axes[1, 0].set(xscale="log", xlim=(.1, 300), xlabel="r [pc]", ylabel="km/s", title="(c) stellar dispersions"); axes[1, 0].legend(fontsize=7)
    axes[1, 1].axhline(0, color="black", lw=.7); axes[1, 1].set(xscale="log", xlim=(.1, 300), ylim=(-1, 1), xlabel="r [pc]", ylabel=r"$\beta$", title="(d) stellar anisotropy"); axes[1, 1].legend(fontsize=7)
    t = [d["t_myr"] for d in diag]
    ax = axes[0, 2]
    ax.plot(t, [d["virial"] for d in diag], label="2K/|W|"); ax.set(xlabel="t [Myr]", ylabel="virial ratio", title="(e) virial ratio"); ax.axhline(1, color="black", lw=.7)
    ax2 = ax.twinx(); ax2.plot(t, [d["energy_error"] for d in diag], color="tab:red", label="dE/E"); ax2.set_ylabel(r"$\Delta E/|E|$", color="tab:red")
    ax = axes[1, 2]
    ax.plot(t, [d["r_half_stars_pc"] for d in diag], label="r_half stars"); ax.plot(t, [d["r_90_stars_pc"] for d in diag], label="r_90 stars")
    ax.set(xlabel="t [Myr]", ylabel="pc", title="(f) Lagrangian radii of bound stars"); ax.legend(fontsize=7)
    fig.suptitle(f"Isolated equilibrium test: {name} (N = {len(mass)}, eps = {json.loads((args.run/'run.json').read_text())['eps_pc']} pc, dt = {json.loads((args.run/'run.json').read_text())['dt_myr']} Myr)")
    plot = ROOT/"plots"/f"nbody_isolation_{name}.png"
    fig.savefig(plot, dpi=140)
    vir = np.array([d["virial"] for d in diag]); de = np.array([d["energy_error"] for d in diag]); rh = np.array([d["r_half_stars_pc"] for d in diag])
    summary = dict(name=name, t_end_myr=t[-1], density=report, virial_range=[float(vir.min()), float(vir.max())], max_energy_error=float(np.max(np.abs(de))),
                   r_half_drift=float(rh[-1]/rh[0]-1), passed=bool((report["stars"]["max_abs_density_dev"] or 0.) < .05
                                                                 and np.max(np.abs(vir-1)) < .02 and np.max(np.abs(de)) < 1e-3 and abs(rh[-1]/rh[0]-1) < .02))
    (args.run/"isolation_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"{name}: t_end {t[-1]:.0f} Myr; max |rho/model-1| (0.3-30 pc): " + ", ".join(f"{k} {v['max_abs_density_dev']:.3f}" for k, v in report.items() if v["max_abs_density_dev"] is not None)
          + f"; virial {vir.min():.4f}-{vir.max():.4f}; max |dE/E| {np.max(np.abs(de)):.1e}; r_half drift {100*(rh[-1]/rh[0]-1):+.2f}%  -> {'PASS' if summary['passed'] else 'FAIL'}")
    print(plot)


if __name__ == "__main__":
    main()
