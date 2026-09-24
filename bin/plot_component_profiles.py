#!/usr/bin/env python3
"""Radial distribution of the mass components (stars, remnants, DM) for selected fits, side by side.

Each model is rebuilt at refined resolution. Top row: density rho(r) of each component
and the total; bottom row: enclosed mass M(<r) with the mass fractions annotated at a
few radii. The data range (0.1-63 pc) is marked; the region inside the innermost
kinematic bin and beyond the outermost is shaded.

Usage: python bin/plot_component_profiles.py --name component_profiles_20260924 \
    "results/df/dftwo_counts_20260923::observed_no_halo_start0::no DM, free stellar mass" \
    "results/df/mpsscan_9_20260923::observed_free_halo_start0::DM halo, photometric stellar mass"
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ocen_dm.kinematics.compact_recovery import mass_profiles, refined_config
from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict

COLORS = dict(stars="tab:blue", remnants="tab:red", halo="tab:green", total="black")
LABELS = dict(stars="stars (follow the light)", remnants="dark remnants (Plummer)", halo="DM halo (cored)", total="total")
MARK_R = (1., 3., 10., 30.)


def read(path):
    return json.loads(Path(path).read_text())


def densities(model, r):
    xyz = np.column_stack((r, 0*r, 0*r))
    stars = model.stellar_potential.density(xyz)
    m = model.config.matter
    remnants = 3*m.M_rem/(4*np.pi*m.a_rem**3)*(1+(r/m.a_rem)**2)**-2.5
    halo = np.asarray(m.halo_density(r), float)
    return dict(stars=np.asarray(stars, float), remnants=remnants, halo=halo, total=np.asarray(stars, float)+remnants+halo)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--name", required=True)
    parser.add_argument("specs", nargs="+", help="BATCH::JOB::LABEL")
    args = parser.parse_args()
    r = np.geomspace(.05, 100., 200)
    n = len(args.specs)
    fig, axes = plt.subplots(2, n, figsize=(6.5*n, 9.5), sharex=True, constrained_layout=True)
    axes = np.atleast_2d(axes.T).T if n > 1 else axes.reshape(2, 1)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), models={})
    for i, spec in enumerate(args.specs):
        batch, job, label = spec.split("::")
        s = read(ROOT/batch/"fits"/job/"summary.json")
        config = refined_config(model_config_from_dict(s["best"]["config"]))
        model = build_df_model(config)
        rho = densities(model, r); enc = mass_profiles(model, r)
        c = config
        ax = axes[0, i]
        for key in ("stars", "remnants", "halo", "total"):
            if key == "halo" and c.matter.rho20 <= 0:
                continue
            ax.plot(r, rho[key], color=COLORS[key], lw=2.2 if key == "total" else 1.5, ls="-" if key != "total" else "--", label=LABELS[key])
        ax.set(xscale="log", yscale="log", ylim=(1e-2, 3e5), ylabel=r"$\rho(r)$ [$M_\odot$ pc$^{-3}$]",
               title=f"{label}\nobjective {s['refined_score']:.1f}; $M_\\star$ {c.M_star:.2e}, $M_{{\\rm rem}}$ {c.matter.M_rem:.2e} ($a$ {c.matter.a_rem:.1f} pc)"
                     + (f", $\\rho_{{20}}$ {c.matter.rho20:.2g}, $r_s$ {c.matter.r_s:.0f} pc" if c.matter.rho20 > 0 else ", no halo"))
        ax.title.set_fontsize(10)
        ax.legend(fontsize=8, loc="lower left")
        ax = axes[1, i]
        for key in ("stars", "remnants", "halo", "total"):
            m = np.array(enc[key])
            if key == "halo" and c.matter.rho20 <= 0:
                continue
            ax.plot(r, np.maximum(m, 1.), color=COLORS[key], lw=2.2 if key == "total" else 1.5, ls="-" if key != "total" else "--", label=LABELS[key])
        tot = np.array(enc["total"])
        text = []
        for rr in MARK_R:
            f = {k: float(np.interp(rr, r, enc[k]))/float(np.interp(rr, r, tot)) for k in ("stars", "remnants", "halo")}
            text.append(f"r < {rr:g} pc: stars {100*f['stars']:.0f}%  remnants {100*f['remnants']:.0f}%" + (f"  DM {100*f['halo']:.0f}%" if c.matter.rho20 > 0 else ""))
            ax.axvline(rr, color="0.8", lw=.6, zorder=0)
        ax.text(.03, .97, "\n".join(text), transform=ax.transAxes, va="top", fontsize=8, family="monospace",
                bbox=dict(facecolor="white", alpha=.8, edgecolor="none"))
        ax.set(xscale="log", yscale="log", ylim=(1e2, 6e6), xlabel="r [pc]", ylabel=r"$M(<r)$ [$M_\odot$]")
        ax.legend(fontsize=8, loc="lower right")
        for a in axes[:, i]:
            a.axvspan(.05, .1, color="0.92", zorder=0); a.axvspan(63., 100., color="0.92", zorder=0)
        record["models"][label] = dict(spec=spec, config=config.to_dict(), r_pc=r.tolist(), density={k: v.tolist() for k, v in rho.items()},
                                       enclosed={k: enc[k] for k in ("stars", "remnants", "halo", "total")}, objective=s["refined_score"])
        print(label, "|", " ; ".join(text), flush=True)
    fig.suptitle("Mass components of the best-fitting models without and with dark matter (grey: outside the kinematic data)")
    plot = ROOT/"plots"/f"{args.name}.png"
    fig.savefig(plot, dpi=150)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/f"{args.name}.json").write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
