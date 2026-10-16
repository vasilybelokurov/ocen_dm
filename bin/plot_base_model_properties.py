#!/usr/bin/env python3
"""Physical properties of a fitted base model (stars + remnants, optional halo and central black hole), by component.

Figure 1 (``<name>_components.png``):
  (a) density rho(r) of stars, remnants (and halo) and the total;
  (b) enclosed mass M(<r) by component;
  (c) circular velocity v_c(r) = sqrt(G M(<r)/r) by component and total;
  (d) remnant share: local rho_rem/rho_tot and cumulative M_rem(<r)/M_tot(<r);
  (e) dynamical-to-stellar mass ratio M_tot(<r)/M_star(<r).
Figure 2 (``<name>_kinematics.png``):
  (a) intrinsic dispersions sigma_r, sigma_t (one tangential component) and the
      local escape speed / sqrt(12) for scale; (b) anisotropy beta(r);
  (c) projected stellar surface density Sigma(R) and the half-mass radii;
  (d) projected dispersions (LOS, PM radial, PM tangential) in km/s from the DF,
      before the published-rotation subtraction used in the likelihood.
The radial range of the kinematic data (0.1-63 pc) is marked. Numbers plotted go to
results/plot_data/<name>.json.

Usage: python bin/plot_base_model_properties.py --name base_model_20261013 \
    results/df/dftwo_counts_area_20261013 observed_no_halo_start3
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ocen_dm.kinematics.compact_recovery import mass_profiles, refined_config
from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict
from ocen_dm.kinematics.positive_df import agama_pc

DATA_RANGE = (0.1, 63.)          # pc, innermost/outermost kinematic bin centres
COLORS = dict(stars="C0", remnants="C3", halo="C2", bh="C1", total="k")


def half_radius(r, m):
    """Radius enclosing half of m[-1] (log interpolation)."""
    return float(np.exp(np.interp(0.5*m[-1], m, np.log(r))))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("batch")
    p.add_argument("job")
    p.add_argument("--name", required=True)
    a = p.parse_args()
    s = json.loads((Path(a.batch)/"fits"/a.job/"summary.json").read_text())
    cfg = refined_config(model_config_from_dict(s["best"]["config"]))
    model = build_df_model(cfg)
    G = agama_pc().G
    mat = cfg.matter

    rmax = cfg.numerics.r_max/100*0.999            # DF moment tables are valid to r_max/100
    r0 = 0.01 if mat.M_bh > 0 else 0.03          # resolve the black hole's sphere of influence
    r = np.geomspace(r0, rmax, 240)
    mp = {k: np.asarray(v) for k, v in mass_profiles(model, r).items() if isinstance(v, list)}
    mom = model.intrinsic_moments(r)
    rho = dict(stars=mom["rho"], remnants=3*mat.M_rem/(4*np.pi*mat.a_rem**3)*(1+(r/mat.a_rem)**2)**-2.5,
               halo=mat.halo_density(r))
    rho["total"] = rho["stars"]+rho["remnants"]+rho["halo"]
    comps = ["stars", "remnants"]+(["halo"] if mat.rho20 > 0 else [])
    mcomps = comps+(["bh"] if mat.M_bh > 0 else [])    # point mass: enclosed mass and v_c only
    vc = {k: np.sqrt(G*mp[k]/r) for k in mcomps+["total"]}
    sig_r = np.sqrt(mom["radial_pressure"]/mom["rho"])
    sig_t = np.sqrt(mom["tangential_pressure"]/mom["rho"])
    xyz = np.column_stack((r, 0*r, 0*r))
    vesc = np.sqrt(-2*model.potential.potential(xyz))
    R = np.geomspace(r0, rmax, 160)
    pm = model.projected_moments(R)
    sig = {k: np.sqrt(pm[k]/pm["Sigma"]) for k in ("los", "pmr", "pmt")}
    r_half = dict(stars=half_radius(r, mp["stars"]), remnants=1.305*mat.a_rem, total=half_radius(r, mp["total"]))
    cum = np.cumsum(pm["Sigma"]*2*np.pi*R*np.gradient(R))
    R_half = float(np.exp(np.interp(0.5*cum[-1], cum, np.log(R))))

    def shade(ax):
        ax.axvspan(r[0], DATA_RANGE[0], color="0.92", lw=0)
        ax.axvspan(DATA_RANGE[1], r[-1], color="0.92", lw=0)
        ax.set_xscale("log")
        ax.set_xlim(r[0], r[-1])

    title = (f"{'Model' if mat.M_bh > 0 else 'Base model'} ({Path(a.batch).name}, {a.job}): M* = {mp['stars'][-1]:.3g}, "
             + (f"M_BH = {mat.M_bh:.3g}, " if mat.M_bh > 0 else "")
             + f"M_rem = {mat.M_rem:.3g} Msun (a = {mat.a_rem:.2f} pc), D = {cfg.distance_kpc:.3f} kpc")
    # ---------------------------------------------------------------- figure 1
    fig, ax = plt.subplots(1, 5, figsize=(22, 4.6), constrained_layout=True)
    for k in comps+["total"]:
        ax[0].loglog(r, rho[k], color=COLORS[k], lw=2 if k == "total" else 1.4, label=k)
    for k in mcomps+["total"]:
        ax[1].loglog(r, mp[k], color=COLORS[k], lw=2 if k == "total" else 1.4, label=k)
        ax[2].plot(r, vc[k], color=COLORS[k], lw=2 if k == "total" else 1.4, label=k)
    for k, v in r_half.items():
        ax[1].axvline(v, color=COLORS[k], ls=":", lw=1)
    ax[0].set(ylabel=r"$\rho$ [M$_\odot$ pc$^{-3}$]", title="(a) density")
    ax[1].set(ylabel=r"$M(<r)$ [M$_\odot$]", title="(b) enclosed mass (dotted: half-mass radii)")
    ax[2].set(ylabel=r"$v_c$ [km/s]", title="(c) circular velocity")
    ax[3].plot(r, rho["remnants"]/rho["total"], color="C3", label=r"local $\rho_{\rm rem}/\rho_{\rm tot}$")
    ax[3].plot(r, mp["remnants"]/mp["total"], color="C3", ls="--", label=r"cumulative $M_{\rm rem}/M_{\rm tot}$")
    ax[3].set(ylabel="fraction", ylim=(0, 1), title="(d) remnant share")
    ax[4].plot(r, mp["total"]/mp["stars"], color="k")
    ax[4].set(ylabel=r"$M_{\rm tot}(<r)/M_\star(<r)$", title="(e) dynamical / stellar mass")
    for x in ax:
        shade(x)
        x.set_xlabel("r [pc]")
    for x in ax[:4]:
        x.legend(fontsize=8)
    fig.suptitle(title, fontsize=10)
    out1 = ROOT/"plots"/f"{a.name}_components.png"
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    # ---------------------------------------------------------------- figure 2
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.6), constrained_layout=True)
    ax[0].plot(r, sig_r, color="C1", label=r"$\sigma_r$")
    ax[0].plot(r, sig_t, color="C4", label=r"$\sigma_t$ (one component)")
    ax[0].plot(r, vesc/np.sqrt(12), color="0.5", ls=":", label=r"$v_{\rm esc}/\sqrt{12}$")
    ax[0].set(ylabel="km/s", title="(a) intrinsic dispersions")
    ax[1].plot(r, mom["beta"], color="k")
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].set(ylabel=r"$\beta = 1-\sigma_t^2/\sigma_r^2$", title="(b) anisotropy")
    ax[2].loglog(R, pm["Sigma"], color="C0")
    ax[2].axvline(R_half, color="C0", ls=":", lw=1, label=f"projected half-mass {R_half:.1f} pc")
    ax[2].axvline(r_half["stars"], color="C0", ls="--", lw=1, label=f"3D half-mass {r_half['stars']:.1f} pc")
    ax[2].set(ylabel=r"$\Sigma_\star$ [M$_\odot$ pc$^{-2}$]", title="(c) stellar surface density")
    for k, lab, c in (("los", "line of sight", "C2"), ("pmr", "PM radial", "C1"), ("pmt", "PM tangential", "C4")):
        ax[3].plot(R, sig[k], color=c, label=lab)
    ax[3].set(ylabel="km/s", xlabel="R [pc]", title="(d) projected dispersions (no rotation term)")
    for x in ax:
        shade(x)
        x.legend(fontsize=8) if x.get_legend_handles_labels()[0] else None
    for x in ax[:2]:
        x.set_xlabel("r [pc]")
    ax[2].set_xlabel("R [pc]")
    fig.suptitle(title, fontsize=10)
    out2 = ROOT/"plots"/f"{a.name}_kinematics.png"
    fig.savefig(out2, dpi=150)
    plt.close(fig)
    # ---------------------------------------------------------------- numbers
    at = np.array([0.3, 1., 3., 10., 30., 100.])
    pick = lambda y: [float(np.exp(np.interp(np.log(x), np.log(r), np.log(y)))) for x in at]
    rec = dict(batch=a.batch, job=a.job, r_pc=at.tolist(),
               M_enclosed={k: pick(mp[k]) for k in mcomps+["total"]}, M_bh=float(mat.M_bh),
               rho={k: pick(rho[k]) for k in comps+["total"]},
               v_c_total=pick(vc["total"]), sigma_r=pick(sig_r), sigma_t=pick(sig_t),
               beta=[float(np.interp(np.log(x), np.log(r), mom["beta"])) for x in at],
               half_mass_radii_pc=dict(r_half, projected_stars=R_half),
               v_c_max=dict(value=float(vc["total"].max()), r_pc=float(r[np.argmax(vc["total"])])),
               central_escape_speed_kms=float(vesc[0]))
    out3 = ROOT/"results"/"plot_data"/f"{a.name}.json"
    out3.parent.mkdir(parents=True, exist_ok=True)
    out3.write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))
    print(out1, out2, out3)


if __name__ == "__main__":
    main()
