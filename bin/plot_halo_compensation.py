#!/usr/bin/env python3
"""Why the DM and no-DM fits coincide: a halo added naively vs a halo with everything refitted.

For each kinematic dataset and the count profiles, plots (model - no-DM best fit)/error
per bin for (i) the no-DM model with the rho20 = 1, r_s = 5 pc halo simply added and
nothing else changed, and (ii) the fit with the same halo where M_star, M_rem, the DF
shape and the distance were refitted. The first is what a halo 'should' do to the
dispersions; the second is what remains after the other components compensate.

Usage: python bin/plot_halo_compensation.py \
    --nodm "results/df/dftwo_counts_20260923::observed_no_halo_start0" \
    --dm "results/df/rho20scan_1_20260923::observed_free_halo_start0"
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
sys.path.insert(0, str(ROOT/"bin"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict
import run_compact_df_recovery as drv

PC = 5.43e3*np.pi/(180*3600)
TITLES = dict(hst_pm_radial_ours="HST PM radial", hst_pm_tangential_ours="HST PM tangential", muse_los_dispersion="MUSE LOS",
              gaia_edr3_ours_radial="Gaia PM radial", gaia_edr3_ours_tangential="Gaia PM tangential",
              ocen_counts_hst_f625w19="HST counts", ocen_counts_gaia_g17="Gaia counts")


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--nodm", default="results/df/dftwo_counts_20260923::observed_no_halo_start0")
    parser.add_argument("--dm", default="results/df/rho20scan_1_20260923::observed_free_halo_start0")
    parser.add_argument("--name", default="halo_compensation_20260924")
    args = parser.parse_args()
    b0, j0 = args.nodm.split("::"); b1, j1 = args.dm.split("::")
    problem = drv.load_problem(ROOT/b0/"observed")
    c0 = model_config_from_dict(read(ROOT/b0/"fits"/j0/"summary.json")["best"]["config"])
    s1 = read(ROOT/b1/"fits"/j1/"summary.json")
    c1 = model_config_from_dict(s1["best"]["config"])
    m_halo_10 = float(np.interp(10., s1["refined_profiles"]["r_pc"], s1["refined_profiles"]["halo"]))
    d0 = c0.to_dict(); d0["matter"]["rho20"] = c1.matter.rho20; d0["matter"]["r_s"] = c1.matter.r_s
    c_naive = model_config_from_dict(d0)

    def ev_of(cfg):
        cfg = refined_config(cfg); return problem.evaluate(cfg, build_df_model(cfg))
    e0, e1, en = ev_of(c0), ev_of(c1), ev_of(c_naive)
    names = [p.name for p in problem.data.profiles]+[c.name for c in problem.counts]
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5), constrained_layout=True)
    axes = axes.ravel()
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), nodm=args.nodm, dm=args.dm,
                  halo=dict(rho20=c1.matter.rho20, r_s=c1.matter.r_s),
                  objective=dict(nodm=e0["objective"], dm_refit=e1["objective"], halo_added=en["objective"]),
                  chi2_kin=dict(nodm=e0["chi2_kinematic"], dm_refit=e1["chi2_kinematic"], halo_added=en["chi2_kinematic"]), datasets={})
    for ax, name in zip(axes, names):
        prof = next((p for p in problem.data.profiles if p.name == name), None)
        if prof is not None:
            r = prof.r*PC; err = 0.5*(prof.err_lo+prof.err_hi)
            dn = (en["predictions"][name]-e0["predictions"][name])/err
            dr = (e1["predictions"][name]-e0["predictions"][name])/err
        else:
            c = next(c for c in problem.counts if c.name == name)
            r = c.r_median*PC; mu0 = e0["counts"][name]["mu"]
            dn = (en["counts"][name]["mu"]-mu0)/np.sqrt(mu0); dr = (e1["counts"][name]["mu"]-mu0)/np.sqrt(mu0)
        ax.plot(r, dn, "o-", color="tab:red", ms=4, label=f"halo added, nothing refitted (max {np.max(np.abs(dn)):.1f}$\\sigma$)")
        ax.plot(r, dr, "s-", color="tab:blue", ms=4, label=f"halo + everything refitted (max {np.max(np.abs(dr)):.2f}$\\sigma$)")
        ax.axhline(0, color="black", lw=.7)
        for y in (-1, 1): ax.axhline(y, color="0.7", lw=.7, ls=":")
        ax.set(xscale="log", title=TITLES.get(name, name), xlabel="projected radius [pc]", ylabel="(model $-$ no-DM fit) / error")
        ax.legend(fontsize=7, loc="best")
        record["datasets"][name] = dict(r_pc=r.tolist(), halo_added=np.asarray(dn).tolist(), refit=np.asarray(dr).tolist())
    axes[-1].axis("off")
    axes[-1].text(0, .95, "\n".join([
        f"halo: cored, rho20 = {c1.matter.rho20:g} Msun/pc^3, r_s = {c1.matter.r_s:.0f} pc",
        f"M_DM(<10 pc) = {m_halo_10:.1e} Msun",
        "",
        f"objective  no-DM fit: {e0['objective']:.1f}",
        f"           halo added, nothing refitted: {en['objective']:.1f}",
        f"           halo + refit: {e1['objective']:.1f}",
        "",
        "what the refit changed to compensate:",
        f"  M_star  {c0.M_star:.3e} -> {c1.M_star:.3e} ({100*(c1.M_star/c0.M_star-1):+.1f}%)",
        f"  M_rem   {c0.matter.M_rem:.2e} -> {c1.matter.M_rem:.2e} ({100*(c1.matter.M_rem/c0.matter.M_rem-1):+.1f}%)",
        f"  a_rem   {c0.matter.a_rem:.2f} -> {c1.matter.a_rem:.2f} pc",
        f"  J_a     {c0.stellar.J_a:.0f} -> {c1.stellar.J_a:.0f};  J_outer {c0.stellar.J_outer:.0f} -> {c1.stellar.J_outer:.0f}",
        f"  b_out   {c0.stellar.b_out:.2f} -> {c1.stellar.b_out:.2f};  b_outer {c0.stellar.b_outer:.2f} -> {c1.stellar.b_outer:.2f}",
        f"  D       {c0.distance_kpc:.3f} -> {c1.distance_kpc:.3f} kpc",
    ]), va="top", fontsize=9, family="monospace", transform=axes[-1].transAxes)
    fig.suptitle("Why the DM and no-DM curves coincide: the halo's kinematic signature (red) is cancelled by refitting the other components (blue)")
    plot = ROOT/"plots"/f"{args.name}.png"
    fig.savefig(plot, dpi=150)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/f"{args.name}.json").write_text(json.dumps(record))
    print(json.dumps(record["objective"]), plot)


if __name__ == "__main__":
    main()
