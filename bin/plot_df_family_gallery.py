#!/usr/bin/env python3
"""Gallery of the regularized exponential DF family: what shapes and observables it can make.

Starts from a fitted model (default: the two-transition, no-halo best fit) and
varies one ingredient at a time, rebuilding each self-consistent equilibrium
(stars + the fit's fixed remnant Plummer sphere; no halo, no point mass):
  (a) envelope exponent alpha, (b) action scale J0 (size), (c) one-transition
  anisotropy b_out (J_a fixed, no second transition), (d) the second transition.
Columns: intrinsic density rho(r) (normalized at 1 pc), intrinsic beta(r),
projected sigma_los(R), and sigma_pmT/sigma_pmR(R), the direct projected
anisotropy signature. These are model illustrations, not fits.

Usage: python bin/plot_df_family_gallery.py [--batch results/df/dftwo_20260923 --job observed_no_halo_start0]
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np

R_INT = np.geomspace(.05, 80., 120)
R_PROJ = np.geomspace(.1, 60., 60)


def variants(base):
    """(row, label, config) for one-at-a-time changes of the base configuration."""
    st = base.stellar
    one = replace(st, b_outer=None, J_outer=None)
    rows = []
    for a in (.7, .85, 1.1, 1.4):
        rows.append(("a", f"alpha = {a:g}", replace(base, stellar=replace(st, alpha=a))))
    for j in (40., 65., 120.):
        rows.append(("b", f"J0 = {j:g}", replace(base, stellar=replace(st, J0=j))))
    for b in (-1., 0., .5, 1.):
        rows.append(("c", f"one transition, b_out = {b:+g}", replace(base, stellar=replace(one, b_out=b, J_a=30.))))
    for bo, jo in ((st.b_outer, st.J_outer), (0., st.J_outer), (st.b_out, st.J_outer), (st.b_outer, 300.)):
        label = f"b_outer = {bo:+.2f}, J_outer = {jo:.0f}"
        rows.append(("d", label, replace(base, stellar=replace(st, b_outer=bo, J_outer=jo))))
    return rows


def compute(item):
    row, label, config_dict = item
    from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict
    try:
        model = build_df_model(model_config_from_dict(config_dict))
        mom = model.intrinsic_moments(R_INT)
        proj = model.projected_moments(R_PROJ)
        s = proj["Sigma"]
        # Mask where the tracer has effectively vanished: moments there are numerical noise.
        live = mom["rho"] > 1e-7*mom["rho"][0]
        live_p = s > 1e-7*s[0]
        nan = lambda v, m: np.where(m, v, np.nan).tolist()
        return dict(row=row, label=label, rho=nan(mom["rho"], live), beta=nan(mom["beta"], live),
                    los=nan(np.sqrt(proj["los"]/s), live_p), pmr=nan(np.sqrt(proj["pmr"]/s), live_p),
                    pmt=nan(np.sqrt(proj["pmt"]/s), live_p), stellar=config_dict["stellar"],
                    half_mass_pc=float(np.interp(.5, np.cumsum(mom["rho"]*R_INT**3)/np.sum(mom["rho"]*R_INT**3), R_INT)))
    except Exception as exc:  # report, do not hide, failed equilibria
        return dict(row=row, label=label, failed=repr(exc)[:200])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--batch", type=Path, default=ROOT/"results/df/dftwo_20260923")
    parser.add_argument("--job", default="observed_no_halo_start0")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    from ocen_dm.kinematics.df_fit import model_config_from_dict
    s = json.loads((args.batch/"fits"/args.job/"summary.json").read_text())
    base = model_config_from_dict(s["best"]["config"])
    items = [(row, label, c.to_dict()) for row, label, c in variants(base)]
    with ProcessPoolExecutor(args.workers, mp_context=multiprocessing.get_context("spawn")) as pool:
        results = list(pool.map(compute, items))
    for r in results:
        if "failed" in r:
            print("FAILED", r["label"], r["failed"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    titles = {"a": r"(a) envelope exponent $\alpha$", "b": r"(b) action scale $J_0$ [pc km/s]",
              "c": r"(c) one transition: $b_{\rm out}$ ($J_a=30$)", "d": r"(d) second transition"}
    fig, axes = plt.subplots(4, 4, figsize=(15, 13), constrained_layout=True)
    for i, row in enumerate("abcd"):
        good = [r for r in results if r["row"] == row and "failed" not in r]
        for r in good:
            rho = np.asarray(r["rho"], float); rho1 = np.interp(np.log(1.), np.log(R_INT), rho)
            lab = r["label"]+f"  (r_half {r['half_mass_pc']:.1f} pc)"
            axes[i, 0].plot(R_INT, rho/rho1, label=lab)
            axes[i, 1].plot(R_INT, r["beta"])
            axes[i, 2].plot(R_PROJ, r["los"])
            axes[i, 3].plot(R_PROJ, np.asarray(r["pmt"], float)/np.asarray(r["pmr"], float))
        axes[i, 0].set(xscale="log", yscale="log", ylim=(1e-7, 20), ylabel=r"$\rho/\rho(1\,{\rm pc})$", title=titles[row])
        axes[i, 0].legend(fontsize=6.5, loc="lower left")
        axes[i, 1].set(xscale="log", ylim=(-3, 1), ylabel=r"$\beta(r)$", title="intrinsic anisotropy")
        axes[i, 1].axhline(0, color="grey", lw=.6)
        axes[i, 2].set(xscale="log", ylabel=r"$\sigma_{\rm los}(R)$ [km/s]", title="projected LOS dispersion")
        axes[i, 3].set(xscale="log", ylim=(.5, 2.), ylabel=r"$\sigma_{\rm pm,T}/\sigma_{\rm pm,R}$", title="projected PM ratio")
        axes[i, 3].axhline(1, color="grey", lw=.6)
        for ax in axes[i]:
            ax.axvspan(63, 80, color="0.92")
    for ax in axes[-1, :2]:
        ax.set_xlabel("r [pc]")
    for ax in axes[-1, 2:]:
        ax.set_xlabel("R [pc]")
    fig.suptitle("The regularized exponential DF family around the two-transition best fit "
                 "(one ingredient varied per row; self-consistent, no halo, fixed remnants;\n"
                 "curves stop where the tracer density falls below 1e-7 of its centre; grey: beyond the data)")
    plot = ROOT/"plots"/"df_family_gallery_20260923.png"
    fig.savefig(plot, dpi=150)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), base=str(args.batch/"fits"/args.job),
                  r_pc=R_INT.tolist(), R_pc=R_PROJ.tolist(), models=results,
                  plot_sha256=hashlib.sha256(plot.read_bytes()).hexdigest())
    (ROOT/"results/plot_data"/"df_family_gallery_20260923.json").write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
