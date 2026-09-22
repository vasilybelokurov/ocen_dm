#!/usr/bin/env python3
"""Bounded compact-DF numerical validation; no data fit or posterior launch.

Fresh directories are required. PNGs go to plots/, arrays to results/plot_data/.
The report records analytical checks, equilibrium closure, changed seed,
doubled grids, actual data-bin numerical shifts, and the independent control.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import gamma

from ocen_dm.kinematics.regularized_df import (
    ContourNumerics, ExponentialParameters, PrescribedMatter, RegularizedDFConfig,
    RegularizedDFModel, RegularizedExponentialDF,
)
from ocen_dm.kinematics.lowered_isothermal import LoweredIsothermalConfig, LoweredIsothermalModel
from ocen_dm.kinematics.positive_df import DFNumerics, agama_pc, spherical_velocity_moments
from ocen_dm.kinematics.likelihood import ProfileLikelihood
from ocen_dm.kinematics.run_io import read_data_snapshot, sha256


def write(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n")


def maximum_relative(a, b):
    return float(np.max(np.abs(np.asarray(a)/np.asarray(b)-1)))


def validation_gates(checks):
    """Recomputable gates from saved numerical results; plain JSON booleans."""
    templates = ("isotropic", "radial", "tangential", "two_transition")
    models = ("stars_only", "cored_halo_remnants", "cusped_halo_remnants")
    values = dict(analytic=checks["harmonic"]["normalization_error"] < 1e-5,
                  contours=all(checks["contours_"+k]["endpoint_error"] < 1e-4 for k in templates),
                  regularity=all(checks["contours_"+k]["low_L_ratio_error"] < 1e-4 for k in templates),
                  isotropy=checks["isotropic_reference"]["exact_max_abs_beta"] < .001,
                  projection=max(v for name in models for v in checks[name+"_independent"]["native_projection_error"].values()) < .005,
                  jeans=max(checks[name+"_independent"]["jeans_relative_error"] for name in models) < .01,
                  refinement=max(checks["resolution"]["prediction_shift_in_observed_sigma"].values()) < .1,
                  seed=max(checks["seed_independence"].values()) < .002)
    return {k: bool(v) for k, v in values.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--data-run", type=Path, default=ROOT/"results/fits/rung2_K1_turnover")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    plotpath = ROOT/"plots"/f"{args.out.name}.png"
    datapath = ROOT/"results/plot_data"/f"{args.out.name}.json"
    if plotpath.exists() or datapath.exists():
        raise FileExistsError("choose a fresh validation name")
    start = time.monotonic()
    ag = agama_pc()
    sources = [Path(__file__), *(ROOT/"src/ocen_dm/kinematics"/s for s in (
        "regularized_df.py", "lowered_isothermal.py", "positive_df.py", "df_fit.py",
        "likelihood.py", "run_io.py"))]
    code = args.out/"code"
    code.mkdir()
    for source in sources:
        shutil.copy2(source, code/source.name)
    result = dict(status="running", started_utc=datetime.now(timezone.utc).isoformat(),
                  source_sha256={p.name: sha256(p) for p in sources},
                  agama_version=getattr(ag, "__version__", None),
                  scope="numerical validation, not model selection or a fit to observations", checks={})
    checks = result["checks"]
    arrays = {}

    def checkpoint(name, value):
        checks[name] = value
        result["elapsed_seconds"] = time.monotonic()-start
        write(args.out/"validation.json", result)
        print(f"{result['elapsed_seconds']:.1f}s: {name}: {value}", flush=True)

    def build(name, cfg, **kw):
        t = time.monotonic()
        def progress(row):
            with (args.out/"iterations.jsonl").open("a") as stream:
                stream.write(json.dumps(dict(model=name, elapsed_seconds=time.monotonic()-t, **row))+"\n")
            print(f"{name}: iteration {row['iteration']}, change {row['max_relative_change']:.3g}", flush=True)
        model = RegularizedDFModel(cfg, progress=progress, **kw)
        checkpoint(name, dict(config=cfg.to_dict(), **model.diagnostics, seconds=time.monotonic()-t))
        return model

    try:
        pars = ExponentialParameters(b_out=0)
        harmonic = RegularizedExponentialDF(3e6, pars, frequency_ratio=lambda jr, l: 2+0*(jr+l))
        exact_norm = pars.J0**3*gamma(3/pars.alpha)/(2*pars.alpha)
        checkpoint("harmonic", dict(normalization_error=abs(harmonic.integral/exact_norm-1),
                                    extended_domain_error=abs(harmonic.normalization_integral(320, 2)/exact_norm-1)))
        potential = ag.Potential(type="Isochrone", mass=3e6, scaleRadius=5.)
        af = ag.ActionFinder(potential)
        r = np.geomspace(.03, 250., 72)
        templates = {
            "isotropic": pars,
            "radial": replace(pars, b_out=.8),
            "tangential": replace(pars, b_out=-.8),
            "two_transition": replace(pars, b_out=.8, J_a=100., b_outer=-.6, J_outer=1200.),
        }
        fixed = {}
        for name, parameters in templates.items():
            df = RegularizedExponentialDF(3e6, parameters, potential=potential)
            rho, pr, pt = spherical_velocity_moments(potential, df, af, r, 64).T
            jr = np.geomspace(.1, 3000, 51)
            step = jr*1e-6
            dr = (df.contour_label(jr+step, 0)-df.contour_label(jr-step, 0))/(2*step)
            dl = (df.contour_label(jr, step)-df.contour_label(jr, 0))/step
            error = maximum_relative(df.contour_label(jr, .3*jr), df.direct_contour_label(jr, .3*jr))
            checkpoint("contours_"+name, dict(endpoint_error=error,
                       low_L_ratio_error=float(np.max(abs(dr/dl-2))),
                       normalization_error=df.normalization_error,
                       mass_error=abs(df.totalMass()/3e6-1), beta_range=[float(min(1-pt/pr)), float(max(1-pt/pr))]))
            fixed[name] = dict(r_pc=r.tolist(), rho=rho.tolist(), beta=(1-pt/pr).tolist(),
                               parameters=vars(parameters))
        approximate = RegularizedExponentialDF(3e6, pars, potential=potential,
                                               numerics=ContourNumerics(frequency_mode="epicycle"))
        _, pr, pt = spherical_velocity_moments(potential, approximate, af, r, 64).T
        fixed["epicycle_reference"] = dict(r_pc=r.tolist(), beta=(1-pt/pr).tolist())
        checkpoint("isotropic_reference", dict(exact_max_abs_beta=float(np.max(abs(np.array(fixed["isotropic"]["beta"])))),
                                               epicycle_max_abs_beta=float(np.max(abs(1-pt/pr)))))
        arrays["fixed_potential"] = fixed
        # Local f(vt) has no velocity-volume Jacobian; compare to the old eta=1
        # expression at identical radius and vr, each divided by its vt=0 value.
        vt = np.linspace(0, 8, 101)
        xv = np.zeros((len(vt), 6))
        xv[:, 0], xv[:, 3], xv[:, 4] = 10., 5., vt
        actions = ag.actions(potential, xv)
        df = RegularizedExponentialDF(3e6, templates["radial"], potential=potential)
        f = df(actions)
        old = np.exp(-((actions[:, 0]+actions[:, 1]+abs(actions[:, 2]))/pars.J0)**pars.alpha)
        arrays["velocity_section"] = dict(vt_kms=vt.tolist(), regularized=(f/f[0]).tolist(),
                                           precursor_eta1=(old/old[0]).tolist())

        n = DFNumerics(r_max=1e4, potential_nodes=128, velocity_nodes=32,
                       moment_nodes=240, projection_nodes=128)
        cfg = RegularizedDFConfig(numerics=n)
        models = {"stars_only": build("stars_only", cfg)}
        for slope, name in ((0., "cored_halo_remnants"), (1., "cusped_halo_remnants")):
            cc = replace(cfg, matter=PrescribedMatter(rho20=2., gamma=slope, M_rem=3e5))
            models[name] = build(name, cc)
        R = np.geomspace(.03, 80., 32)
        for name, model in models.items():
            t = model.moment_table
            moments = model.intrinsic_moments(r[r <= 100])
            projected = model.projected_moments(R)
            direct = model.direct_projected_moments(R[::3])
            fast = model.projected_moments(R[::3])
            native_error = {k: maximum_relative(fast[k], direct[k]) for k in fast}
            rr = np.geomspace(.05, 80., 60)
            m = model.intrinsic_moments(rr)
            slope = CubicSpline(np.log(t["r"]), np.log(t["radial_pressure"]))(np.log(rr), 1)
            lhs = m["radial_pressure"]*(slope+2*m["beta"])/rr
            rhs = m["rho"]*model.potential.force(np.column_stack((rr, rr*0, rr*0)))[:, 0]
            checkpoint(name+"_independent", dict(native_projection_error=native_error,
                                                 jeans_relative_error=maximum_relative(lhs, rhs)))
            arrays[name] = dict(r_pc=r[r <= 100].tolist(), **{k: v.tolist() for k, v in moments.items()},
                                R_pc=R.tolist(), projected={k: v.tolist() for k, v in projected.items()})
        base = models["cusped_halo_remnants"]
        fine_n = replace(n, potential_nodes=256, velocity_nodes=64, moment_nodes=480,
                         projection_nodes=256, iteration_tolerance=n.iteration_tolerance/2)
        fine_c = replace(base.config.contours, action_nodes=384, circularity_nodes=129,
                         normalization_order=320, ode_rtol=base.config.contours.ode_rtol/2)
        fine = build("refined_cusp", replace(base.config, numerics=fine_n, contours=fine_c),
                     initial_stellar_potential=base.stellar_potential)
        coarse_moments, fine_moments = base.projected_moments(R), fine.projected_moments(R)
        refined_error = {k: maximum_relative(coarse_moments[k], fine_moments[k]) for k in coarse_moments}
        metadata = json.loads((args.data_run/"summary.json").read_text())
        data = read_data_snapshot(args.data_run, metadata["data_snapshot"])
        # Zero streaming terms here: this checks raw-moment numerical accuracy,
        # without assessing the adopted streaming subtraction or fitting data.
        from ocen_dm.kinematics.likelihood import KinematicData
        data = KinematicData(tuple(replace(p, streaming2=None) for p in data.profiles))
        likelihood = ProfileLikelihood(data)
        low, high = (likelihood.predict(m, cfg.distance_kpc) for m in (base, fine))
        shifts = {p.name: float(np.max(abs(low[p.name]-high[p.name])/np.minimum(p.err_lo, p.err_hi)))
                  for p in data.profiles}
        checkpoint("resolution", dict(projected_relative_error=refined_error,
                                       prediction_shift_in_observed_sigma=shifts,
                                       data_snapshot=metadata["data_snapshot"]))
        alternate = build("alternate_seed", replace(cfg, seed_radius=15.))
        original = models["stars_only"].projected_moments(R)
        other = alternate.projected_moments(R)
        checkpoint("seed_independence", {k: maximum_relative(original[k], other[k]) for k in original})
        controls = {}
        for slope in (0., 1.):
            li = LoweredIsothermalModel(LoweredIsothermalConfig(
                matter=PrescribedMatter(rho20=2., gamma=slope, M_rem=3e5)))
            mm = li.intrinsic_moments(r)
            controls[str(slope)] = dict(r_pc=r.tolist(), beta=mm["beta"].tolist(), rho=mm["rho"].tolist())
            checkpoint("lowered_isothermal_"+str(slope), li.diagnostics)
        arrays["lowered_isothermal"] = controls

        gates = validation_gates(checks)
        result.update(status="passed" if all(gates.values()) else "failed_validation", gates=gates,
                      elapsed_seconds=time.monotonic()-start)
        datapath.parent.mkdir(parents=True, exist_ok=True)
        write(datapath, dict(validation=str(args.out.resolve()), source_sha256=result["source_sha256"], **arrays))
        fig, ax = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
        for name, row in fixed.items():
            ax[0, 0].semilogx(row["r_pc"], row["beta"], label=name.replace("_", " "),
                              ls="--" if name == "epicycle_reference" else "-")
        ax[0, 0].set(title="Computed anisotropy in a fixed isochrone potential", xlabel="r [pc]", ylabel="β(r)")
        ax[0, 0].axhline(0, color="k", lw=.5)
        ax[0, 0].legend(fontsize=7)
        row = arrays["velocity_section"]
        ax[0, 1].plot(row["vt_kms"], row["regularized"], label="Regularized")
        ax[0, 1].plot(row["vt_kms"], row["precursor_eta1"], label="Precursor, η=1", ls="--")
        ax[0, 1].set(title="Local DF at r=10 pc, vr=5 km/s", xlabel="vt [km/s]", ylabel="f(vt) / f(0)")
        ax[0, 1].legend(fontsize=8)
        for name in models:
            row = arrays[name]
            label = name.replace("_", " ")
            ax[0, 2].loglog(row["r_pc"], row["rho"], label=label)
            ax[1, 0].semilogx(row["r_pc"], row["beta"], label=label)
        ax[0, 2].set(title="Self-consistent stellar densities", xlabel="r [pc]", ylabel="ρstar [Msun/pc³]")
        ax[0, 2].legend(fontsize=7)
        ax[1, 0].set(title="Self-consistent anisotropy (same stellar parameters)", xlabel="r [pc]", ylabel="β(r)")
        for k in ("los", "pmr", "pmt"):
            lo = np.sqrt(coarse_moments[k]/coarse_moments["Sigma"])
            hi = np.sqrt(fine_moments[k]/fine_moments["Sigma"])
            ax[1, 1].semilogx(R, 100*(lo/hi-1), label=k)
        ax[1, 1].set(title="Doubling numerical grids: cusped halo case", xlabel="R [pc]", ylabel="Dispersion difference [%]")
        ax[1, 1].legend(fontsize=8)
        for slope, row in controls.items():
            mask = np.array(row["rho"]) > 0
            ax[1, 2].semilogx(np.array(row["r_pc"])[mask], np.array(row["beta"])[mask], label=f"halo γ={slope}")
        ax[1, 2].set(title="Independent lowered-isothermal control", xlabel="r [pc]", ylabel="β(r)")
        ax[1, 2].legend(fontsize=8)
        fig.suptitle("Compact DF infrastructure: numerical validation, not fits to Omega Cen")
        fig.savefig(plotpath, dpi=170)
        plt.close(fig)
        result.update(plot=str(plotpath), plot_data=str(datapath))
        write(args.out/"validation.json", result)
        print(json.dumps(dict(status=result["status"], gates=gates, seconds=result["elapsed_seconds"]), indent=2))
        if not all(gates.values()):
            raise SystemExit(1)
    except Exception as exc:
        result.update(status="error", error=repr(exc), elapsed_seconds=time.monotonic()-start)
        write(args.out/"validation.json", result)
        raise


if __name__ == "__main__":
    main()
