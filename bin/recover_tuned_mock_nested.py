#!/usr/bin/env python3
"""Step 5 with the production sampler (2026-11-13): UltraNest fit of one spherical Jeans family to a noiseless tuned
flattened JAM mock.

The local optimiser in recover_tuned_mock.py gave start-dependent likelihoods on the flattened mocks, so the fits
are redone with the same nested sampler and settings as the production runs (fit.run_nested: 400 live points,
dlogz 0.5, slice sampler with the production walk length). Mock data: the saved tuned prediction (noiseless) with
the corrected real errors and streaming terms, built exactly as in recover_tuned_mock.py.

Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/recover_tuned_mock_nested.py \
           --truth q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1 --family K1
       --truth real fits the real data with the corrected rotation terms (same data as the JAM fits).
       --per-bin N: stars per HST/MUSE bin used for the rotation terms (default 10, as the 2026-11-13 runs; 0 = all
       stars, as the JAM fits since 2026-11-15). --label L: output suffix.
       --anisotropy spline (2026-11-16): beta linear in ln r through the JAM spline nodes (anisotropy.SPLINE_NODES_ARCSEC)
       in place of the turnover profile; all other parameters and priors as the production run.
       Under MPI (mpiexec -n N) the output check is made before the slow Setup, so all ranks see the same answer.
Output: results/fits/recover_<truth>_<family>[_<label>]/ (summary.json, posterior.ecsv, ...)
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np
import yaml

from ocen_dm.kinematics import FitProblem, KinematicData
from ocen_dm.kinematics.fit import run_nested
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config
from tune_jam_to_data import Setup

FAMILIES = {"K1": "jeans_K1_turnover_counts_c", "NFW": "jeans_K2_nfw_turnover_counts_c"}


TURN = ("beta_0", "beta_mid", "beta_inf", "r_beta", "delta_r_beta")


def spline_config(cfg):
    """The production family config with the turnover anisotropy replaced by the spline (priors from fit.py)."""
    from ocen_dm.kinematics.fit import NoDarkMatterModel
    from ocen_dm.kinematics.run_io import _prior_config
    tmp = NoDarkMatterModel(mge_fit=family_from_config(cfg).mge_fit, anisotropy_profile="spline",
                            beta0_max=cfg["beta0_max"], instruments=cfg["instruments"], fixed=cfg["fixed"])
    spl = [{"name": p.name, "prior": _prior_config(p.prior), "unit": p.unit, "latex": p.latex}
           for p in tmp.parameters if p.name.startswith("b") and p.name[1:].isdigit()]
    params, done = [], False
    for p in cfg["parameters"]:
        if p["name"] in TURN:
            if not done:
                params += spl; done = True
            continue
        params.append(p)
    return {**cfg, "anisotropy_profile": "spline", "parameters": params}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="jam_tuned tag, or 'real' for the corrected real data")
    ap.add_argument("--family", choices=sorted(FAMILIES), required=True)
    ap.add_argument("--per-bin", type=int, default=10, help="stars per HST/MUSE bin for the rotation terms; 0 = all")
    ap.add_argument("--label", default="")
    ap.add_argument("--anisotropy", choices=("turnover", "spline"), default="turnover")
    args = ap.parse_args()
    out = ROOT/"results/fits"/(f"recover_{args.truth}_{args.family}"+(f"_{args.label}" if args.label else ""))
    if out.exists():
        sys.exit(f"{out} exists; use --label")
    S = Setup(per_bin=None if args.per_bin == 0 else args.per_bin)   # corrected errors and streaming terms
    if args.truth == "real":                  # the real data with the corrected HST rotation terms (no mock)
        mdata = S.data
        prov = {"kind": "real", "rotation": "results/rotation/rotation_field_cylindrical_20261113.json (MAP) for HST, "
                "MUSE and Gaia, errors with the field-chain spread (tune_jam_to_data.Setup._corrected_data)",
                "per_bin": args.per_bin, "script": "bin/recover_tuned_mock_nested.py --truth real"}
    else:
        T = json.loads((ROOT/f"results/rotation/jam_tuned_{args.truth}.json").read_text())
        mdata = KinematicData(tuple(replace(p, value=np.array(T["pred"][p.name])) for p in S.data.profiles))
        prov = {"kind": "mock", "truth": f"results/rotation/jam_tuned_{args.truth}.json",
                "noise": "none", "script": "bin/recover_tuned_mock_nested.py"}
    prod = ROOT/"results/fits"/FAMILIES[args.family]
    cfg = family_config(problem_for(load_run_metadata(prod)).family)
    if args.anisotropy == "spline":
        cfg = spline_config(cfg)
        prov["anisotropy"] = "spline (anisotropy.SPLINE_NODES_ARCSEC)"
    fam = family_from_config(cfg)
    smp = yaml.safe_load((prod/"run.yaml").read_text())["sampler"]
    s = run_nested(FitProblem(fam, mdata), out, n_live=smp["n_live"], dlogz=smp["dlogz"], seed=smp["seed"],
                   step_sampler=smp["step_sampler"], slice_nsteps=smp["slice_nsteps"],
                   data_provenance=prov)
    print(f"{args.truth} {args.family}: lnL_max {s['lnL_max']:.2f}  logZ {s['logz']:.2f}+-{s['logzerr']:.2f}  "
          f"chi2 {s['chi2_ml_total']}  ({s['elapsed_s']/60:.0f} min)", flush=True)


if __name__ == "__main__":
    main()
