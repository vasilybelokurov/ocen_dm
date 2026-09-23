#!/usr/bin/env python3
"""Compare the intrinsic anisotropy beta(r) of the best compact-DF fits with Jeans fits.

beta = 1 - sigma_t^2/sigma_r^2 (sigma_t one tangential component). DF curves come
from the refined model of each selected fit; Jeans curves use the maximum-
likelihood TurnoverAnisotropy of the rung-2 fits (same 89-bin data snapshot).
The shaded bands show the radial coverage of HST/MUSE (<= 9.5 pc) and Gaia (9.4-63 pc).

Usage: python bin/plot_beta_profiles.py results/df/observed_df_wide_20260923 \
           observed_free_halo_start0 observed_no_halo_start0
       python bin/plot_beta_profiles.py --name NAME BATCH::JOB::LABEL [...]
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

from ocen_dm.kinematics.anisotropy import TurnoverAnisotropy
from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict

JEANS = {"rung2_K1_turnover": "Jeans turnover, no DM (chi2 190.0)",
         "rung2_K2_cored_turnover": "Jeans turnover, cored DM (chi2 186.8)"}


def read(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("specs", nargs="+",
                        help="BATCH JOB [JOB...] or, with --name, BATCH::JOB::LABEL entries")
    parser.add_argument("--name", help="output name; enables BATCH::JOB::LABEL specs across batches")
    args = parser.parse_args()
    if args.name:
        entries = [tuple(e.split("::")) for e in args.specs]
        out_name = args.name
    else:
        batch = Path(args.specs[0])
        entries = [(str(batch), job, None) for job in args.specs[1:]]
        out_name = batch.resolve().name+"_beta"
    out = Path(entries[0][0]).resolve()
    r = np.geomspace(.05, 80., 200)
    record = dict(batch=str(out), created_utc=datetime.now(timezone.utc).isoformat(), r_pc=r.tolist(), curves={})
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.axvspan(.1, 9.5, color="tab:blue", alpha=.07, label="HST/MUSE coverage")
    ax.axvspan(9.4, 63, color="tab:green", alpha=.07, label="Gaia coverage")
    for batch, job, name in entries:
        s = read(Path(batch)/"fits"/job/"summary.json")
        model = build_df_model(refined_config(model_config_from_dict(s["best"]["config"])))
        beta = model.intrinsic_moments(r)["beta"]
        g = s["gates"]
        label = (f"DF {name or job.removeprefix('observed_')} "
                 f"(chi2_kin/N {g['chi2_kinematic']/g['n_kinematic']:.2f}, N={g['n_kinematic']})")
        ax.plot(r, beta, lw=2, label=label)
        record["curves"][label] = beta.tolist()
    for run, label in JEANS.items():
        p = read(ROOT/"results/fits"/run/"summary.json")["parameters"]
        a = TurnoverAnisotropy(**{k: p[k]["ml"] for k in ("beta_0", "beta_mid", "beta_inf", "r_beta", "delta_r_beta")})
        beta = a.beta(r)
        ax.plot(r, beta, ls="--", label=label)
        record["curves"][label] = beta.tolist()
    ax.axhline(0, color="black", lw=.7)
    ax.set(xscale="log", xlabel="Radius r [pc]", ylabel=r"$\beta = 1-\sigma_t^2/\sigma_r^2$",
           title="Intrinsic anisotropy: compact DF fits vs Jeans turnover fits")
    ax.legend(fontsize=7)
    plot = ROOT/"plots"/(out_name+".png")
    fig.savefig(plot, dpi=160)
    record["plot_sha256"] = hashlib.sha256(plot.read_bytes()).hexdigest()
    (ROOT/"results/plot_data"/(out_name+".json")).write_text(json.dumps(record))
    print(plot)


if __name__ == "__main__":
    main()
