#!/usr/bin/env python3
"""Locate fast-vs-direct projected-moment discrepancies for a saved compact-DF fit.

Rebuilds the fit's best configuration at nominal and refined resolution and
reports |fast/direct - 1| for Sigma, los, pmr and pmt at each projected radius.
The fitter's validation gate uses the maximum over R = geomspace(0.05, 80, 12) pc.

Usage: python bin/probe_projection_accuracy.py results/df/compact_recovery_20260923_v2 cored_dm_free_halo_start1
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np

from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.df_fit import build_df_model, model_config_from_dict


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batch", type=Path)
    parser.add_argument("job")
    args = parser.parse_args()
    out = args.batch.resolve()
    best = json.loads((out/"fits"/args.job/"best.json").read_text())
    config = model_config_from_dict(best["config"])
    R = np.geomspace(.05, 80., 24)
    result = dict(batch=str(out), job=args.job, R_pc=R.tolist(), resolutions={})
    for name, c in (("nominal", config), ("refined", refined_config(config))):
        model = build_df_model(c)
        fast, direct = model.projected_moments(R), model.direct_projected_moments(R)
        err = {k: np.abs(fast[k]/direct[k]-1).tolist() for k in fast}
        result["resolutions"][name] = dict(relative_error=err, diagnostics=model.diagnostics)
        print(name, {k: "max %.4f at R=%.3g pc" % (max(v), R[int(np.argmax(v))]) for k, v in err.items()}, flush=True)
    dest = ROOT/"results/maintenance"/f"projection_accuracy_{out.name}_{args.job}.json"
    dest.write_text(json.dumps(result, default=float))
    print(dest)


if __name__ == "__main__":
    main()
