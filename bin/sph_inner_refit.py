#!/usr/bin/env python3
"""Step 1 (docs/PLAN_2026-11-19_limits.md): local MAP refits of the spherical spline K1 / NFW with the fixed inner node
beta(2") = -0.5 and b10 in [-1, 0]; legacy (no inner node, b10 <= -0.5) as control. Several starts per case.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/sph_inner_refit.py KIND {inner,legacy}
Output: results/limits/sph_local_<KIND>_<mode>.json
"""
import json, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"bin"))
from sph_local import Local

kind, mode = sys.argv[1], sys.argv[2]
L = Local(kind, inner_node=2.0 if mode == "inner" else None)
x0 = np.load(ROOT/f"results/fits/recover_real_{kind}_spline_nw/ml_x.npy")
ib = L.names.index("b10")
starts = {"nested_ML": x0}
if mode == "inner":
    for b in (0.0, -0.25, -0.75):
        s = x0.copy(); s[ib] = b; starts[f"b10={b}"] = s
out = {}
for lab, s in starts.items():
    t = time.time(); r = L.fit(s)
    print(f"{kind} {mode} {lab}: lnpost {r['lnpost']:.3f} lnL {r['lnL']:.3f} lnpriorD {r['lnprior_D']:.3f} "
          f"chi2 {r['chi2_total']:.2f} b10 {r['theta']['b10']:.3f} D {r['theta']['distance']:.3f} "
          f"Mrem/M* {r['theta']['M_rem']/r['theta']['M_star']:.2f} ({time.time()-t:.0f}s)", flush=True)
    out[lab] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()}
best = max(out, key=lambda k: out[k]["lnpost"])
res = dict(kind=kind, mode=mode, best_start=best, best=out[best], all=out, names=L.names,
           spread_lnpost=float(max(o["lnpost"] for o in out.values()) - min(o["lnpost"] for o in out.values())))
(ROOT/f"results/limits/sph_local_{kind}_{mode}.json").write_text(json.dumps(res, indent=1))
print("BEST", kind, mode, best, json.dumps({k: round(v, 4) for k, v in out[best]["theta"].items()}))
