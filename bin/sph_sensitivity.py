#!/usr/bin/env python3
"""Plan steps 5-6 (docs/PLAN_2026-11-19_limits.md): spherical spline (inner node) K1 vs NFW local MAP fits under
prior / likelihood variants. The halo gain is reported as 2*Delta lnpost (NFW - K1) for each case.

Cases
  base      : production priors
  floor     : M_star >= 1.3e6 Msun (M/L_V >= 1.2 for L_V ~ 1.1e6 Lsun)
  mremHN    : extra prior -0.5 (f/0.10)^2 on f = M_rem/M_star, f <= 0.3
  floor+HN  : both
  rotmarg   : rotation field marginalised over 200 chain draws (bin/rotation_draws.py) instead of per-bin errors
Starts: the inner-node K1 and NFW optima (results/limits/sph_local_{K1,NFW}_inner.json), projected into the
constrained box (M_star to the floor, M_rem to 0.2 M_star), plus the K1 optimum with a small halo for NFW.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/sph_sensitivity.py --workers 4
Output: results/limits/sph_sensitivity.json
"""
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"bin"))

CASES = {"base": {}, "floor": dict(floor=1.3e6), "mremHN": dict(hn=0.10), "floor+HN": dict(floor=1.3e6, hn=0.10),
         "rotmarg": dict(rotmarg=True)}


def m50(th):
    m = lambda x: np.log1p(x) - x/(1+x)
    return th["M_dm_100"]*m(50/th["r_s"])/m(100/th["r_s"])


def run(job):
    from sph_local import Local
    case, kind = job
    c = CASES[case]
    po = {"M_star": (c["floor"], 1e7)} if "floor" in c else None
    L = Local(kind, inner_node=2.0, prior_overrides=po, rotmarg=c.get("rotmarg", False), mrem_hn=c.get("hn"))
    src = {k: json.loads((ROOT/f"results/limits/sph_local_{k}_inner.json").read_text()) for k in ("K1", "NFW")}
    starts = []
    for k in (("K1",) if kind == "K1" else ("NFW", "K1")):
        th = dict(src[k]["best"]["theta"])
        if kind == "NFW" and k == "K1":
            th.update(M_dm_100=3e4, r_s=300.0)
        if "floor" in c:
            th["M_star"] = max(th["M_star"], c["floor"]*1.01)
        if "hn" in c and th["M_rem"]/th["M_star"] > 0.25:
            th["M_rem"] = 0.2*th["M_star"]
        starts.append((k, np.array([th[n] for n in L.names])))
    res = []
    for lab, x0 in starts:
        t = time.time(); r = L.fit(x0)
        r = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in r.items()}
        r["start"], r["elapsed_s"] = lab, round(time.time()-t, 1)
        res.append(r)
        print(f"{case:9s} {kind:3s} start {lab}: lnpost {r['lnpost']:.3f} lnL {r['lnL']:.3f} prior {r['lnprior_D']:.3f} "
              f"chi2 {r['chi2_total']:.2f} M*/1e6 {r['theta']['M_star']/1e6:.2f} Mrem/M* "
              f"{r['theta']['M_rem']/r['theta']['M_star']:.2f} ({r['elapsed_s']:.0f}s)", flush=True)
    best = max(res, key=lambda r: r["lnpost"])
    return case, kind, dict(best=best, spread=max(r["lnpost"] for r in res) - min(r["lnpost"] for r in res))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--cases", nargs="*", default=list(CASES)); a = ap.parse_args()
    jobs = [(c, k) for c in a.cases for k in ("NFW", "K1")]
    with Pool(a.workers) as p:
        out = p.map(run, jobs, chunksize=1)
    R = {}
    for case, kind, r in out:
        R.setdefault(case, {})[kind] = r
    for case, d in R.items():
        k1, nf = d["K1"]["best"], d["NFW"]["best"]
        d["two_dlnpost"] = 2*(nf["lnpost"] - k1["lnpost"])
        d["two_dlnL"] = 2*(nf["lnL"] - k1["lnL"])
        d["NFW_M50"] = float(m50(nf["theta"]))
        print(f"{case:9s} 2dlnpost {d['two_dlnpost']:.2f}  2dlnL {d['two_dlnL']:.2f}  chi2 K1 {k1['chi2_total']:.1f} "
              f"NFW {nf['chi2_total']:.1f}  NFW M100 {nf['theta']['M_dm_100']:.3g} r_s {nf['theta']['r_s']:.0f} "
              f"M(<50pc) {d['NFW_M50']:.3g}  Mrem/M* K1 {k1['theta']['M_rem']/k1['theta']['M_star']:.2f} "
              f"NFW {nf['theta']['M_rem']/nf['theta']['M_star']:.2f}")
    p = ROOT/"results/limits/sph_sensitivity.json"
    old = json.loads(p.read_text()) if p.exists() else {}
    old.update(R); p.write_text(json.dumps(old, indent=1))


if __name__ == "__main__":
    main()
