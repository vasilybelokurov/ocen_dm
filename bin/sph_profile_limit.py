#!/usr/bin/env python3
"""Plan step 4 (docs/PLAN_2026-11-19_limits.md): profile likelihood of the spherical NFW (inner-node spline) in the halo
mass at fixed r_s, and the dark mass inside 50 pc.

Grid: r_s in {30, 100, 300, 1000} pc x log10 M_dm(<100 pc) in [3, 7] (step 0.25); at each point every other parameter
(incl. D) is refitted (sph_local.Local, restarted Nelder-Mead), warm-started from the neighbouring grid point and from
the K1 optimum with that halo. M(<50 pc) = M_dm_100 m(50/r_s)/m(100/r_s), m(x) = ln(1+x) - x/(1+x) (the NFW
truncation r_t = 1000 pc lies far outside). Upper limit: largest M with q(M) = 2[max lnpost - lnpost(M)] <= threshold
on the high-M side (threshold: --thr, default 2.71; the calibrated value from the injections replaces it).
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/sph_profile_limit.py [--rotmarg] [--workers 4]
Output: results/limits/sph_profile[_rotmarg].json
"""
import argparse, json, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
from sph_local import Local

RS = (30.0, 100.0, 300.0, 1000.0)
LOGM = np.arange(3.0, 7.0001, 0.25)


def mnfw(x):
    return np.log1p(x) - x/(1+x)


def m50(m100, rs):
    return m100*mnfw(50/rs)/mnfw(100/rs)


def profile_rs(job):
    rs, rotmarg = job
    L = Local("NFW", rotmarg=rotmarg)
    K = json.loads((ROOT/"results/limits/sph_local_K1_inner.json").read_text())["best"]["theta"]
    N = json.loads((ROOT/"results/limits/sph_local_NFW_inner.json").read_text())["best"]["theta"]
    rows, prev = [], None
    for lm in LOGM:
        fixed = {"M_dm_100": 10**lm, "r_s": rs}
        starts = [np.array([{**K, **fixed}.get(n) for n in L.names]), np.array([{**N, **fixed}[n] for n in L.names])]
        if prev is not None:
            starts.append(np.array([{**prev, **fixed}[n] for n in L.names]))
        r = max((L.fit(s, fixed=fixed) for s in starts), key=lambda r: r["lnpost"])
        prev = r["theta"]
        rows.append(dict(logM100=float(lm), M50=float(m50(10**lm, rs)), lnpost=r["lnpost"], lnL=r["lnL"],
                         lnprior_D=r["lnprior_D"], chi2=r["chi2_total"], theta=r["theta"]))
        print(f"r_s {rs:g} logM100 {lm:.2f}: lnpost {r['lnpost']:.3f} M_star {r['theta']['M_star']:.3g} "
              f"D {r['theta']['distance']:.3f}", flush=True)
    return rs, rows


def upper_limit(rows, lnmax, thr):
    lm = np.array([r["logM100"] for r in rows]); q = 2*(lnmax - np.array([r["lnpost"] for r in rows]))
    m5 = np.log10([r["M50"] for r in rows])
    i0 = int(np.argmin(q))
    for i in range(i0, len(q)-1):
        if q[i] <= thr < q[i+1]:
            f = (thr - q[i])/(q[i+1]-q[i])
            return float(10**(lm[i] + f*(lm[i+1]-lm[i]))), float(10**(m5[i] + f*(m5[i+1]-m5[i])))
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rotmarg", action="store_true"); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--thr", type=float, default=2.71)
    a = ap.parse_args()
    with Pool(a.workers) as pool:
        res = dict(pool.map(profile_rs, [(rs, a.rotmarg) for rs in RS]))
    K1 = Local("K1", rotmarg=a.rotmarg).fit(np.array(list(json.loads(
        (ROOT/"results/limits/sph_local_K1_inner.json").read_text())["best"]["theta"].values())))
    LN = Local("NFW", rotmarg=a.rotmarg)
    NF = LN.fit(np.array([json.loads((ROOT/"results/limits/sph_local_NFW_inner.json").read_text())["best"]["theta"][n]
                          for n in LN.names]))  # free NFW optimum (may lie between grid points / at the r_s bound)
    lnmax = max(NF["lnpost"], max(max(r["lnpost"] for r in rows) for rows in res.values()))
    lim = {str(rs): dict(zip(("M100_lim", "M50_lim"), upper_limit(rows, lnmax, a.thr))) for rs, rows in res.items()}
    lim_rs = {str(rs): dict(zip(("M100_lim", "M50_lim"), upper_limit(rows, max(r["lnpost"] for r in rows), a.thr)))
              for rs, rows in res.items()}
    out = dict(threshold=a.thr, rotmarg=a.rotmarg, lnpost_K1=K1["lnpost"], lnpost_max=lnmax, NFW_free=NF["theta"],
               lnpost_NFW_free=NF["lnpost"],
               gain_2dlnL=2*(lnmax-K1["lnpost"]), limits_global_max=lim, limits_per_rs_max=lim_rs,
               grid={str(k): v for k, v in res.items()})
    f = ROOT/f"results/limits/sph_profile{'_rotmarg' if a.rotmarg else ''}.json"
    f.write_text(json.dumps(out, indent=1))
    print("K1 lnpost", round(K1["lnpost"], 3), "max", round(lnmax, 3), "2dlnL", round(out["gain_2dlnL"], 2))
    print("limits (global max):", json.dumps(lim)); print("limits (per-r_s max):", json.dumps(lim_rs)); print(f)


if __name__ == "__main__":
    main()
