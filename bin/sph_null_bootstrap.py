#!/usr/bin/env python3
"""Plan step 3 (docs/PLAN_2026-11-19_limits.md): parametric-bootstrap calibration of the spherical halo test.

Mocks: the model prediction at a truth point plus split-normal noise drawn from the data error bars (the exact
generative model of the likelihood: below the model with scale err_hi, above with err_lo, probabilities
err_hi/(err_lo+err_hi) and err_lo/(...)).
--mode null : truth = the local K1 optimum with the inner node (results/limits/sph_local_K1_inner.json). Each mock is
              fitted with K1 (start: truth) and NFW (r_s free; starts: real-data NFW optimum, truth + tiny halo);
              statistic t = 2 [lnpost_NFW - lnpost_K1] (the D-prior term is common).
--mode inject --mdm M : truth = NFW with M_dm(<100 pc) = M, r_s = 300 pc, other parameters refitted to the real
              data at that fixed halo. Each mock: profile at r_s = 300 fixed, q = 2 [max_M lnpost - lnpost(M_true)]
              (one-sided: 0 if M_hat > M_true) and the free-r_s detection statistic t.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src:bin python bin/sph_null_bootstrap.py --mode null --n 200 --workers 4
Output: results/limits/boot_<mode>[_<M>].jsonl (one line per mock; resumable).
"""
import argparse, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
from ocen_dm.kinematics import KinematicData
from ocen_dm.kinematics.fit import BinnedProfile_replace
from sph_local import Local, real_data

RS_INJ = 300.0


def draw(pred, data, rng):
    out = []
    for p in data.profiles:
        m = np.asarray(pred[p.name], float)
        lo, hi = np.asarray(p.err_lo, float), np.asarray(p.err_hi, float)
        below = rng.random(m.size) < hi/(lo+hi)
        z = np.abs(rng.normal(size=m.size))
        out.append(BinnedProfile_replace(p, value=np.where(below, m - z*hi, m + z*lo)))
    return KinematicData(tuple(out))


def truths(mode, mdm):
    K = json.loads((ROOT/"results/limits/sph_local_K1_inner.json").read_text())["best"]["theta"]
    N = json.loads((ROOT/"results/limits/sph_local_NFW_inner.json").read_text())["best"]["theta"]
    if mode == "null":
        return "K1", K, N
    f = ROOT/f"results/limits/inject_truth_{mdm:.3g}.json"
    if f.exists():
        return "NFW", json.loads(f.read_text()), N
    L = Local("NFW")
    x0 = np.array([N[n] for n in L.names])
    r = L.fit(x0, fixed={"M_dm_100": mdm, "r_s": RS_INJ})
    f.write_text(json.dumps(r["theta"], indent=1))
    return "NFW", r["theta"], N


def one(job):
    i, mode, mdm, seed = job
    kind, T, Nreal = truths(mode, mdm)
    data = real_data()
    Lt = Local(kind, data=data)
    pred = Lt.P.predict(np.array([T[n] for n in Lt.names]))
    mock = draw(pred, data, np.random.default_rng(seed))
    LK, LN = Local("K1", data=mock), Local("NFW", data=mock)
    t0 = time.time()
    xk = np.array([T[n] for n in LK.names])
    rk = LK.fit(xk)
    starts = [np.array([Nreal[n] for n in LN.names])]
    s2 = dict(T); s2.setdefault("M_dm_100", 1e4); s2.setdefault("r_s", 300.0)
    starts.append(np.array([s2[n] for n in LN.names]))
    rn = max((LN.fit(s) for s in starts), key=lambda r: r["lnpost"])
    out = dict(i=i, seed=seed, lnpost_K1=rk["lnpost"], lnpost_NFW=rn["lnpost"], t=2*(rn["lnpost"]-rk["lnpost"]),
               NFW_theta=rn["theta"], K1_theta=rk["theta"], chi2_K1=rk["chi2_total"])
    if mode == "inject":
        xs = np.array([s2[n] for n in LN.names])
        prof = LN.fit(xs, fixed={"r_s": RS_INJ})                       # max over M at r_s = 300
        ptrue = LN.fit(xs, fixed={"r_s": RS_INJ, "M_dm_100": mdm})     # at the true M
        mhat = prof["theta"]["M_dm_100"]
        out.update(Mhat_rs300=mhat, lnpost_rs300=prof["lnpost"], lnpost_true=ptrue["lnpost"],
                   q=0.0 if mhat > mdm else 2*max(prof["lnpost"]-ptrue["lnpost"], 0.0))
    out["elapsed_s"] = time.time()-t0
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("null", "inject"), required=True)
    ap.add_argument("--mdm", type=float, default=None)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed0", type=int, default=1000)
    a = ap.parse_args()
    f = ROOT/("results/limits/boot_null.jsonl" if a.mode == "null" else f"results/limits/boot_inject_{a.mdm:.3g}.jsonl")
    done = {json.loads(l)["i"] for l in f.read_text().splitlines()} if f.exists() else set()
    if a.mode == "inject":
        truths(a.mode, a.mdm)                       # build the truth once before forking
    jobs = [(i, a.mode, a.mdm, a.seed0 + i + (0 if a.mode == "null" else int(np.log10(a.mdm)*1e4)))
            for i in range(a.n) if i not in done]
    with Pool(a.workers) as pool, open(f, "a") as fh:
        for r in pool.imap_unordered(one, jobs):
            fh.write(json.dumps(r)+"\n"); fh.flush()
            print(f"mock {r['i']}: t {r['t']:.2f}" + (f" q {r['q']:.2f} Mhat {r['Mhat_rs300']:.3g}" if "q" in r else "")
                  + f" ({r['elapsed_s']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
