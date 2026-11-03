#!/usr/bin/env python3
"""Centre-marginalised likelihood of the saved fixed-BH DF fits (docs/CENTRE_MARGINALISATION.md).

For each centre of the Gauss-Hermite grids (sigma_c = 1, 2, 3 arcsec with n = 6; sigma_c = 2 with
n = 10) the centre-dependent data are rebuilt (ocen_dm.centre_marg) and the saved best fits (no BH,
8,200 and 4e4 Msun; not refitted) are scored with the full log-likelihood. Also runs:
T2/T3 (zero offset = frozen problem), T4 (n = 10 vs 6), T7 (HST PM bins 13.8-58.5" at the
+-2 sigma_c centres), T8 (unbinned counts-only cross-check in a fixed 2-58.8" window), and the
centre-marginalised Poisson probability of the two innermost HST count bins.

Outputs: results/plot_data/centre_marginal_<label>.json, plots/centre_marginal_<label>.png,
plots/centre_marginal_<label>_data.png.
Usage: python bin/centre_marginal.py [--label 20261103]   (1 core; ~30-60 min, mostly model builds)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.stats import poisson

from ocen_dm import centre_marg as cm
from ocen_dm.kinematics.compact_recovery import refined_config
from ocen_dm.kinematics.counts import CountProfile
from ocen_dm.kinematics.df_fit import DFJointProblem, build_df_model, model_config_from_dict
from ocen_dm.kinematics.likelihood import ARCSEC_PER_RAD
from ocen_dm.kinematics.run_io import read_data_snapshot

RUNS = {"no BH": "results/df/dftwo_counts_area_20261013",
        "BH 8200": "results/df/dfbh8200_restart_20261015",
        "BH 40000": "results/df/dfbh40000_restart_20261015"}
GRIDS = {"s1_n6": (1., 6), "s2_n6": (2., 6), "s3_n6": (3., 6), "s2_n10": (2., 10)}
JOB = "observed_no_halo_start0"


def load_run(path):
    path = ROOT/path
    s = json.loads((path/"fits"/JOB/"summary.json").read_text())
    cfg = refined_config(model_config_from_dict(s["best"]["config"]))
    od = path/"observed"
    rec = json.loads((od/"problem.json").read_text())
    counts = [CountProfile.from_dict(json.loads((od/f"{n}.json").read_text())) for n in rec["counts"]]
    return cfg, read_data_snapshot(od, rec["data_snapshot"]), counts, rec


def scores(ev):
    out = dict(lnl=ev["loglike_profiled"], lnl_kin=ev["loglike_kinematic"], chi2_kin=ev["chi2_kinematic"],
               dev_counts=ev["deviance_counts"])
    for k, t in ev["terms"].items():
        out[f"lnl:{k}"] = t["loglike"]
    for k, t in ev["counts"].items():
        out[f"lnl:{k}"] = t["loglike"]
    hst = next(v for k, v in ev["counts"].items() if k.startswith(cm.HST_NAME))
    out["mu_inner2"] = float(np.sum(hst["mu"][:2]))
    out["mu_hst"] = [float(v) for v in hst["mu"][:14]]
    for k in ("hst_pm_radial_ours", "hst_pm_tangential_ours", "muse_los_dispersion"):
        out[f"pred:{k}"] = [float(v) for v in ev["predictions"][k][:8]]
    return out


class Unbinned:
    """Counts-only unbinned Poisson likelihood in a fixed window 2 <= R0 < r_out about AvdM10 (T8)."""

    def __init__(self, r_in, r_out):
        c = cm._hst_catalogue()
        R0 = np.hypot(c["x"], c["y"])
        m = (R0 >= r_in) & (R0 < r_out)
        self.x, self.y, self.n = c["x"][m], c["y"][m], int(m.sum())
        u, wu = np.polynomial.legendre.leggauss(160)
        lr = np.log(r_in)+(u+1)/2*np.log(r_out/r_in)
        rr = np.exp(lr)
        wr = wu/2*np.log(r_out/r_in)*rr*rr                       # dA = R dR dphi = R^2 dlnR dphi
        phi = (np.arange(512)+.5)/512*2*np.pi
        self.qx = (rr[:, None]*np.sin(phi)[None, :]).ravel()
        self.qy = (rr[:, None]*np.cos(phi)[None, :]).ravel()
        self.qw = (wr[:, None]*np.full(512, 2*np.pi/512)[None, :]).ravel()

    def lnl(self, lnsig, c):
        """Profile (amplitude) log-likelihood up to an N-only constant."""
        Ri = np.hypot(self.x-c[0], self.y-c[1])
        Rq = np.hypot(self.qx-c[0], self.qy-c[1])
        return float(np.sum(lnsig(Ri))-self.n*np.log(np.sum(self.qw*np.exp(lnsig(Rq)))))


def model_lnsigma(model, cfg):
    pc = cfg.distance_kpc*1e3/ARCSEC_PER_RAD
    rmin = cfg.numerics.r_min*1.0001/pc if hasattr(cfg, "numerics") else 1e-3
    R = np.geomspace(max(rmin, 1e-3), 80., 600)
    S = model.projected_moments(R*pc)["Sigma"]
    f = PchipInterpolator(np.log(R), np.log(S), extrapolate=True)
    lo = np.log(R[0])
    return lambda r: f(np.maximum(np.log(np.maximum(r, 1e-12)), lo))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="20261103")
    ap.add_argument("--smoke", action="store_true", help="no-BH model only, 3 centres, no T7 (code check)")
    args = ap.parse_args()
    t_start = time.time()
    runs = {k: load_run(v) for k, v in RUNS.items() if not args.smoke or k == "no BH"}
    shas = {json.dumps(r[3]["data_snapshot"]["sha256"])+json.dumps(r[3]["counts_sha256"]) for r in runs.values()}
    if len(shas) != 1:
        raise SystemExit("the three runs do not share the same frozen data")
    _, data, counts, _ = runs["no BH"]
    cd = cm.CentreData(data, counts)

    models, frozen, lnsig = {}, {}, {}
    for k, (cfg, d, c, _) in runs.items():
        t0 = time.time()
        models[k] = build_df_model(cfg)
        ev = DFJointProblem(d, None, c).evaluate(cfg, models[k])
        frozen[k] = scores(ev)
        lnsig[k] = model_lnsigma(models[k], cfg)
        print(f"{k}: built and scored in {time.time()-t0:.0f}s; lnL {frozen[k]['lnl']:.4f} objective {ev['objective']:.3f}",
              flush=True)

    grids = {g: cm.centre_grid(*v) for g, v in GRIDS.items()}
    centres = {(0., 0.)}
    for c, _ in grids.values():
        centres |= {tuple(np.round(x, 12)) for x in c}
    centres = sorted(centres)
    if args.smoke:
        centres = [(0., 0.), (2., -1.), (-3., 1.5)]
        grids = {"smoke": (np.array([[2., -1.], [-3., 1.5]]), np.array([.5, .5]))}
    hst_tmpl = next(c for c in counts if c.name.startswith(cm.HST_NAME))
    edges = np.r_[hst_tmpl.r_lower, hst_tmpl.r_upper[-1]]
    n_win = 14                                                      # bins 2-58.8"
    ub = Unbinned(edges[0], edges[n_win])
    sub = lambda p: CountProfile(p.name+"_win", p.r_lower[:n_win], p.r_upper[:n_win], p.r_nodes[:n_win],
                                 p.counts[:n_win], p.area_arcsec2[:n_win], p.source, p.selection, False, p.note,
                                 {}, node_weights=p.node_weights[:n_win])
    rows = []
    for i, c in enumerate(centres):
        t0 = time.time()
        kin, cps = cd.kinematics(c), cd.count_products(c)
        hst = next(p for p in cps if p.name.startswith(cm.HST_NAME))
        row = dict(c=list(c), n_inner2=int(hst.counts[:2].sum()), counts_hst=hst.counts.tolist(),
                   hst_pmr=[float(v) for v in next(p for p in kin.profiles if p.name == "hst_pm_radial_ours").value[:cd.n_hst]],
                   hst_pmt=[float(v) for v in next(p for p in kin.profiles if p.name == "hst_pm_tangential_ours").value[:cd.n_hst]],
                   muse=[float(v) for v in next(p for p in kin.profiles if p.name == "muse_los_dispersion").value[:cd.n_muse]],
                   models={})
        w = sub(hst)
        for k, (cfg, *_rest) in runs.items():
            ev = DFJointProblem(kin, None, cps).evaluate(cfg, models[k])
            s = scores(ev)
            pc = cfg.distance_kpc*1e3/ARCSEC_PER_RAD
            sig = models[k].projected_moments(np.asarray(w.r_nodes).ravel()*pc)["Sigma"]
            s["lnl_binned_win"] = w.compare(sig)["loglike"]
            s["lnl_unbinned_win"] = ub.lnl(lnsig[k], c)
            s["P_inner2"] = float(poisson.cdf(row["n_inner2"], s["mu_inner2"]))
            row["models"][k] = s
        rows.append(row)
        print(f"[{i+1}/{len(centres)}] c=({c[0]:+.2f},{c[1]:+.2f}) {time.time()-t0:.1f}s "
              + " ".join(f"{k}:{row['models'][k]['lnl']:.2f}" for k in runs), flush=True)

    index = {tuple(r["c"]): r for r in rows}
    r0 = index[(0., 0.)]
    tests = {"T2_max_abs_diff": max(abs(r0["models"][k]["lnl"]-frozen[k]["lnl"]) for k in runs)}
    out = dict(runs=RUNS, frozen=frozen, grids={}, tests=tests, n_unbinned=ub.n, window_arcsec=[edges[0], edges[n_win]])
    for g, (C, W) in grids.items():
        res = {}
        for key in ("lnl", "lnl_binned_win", "lnl_unbinned_win"):
            L = np.array([[index[tuple(np.round(x, 12))]["models"][k][key] for x in C] for k in runs])
            res[key] = {k: float(v) for k, v in zip(runs, cm.marginal_loglike(L, W))}
            d = L-L[0]                                                   # per-centre difference to no BH
            res[key+"_dm2"] = {k: dict(marginal=float(-2*(cm.marginal_loglike(L[j], W)-cm.marginal_loglike(L[0], W))),
                                       mean=float(-2*W@d[j]), sd=float(2*np.sqrt(W@(d[j]-W@d[j])**2)),
                                       min=float(-2*d[j].max()), max=float(-2*d[j].min()))
                               for j, k in enumerate(runs)}
        res["P_inner2"] = {k: float(W@np.array([index[tuple(np.round(x, 12))]["models"][k]["P_inner2"] for x in C]))
                           for k in runs}
        out["grids"][g] = res
    if args.smoke:
        out["rows"] = rows
        print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1)); return
    out["tests"]["T3_delta_prior"] = max(abs(cm.marginal_loglike([r0["models"][k]["lnl"]], [1.])-frozen[k]["lnl"]) for k in runs)
    out["tests"]["T4_n10_minus_n6"] = {k: out["grids"]["s2_n10"]["lnl_dm2"][k]["marginal"]-out["grids"]["s2_n6"]["lnl_dm2"][k]["marginal"]
                                       for k in runs}
    # T7: HST PM bins 13.8-58.8" at the +-2 sigma_c (4") centres
    p0 = {p.name: p for p in data.profiles}
    nb = int(np.count_nonzero(p0["hst_pm_radial_ours"].r_upper <= 58.9))
    e7 = np.r_[p0["hst_pm_radial_ours"].r_lower[cd.n_hst:nb], p0["hst_pm_radial_ours"].r_upper[nb-1]]
    t7 = {}
    for c in ((4., 0.), (-4., 0.), (0., 4.), (0., -4.)):
        t = cm.hst_pm_table(e7, c)
        t7[str(c)] = {comp: ((np.asarray(t[col])-p0[name].value[cd.n_hst:nb])/p0[name].err_lo[cd.n_hst:nb]).round(3).tolist()
                      for comp, col, name in (("pmr", "sigma_pmr", "hst_pm_radial_ours"),
                                              ("pmt", "sigma_pmt", "hst_pm_tangential_ours"))}
    out["tests"]["T7_shift_over_err"] = t7
    out["tests"]["T7_max_abs"] = float(max(abs(v) for d in t7.values() for a in d.values() for v in a))
    out["rows"] = rows
    out["elapsed_s"] = time.time()-t_start
    path = ROOT/f"results/plot_data/centre_marginal_{args.label}.json"
    path.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k not in ("rows",)}, indent=1))
    print("written", path)


if __name__ == "__main__":
    main()
