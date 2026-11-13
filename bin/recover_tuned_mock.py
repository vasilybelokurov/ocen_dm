#!/usr/bin/env python3
"""Step 5 of the flattened-mock plan (2026-11-13): fit the SPHERICAL Jeans model to a tuned flattened JAM mock.

Truth: a tuned no-DM flattened, rotating JAM model (bin/tune_jam_to_data.py output, results/rotation/jam_tuned_<tag>.json).
Mock data (noiseless): the truth's binned dispersions as the real estimators would measure them. With the measured
rotation field as the mock's mean field, each estimator removes exactly the data's rotation term, so the mock value is
sqrt(<m2>_bin - streaming2) = the saved tuned prediction (all stars, Sigma-weighted Gaia rings); errors and streaming2
are those of the corrected real data (tune_jam_to_data.Setup). Noisy mocks (--noise N): value + Gaussian noise with the
data errors (split normal sides chosen at random with probability err_hi/(err_lo+err_hi)), seeds 0..N-1.
Fits: production spherical Jeans families K1 (no DM) and NFW (jeans_K1_turnover_counts_c, jeans_K2_nfw_turnover_counts_c),
max likelihood (Powell + Nelder-Mead polish, 3 starts: truth and two perturbations), the same as jam_mock_rotation_test.fit.
Compared with truth: M(<r) at 1/5/20/50 pc -- the truth's is the mass of the OBLATE model inside spheres (deprojected
MGE + BH), not the round parent -- total mass, M_bh, D, and spherical-equivalent beta at 4/15/25/40 pc (truth: shell
average of JamPy's intrinsic moments, beta_sph = 1-(s_th^2+<vphi^2>)/(2 s_r^2)). Bias is quoted also in units of the
production K1 posterior s.d. (weighted_post.txt, old rotation terms: a representative width only).

Usage: PYTHONPATH=src:bin python bin/recover_tuned_mock.py --truth q0.90_r0.90_i44_from_q0.85_r0.85_i44_v1 [--noise 0]
Output: results/rotation/recover_<truth>[_noise<N>].json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src")); sys.path.insert(0, str(ROOT/"bin"))
import numpy as np
from jampy.axi.jam_axi_intr import intrinsic_moments

import jam_mock_rotation_test as J
from ocen_dm.kinematics import FitProblem, KinematicData
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config
from tune_jam_to_data import Setup

R_MASS = np.array([1., 5., 20., 50.])          # pc
R_BETA = np.array([4., 15., 25., 40.])         # pc
FAMILIES = (("K1", "jeans_K1_turnover_counts_c"), ("NFW", "jeans_K2_nfw_turnover_counts_c"))


def oblate_enclosed(pot, mbh, D, inc_deg, r_pc):
    """Mass inside spheres of radius r [pc] for a flattened mass MGE (JamPy convention: surf [Msun/pc^2], sigma [arcsec]
    major axis, observed q') at inclination i, plus the BH. Each Gaussian: M = 2 pi surf sigma^2 q' (pc units), intrinsic
    q from q'^2 = cos^2 i + q^2 sin^2 i; the shell average of exp(-s^2 (1-mu^2 + mu^2/q^2)/2 sigma^2) over mu = cos theta."""
    surf, sig, qp = map(np.asarray, pot)
    pc = D*1e3*np.pi/648000.
    i = np.radians(inc_deg)
    sp = sig*pc
    qi = np.sqrt(np.clip(qp**2-np.cos(i)**2, 1e-12, None))/np.sin(i) if inc_deg < 89.999 else qp
    M = 2*np.pi*surf*sp**2*qp
    mu = np.linspace(0, 1, 201)
    out = []
    for r in np.atleast_1d(r_pc):
        s = np.r_[0., np.geomspace(1e-4, r, 3000)]
        tot = 0.
        for Mk, sk, qk in zip(M, sp, qi):
            rho0 = Mk/((2*np.pi)**1.5*sk**3*qk)
            e = np.exp(-0.5*s[:, None]**2*(1-mu[None]**2+mu[None]**2/qk**2)/sk**2)
            shell = np.trapezoid(e, mu, axis=1)                    # <...>_mu over [0,1] = full-sphere average
            tot += np.trapezoid(4*np.pi*s**2*rho0*shell, s)
        out.append(tot+mbh)
    return np.array(out)


def beta_sph_truth(lum, pot, mbh, D, inc_deg, anis, r_pc):
    i = np.radians(inc_deg); pc = D*1e-3*np.pi/0.648
    def deproj(m):
        s, sg, q = map(np.asarray, m); qi = np.sqrt(q**2-np.cos(i)**2)/np.sin(i); sgp = sg*pc
        return (s*q/(sgp*qi*np.sqrt(2*np.pi)), sgp, qi)
    L, P = deproj(lum), deproj(pot)
    th = np.linspace(0.02, np.pi/2-0.02, 40)
    R = (r_pc[:, None]*np.sin(th)[None]).ravel(); z = (r_pc[:, None]*np.cos(th)[None]).ravel()
    s2r, s2t, _, v2p, nu = [x.reshape(len(r_pc), len(th)) for x in intrinsic_moments(
        R, z, L, P, mbh, lambda r, t: (anis.beta(r), 0.), 0., None, 31)]
    w = nu*np.sin(th)[None]; av = lambda x: np.sum(w*x, 1)/np.sum(w, 1)
    return 1-(av(s2t)+av(v2p))/(2*av(s2r))


def summarise(family, x):
    jeans, D, _ = family.build(dict(zip(family.names, x)))
    th = dict(zip(family.names, x))
    return dict(M_r=jeans.mass.enclosed_mass(R_MASS).tolist(), M_tot=float(np.ravel(jeans.mass.enclosed_mass(1e4))[0]),
                M_bh=float(th["M_bh"]), D=float(D), beta=jeans.anisotropy.beta(R_BETA).tolist(),
                M_dm_100=float(th.get("M_dm_100", 0.)))


def posterior_sd(n=300):
    """Representative widths: s.d. over the production K1 posterior of the summary quantities."""
    run = load_run_metadata(ROOT/"results/fits/jeans_K1_turnover_counts_c")
    fam = problem_for(run).family
    d = ROOT/"results/fits/jeans_K1_turnover_counts_c/ultranest/chains/weighted_post.txt"
    names = open(d).readline().split(); rows = np.loadtxt(d, skiprows=1)
    w = rows[:, names.index("weight")]
    k = np.random.default_rng(3).choice(len(rows), n, p=w/w.sum())
    ss = [summarise(fam, np.array([rows[j, names.index(nm)] for nm in fam.names])) for j in k]
    f = lambda key: np.std([s[key] for s in ss], axis=0).tolist()
    return dict(M_r=f("M_r"), M_tot=f("M_tot"), M_bh=f("M_bh"), D=f("D"), beta=f("beta"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--noise", type=int, default=0, help="number of noisy realisations (0 = noiseless only)")
    args = ap.parse_args()
    t0 = time.perf_counter()
    T = json.loads((ROOT/f"results/rotation/jam_tuned_{args.truth}.json").read_text())
    S = Setup(per_bin=10)                       # data (corrected terms and errors) only; predictions are saved
    theta = dict(zip(T["names"], T["x"]))
    jeans, D, _ = S.fam.build(theta)
    lum, pot, mbh = J.jam_models(jeans, D, T["qobs"], T["qrem"])
    truth = dict(M_r=oblate_enclosed(pot, mbh, D, T["inc"], R_MASS).tolist(),
                 M_tot=float(J.mge_total(*pot, D*1e3*np.pi/648000.)+mbh), M_bh=float(mbh), D=float(D),
                 beta=beta_sph_truth(lum, pot, mbh, D, T["inc"], jeans.anisotropy, R_BETA).tolist(),
                 round_parent_M_r=jeans.mass.enclosed_mass(R_MASS).tolist())
    sd = posterior_sd()
    print(f"truth {args.truth}: M(<1,5,20,50 pc) {np.round(np.array(truth['M_r'])/1e6, 4)} (round parent "
          f"{np.round(np.array(truth['round_parent_M_r'])/1e6, 4)}) 1e6; M_bh {mbh:.0f}; D {D:.3f}; "
          f"beta_sph(4,15,25,40 pc) {np.round(truth['beta'], 2)}", flush=True)

    rng = np.random.default_rng(0)
    real = (["noiseless"]+[f"seed{k}" for k in range(args.noise)])
    out = dict(truth_tag=args.truth, truth=truth, posterior_sd=sd, r_mass_pc=R_MASS.tolist(), r_beta_pc=R_BETA.tolist(),
               created=time.strftime("%Y-%m-%dT%H:%M:%S"), fits={})
    fams = {}
    for label, cfg in FAMILIES:
        r = load_run_metadata(ROOT/"results/fits"/cfg)
        fams[label] = family_from_config(family_config(problem_for(r).family))
    for name in real:
        profs = []
        for p in S.data.profiles:
            v = np.array(T["pred"][p.name])
            if name != "noiseless":
                k = int(name[4:]); g = np.random.default_rng(1000+k)
                z = np.abs(g.standard_normal(p.n))
                up = g.random(p.n) < p.err_hi/(p.err_lo+p.err_hi)
                v = v+np.where(up, z*p.err_hi, -z*p.err_lo)
            profs.append(replace(p, value=v))
        mdata = KinematicData(tuple(profs))
        res = {}
        for label, family in fams.items():
            pb = FitProblem(family, mdata)
            x0 = np.array([theta.get(n, np.nan) for n in family.names])
            if label == "NFW":
                x0[family.names.index("M_dm_100")] = 2e3
                x0[family.names.index("r_s")] = 100.
            xb, lnl = J.fit(pb, x0, rng=rng)
            chi2 = pb.chi2(xb)
            res[label] = dict(names=list(family.names), x=xb.tolist(), lnL=float(lnl),
                              chi2={k: v[0] for k, v in chi2.items()}, chi2_total=sum(v[0] for v in chi2.values()),
                              summary=summarise(family, xb))
        res["dlnL_NFW_K1"] = res["NFW"]["lnL"]-res["K1"]["lnL"]
        out["fits"][name] = res
        s = res["K1"]["summary"]
        bias = lambda key: ((np.array(s[key])-np.array(truth[key]))/np.array(sd[key]))
        print(f"{name}: K1 chi2 {res['K1']['chi2_total']:.1f}  NFW chi2 {res['NFW']['chi2_total']:.1f}  "
              f"dlnL(NFW-K1) {res['dlnL_NFW_K1']:+.2f}  M_dm_100 {res['NFW']['summary']['M_dm_100']:.3g}  "
              f"({time.perf_counter()-t0:.0f} s)", flush=True)
        print(f"   K1 M(<1,5,20,50) / truth - 1: {np.round(np.array(s['M_r'])/np.array(truth['M_r'])-1, 3)}  "
              f"[in post. s.d.: {np.round(bias('M_r'), 2)}]")
        print(f"   K1 M_tot/truth-1 {s['M_tot']/truth['M_tot']-1:+.3f} [{bias('M_tot'):+.2f} sd]; "
              f"M_bh/truth-1 {s['M_bh']/truth['M_bh']-1:+.3f} [{bias('M_bh'):+.2f} sd]; "
              f"D-truth {s['D']-truth['D']:+.3f} kpc [{bias('D'):+.2f} sd]")
        print(f"   K1 beta(4,15,25,40) {np.round(s['beta'], 2)} vs truth beta_sph {np.round(truth['beta'], 2)}",
              flush=True)
    tag = args.truth+(f"_noise{args.noise}" if args.noise else "")
    (ROOT/f"results/rotation/recover_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
