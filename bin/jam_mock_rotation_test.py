#!/usr/bin/env python3
"""Step 4 of the rotation plan (2026-11-13): can the spherical Jeans pipeline recover a rotating, flattened cluster?

Truth: the no-DM count-tracer Jeans model (jeans_K1_turnover_counts_c, MAP sample) made axisymmetric with Cappellari's
JamPy 9 (jampy.axi.proj; velocity ellipsoid aligned in spherical coordinates, Cappellari 2026): the stellar/tracer
MGE and the remnant Plummer (as an NNLS MGE, backends.mge_from_component) flattened to observed q' (default 0.90),
the black hole as a point mass, the turnover beta(r) of the MAP sample as JAM's beta = 1 - <v_th^2>/<v_r^2>.
Inclination i = 44 deg, distance = MAP distance. JAM gives the projected second-moment tensor, which is independent
of how <v_phi^2> is split into mean rotation and dispersion (JamPy docs: gamma/kappa affect first moments only).
JAM's own (gamma, kappa) cannot reproduce the measured rotation-curve shape (tested: rising outwards), so the mean
field is the measured one: the MAP cylindrical rotation field of bin/fit_rotation_field.py (same tracer MGE,
q' = 0.90, i = 44 deg), projected with that script's code. The mean field only sets what each estimator subtracts.

Mock data: JAM first moments and full second-moment tensors at the real star positions (HST: oMEGACat flagged PM
stars; MUSE: selection_hq_los stars; Gaia: rings at the profile's radial nodes x 64 azimuths), averaged per real bin
with the estimator of each dataset:
  HST PM (local mean removed):       sigma^2 = <<v_k^2> - vbar_k^2>;  streaming2 (pmt) = <vbar_T^2>; pmr none
  MUSE LOS (ring sinusoid removed):  sigma^2 = <<v_z^2> - vbar_z^2> + <(vbar_z - s)^2>; streaming2 = <s^2>
  Gaia PM (one mean per ring):       sigma^2 = <v_k^2> - <vbar_k>^2;  streaming2 (pmt) = <vbar_T>^2; pmr none
with the real error bars (noiseless mock: any offset is bias). Then the spherical pipeline (FitProblem, same family,
bins, nodes, priors) is maximised for the no-DM (K1) and NFW families and compared with the truth.
Control: q' = 1, no rotation (kappa = 0): must reproduce the truth (tests JAM vs our Jeans solver).

Usage: PYTHONPATH=src:bin python bin/jam_mock_rotation_test.py [--qobs 0.9] [--qrem Q] [--rot-scale 1] [--control] [--tag NAME] [--tune-only]
Output: results/rotation/jam_mock_<tag>.json (git-ignored), printed summary.
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
from astropy.table import Table
from scipy.optimize import minimize

import jampy as jam
from ocen_dm.cluster import OCEN_DEC, OCEN_RA
from ocen_dm.kinematics import FitProblem, KinematicData
from ocen_dm.kinematics.backends import _split_mass_model, mge_from_component
from ocen_dm.kinematics.hst_profile import load_hst_sample
from ocen_dm.kinematics.report import load_run_metadata, problem_for
from ocen_dm.kinematics.run_io import family_config, family_from_config
from ocen_dm.paths import raw_dir
from plot_jeans_centre2 import map_point

ARCSEC_PER_RAD = 206264.806
PA = 104.2
INC = 44.0


def star_positions(data, rng, per_bin=250):
    """Per dataset: (x_east, y_north, bin index, node id) of sample points. HST/MUSE: real star positions (node -1,
    equal weights). Gaia: 64 azimuths at each radial node (equal-count quantiles of the stars); these are weighted by
    the flattened tracer surface density within each node (mock weights), as stars of a flattened cluster are."""
    hst = load_hst_sample(require_flag=True)
    hx, hy = hst.r_arcsec*np.sin(hst.phi), hst.r_arcsec*np.cos(hst.phi)
    cat = Table.read(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")
    s = np.asarray(cat["selection_hq_los"], bool) & np.isfinite(np.asarray(cat["vlos"], float))
    cosd = np.cos(np.radians(OCEN_DEC))
    mx = (np.asarray(cat["RA"], float)[s]-OCEN_RA)*cosd*3600; my = (np.asarray(cat["DEC"], float)[s]-OCEN_DEC)*3600
    out = {}
    for p in data.profiles:
        xs, ys, bins, nodes = [], [], [], []
        if p.instrument in ("HST", "MUSE"):
            X, Y = (hx, hy) if p.instrument == "HST" else (mx, my)
            R = np.hypot(X, Y)
            for k, (lo, hi) in enumerate(zip(p.r_lower, p.r_upper)):
                idx = np.flatnonzero((R >= lo) & (R < hi))
                if len(idx) > per_bin:
                    idx = rng.choice(idx, per_bin, replace=False)
                xs.append(X[idx]); ys.append(Y[idx]); bins.append(np.full(len(idx), k)); nodes.append(np.full(len(idx), -1))
        else:                                   # Gaia: complete rings at the bin's radial nodes
            ph = (np.arange(64)+.5)/64*2*np.pi
            for k in range(p.n):
                rr = np.asarray(p.r_nodes[k]) if p.r_nodes is not None else np.array([p.r[k]])
                xs.append((rr[:, None]*np.sin(ph)).ravel()); ys.append((rr[:, None]*np.cos(ph)).ravel())
                bins.append(np.full(xs[-1].size, k)); nodes.append(np.repeat(1000*k+np.arange(len(rr)), len(ph)))
        out[p.name] = (np.concatenate(xs), np.concatenate(ys), np.concatenate(bins), np.concatenate(nodes))
    return out


def flatten(m, s_pc, q, pc):
    """Flattened MGE for JamPy from a circular one: masses m [Msun or L], circular sigmas s_pc [pc], observed q'
    (scalar or per Gaussian). sigma_major = sigma/sqrt(q') keeps each Gaussian's total (2 pi surf q' sigma_maj^2 =
    2 pi surf sigma^2) and its circularised profile; surf is unchanged. (Fixed 2026-11-13: the first version kept
    sigma_major = sigma, which removed a fraction 1-q' of every Gaussian's mass and shrank the profile by sqrt(q').)"""
    q = np.broadcast_to(np.asarray(q, float), np.shape(s_pc)).copy()
    return m/(2*np.pi*s_pc**2), s_pc/np.sqrt(q)/pc, q


def mge_total(surf, sig_as, q, pc):
    return float(np.sum(2*np.pi*surf*q*(sig_as*pc)**2))


def jam_models(jeans, D, qobs, qrem):
    """(lum, pot, mbh) for JamPy. qobs: tracer/stellar q' (scalar or one per stellar Gaussian); qrem: remnant q'."""
    pc = D*1e3/ARCSEC_PER_RAD
    gaussians, mbh, other = _split_mass_model(jeans.mass)
    stars = gaussians[0]
    rem = [mge_from_component(c, np.geomspace(0.05, 500., 40)) for c in other]
    ms, ss = np.asarray(stars.masses)*stars.mass_to_light, np.asarray(stars.sigmas)
    parts = [flatten(ms, ss, qobs, pc)]
    for g in rem:
        parts.append(flatten(np.asarray(g.masses)*g.mass_to_light, np.asarray(g.sigmas), qrem, pc))
    pot = tuple(np.concatenate([pp[i] for pp in parts]) for i in range(3))
    want = float(np.sum(ms))+sum(float(np.sum(np.asarray(g.masses)*g.mass_to_light)) for g in rem)
    got = mge_total(*pot, pc)
    assert abs(got/want-1) < 1e-10, (got, want)
    # JamPy 9 limits its LOS quadrature and intrinsic grid to the largest tracer sigma (jam_axi_proj.vmom_proj:
    # rmax = max(sigma_lum)); beyond it the moments are wrong (tested 2026-11-13 against an exact Gaussian-tracer /
    # point-mass case: flat 14.7 km/s instead of 7.7 -> 4.4). A negligible tracer Gaussian at 5x the largest sigma,
    # 1e-10 of the peak, extends the grid; with it JAM matches the exact case to 0.001 km/s at 300-2150".
    sl, sg, qq = parts[0]
    lum = (np.r_[sl, 1e-10*sl.max()], np.r_[sg, 5*sg.max()], np.r_[qq, qq[-1]])
    return lum, pot, mbh


def surface_density(lum, xp, yp):
    """Projected tracer density of a flattened MGE at (x' major, y' minor) [arcsec]."""
    surf, sig, q = lum
    return np.sum(surf[:, None]*np.exp(-0.5*(xp[None, :]**2+(yp[None, :]/q[:, None])**2)/sig[:, None]**2), axis=0)


class MeanField:
    """Projected mean velocity (v_x', v_y', v_z) [km/s] from the fitted rotation field, in the (x', y') frame."""
    def __init__(self, scale=1.0, params=None):
        """params: a rotation-field parameter vector (nodes, incl, pa, vsys); default the fit's MAP."""
        import fit_rotation_field as F
        rot = json.loads((ROOT/"results/rotation/rotation_field_cylindrical_20261113.json").read_text())
        p = np.array(rot["map"] if params is None else params); nn = len(F.NODES)
        mlos, vt, vr = F.Field("cylindrical").project(scale*p[:nn], p[-3])
        self.F = F
        g = lambda a: np.c_[a, a[:, :1]]
        from scipy.interpolate import RegularGridInterpolator as RGI
        ax = (np.log(F.RP), np.r_[F.PHI, 2*np.pi])
        self.ilos, self.it, self.ir = (RGI(ax, g(a)) for a in (mlos, vt, vr))

    def __call__(self, xp, yp):
        F = self.F
        r = np.clip(np.hypot(xp, yp), F.RP[0], F.RP[-1]); ph = np.mod(np.arctan2(yp, xp), 2*np.pi)
        q = np.column_stack([np.log(r), ph])
        vlos, vt, vr = self.ilos(q), self.it(q), self.ir(q)
        return np.array([vr*np.cos(ph)-vt*np.sin(ph), vr*np.sin(ph)+vt*np.cos(ph), vlos])


def run_jam(lum, pot, mbh, D, anis, x, y, gamma, kappa, inc):
    a = np.radians(PA)
    xp = x*np.sin(a)+y*np.cos(a)
    yp = -x*np.cos(a)+y*np.sin(a)
    beta = lambda r, th: (anis.beta(r), 0.)
    out = jam.axi.proj(lum, pot, inc, mbh, D/1e3, xp, yp, beta=beta, gamma=gamma, kappa=kappa, ml=1, moment="zz",
                       interp=False, plot=False, quiet=True, nang=41, nlos=400)
    vel = MEAN(xp, yp) if MEAN is not None else 0*out.vel
    return xp, yp, vel, out.vel2


MEAN = None


def weights(lum, xp, yp, nodes):
    w = np.ones(len(xp))
    g = nodes >= 0
    if np.any(g):
        sd = surface_density(lum, xp[g], yp[g])
        u, inv = np.unique(nodes[g], return_inverse=True)
        w[g] = sd/np.bincount(inv, sd)[inv]*np.bincount(inv)[inv]
    return w


def estimator(p, xp, yp, bins, vel, vel2, D, w=None):
    """Binned mock dispersion and streaming2 for one dataset, in the dataset's units (weighted means; w = None: equal).
    HST (local mean removed): sigma^2 = <<v^2> - vbar^2>, streaming2 = <vbar^2> for BOTH components (the radial local
    mean is small but removed too). MUSE: ring sinusoid removed. Gaia: one mean per ring removed."""
    conv = 1.0 if p.kind == "los" else 1/(4.74047*D)
    R = np.hypot(xp, yp)
    eR = np.array([xp/R, yp/R]); eT = np.array([-yp/R, xp/R])
    if p.kind == "los":
        m1, m2 = vel[2], vel2[2, 2]
    else:
        e = eR if p.kind == "pmr" else eT
        m1 = e[0]*vel[0]+e[1]*vel[1]
        m2 = e[0]**2*vel2[0, 0]+2*e[0]*e[1]*vel2[0, 1]+e[1]**2*vel2[1, 1]
    w = np.ones(len(xp)) if w is None else w
    sig2, st2 = np.empty(p.n), np.zeros(p.n)
    for k in range(p.n):
        j = bins == k
        mean = lambda a: np.sum(w[j]*a)/np.sum(w[j])
        if p.instrument == "HST":
            sig2[k] = mean(m2[j]-m1[j]**2); st2[k] = mean(m1[j]**2)
        elif p.kind == "los":
            ph = np.arctan2(yp[j], xp[j])
            A = np.column_stack([np.cos(ph), np.sin(ph)])
            sw = np.sqrt(w[j])
            c, *_ = np.linalg.lstsq(A*sw[:, None], m1[j]*sw, rcond=None)
            s = A @ c
            sig2[k] = mean(m2[j]-m1[j]**2)+mean((m1[j]-s)**2); st2[k] = mean(s**2)
        else:
            sig2[k] = mean(m2[j])-mean(m1[j])**2; st2[k] = mean(m1[j])**2
    sig = np.sqrt(sig2)*conv
    st2 = st2*conv**2
    keep_streaming = p.kind in ("pmt", "los") or p.instrument == "HST"
    return sig, (st2 if keep_streaming else None)


def fit(pb, x0, n_restart=3, rng=None):
    best = None
    lo = np.array([pp.prior.transform(1e-6) for pp in pb.family.parameters])
    hi = np.array([pp.prior.transform(1-1e-6) for pp in pb.family.parameters])
    f = lambda x: -pb.loglike_vector(np.clip(x, lo, hi)) if np.all(np.isfinite(x)) else 1e100
    starts = [x0]+[np.clip(x0*(1+0.05*rng.standard_normal(len(x0))), lo, hi) for _ in range(n_restart-1)]
    for s in starts:
        r = minimize(f, s, method="Powell", options=dict(maxiter=40000, xtol=1e-5, ftol=1e-9))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-6, fatol=1e-9))
        if best is None or r.fun < best.fun:
            best = r
    return np.clip(best.x, lo, hi), -best.fun


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qobs", type=float, default=0.90)
    ap.add_argument("--qrem", type=float, default=None, help="remnant q' (default = qobs)")
    ap.add_argument("--control", action="store_true")
    ap.add_argument("--rot-scale", type=float, default=1.0, help="scale the measured rotation (sensitivity test)")
    ap.add_argument("--tag", default="default")
    ap.add_argument("--tune-only", action="store_true")
    ap.add_argument("--mock-only", action="store_true", help="build and compare the mock, no spherical refits")
    args = ap.parse_args()
    rng = np.random.default_rng(1)
    run = load_run_metadata(ROOT/"results/fits/jeans_K1_turnover_counts_c")
    pb0 = problem_for(run)
    fam = pb0.family
    x_true, lnl_true, _ = map_point("jeans_K1_turnover_counts_c", fam)
    theta = dict(zip(fam.names, x_true))
    jeans, D, _ = fam.build(theta)
    qobs = 1.0 if args.control else args.qobs     # --qobs 1 --rot-scale 1: rotation without flattening
    qrem = qobs if args.qrem is None else args.qrem
    inc = 90.0 if args.control else INC
    lum, pot, mbh = jam_models(jeans, D, qobs, qrem)
    pos = star_positions(pb0.data, rng)
    rot = json.loads((ROOT/"results/rotation/rotation_field_cylindrical_20261113.json").read_text())

    def mock(gamma, kappa):
        res = {}
        for p in pb0.data.profiles:
            x, y, b, nd = pos[p.name]
            xp, yp, vel, vel2 = run_jam(lum, pot, mbh, D, jeans.anisotropy, x, y, gamma, kappa, inc)
            res[p.name] = (p, xp, yp, b, vel, vel2, weights(lum, xp, yp, nd))
        return res

    def rotation_curves(res):
        """MUSE ring LOS amplitude and ring-mean |PM_T| (km/s) at the radii of the rotation-field data."""
        p, xp, yp, b, vel, _, _ = res["muse_los_dispersion"]
        A = []
        for k in range(p.n):
            j = b == k; ph = np.arctan2(yp[j], xp[j])
            c, *_ = np.linalg.lstsq(np.column_stack([np.cos(ph), np.sin(ph)]), vel[2][j], rcond=None)
            A.append(np.hypot(*c))
        ring_r = np.array(rot["pm_r"])
        ph = (np.arange(64)+.5)/64*2*np.pi
        xx = (ring_r[:, None]*np.sin(ph)).ravel(); yy = (ring_r[:, None]*np.cos(ph)).ravel()   # (diagnostic only)
        xp2, yp2, v, _ = run_jam(lum, pot, mbh, D, jeans.anisotropy, xx, yy, gam, kap, inc)
        R = np.hypot(xp2, yp2)
        vt = (-yp2*v[0]+xp2*v[1])/R
        return np.array(A), np.abs(vt.reshape(len(ring_r), -1).mean(axis=1))

    t0 = time.perf_counter()
    global MEAN
    gam, kap = 0.0, 0.0
    MEAN = None if (args.control or args.rot_scale == 0) else MeanField(args.rot_scale)
    res = mock(gam, kap)
    print(f"JAM mock built in {time.perf_counter()-t0:.1f} s (q' {qobs}, i {inc}, rotation {'none' if MEAN is None else args.rot_scale})")
    if not args.control:
        A, ring = rotation_curves(res)
        pm_v = np.array(rot["pm_v"])
        muse = {p.name: p for p in pb0.data.profiles}["muse_los_dispersion"]
        print("PM ring rotation  data:", np.round(pm_v, 2).tolist())
        print("PM ring rotation  mock:", np.round(ring, 2).tolist())
        print("MUSE ring amplitude (published v_rot):", np.round(np.sqrt(2*muse.streaming2), 2).tolist())
        print("MUSE ring amplitude mock             :", np.round(A, 2).tolist())
        if args.tune_only:
            return

    profiles = []
    for p in pb0.data.profiles:
        _, xp, yp, b, vel, vel2, w = res[p.name]
        sig, st2 = estimator(p, xp, yp, b, vel, vel2, D, w)
        profiles.append(replace(p, value=sig, streaming2=st2))
        tot = sig**2+(0 if st2 is None else st2); real = p.value**2+(0 if p.streaming2 is None else p.streaming2)
        d = (np.sqrt(tot)-np.sqrt(real))/p.err_lo
        print(f"  mock vs real data {p.name:26s} (mock-data)/err mean {d.mean():+6.2f} rms {np.sqrt(np.mean(d**2)):5.2f}")
    mdata = KinematicData(tuple(profiles))
    out = dict(tag=args.tag, created=time.strftime("%Y-%m-%dT%H:%M:%S"), qobs=qobs, qrem=qrem, inc=inc,
               rot_scale=None if MEAN is None else args.rot_scale, distance=D,
               truth=dict(zip(fam.names, map(float, x_true))), fits={},
               mock={p.name: dict(value=p.value.tolist(), streaming2=None if p.streaming2 is None else
                                  np.asarray(p.streaming2).tolist()) for p in profiles})
    if args.mock_only:
        (ROOT/f"results/rotation/jam_mock_{args.tag}.json").write_text(json.dumps(out, indent=1))
        return
    for label, cfg_run in (("K1", "jeans_K1_turnover_counts_c"), ("NFW", "jeans_K2_nfw_turnover_counts_c")):
        r = load_run_metadata(ROOT/"results/fits"/cfg_run)
        family = family_from_config(family_config(problem_for(r).family))
        pb = FitProblem(family, mdata)
        x0 = np.array([theta.get(n, np.nan) for n in family.names])
        if label == "NFW":
            x0[family.names.index("M_dm_100")] = 2e3
            x0[family.names.index("r_s")] = 100.
        lnl0 = pb.loglike_vector(x0)
        xb, lnlb = fit(pb, x0, rng=rng)
        chi2 = pb.chi2(xb)
        out["fits"][label] = dict(names=list(family.names), x=xb.tolist(), lnL=lnlb, lnL_truth=lnl0,
                                  chi2={k: v[0] for k, v in chi2.items()}, chi2_total=sum(v[0] for v in chi2.values()))
        print(f"\n{label}: lnL at truth {lnl0:.2f}, best {lnlb:.2f}; chi2 {out['fits'][label]['chi2_total']:.2f}",
              {k: round(v[0], 2) for k, v in chi2.items()})
        for n, v in zip(family.names, xb):
            t = theta.get(n)
            print(f"   {n:14s} {v:12.5g}   truth {'' if t is None else f'{t:.5g}'}")
    if "NFW" in out["fits"]:
        print(f"\nlnL(NFW) - lnL(K1) = {out['fits']['NFW']['lnL']-out['fits']['K1']['lnL']:.2f}")
    (ROOT/"results/rotation").mkdir(parents=True, exist_ok=True)
    (ROOT/f"results/rotation/jam_mock_{args.tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
