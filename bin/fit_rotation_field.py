#!/usr/bin/env python3
"""Step 2 of the rotation plan (2026-11-13): one axisymmetric rotation field fitted to LOS and PM data together.

Model. Oblate tracer: the count MGE (data/processed/kinematics/tracer_mge_counts_c.json), every Gaussian with the same
intrinsic axis ratio q, set by the observed axis ratio q' = 0.90 and the inclination i: q^2 = (q'^2 - cos^2 i)/sin^2 i.
Rotation: cylindrical, v_phi(R_cyl) independent of z, a monotone cubic (PCHIP) through nodes at R = 30, 90, 200, 450,
1000, 2200" (v = 0 at R = 0, linear inside the first node, constant beyond the last). Projection: tracer-weighted
mean along each line of sight of v_los = v_phi (X/R) sin i and of the two sky components, then the tangential sky
component about the centre. Parameters: 6 nodes, i, PA (position angle of the major axis / maximum LOS velocity,
N through E), v_sys. Distance fixed at 5.43 kpc (PM <-> km/s).

Data.
- MUSE per-star LOS velocities (oMEGACat VI catalog_and_selections.fits, selection_hq_los): v_i ~ N(v_sys +
  vbar_los(x_i, y_i), sigma_los(R_i)^2 + e_i^2), sigma_los(R) interpolated from the published oMEGACat VI profile
  (dispersion about the rotating mean, consistent with this likelihood).
- Inner PM rotation: oMEGACat II Table 10 (Haeberle et al. 2024, arXiv:2404.03722; Gaia FPR minus HST), ring-averaged
  mean tangential PM, 12 rings to 360" (reproduced by bin/measure_inner_pm_rotation.py, 12/12 within 2 sigma).
- Outer PM rotation: our Gaia EDR3 mixture fits, |mean_pmt| per annulus (outer_gaia.load_edr3_profile), 9 rings
  357-2149", error sigma_pmt / sqrt((1-f_field) N) (as in likelihood.py, rotation="ours").
Model ring averages: azimuthal average of the projected mean tangential PM at the ring's median radius.
Fit: maximum likelihood (Nelder-Mead / Powell from several starts), then emcee for uncertainties.

Usage: PYTHONPATH=src python bin/fit_rotation_field.py [--steps N] [--law cylindrical|spherical]
Output: results/rotation/rotation_field_<law>_20261113.json, plots/rotation_field_<law>_20261113.png
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import numpy as np
from astropy.table import Table
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator
from scipy.optimize import minimize

from ocen_dm.cluster import OCEN_DEC, OCEN_RA
from ocen_dm.paths import processed_dir, raw_dir

KMS = 4.74047*5.43
QPROJ = 0.90
NODES = np.array([30., 90., 200., 450., 1000., 2200.])
T10 = np.array([[20.66, .0548, .0324], [47.6, .086, .0181], [76.65, .1478, .0131], [105.91, .2088, .0105],
                [135.99, .2556, .0093], [165.6, .2562, .0083], [195.11, .2678, .008], [225.38, .2648, .0078],
                [254.06, .2797, .0083], [282.86, .2619, .0103], [309.92, .2658, .0176], [340.91, .2455, .0421]])
NAMES = [f"v{int(n)}" for n in NODES]+["incl", "pa", "vsys"]

# projection grid: sky polar (Rp, phi) and line-of-sight coordinate s
RP = np.geomspace(2., 3000., 56)
PHI = np.linspace(0, 2*np.pi, 49)[:-1]
U = np.linspace(-1, 1, 161)
S_OF_U = np.sinh(9.*U)                    # s = scale * sinh(9 u): dense near s = 0
DS_DU = 9.*np.cosh(9.*U)


class Field:
    def __init__(self, law="cylindrical"):
        mge = json.loads((processed_dir()/"kinematics"/"tracer_mge_counts_c.json").read_text())
        keep = np.asarray(mge["fractions"]) > 1e-6
        self.sig = np.asarray(mge["sigmas_arcsec"])[keep]
        self.frac = np.asarray(mge["fractions"])[keep]
        self.law = law

    def project(self, vnodes, incl):
        """Projected mean fields on the (RP, PHI) grid in the frame x' = major axis, y' = minor axis.
        Returns vlos, v_tan (tangential sky component about the centre), in km/s."""
        i = np.radians(incl)
        q2 = (QPROJ**2-np.cos(i)**2)/np.sin(i)**2
        if not 0.05 < q2 <= 1.0:
            return None
        q = np.sqrt(q2)
        xp = (RP[:, None]*np.cos(PHI)[None, :])[..., None]
        yp = (RP[:, None]*np.sin(PHI)[None, :])[..., None]
        scale = max(RP.max(), 3*self.sig.max())/np.sinh(9.)
        s = scale*S_OF_U[None, None, :]
        w_s = scale*DS_DU[None, None, :]
        # LOS l = (0, sin i, cos i); e1 = (1,0,0); e2 = (0, cos i, -sin i); point = x' e1 + y' e2 + s l
        X = xp
        Y = yp*np.cos(i)+s*np.sin(i)
        Z = -yp*np.sin(i)+s*np.cos(i)
        Rc = np.sqrt(X*X+Y*Y)
        m2 = Rc*Rc+Z*Z/q2
        rho = np.zeros(np.broadcast_shapes(X.shape, s.shape))
        for sg, f in zip(self.sig, self.frac):
            rho = rho+f/(sg**3*q)*np.exp(-0.5*m2/sg**2)
        rr = Rc if self.law == "cylindrical" else np.sqrt(Rc*Rc+Z*Z)
        v = vcurve(vnodes, rr)
        with np.errstate(invalid="ignore", divide="ignore"):
            cx, cy = np.where(Rc > 0, -Y/Rc, 0.), np.where(Rc > 0, X/Rc, 0.)   # v_phi direction in (X, Y)
        vX, vY = v*cx, v*cy
        vlos = vY*np.sin(i)
        v1 = vX                                   # along e1 (x')
        v2 = vY*np.cos(i)                         # along e2 (y'); Z-velocity is zero
        wsum = np.sum(rho*w_s, axis=2)
        mean = lambda a: np.sum(a*rho*w_s, axis=2)/wsum
        mlos, m1, m2_ = mean(vlos), mean(v1), mean(v2)
        ph = PHI[None, :]
        vtan = -m1*np.sin(ph)+m2_*np.cos(ph)       # tangential (counter-clockwise in x', y')
        vrad = m1*np.cos(ph)+m2_*np.sin(ph)
        return mlos, vtan, vrad


def vcurve(vnodes, r):
    r = np.asarray(r, float)
    p = PchipInterpolator(np.r_[0., NODES], np.r_[0., vnodes], extrapolate=False)
    out = p(np.clip(r, 0, NODES[-1]))
    return np.where(r > NODES[-1], vnodes[-1], out)


def load_data():
    hst = Table.read(raw_dir()/"omegacat_vi_kinematics"/"catalog_and_selections.fits")
    sel = np.asarray(hst["selection_hq_los"], bool) & np.isfinite(np.asarray(hst["vlos"], float)) \
        & np.isfinite(np.asarray(hst["vlos_err"], float))
    cosd = np.cos(np.radians(OCEN_DEC))
    x = (np.asarray(hst["RA"], float)[sel]-OCEN_RA)*cosd*3600      # East
    y = (np.asarray(hst["DEC"], float)[sel]-OCEN_DEC)*3600          # North
    v, e = np.asarray(hst["vlos"], float)[sel], np.asarray(hst["vlos_err"], float)[sel]
    disp = Table.read(processed_dir()/"kinematics"/"omegacat_vi_los_dispersion.ecsv")
    R = np.hypot(x, y)
    sig = np.interp(R, np.asarray(disp["r_median"], float), np.asarray(disp["sigma_los"], float))
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    g = load_edr3_profile()
    gerr = np.asarray(g["sigma_pmt"])/np.sqrt(np.maximum((1-np.asarray(g["f_field"]))*np.asarray(g["n_stars"]), 1))
    return dict(x=x, y=y, v=v, var=sig**2+e**2, R=R,
                pm_r=np.r_[T10[:, 0], np.asarray(g["r_median"], float)],
                pm_v=np.r_[T10[:, 1], np.abs(np.asarray(g["mean_pmt"], float))]*KMS,
                pm_e=np.r_[T10[:, 2], gerr]*KMS,
                pm_src=np.r_[np.zeros(len(T10)), np.ones(len(g))])


class Problem:
    def __init__(self, data, field):
        self.d, self.f = data, field
        self.lo = np.r_[np.zeros(len(NODES)), 20., 60., 200.]
        self.hi = np.r_[np.full(len(NODES), 20.), 89., 150., 260.]

    def model(self, p):
        vn, incl, pa, vsys = p[:len(NODES)], p[-3], p[-2], p[-1]
        out = self.f.project(vn, incl)
        if out is None:
            return None
        mlos, vtan, vrad = out
        # star positions -> (x', y') frame: x' along PA (N through E), y' = PA + 90
        a = np.radians(pa)
        xp = self.d["x"]*np.sin(a)+self.d["y"]*np.cos(a)
        yp = -self.d["x"]*np.cos(a)+self.d["y"]*np.sin(a)
        rp = np.clip(np.hypot(xp, yp), RP[0], RP[-1])
        ph = np.mod(np.arctan2(yp, xp), 2*np.pi)
        phig = np.r_[PHI, 2*np.pi]
        grid = np.c_[mlos, mlos[:, :1]]
        vl = RegularGridInterpolator((np.log(RP), phig), grid)(np.column_stack([np.log(rp), ph]))
        ring = np.interp(np.log(self.d["pm_r"]), np.log(RP), np.mean(vtan, axis=1))
        ringr = np.interp(np.log(self.d["pm_r"]), np.log(RP), np.mean(vrad, axis=1))
        return vsys+vl, np.abs(ring), ringr

    def lnpost(self, p):
        if np.any(p < self.lo) or np.any(p > self.hi):
            return -np.inf
        m = self.model(p)
        if m is None:
            return -np.inf
        vl, ring, _ = m
        l1 = -0.5*np.sum((self.d["v"]-vl)**2/self.d["var"])
        l2 = -0.5*np.sum((self.d["pm_v"]-ring)**2/self.d["pm_e"]**2)
        return l1+l2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=1500)
    ap.add_argument("--law", default="cylindrical", choices=("cylindrical", "spherical"))
    args = ap.parse_args()
    data = load_data()
    pb = Problem(data, Field(args.law))
    print(f"MUSE stars {len(data['v'])}; PM rings {len(data['pm_r'])}")
    t0 = time.perf_counter(); pb.lnpost(np.r_[2, 4, 7, 6, 3, 1, 45, 103, 233]); print(f"{time.perf_counter()-t0:.3f} s/call")
    best = None
    for start in ([2, 4, 7, 6, 3, 1, 45, 103, 233], [1, 3, 6, 7, 4, 1, 55, 100, 233], [2, 5, 7, 5, 2, .5, 38, 105, 233]):
        r = minimize(lambda p: -pb.lnpost(p), np.array(start, float), method="Powell",
                     options=dict(maxiter=20000, xtol=1e-4, ftol=1e-6))
        print("start", start, "->", np.round(r.x, 3), -r.fun)
        if best is None or r.fun < best.fun:
            best = r
    import emcee
    ndim, nw = len(NAMES), 32
    rng = np.random.default_rng(42)
    p0 = best.x+1e-3*rng.standard_normal((nw, ndim))*np.maximum(np.abs(best.x), 1)
    sampler = emcee.EnsembleSampler(nw, ndim, pb.lnpost)
    sampler.run_mcmc(p0, args.steps, progress=False)
    chain = sampler.get_chain(discard=args.steps//3, thin=5, flat=True)
    lnp = sampler.get_log_prob(discard=args.steps//3, thin=5, flat=True)
    try:
        tau = sampler.get_autocorr_time(discard=args.steps//3, quiet=True).tolist()
    except Exception:
        tau = None
    q = np.percentile(chain, [16, 50, 84], axis=0)
    pbest = chain[np.argmax(lnp)]
    vl, ring, ringr = pb.model(pbest)
    chi2_pm = float(np.sum((data["pm_v"]-ring)**2/data["pm_e"]**2))
    resid = (data["v"]-vl)/np.sqrt(data["var"])
    out = dict(law=args.law, q_proj=QPROJ, nodes_arcsec=NODES.tolist(), names=NAMES,
               ml=best.x.tolist(), ml_lnpost=float(-best.fun), map=pbest.tolist(),
               p16=q[0].tolist(), p50=q[1].tolist(), p84=q[2].tolist(), autocorr=tau,
               acceptance=float(np.mean(sampler.acceptance_fraction)), n_muse=int(len(data["v"])),
               chi2_pm=chi2_pm, n_pm=int(len(ring)), muse_resid_rms=float(np.std(resid)),
               pm_r=data["pm_r"].tolist(), pm_v=data["pm_v"].tolist(), pm_e=data["pm_e"].tolist(),
               pm_model=ring.tolist(), pm_model_radial=ringr.tolist())
    (ROOT/"results/rotation").mkdir(parents=True, exist_ok=True)
    np.save(ROOT/f"results/rotation/rotation_field_{args.law}_chain_20261113.npy", chain)
    (ROOT/f"results/rotation/rotation_field_{args.law}_20261113.json").write_text(json.dumps(out, indent=1))
    for n, a, b, c in zip(NAMES, *q):
        print(f"{n:6s} {b:8.3f}  [{a:8.3f}, {c:8.3f}]")
    print(f"chi2 PM rings {chi2_pm:.1f}/{len(ring)}; MUSE normalised residual rms {np.std(resid):.3f}; "
          f"acceptance {out['acceptance']:.2f}; tau {None if tau is None else np.round(tau, 0).tolist()}")


if __name__ == "__main__":
    main()
