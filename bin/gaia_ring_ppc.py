#!/usr/bin/env python3
"""Posterior-predictive check of the outermost Gaia ring mixture fit (2026-11-19).

For each star i the fitted model density in its residual (mu_R, mu_T) plane is
  f * F(absolute PM) + (1-f) * N2(mean, diag(sigma_R^2 + e_R,i^2, sigma_T^2 + e_T,i^2)),
F = the 2D field template (per (mas/yr)^2), as in outer_profile.dispersion_2d (depth term < 1%, omitted).
Summing over stars gives the expected counts on a (mu_R, mu_T) grid; compared with the observed counts as
marginals and as a function of distance from the cluster centre in PM space.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src python bin/gaia_ring_ppc.py [--ring -1] -> plots/gaia_ring<k>_ppc.png
"""
import argparse, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/"src"))
from astropy.table import Table
from ocen_dm.cluster import OCEN_DISTANCE_KPC
from ocen_dm.paths import processed_dir
from ocen_dm.kinematics import outer_gaia as G
from ocen_dm.kinematics.outer_profile import dispersion_2d, load_members
from ocen_dm.kinematics.vb2021_replication import published_profile
from ocen_dm.selection.field_template import field_density_2d
from ocen_dm.kinematics.perspective import depth_dispersion


def setup(ring):
    prof = G.load_edr3_profile(); row = prof[ring]; lo, hi = float(row["r_lower"]), float(row["r_upper"])
    D = OCEN_DISTANCE_KPC
    cat = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv"); qf = np.asarray(cat["quality_flag"], int)
    rp, sp = published_profile(); base = load_members(exact=True, distance_kpc=D)
    err = 0.5*(base.err_r + base.err_t)
    keep = ((qf & G.QUALITY_BIT) > 0) & (err < 0.4*np.interp(base.r_arcsec, rp, sp))
    m = keep & (base.r_arcsec >= lo) & (base.r_arcsec < hi)
    gm = base.g_mag[m]
    dens = field_density_2d(err_max=float(np.max(err[m])), g_range=(float(gm.min())-0.5, float(gm.max())+0.5))
    sd = depth_dispersion(G._composite_tracer(D), base.r_arcsec[m]*D*1e3/206264.806, float(np.hypot(*base.mu_sys)), D)
    return base, m, dens, sd, lo, hi


def expected(base, m, dens, o, g):
    """Expected counts per cell of the (mu_R, mu_T) residual grid g (cell centres), summed over stars."""
    s = base.select(m); phi = np.asarray(s.phi, float); cp, sp_ = np.sin(phi), np.cos(phi)
    R, T = np.meshgrid(g, g, indexing="ij"); dA = (g[1]-g[0])**2
    tot_f, tot_c = np.zeros_like(R), np.zeros_like(R)
    for i in range(len(phi)):
        a = R*cp[i] - T*sp_[i] + s.sys_a[i]; d = R*sp_[i] + T*cp[i] + s.sys_d[i]
        tot_f += np.asarray(dens(a.ravel(), d.ravel())).reshape(R.shape)
        vr, vt = o["sigma_r"]**2 + s.err_r[i]**2, o["sigma_t"]**2 + s.err_t[i]**2
        tot_c += np.exp(-0.5*((R-o["mean_r"])**2/vr + (T-o["mean_t"])**2/vt))/(2*np.pi*np.sqrt(vr*vt))
    return o["f"]*tot_f*dA, (1-o["f"])*tot_c*dA


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ring", type=int, default=-1); a = ap.parse_args()
    base, m, dens, sd, lo, hi = setup(a.ring)
    o = dispersion_2d(base, m, dens, depth_var=sd**2, field_at=base.absolute_pm)
    w = 0.1; e = np.arange(-3, 3+1e-9, w); g = 0.5*(e[1:]+e[:-1])
    Ef, Ec = expected(base, m, dens, o, g)
    mr, mt = base.mu_r[m], base.mu_t[m]
    H, _, _ = np.histogram2d(mr, mt, [e, e])
    print(f"ring {lo:.0f}-{hi:.0f}: N {m.sum()}, in window obs {H.sum():.0f} exp {Ef.sum()+Ec.sum():.1f} "
          f"(field {Ef.sum():.1f}, cluster {Ec.sum():.1f}); f {o['f']:.4f} sigma_R {o['sigma_r']:.4f} sigma_T {o['sigma_t']:.4f}")
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
    for k, (lab, axis) in enumerate((("R", 1), ("T", 0))):
        ob, ef, ec = H.sum(axis), Ef.sum(axis), Ec.sum(axis)
        A = ax[k]; A.errorbar(g, ob, np.sqrt(np.maximum(ob, 1)), fmt="ko", ms=3, label="observed (all stars)")
        A.plot(g, ef, "C7-", label="model: field template"); A.plot(g, ec, "C0-", label="model: cluster")
        A.plot(g, ef+ec, "r-", lw=2, label="model: total")
        chi = (ob-ef-ec)/np.sqrt(np.maximum(ef+ec, 1)); A.set_title(f"mu_{lab} marginal: chi2 {np.sum(chi**2):.0f} / {len(g)}")
        A.set(xlabel=fr"$\mu_{lab}$ residual (mas/yr)", ylabel="stars / 0.1 mas/yr"); A.legend(fontsize=7)
    # radial profile in PM space around the cluster mean
    gr, gt = np.meshgrid(g, g, indexing="ij"); rr = np.hypot(gr-o["mean_r"], gt-o["mean_t"])
    be = np.arange(0, 2.01, 0.1); bc = 0.5*(be[1:]+be[:-1])
    ob = np.histogram(np.hypot(mr-o["mean_r"], mt-o["mean_t"]), be)[0]
    ef = np.histogram(rr, be, weights=Ef)[0]; ec = np.histogram(rr, be, weights=Ec)[0]
    area = np.pi*(be[1:]**2-be[:-1]**2)
    A = ax[2]; A.errorbar(bc, ob/area, np.sqrt(np.maximum(ob, 1))/area, fmt="ko", ms=3, label="observed")
    A.plot(bc, ef/area, "C7-", label="field"); A.plot(bc, ec/area, "C0-", label="cluster"); A.plot(bc, (ef+ec)/area, "r-", lw=2, label="total")
    A.set(yscale="log", xlabel=r"$|\Delta\mu|$ from cluster mean (mas/yr)", ylabel="stars per (mas/yr)$^2$",
          title="PM-space radial profile"); A.legend(fontsize=7)
    A = ax[3]; A.plot(bc, (ob-ef-ec)/np.sqrt(np.maximum(ef+ec, 1)), "ko-"); A.axhline(0, c="k", lw=.5)
    A.axhspan(-1, 1, color="0.9"); A.set(xlabel=r"$|\Delta\mu|$ (mas/yr)", ylabel="(obs - model)/sqrt(model)", title="residuals")
    fig.suptitle(f"Gaia ring {lo:.0f}-{hi:.0f}\": posterior-predictive check of the field + cluster mixture", fontsize=10)
    fig.tight_layout(); f = ROOT/f"plots/gaia_ring{a.ring % len(G.load_edr3_profile())}_ppc.png"; fig.savefig(f, dpi=100); print(f)


if __name__ == "__main__":
    main()
