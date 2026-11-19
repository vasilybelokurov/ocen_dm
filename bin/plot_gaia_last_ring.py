#!/usr/bin/env python3
"""How the outermost Gaia EDR3 dispersion point is measured (2026-11-19).

Re-runs the production per-ring measurement (outer_gaia.build_edr3_profile, raw errors) for one ring only and shows:
(a) radii of the selected stars: all (cluster + field) vs membership-weighted, with the 8 model-averaging nodes;
(b) the residual proper motions (radial vs tangential) coloured by cluster membership probability;
(c, d) mu_R and mu_T histograms: all stars, split into membership-weighted cluster and field parts, with the fitted
cluster Gaussian (fitted sigma + median per-star error, scaled by the fitted cluster count);
(e) the fitted sigma_R, sigma_T of all Gaia rings vs the JAM C6 (inner node) K1 and NFW models.
Checks that the re-measured values equal the stored product.
Usage: OMP_NUM_THREADS=1 PYTHONPATH=src python bin/plot_gaia_last_ring.py [--ring -1] -> plots/gaia_ring<k>_measurement.png
"""
import argparse, json, sys
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


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ring", type=int, default=-1); a = ap.parse_args()
    prof = G.load_edr3_profile()
    row = prof[a.ring]; lo, hi = float(row["r_lower"]), float(row["r_upper"])
    D = OCEN_DISTANCE_KPC
    cat = Table.read(processed_dir()/"tails"/"vasiliev2021_ocen_members.ecsv")
    qf = np.asarray(cat["quality_flag"], int)
    rp, sp = published_profile()
    base = load_members(exact=True, distance_kpc=D)
    err = 0.5*(base.err_r + base.err_t)
    keep = ((qf & G.QUALITY_BIT) > 0) & (err < 0.4*np.interp(base.r_arcsec, rp, sp))
    m = keep & (base.r_arcsec >= lo) & (base.r_arcsec < hi)
    gm = base.g_mag[m]
    dens2d = field_density_2d(err_max=float(np.max(err[m])), g_range=(float(gm.min())-0.5, float(gm.max())+0.5))
    r_pc = base.r_arcsec[m]*D*1e3/206264.806
    sd = depth_dispersion(G._composite_tracer(D), r_pc, float(np.hypot(*base.mu_sys)), D)
    o = dispersion_2d(base, m, dens2d, depth_var=sd**2, field_at=base.absolute_pm)
    p = np.asarray(o["p_member"]); r = base.r_arcsec[m]
    print(f"ring {lo:.0f}-{hi:.0f}\": N {m.sum()}, f_field {o['f']:.3f}, cluster N {np.sum(p):.1f}; "
          f"sigma_R {o['sigma_r']:.4f} +- {o['sigma_r_err']:.4f} (product {row['sigma_pmr']:.4f}), "
          f"sigma_T {o['sigma_t']:.4f} +- {o['sigma_t_err']:.4f} (product {row['sigma_pmt']:.4f}); "
          f"mean_R {o['mean_r']:.4f} mean_T {o['mean_t']:.4f}")
    assert abs(o["sigma_r"] - row["sigma_pmr"]) < 2e-3 and abs(o["sigma_t"] - row["sigma_pmt"]) < 2e-3

    murr, mutt = base.mu_r[m], base.mu_t[m]          # residual PMs in the radial / tangential frame
    er, et = base.err_r[m], base.err_t[m]
    fig, ax = plt.subplots(1, 5, figsize=(24, 4.8))
    b = np.linspace(lo, hi, 25)
    ax[0].hist(r, b, color="0.75", label=f"all selected ({m.sum()})")
    ax[0].hist(r, b, weights=p, color="C0", alpha=.8, label=f"membership-weighted ({p.sum():.0f})")
    nodes_w = G._radial_nodes(r, p); nodes_u = G._radial_nodes(r)
    for k, (nn, c, l) in enumerate(((nodes_w, "C0", "model nodes (member-weighted, used)"),
                                    (nodes_u, "C3", "nodes before 2026-11-16 fix"))):
        for j, x in enumerate(nn):
            ax[0].axvline(x, color=c, ls="-" if k == 0 else ":", lw=1, label=l if j == 0 else None)
    ax[0].set(xlabel="R (arcsec)", ylabel="stars", yscale="log", title=f"ring {lo:.0f}-{hi:.0f}\" ({lo*D*1e3/206265:.0f}-{hi*D*1e3/206265:.0f} pc)")
    ax[0].legend(fontsize=7)
    o_ = np.argsort(p)
    sc = ax[1].scatter(murr[o_], mutt[o_], c=p[o_], s=4, cmap="viridis", vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax[1], label="P(member)")
    t = np.linspace(0, 2*np.pi, 200)
    for k in (1, 2):
        ax[1].plot(o["mean_r"] + k*np.hypot(o["sigma_r"], np.median(er))*np.cos(t),
                   o["mean_t"] + k*np.hypot(o["sigma_t"], np.median(et))*np.sin(t), "r-", lw=1)
    ax[1].set(xlim=(-3, 3), ylim=(-3, 3), xlabel=r"$\mu_R$ residual (mas/yr)", ylabel=r"$\mu_T$ residual (mas/yr)",
              title="cluster 1,2 sigma (fit + median error)")
    for i, (v, s, mu, e, lab) in enumerate(((murr, o["sigma_r"], o["mean_r"], er, "R"),
                                           (mutt, o["sigma_t"], o["mean_t"], et, "T"))):
        A = ax[2+i]; bb = np.linspace(-2.5, 2.5, 51); w = bb[1]-bb[0]
        A.hist(v, bb, color="0.8", label="all")
        A.hist(v, bb, weights=1-p, histtype="step", color="k", lw=1.2, label="field-weighted")
        A.hist(v, bb, weights=p, histtype="step", color="C0", lw=2, label="member-weighted")
        x = np.linspace(-2.5, 2.5, 400)
        # cluster model = sum over members of N(mu, sigma^2 + e_i^2) (depth term omitted, < 1%)
        g = sum(pi*np.exp(-0.5*(x-mu)**2/(s**2+ei**2))/np.sqrt(2*np.pi*(s**2+ei**2)) for pi, ei in zip(p, e))
        A.plot(x, g*w, "r-", lw=1.5, label=f"fit: sigma_{lab} = {s:.3f} +- {o[f'sigma_{lab.lower()}_err']:.3f}")
        A.set(xlabel=fr"$\mu_{lab}$ residual (mas/yr)", ylabel="stars / bin", title=f"{lab} component (median err {np.median(e):.3f})")
        A.legend(fontsize=7)
    # all rings vs JAM C6 models
    for comp, c in (("pmr", "C0"), ("pmt", "C1")):
        ax[4].errorbar(prof["r_median"], prof[f"sigma_{comp}"], prof[f"sigma_{comp}_err"], fmt="o", color=c, ms=4,
                       label=f"Gaia sigma_{comp[-1].upper()} (data)")
    for fam, ls in (("K1", "-"), ("NFW", "--")):
        J = json.loads((ROOT/f"results/rotation/jamfit_real_{fam}_lsq_qhyb600_remstars_bspline6_in2_nw.json").read_text())
        for nm, c in (("gaia_edr3_ours_radial", "C0"), ("gaia_edr3_ours_tangential", "C1")):
            d = J["data_used"][nm]
            ax[4].plot(d["r"] if "r" in d else prof["r_median"], J["pred"][nm], ls, color=c, lw=1, label=f"JAM C6 {fam}" if c == "C0" else None)
    ax[4].axvspan(lo, hi, color="0.9", zorder=0)
    ax[4].set(xscale="log", xlabel="R (arcsec)", ylabel="dispersion (mas/yr)", title="all Gaia rings (shaded = this ring)")
    ax[4].legend(fontsize=7)
    fig.tight_layout(); k = a.ring % len(prof)
    f = ROOT/f"plots/gaia_ring{k}_measurement.png"; fig.savefig(f, dpi=100); print(f)


if __name__ == "__main__":
    main()
