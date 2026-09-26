#!/usr/bin/env python3
"""The N-body tails of models A (no DM) and B (DM) against the observed omega Cen tails.

Observed quantities (values quoted in the papers; nothing digitised from figures):
  K21 = Kuzma, Ferguson & Penarrubia 2021, MNRAS 507, 1127 (arXiv:2108.02531), Gaia EDR3:
        slopes of the star-count profile beyond the King r_t = 46.4', on-axis -3.40 +- 0.20,
        off-axis -5.22 +- 0.26 (45-deg wedges, out to 5 deg); tail position angle
        122.88 +- 2.14 deg E of N; >= 0.1 per cent of the stellar mass between the Wilson
        radius 70.6' and 5 deg (lower limit, mass-follows-light LIMEPY).
  DC12 = Da Costa 2012: LOS dispersion 6.7 km/s at the tidal radius (as used by K26).
  K26 = Kuzma et al. 2026 (arXiv:2605.23474): 5 members beyond r_J = 106.9', out to 3.2 deg;
        outermost P_mem >= 0.3 star at x'' = 3 deg offset by ~7 km/s (sign not stated);
        "no evidence of strong gradients".
  I19 = Ibata et al. 2019, Nat. Astron. (arXiv:1902.09544): Fimbulthul at b ~ 35 deg along
        l = -52..-32 deg, ~1.5 kpc closer than omega Cen (STREAMFINDER distance 4.1 kpc);
        5 CFHT stars, <v_helio> = 199.7 km/s, rms 5.4 km/s (their Table 2).
Model quantities are computed for the stars (all, bound + unbound) at the present-day
snapshot, relative to each model's own centre for the cluster-centric measures and in
absolute sky coordinates for Fimbulthul (model B's cluster is 0.6 kpc / 6 deg in b off
omega Cen's position, see JOURNAL 2026-09-25).

Usage: python bin/nbody/plot_observed_comparison.py [--name nbody_vs_observed]
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(key, "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"bin"/"nbody"))
from run_nbody import cluster_centre, bound_mask  # noqa: E402
from analyse_tails import to_sky  # noqa: E402

MODELS = (("A_nodm", "A: no DM", "tab:blue"), ("B_dm_phot", "B: DM halo", "tab:orange"))
R_T, R_W, R_J = 46.4/60, 70.6/60, 106.9/60           # deg
K21 = dict(on=(-3.40, 0.20), off=(-5.22, 0.26), pa=(122.88, 2.14), frac_min=0.001)
DC12_SIG = 6.7
I19_V = np.array([205.52, 191.87, 202.72, 203.74, 194.74]); I19_POS = np.array([[200.291401, -26.516432], [202.470324, -24.439468],
                   [202.470611, -25.084004], [205.597014, -25.009721], [208.831549, -22.971256]])
FIMB_L, FIMB_B = (-52., -32.), (30., 40.)


def load(m):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    ics = np.load(ROOT/f"results/nbody/{m}/ics.npz"); sp = ics["species"]; mass = ics["mass"].astype(float)
    s = np.load(ROOT/f"results/nbody/{m}/orbit/snap_today.npz"); pos, vel = s["pos"].astype(float), s["vel"].astype(float)
    lum = sp <= 1
    c, vc = cluster_centre(pos, vel, mass, lum, np.median(pos[lum], axis=0))
    bound = bound_mask(pos, vel, mass, 0.3, c, vc)
    st = sp == 0
    l, b, d, vr, *_ = to_sky(pos[st], vel[st]); lc, bc, dc, vrc = [x[0] for x in to_sky(c[None], vc[None])[:4]]
    eq = SkyCoord(l=l*u.deg, b=b*u.deg, frame="galactic").icrs; ceq = SkyCoord(l=lc*u.deg, b=bc*u.deg, frame="galactic").icrs
    ra, dec = eq.ra.rad, eq.dec.rad; a0, d0 = ceq.ra.rad, ceq.dec.rad
    cosc = np.sin(d0)*np.sin(dec)+np.cos(d0)*np.cos(dec)*np.cos(ra-a0)
    xi = np.degrees(np.cos(dec)*np.sin(ra-a0)/cosc); eta = np.degrees((np.cos(d0)*np.sin(dec)-np.sin(d0)*np.cos(dec)*np.cos(ra-a0))/cosc)
    return dict(l=np.where(l > 180, l-360, l), b=b, d=d, vr=vr, xi=xi, eta=eta, R=np.hypot(xi, eta), pa=np.degrees(np.arctan2(xi, eta)) % 360,
                dv=vr-vrc, dc=dc, vrc=vrc, lc=lc if lc < 180 else lc-360, bc=bc, ub=~bound[st], mstar=mass[st][0],
                ra=np.degrees(ra), dec=np.degrees(dec))


def wedge_profile(d, axis_pa, edges, on=True):
    dpa = np.abs(((d["pa"]-axis_pa+90) % 180)-90)          # angle from the tail axis, folded to 0-90
    sel = dpa < 22.5 if on else dpa > 67.5
    n, _ = np.histogram(d["R"][sel], edges)
    area = np.pi*(edges[1:]**2-edges[:-1]**2)*(45*2/360)    # two 45-deg wedges
    return n/area, np.sqrt(np.maximum(n, 1))/area, n


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--name", default="nbody_vs_observed")
    args = parser.parse_args()
    D = {m: load(m) for m, _, _ in MODELS}
    fig, axes = plt.subplots(2, 3, figsize=(18, 11), constrained_layout=True)
    ax = axes.ravel()
    edges = np.geomspace(0.2, 5., 22); mid = np.sqrt(edges[1:]*edges[:-1])
    fit = (mid > R_T) & (mid < 5.)
    rec = dict(models={})
    for m, lab, col in MODELS:
        d = D[m]
        # tail PA of the model: mean axis of unbound stars at 1.5-5 deg (axial average)
        k = d["ub"] & (d["R"] > 1.5) & (d["R"] < 5.)
        ang = np.radians(2*d["pa"][k]); pa_model = (np.degrees(np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()))/2) % 180
        out = dict(pa_tail=float(pa_model))
        for on, ls in ((True, "-"), (False, "--")):
            s, e, n = wedge_profile(d, pa_model, edges, on)
            good = fit & (n >= 5)
            p = np.polyfit(np.log10(mid[good]), np.log10(s[good]), 1, w=1/np.maximum(e[good]/s[good]/np.log(10), 1e-3))
            ax[0].errorbar(mid, s, e, color=col, ls=ls, marker="o" if on else "s", ms=3, capsize=0,
                           label=f"{lab}, {'on' if on else 'off'}-axis: slope {p[0]:.2f} (r$_t$-5$^\\circ$)")
            out[f"slope_{'on' if on else 'off'}"] = float(p[0])
        # mass fraction between the Wilson radius and 5 deg (all stars, bound + unbound)
        out["frac_wilson_5deg"] = float(np.sum((d["R"] > R_W) & (d["R"] < 5.))/len(d["R"]))
        rec["models"][m] = out
    # K21 slopes, anchored at the model A on-axis density at r_t
    sA, _, _ = wedge_profile(D["A_nodm"], rec["models"]["A_nodm"]["pa_tail"], edges, True)
    s0 = np.interp(R_T, mid, sA); rr = np.geomspace(R_T, 5, 20)
    for key, ls in (("on", "-"), ("off", "--")):
        g, ge = K21[key]
        ax[0].plot(rr, s0*(rr/R_T)**g, color="black", ls=ls, lw=2.5, alpha=.6, label=f"K21 {key}-axis slope {g} $\\pm$ {ge} (normalised at r$_t$)")
    for R, t in ((R_T, "r$_t$"), (R_W, "r$_W$"), (R_J, "r$_J$")):
        ax[0].axvline(R, color="0.7", lw=.8); ax[0].text(R, 1e5, t, fontsize=8, ha="center")
    ax[0].set(xscale="log", yscale="log", xlabel="projected radius [deg]", ylabel="star particles per deg$^2$",
              title="(a) star-count profile along (solid) and across (dashed) the tails vs Kuzma+2021")
    ax[0].legend(fontsize=7)
    # (b) summary numbers vs K21
    names = ["slope on", "slope off"]
    for j, (m, lab, col) in enumerate(MODELS):
        o = rec["models"][m]
        ax[1].errorbar([0+(j-.5)*.2, 1+(j-.5)*.2], [o["slope_on"], o["slope_off"]], fmt="o", color=col, ms=8, label=lab)
    ax[1].errorbar([0.3, 1.3], [K21["on"][0], K21["off"][0]], [K21["on"][1], K21["off"][1]], fmt="*", color="black", ms=12, capsize=3, label="Kuzma+2021")
    ax[1].set(xticks=[0, 1], xticklabels=names, xlim=(-.5, 1.7), ylabel="power-law slope beyond r$_t$", title="(b) profile slopes, models vs observed")
    ax[1].legend(fontsize=8)
    bx = ax[1].inset_axes([.08, .08, .25, .38])
    for j, (m, lab, col) in enumerate(MODELS):
        bx.bar(j, 100*rec["models"][m]["frac_wilson_5deg"], color=col)
    bx.axhline(100*K21["frac_min"], color="black", ls="--"); bx.text(1.5, 100*K21["frac_min"], " K21 lower limit", va="bottom", fontsize=7)
    bx.set(xticks=[0, 1], xticklabels=["A", "B"], ylabel="% of stars", title="stars between r$_W$ and 5$^\\circ$", xlim=(-.5, 2.5)); bx.title.set_fontsize(8); bx.tick_params(labelsize=7)
    # (c) sky view in the tangent plane with PAs
    for m, lab, col in MODELS:
        d = D[m]; k = d["ub"] & (d["R"] < 6)
        ax[2].scatter(d["xi"][k], d["eta"][k], s=2, color=col, alpha=.5, label=f"{lab}: tail PA {rec['models'][m]['pa_tail']:.1f}$^\\circ$")
    t = np.radians(K21["pa"][0]); ax[2].plot([-6*np.sin(t), 6*np.sin(t)], [-6*np.cos(t), 6*np.cos(t)], color="black", lw=1.5, label=f"K21 PA {K21['pa'][0]} $\\pm$ {K21['pa'][1]}$^\\circ$")
    for R in (R_T, R_J, 3.2):
        ax[2].add_patch(plt.Circle((0, 0), R, fill=False, color="0.6", ls=":"))
    ax[2].set(xlim=(6, -6), ylim=(-6, 6), aspect="equal", xlabel=r"$\xi$ [deg] (East to the left)", ylabel=r"$\eta$ [deg]",
              title="(c) unbound stars within 6$^\\circ$ (circles: r$_t$, r$_J$, 3.2$^\\circ$ = K26 outermost field)")
    ax[2].legend(fontsize=7, markerscale=4)
    # (d) v_los along the tails vs DC12 / K26
    for m, lab, col in MODELS:
        d = D[m]; pa = np.radians(rec["models"][m]["pa_tail"])
        x2 = d["xi"]*np.sin(pa)+d["eta"]*np.cos(pa); y2 = -d["xi"]*np.cos(pa)+d["eta"]*np.sin(pa)
        k = (np.abs(y2) < 0.5) & (np.abs(x2) < 4) & (d["R"] > R_T) & (np.abs(d["d"]-d["dc"]) < 0.5)
        ax[3].scatter(x2[k], d["dv"][k], s=3, color=col, alpha=.5, label=f"{lab}: {k.sum()} particles beyond r$_t$")
        bins = np.linspace(-4, 4, 17); ib = np.digitize(x2[k], bins)-1
        sg = [np.std(d["dv"][k][ib == j]) if np.sum(ib == j) > 10 else np.nan for j in range(16)]
        ax[3].plot(0.5*(bins[1:]+bins[:-1]), sg, color=col, lw=2, ls="--")
        rec["models"][m]["sigma_los_beyond_rt_within_4deg"] = float(np.std(d["dv"][k]))
    ax[3].axhspan(-DC12_SIG, DC12_SIG, color="0.85", zorder=0, label=f"DC12 $\\sigma$ = {DC12_SIG} km/s at r$_t$ ($\\pm$1$\\sigma$ band)")
    ax[3].plot([3, -3], [7, 7], "k^", ms=9, mfc="none", label="K26 outermost member: |$\\Delta v$| ~ 7 km/s at |x''| = 3$^\\circ$ (sign not given)")
    ax[3].plot([3, -3], [-7, -7], "kv", ms=9, mfc="none")
    ax[3].set(xlim=(-4, 4), ylim=(-40, 40), xlabel="x'' along the tail [deg]", ylabel=r"$v_{\rm los}-v_{\rm sys}$ [km/s]",
              title="(d) velocities beyond r$_t$, |y''| < 0.5$^\\circ$ (dashed: model $\\sigma$ per bin)")
    ax[3].legend(fontsize=7, loc="lower left")
    # (e,f) Fimbulthul
    for m, lab, col in MODELS:
        d = D[m]
        k = d["ub"] & (d["l"] > FIMB_L[0]) & (d["l"] < FIMB_L[1]) & (d["b"] > FIMB_B[0]) & (d["b"] < FIMB_B[1])
        near = k & (np.abs(d["d"]-(d["dc"]-1.5)) < 1.)                    # the stream is ~1.5 kpc closer
        ax[4].hist(d["vr"][k], np.linspace(100, 320, 56), histtype="step", color=col, lw=1.5, label=f"{lab}: all debris in the box ({k.sum()})")
        ax[4].hist(d["vr"][near], np.linspace(100, 320, 56), histtype="stepfilled", color=col, alpha=.35, label=f"{lab}: at D$_\\odot$ - 1.5 $\\pm$ 1 kpc ({near.sum()})")
        rec["models"][m]["fimbulthul"] = dict(n_box=int(k.sum()), frac_box=float(k.sum()/np.sum(~d["ub"])), n_near=int(near.sum()),
                                              v_median=float(np.median(d["vr"][near])) if near.sum() else None,
                                              v_std=float(np.std(d["vr"][near])) if near.sum() > 2 else None,
                                              d_median=float(np.median(d["d"][k])) if k.sum() else None)
        ax[5].scatter(d["l"][d["ub"]], d["b"][d["ub"]], s=1, color=col, alpha=.3, rasterized=True)
        ax[5].plot(d["lc"], d["bc"], "+", color=col, ms=12, mew=2, label=f"{lab}: model cluster (l, b) = ({d['lc']:.1f}, {d['bc']:.1f})")
    for v in I19_V:
        ax[4].axvline(v, color="black", lw=1)
    ax[4].axvspan(199.7-5.4, 199.7+5.4, color="0.8", zorder=0, label="I19 Fimbulthul: 199.7, rms 5.4 km/s (5 stars, lines)")
    ax[4].set(xlabel=r"$v_{\rm helio}$ [km/s]", ylabel="star particles", title="(e) Fimbulthul box: l -52..-32, b 30..40 deg")
    ax[4].legend(fontsize=7)
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    g = SkyCoord(ra=I19_POS[:, 0]*u.deg, dec=I19_POS[:, 1]*u.deg).galactic
    ax[5].plot(np.where(g.l.deg > 180, g.l.deg-360, g.l.deg), g.b.deg, "k*", ms=12, label="I19 CFHT Fimbulthul stars")
    ax[5].add_patch(plt.Rectangle((FIMB_L[0], FIMB_B[0]), FIMB_L[1]-FIMB_L[0], FIMB_B[1]-FIMB_B[0], fill=False, color="black", ls="--"))
    ax[5].plot(-50.898, 14.968, "kx", ms=10, mew=2, label="omega Cen observed")
    ax[5].set(xlim=(-10, -80), ylim=(-10, 55), xlabel="l [deg]", ylabel="b [deg]", title="(f) unbound stars on the sky around Fimbulthul (dashed box)")
    ax[5].legend(fontsize=7, markerscale=1)
    fig.suptitle("N-body tails of the no-DM (A) and DM (B) models vs observed omega Cen tails (K21 = Kuzma+2021, DC12 = Da Costa 2012, K26 = Kuzma+2026, I19 = Ibata+2019)")
    plot = ROOT/"plots"/f"{args.name}_today.png"
    fig.savefig(plot, dpi=130)
    (ROOT/"results/plot_data"/f"{args.name}_today.json").write_text(json.dumps(rec, indent=1, default=float))
    print(json.dumps(rec, indent=1, default=float)); print(plot)


if __name__ == "__main__":
    main()
