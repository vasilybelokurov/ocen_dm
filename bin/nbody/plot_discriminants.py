#!/usr/bin/env python3
"""Observable differences between the no-DM (A) and DM (B) disruptions at the present day.

Panels (stars only unless stated; counts normalised to each model's bound stellar count;
sky offsets relative to each model's own cluster; 'same wrap' = within 0.5 kpc of the
cluster's heliocentric distance):
 (a) cumulative unbound/bound star ratio versus angular distance from the cluster;
 (b) ratio A/B of (a);
 (c) v_los - v_los(cluster) of same-wrap tail stars within 10 deg, histograms with
     std and robust (MAD) dispersions and their 1-sigma errors sigma/sqrt(2N);
 (d) projected LOS dispersion of all stars versus projected radius (bound cluster and
     escapers), with errors;
 (e) stellar mass-loss rate versus time (per Gyr, fraction of the bound stellar mass,
     smoothed over one radial period) with pericentres;
 (f) tail width (rms offset across the tail) versus angular distance;
 (g) the dark budget of the debris: unbound stellar, remnant and DM mass.

Usage: python bin/nbody/plot_discriminants.py [--name nbody_discriminants]
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
from run_nbody import bound_mask, cluster_centre  # noqa: E402
from analyse_tails import to_sky  # noqa: E402

MODELS = (("A_nodm", "A: no DM", "tab:blue"), ("B_dm_phot", "B: DM halo", "tab:orange"))
SUN = np.array([-8178., 0., 0.])      # simulation (Baumgardt) frame, pc
T_RAD = 89.                           # radial period [Myr]


def load(m):
    ics = np.load(ROOT/f"results/nbody/{m}/ics.npz"); sp = ics["species"]; mass = ics["mass"].astype(float)
    s = np.load(ROOT/f"results/nbody/{m}/orbit/snap_today.npz"); pos, vel = s["pos"].astype(float), s["vel"].astype(float)
    lum = sp <= 1
    c, vc = cluster_centre(pos, vel, mass, lum, np.median(pos[lum], axis=0))
    bound = bound_mask(pos, vel, mass, 0.3, c, vc)
    st = sp == 0
    l, b, d, vr, *_ = to_sky(pos[st], vel[st]); lc, bc, dc, vrc = [x[0] for x in to_sky(c[None], vc[None])[:4]]
    dl = ((l-lc+180) % 360-180)*np.cos(np.radians(bc)); db = b-bc
    return dict(sp=sp, mass=mass, pos=pos, vel=vel, c=c, vc=vc, bound=bound, st=st, dl=dl, db=db, sep=np.hypot(dl, db),
                dv=vr-vrc, dd=d-dc, ub=~bound[st], nb=int(np.sum(bound & st)))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--name", default="nbody_discriminants")
    args = parser.parse_args()
    D = {m: load(m) for m, _, _ in MODELS}
    fig, axes = plt.subplots(3, 3, figsize=(17, 14), constrained_layout=True)
    ax = axes.ravel()
    radii = np.geomspace(1, 60, 40)
    cum = {}
    for m, lab, col in MODELS:
        d = D[m]
        n = np.array([np.sum(d["ub"] & (d["sep"] < R)) for R in radii])
        cum[m] = n/d["nb"]
        ax[0].plot(radii, cum[m], color=col, lw=2, label=lab)
        ax[0].fill_between(radii, (n-np.sqrt(n))/d["nb"], (n+np.sqrt(n))/d["nb"], color=col, alpha=.2)
    ax[0].set(xscale="log", yscale="log", xlabel="angular distance from the cluster [deg]", ylabel="unbound / bound stars (cumulative)",
              title="(a) tail stars relative to the cluster (band: Poisson)"); ax[0].legend()
    nA = cum["A_nodm"]*D["A_nodm"]["nb"]; nB = cum["B_dm_phot"]*D["B_dm_phot"]["nb"]
    ratio = cum["A_nodm"]/cum["B_dm_phot"]; err = ratio*np.sqrt(1/np.maximum(nA, 1)+1/np.maximum(nB, 1))
    ax[1].plot(radii, ratio, color="black", lw=2); ax[1].fill_between(radii, ratio-err, ratio+err, color="0.7")
    ax[1].axhline(1, color="0.5", ls=":"); ax[1].set(xscale="log", ylim=(0, 3), xlabel="angular distance [deg]", ylabel="A / B",
                                                   title="(b) ratio of (a): no-DM over DM (band: Poisson)")
    bins = np.linspace(-40, 40, 41)
    for m, lab, col in MODELS:
        d = D[m]; k = d["ub"] & (d["sep"] < 10) & (np.abs(d["dd"]) < 0.5); v = d["dv"][k]
        mad = 1.4826*np.median(np.abs(v-np.median(v))); N = k.sum()
        ax[2].hist(v, bins, histtype="step", color=col, lw=2, density=True,
                   label=f"{lab}: N {N}, std {np.std(v):.1f}$\\pm${np.std(v)/np.sqrt(2*N):.1f}, robust {mad:.1f}$\\pm${mad/np.sqrt(2*N):.1f} km/s")
    ax[2].set(xlabel=r"$v_{\rm los}-v_{\rm los}$(cluster) [km/s]", ylabel="density", title="(c) near-tail velocities (< 10 deg, same wrap)"); ax[2].legend(fontsize=8)
    edges = np.array([1, 3, 5, 10, 15, 20, 30, 45, 60, 80, 100, 150, 200, 300])
    for m, lab, col in MODELS:
        d = D[m]; rel = d["pos"]-d["c"]; los = (d["c"]-SUN)/np.linalg.norm(d["c"]-SUN)
        Rp = np.linalg.norm(rel-np.outer(rel@los, los), axis=1); vl = (d["vel"]-d["vc"])@los
        sel = d["st"] & (np.abs(rel@los) < 500)
        sig, er, mid = [], [], []
        for a, b in zip(edges[:-1], edges[1:]):
            k = sel & (Rp >= a) & (Rp < b)
            if k.sum() > 30:
                sig.append(np.std(vl[k])); er.append(np.std(vl[k])/np.sqrt(2*k.sum())); mid.append(np.sqrt(a*b))
        ax[3].errorbar(mid, sig, er, color=col, marker="o", ms=4, label=lab, capsize=2)
    ax[3].set(xscale="log", xlabel="projected radius [pc]", ylabel=r"$\sigma_{\rm los}$ of all stars [km/s]",
              title="(d) projected dispersion, cluster to escapers (|depth| < 0.5 kpc)"); ax[3].axvspan(60, 300, color="0.93", zorder=0); ax[3].legend()
    ins = ax[3].inset_axes([.55, .5, .42, .45])
    for m, lab, col in MODELS:
        d = D[m]; rel = d["pos"]-d["c"]; los = (d["c"]-SUN)/np.linalg.norm(d["c"]-SUN)
        Rp = np.linalg.norm(rel-np.outer(rel@los, los), axis=1); vl = (d["vel"]-d["vc"])@los; sel = d["st"] & (np.abs(rel@los) < 500)
        s2, e2, m2 = [], [], []
        for a, b in zip(edges[7:-1], edges[8:]):
            k = sel & (Rp >= a) & (Rp < b)
            if k.sum() > 30:
                s2.append(np.std(vl[k])); e2.append(np.std(vl[k])/np.sqrt(2*k.sum())); m2.append(np.sqrt(a*b))
        ins.errorbar(m2, s2, e2, color=col, marker="o", ms=3, capsize=2)
    ins.set(xscale="log", title="zoom 45-300 pc", xlabel="pc"); ins.tick_params(labelsize=7); ins.title.set_fontsize(8)
    for m, lab, col in MODELS:
        rows = [json.loads(x) for x in open(ROOT/f"results/nbody/{m}/orbit/diagnostics.jsonl")]
        t = np.array([r["t_myr"] for r in rows]); ms = np.array([r["bound_mass"]["stars"] for r in rows])
        k = t <= 1956.; t, ms = t[k], ms[k]
        w = int(T_RAD)
        rate = -(ms[w:]-ms[:-w])/ms[w:]/(T_RAD*1e-3)*100
        ax[4].plot(t[w:]-T_RAD/2, rate, color=col, lw=1.5, label=lab)
    ax[4].set(xlabel="t [Myr] (present day = 1956)", ylabel="stellar mass-loss rate [% / Gyr]", yscale="log",
              title="(e) stripping rate, averaged over one radial period"); ax[4].legend()
    seps = np.array([1, 2, 4, 6, 8, 10, 14, 20, 30, 45])
    for m, lab, col in MODELS:
        d = D[m]; k0 = d["ub"] & (np.abs(d["dd"]) < 0.5)
        # direction of the tail near the cluster: principal axis of same-wrap stars within 10 deg
        k10 = k0 & (d["sep"] < 10); X = np.column_stack((d["dl"][k10], d["db"][k10]))
        _, _, vt = np.linalg.svd(X-X.mean(0), full_matrices=False); e1 = vt[0]; e2 = np.array([-e1[1], e1[0]])
        across = d["dl"]*e2[0]+d["db"]*e2[1]
        wv, we, wm = [], [], []
        for a, b in zip(seps[:-1], seps[1:]):
            k = d["ub"] & (d["sep"] >= a) & (d["sep"] < b)
            if k.sum() > 20:
                wv.append(np.std(across[k])); we.append(np.std(across[k])/np.sqrt(2*k.sum())); wm.append(np.sqrt(a*b))
        ax[5].errorbar(wm, wv, we, color=col, marker="o", ms=4, capsize=2, label=lab)
    ax[5].set(xscale="log", xlabel="angular distance from the cluster [deg]", ylabel="rms offset across the tail [deg]", title="(f) tail width"); ax[5].legend()
    cats = ["stars", "remnants", "DM"]
    x = np.arange(3)
    for j, (m, lab, col) in enumerate(MODELS):
        d = D[m]; ub = ~d["bound"]
        vals = [d["mass"][ub & (d["sp"] == s)].sum() for s in (0, 1, 2)]
        ax[6].bar(x+(j-.5)*.4, np.maximum(vals, 1.), .4, color=col, label=lab)
        for xx, v in zip(x+(j-.5)*.4, vals):
            if v > 0: ax[6].text(xx, v*1.2, f"{v:.1e}", ha="center", fontsize=8)
    ax[6].set(yscale="log", ylim=(1e2, 1e7), xticks=x, xticklabels=cats, ylabel=r"unbound mass [$M_\odot$]", title="(g) what the debris is made of (only stars are visible)")
    ax[6].legend()
    # (h) sky map of near tails, same wrap, both models
    for m, lab, col in MODELS:
        d = D[m]; k = d["ub"] & (d["sep"] < 20)
        ax[7].scatter(d["dl"][k], d["db"][k], s=2 if m == "A_nodm" else 3, color=col, alpha=.5, label=f"{lab}: {k.sum()} particles")
    ax[7].set(xlim=(20, -20), ylim=(-20, 20), aspect="equal", xlabel=r"$\Delta l\,\cos b$ [deg]", ylabel=r"$\Delta b$ [deg]",
              title="(h) unbound stars within 20 deg (particle numbers differ: 6.5 vs 9.4 Msun each)"); ax[7].legend(fontsize=8, markerscale=4)
    # (i) summary text
    ax[8].axis("off")
    rows = [("tail/cluster, all", cum and D["A_nodm"]["ub"].sum()/D["A_nodm"]["nb"], D["B_dm_phot"]["ub"].sum()/D["B_dm_phot"]["nb"])]
    txt = ["A/B at the present day:",
           f"  tail/cluster stars, whole tail   {rows[0][1]/rows[0][2]:.2f}",
           f"  ... within 5 / 10 / 20 deg       {np.interp(5, radii, ratio):.2f} / {np.interp(10, radii, ratio):.2f} / {np.interp(20, radii, ratio):.2f}",
           "  stripping rate now (panel e)     ~1.1",
           "  near-tail sigma_los (panel c)    ~1.2",
           "  sigma_los at 100-200 pc (d)      ~1.15",
           "  tail width (f)                   ~1.0",
           "  tail dE, dLz spreads             1.0 (within 2%)",
           "",
           "B's debris is 93% dark matter (g):",
           "invisible, but it perturbs B's orbit",
           "(0.6 kpc off omega Cen's position today)."]
    ax[8].text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=10, transform=ax[8].transAxes)
    fig.suptitle("Observable differences between the no-DM (A) and DM (B) disruptions at the present day (t = 1956 Myr)")
    plot = ROOT/"plots"/f"{args.name}_today.png"
    fig.savefig(plot, dpi=130)
    print(plot)


if __name__ == "__main__":
    main()
