#!/usr/bin/env python3
"""oMEGACat field stars selected as EVERYTHING outside the cluster's PM clump, vs radius.

Parent: stars with PMs and PM errors < 0.5 mas/yr (as bin/plot_field_counts.py). "Outside": |mu_rel| > --rcl
mas/yr from the cluster AND chi2 > 25 against the cluster (intrinsic sigma(R) interpolated from the hq values:
0.81, 0.78, 0.73, 0.66, 0.58 mas/yr at 15, 45, 90, 160, 270'').
Leakage of cluster stars (PM-error tails) is estimated from the half of the PM plane facing AWAY from the field
clump (position angle of mu more than 90 deg from the field-centroid direction, (-2.58, +5.37) mas/yr): the
tails are taken as isotropic about the cluster, so leakage over all directions = 2 x the far-half count; the field
clump (6.0 mas/yr from the cluster, sd ~1.3) contributes negligibly to the far half. Net = all - 2 x far half.
Shown for --rcl = 4 (raw, leakage, net) and the nets for 3 and 5 mas/yr as a robustness check. Areas from 3''
cell occupancy; Poisson errors. A second panel repeats the inner 10'' as counts per 2'' annulus.

Usage: python bin/plot_field_counts_outside.py
Output: plots/omegacat_field_counts_outside.png, results/plot_data/omegacat_field_counts_outside.json
"""
from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
from astropy.io import fits

SIG_R = np.array([15., 45., 90., 160., 270.])
SIG_V = np.array([0.81, 0.78, 0.73, 0.66, 0.58])
FC = np.array([-2.58, 5.37])
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
RAMP = LinearSegmentedColormap.from_list("blue", ["#e8f0fb", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])


def main():
    t = fits.open(ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits", memmap=True)[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.); y = 0.04*(np.asarray(t["y"], float)-15000.)
    R = np.hypot(x, y)
    ma, md = np.asarray(t["pmra_corrected"], float), np.asarray(t["pmdec_corrected"], float)
    ea, ed = np.asarray(t["pmra_corrected_err"], float), np.asarray(t["pmdec_corrected_err"], float)
    pos = np.isfinite(x) & np.isfinite(y)
    par = pos & np.isfinite(ma) & np.isfinite(md) & (ea < 0.5) & (ed < 0.5)
    s = np.interp(R, SIG_R, SIG_V)
    chi2 = ma**2/(s**2+ea**2) + md**2/(s**2+ed**2)
    mu = np.hypot(ma, md)
    cosang = (ma*FC[0]+md*FC[1])/np.maximum(mu*np.hypot(*FC), 1e-12)
    far = cosang < 0                                       # more than 90 deg from the field direction
    c = 3.
    ix = np.floor(x[pos]/c).astype(np.int64)+10**5; iy = np.floor(y[pos]/c).astype(np.int64)+10**5
    cells = np.unique(ix*10**6+iy)
    Rc = np.hypot((cells//10**6-10**5+.5)*c, (cells % 10**6-10**5+.5)*c)

    def profile(edges, sel, area=None):
        n = np.histogram(R[sel], edges)[0].astype(float)
        return n
    edges = np.arange(0, 480, 30.); mid = .5*(edges[1:]+edges[:-1])
    area = np.array([np.sum((Rc >= a) & (Rc < b))*c**2/3600. for a, b in zip(edges[:-1], edges[1:])])
    good = area > 0.05
    out = dict(edges_arcsec=edges.tolist(), area_arcmin2=area.tolist(), cuts={})
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.6), gridspec_kw=dict(width_ratios=[1.5, 1, 1]))
    a = ax[0]
    for rcl, col_net in ((3., "#1baf7a"), (4., "#0d366b"), (5., "#eda100")):
        sel = par & (chi2 > 25) & (mu > rcl)
        n_all = profile(edges, sel); n_far = profile(edges, sel & far)
        net = n_all-2*n_far; err = np.sqrt(n_all+4*n_far)
        dn = lambda v: np.where(good, v/np.maximum(area, 1e-9), np.nan)
        out["cuts"][str(rcl)] = dict(n_all=n_all.tolist(), n_far_half=n_far.tolist(), net_density=dn(net).tolist(),
                                     net_err=dn(err).tolist(), raw_density=dn(n_all).tolist(), leak_density=dn(2*n_far).tolist())
        if rcl == 4.:
            a.errorbar(mid, dn(n_all), dn(np.sqrt(n_all)), fmt="o-", color="#2a78d6", ms=4, capsize=0,
                       label=f"all outside |mu| > 4 (N = {int(n_all.sum()):,})")
            a.errorbar(mid, dn(2*n_far), dn(2*np.sqrt(n_far)), fmt="s--", color=MUTED, ms=4, capsize=0,
                       label=f"leakage = 2 x far half (N = {int(2*n_far.sum()):,})")
        a.errorbar(mid, dn(net), dn(err), fmt="D-", color=col_net, ms=5 if rcl == 4. else 4, lw=2.2 if rcl == 4. else 1.3,
                   capsize=0, label=f"net, |mu| > {rcl:g} mas/yr")
    a.set_yscale("log"); a.set_xlabel("R [arcsec]", color=INK); a.set_ylabel("field-star surface density [arcmin$^{-2}$]", color=INK)
    a.set_title("everything outside the cluster PM clump (chi2 > 25), PM err < 0.5 mas/yr", fontsize=10, color=INK, loc="left")
    a.grid(True, color=GRID, lw=0.6); a.legend(fontsize=8, frameon=False)
    # inner 10''
    a = ax[1]
    e2 = np.arange(0, 10.01, 2.); m2 = .5*(e2[1:]+e2[:-1]); ar2 = np.pi*(e2[1:]**2-e2[:-1]**2)/3600.
    sel = par & (chi2 > 25) & (mu > 4)
    n_all = np.histogram(R[sel], e2)[0]; n_far = np.histogram(R[sel & far], e2)[0]
    net_out = np.array(out["cuts"]["4.0"]["net_density"])
    dens_ref = float(np.nanmedian(net_out[(mid > 30) & (mid < 330)]))
    w = 0.5
    a.bar(m2-w/2, n_all, w, color="#2a78d6", label=f"all outside (N = {n_all.sum()})")
    a.bar(m2+w/2, 2*n_far, w, color=MUTED, label=f"2 x far half (N = {2*n_far.sum()})")
    a.step(e2, np.r_[dens_ref*ar2, dens_ref*ar2[-1]], where="post", color="#eb6834", lw=1.8,
           label=f"uniform field at the net {dens_ref:.0f} arcmin$^{{-2}}$ (30-330'')")
    a.set_xlabel("R [arcsec]", color=INK); a.set_ylabel("stars per 2'' annulus", color=INK)
    a.set_title("inner 10'', |mu| > 4 mas/yr", fontsize=10, color=INK, loc="left"); a.legend(fontsize=8, frameon=False)
    a.grid(True, color=GRID, lw=0.6, axis="y")
    out["inner10"] = dict(edges=e2.tolist(), n_all=n_all.tolist(), n_far=n_far.tolist(), reference_density=dens_ref)
    # VPD with the regions
    a = ax[2]
    H, xe, ye = np.histogram2d(ma[par], md[par], bins=(np.linspace(-10, 10, 201), np.linspace(-10, 10, 201)))
    a.pcolormesh(xe, ye, H.T, cmap=RAMP, norm=LogNorm(vmin=1, vmax=H.max()), rasterized=True)
    th = np.linspace(0, 2*np.pi, 300)
    for rcl, ls in ((3., ":"), (4., "-"), (5., "--")):
        a.plot(rcl*np.cos(th), rcl*np.sin(th), ls, color="#eb6834", lw=1.2, label=f"|mu| = {rcl:g}")
    u = -FC/np.hypot(*FC); perp = np.array([-u[1], u[0]])
    a.plot([perp[0]*-10, perp[0]*10], [perp[1]*-10, perp[1]*10], "-", color=MUTED, lw=1)
    a.text(*(u*7), "far half\n(leakage)", color=MUTED, ha="center", fontsize=9)
    a.plot(*FC, "x", color="#2a78d6", ms=10, mew=2, label="field centroid")
    a.set_xlim(-10, 10); a.set_ylim(-10, 10); a.set_aspect("equal")
    a.set_xlabel(r"$\mu_{\alpha*}$ rel. to cluster [mas/yr]", color=INK); a.set_ylabel(r"$\mu_\delta$ rel. [mas/yr]", color=INK)
    a.set_title("selection in the vector-point diagram", fontsize=10, color=INK, loc="left"); a.legend(fontsize=8, frameon=False, loc="upper right")
    for a in ax:
        a.tick_params(colors=MUTED, labelsize=9)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.suptitle("oMEGACat: field stars as everything outside the cluster PM clump; leakage from the half-plane away from the field",
                 fontsize=11.5, color=INK)
    fig.tight_layout()
    fig.savefig(ROOT/"plots/omegacat_field_counts_outside.png", dpi=130)
    (ROOT/"results/plot_data/omegacat_field_counts_outside.json").write_text(json.dumps(out, indent=1))
    for k, v in out["cuts"].items():
        print("|mu| >", k, "net:", np.round(v["net_density"], 1).tolist())
        print("        raw:", np.round(v["raw_density"], 1).tolist(), " leak:", np.round(v["leak_density"], 1).tolist())
    print("inner10", out["inner10"])


if __name__ == "__main__":
    main()
