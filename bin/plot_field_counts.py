#!/usr/bin/env python3
"""Surface density of proper-motion NON-members (field stars) in the oMEGACat HST catalogue vs radius.

Data: data/raw/omegacat_vi_kinematics/catalog_and_selections.fits (oMEGACat, Zenodo 10.5281/zenodo.14978551).
PMs are relative to the cluster, so members cluster at the origin and field stars sit ~7.5 mas/yr away.
Parent: stars with PMs and PM errors < --emax mas/yr in both components (no hq flag, so the sample reaches
the catalogue edge, ~466 arcsec). Field star: within --rwin mas/yr of the field centroid in relative PM
(measured: (-2.58, +5.37) mas/yr, 6.0 mas/yr from the cluster; iterated median of stars > 4 mas/yr from the
cluster), more than --rcl mas/yr from the cluster, AND chi2_cl = mu_a^2/(s(R)^2+e_a^2) + mu_d^2/(s(R)^2+e_d^2) > --chi2,
where s(R) is the cluster's intrinsic PM dispersion, interpolated from the hq-astrometry values measured in
bin/plot_cmd_radial.py (0.81, 0.78, 0.73, 0.66, 0.58 mas/yr at 15, 45, 90, 160, 270 arcsec; held constant
outside). Radius about pixel (15000, 15000), 0.04 arcsec/px. Area per annulus from the occupancy of 3-arcsec
cells by any catalogue star (as bin/build_count_profiles.py). Poisson errors.
CONTROL: the same window mirrored through the cluster (+2.58, -5.37) contains no field population, so its
counts measure leakage of cluster stars (PM-error tails) into a window at that distance, per radius.
NB: a chi2 > 25 cut alone does not isolate the field: it is dominated by cluster stars in the non-Gaussian
tails of the PM errors and is strongly centrally concentrated (tested 2026-10-08).
A uniform field population should give a flat profile; departures measure radius-dependent incompleteness
(crowding) of the selection, or residual member contamination (see the control).

Usage: python bin/plot_field_counts.py [--emax 0.5] [--chi2 25]
Output: plots/omegacat_field_counts.png, results/plot_data/omegacat_field_counts.json
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
import numpy as np
from astropy.io import fits

SIG_R = np.array([15., 45., 90., 160., 270.])
SIG_V = np.array([0.81, 0.78, 0.73, 0.66, 0.58])
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
COLORS = ["#2a78d6", "#6b6a64", "#eb6834", "#1baf7a", "#eda100"]
RAMP = LinearSegmentedColormap.from_list("blue", ["#e8f0fb", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--emax", type=float, default=0.5)
    p.add_argument("--chi2", type=float, default=25.)
    p.add_argument("--cell", type=float, default=3.)
    p.add_argument("--rwin", type=float, default=3.0, help="radius of the field PM window [mas/yr]")
    p.add_argument("--rcl", type=float, default=4.0, help="minimum distance from the cluster PM [mas/yr]")
    args = p.parse_args()
    t = fits.open(ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits", memmap=True)[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.); y = 0.04*(np.asarray(t["y"], float)-15000.)
    R = np.hypot(x, y)
    ma, md = np.asarray(t["pmra_corrected"], float), np.asarray(t["pmdec_corrected"], float)
    ea, ed = np.asarray(t["pmra_corrected_err"], float), np.asarray(t["pmdec_corrected_err"], float)
    f6 = np.asarray(t["f625w"], float)
    pos = np.isfinite(x) & np.isfinite(y)
    parent = pos & np.isfinite(ma) & np.isfinite(md) & (ea < args.emax) & (ed < args.emax)
    s = np.interp(R, SIG_R, SIG_V)
    chi2 = ma**2/(s**2+ea**2) + md**2/(s**2+ed**2)
    FC = np.array([-2.58, 5.37])
    away = parent & (chi2 > args.chi2) & (np.hypot(ma, md) > args.rcl)
    field = away & (np.hypot(ma-FC[0], md-FC[1]) < args.rwin)
    control = away & (np.hypot(ma+FC[0], md+FC[1]) < args.rwin)
    # footprint
    c = args.cell
    ix = np.floor(x[pos]/c).astype(np.int64)+10**5; iy = np.floor(y[pos]/c).astype(np.int64)+10**5
    cells = np.unique(ix*10**6+iy)
    Rc = np.hypot((cells//10**6-10**5+.5)*c, (cells % 10**6-10**5+.5)*c)
    edges = np.arange(0, 480, 30.)
    mid = 0.5*(edges[1:]+edges[:-1])
    area = np.array([np.sum((Rc >= a) & (Rc < b))*c**2/3600. for a, b in zip(edges[:-1], edges[1:])])   # arcmin^2
    sets = [("all field stars", field), ("CONTROL: mirrored window", control), ("F625W < 20", field & (f6 < 20)),
            ("20 <= F625W < 22", field & (f6 >= 20) & (f6 < 22)), ("22 <= F625W < 24", field & (f6 >= 22) & (f6 < 24))]
    out = dict(selection=dict(emax=args.emax, chi2=args.chi2, sigma_R=SIG_R.tolist(), sigma=SIG_V.tolist()),
               edges_arcsec=edges.tolist(), area_arcmin2=area.tolist(), n_parent=int(parent.sum()), n_field=int(field.sum()),
               profiles={})
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
    a = ax[0]
    for col, (lab, sel) in zip(COLORS, sets):
        n = np.histogram(R[sel], edges)[0].astype(float)
        good = area > 0.05
        dens = np.where(good, n/np.maximum(area, 1e-9), np.nan)
        err = np.where(good, np.sqrt(n)/np.maximum(area, 1e-9), np.nan)
        out["profiles"][lab] = dict(n=n.tolist(), density_per_arcmin2=dens.tolist(), err=err.tolist())
        a.errorbar(mid, dens, err, fmt="s--" if lab.startswith("CONTROL") else "o-", color=col, ms=4, lw=1.5, capsize=0,
                   label=f"{lab} (N = {int(n.sum()):,})")
    pf, pc = out["profiles"]["all field stars"], out["profiles"]["CONTROL: mirrored window"]
    net = np.array(pf["density_per_arcmin2"])-np.array(pc["density_per_arcmin2"])
    nerr = np.hypot(np.array(pf["err"]), np.array(pc["err"]))
    out["profiles"]["field minus control"] = dict(density_per_arcmin2=net.tolist(), err=nerr.tolist())
    a.errorbar(mid, net, nerr, fmt="D-", color="#0d366b", ms=5, lw=2.2, capsize=0, label="field minus control (leakage removed)")
    a.set_yscale("log"); a.set_xlabel("R [arcsec]", color=INK); a.set_ylabel("field-star surface density [arcmin$^{-2}$]", color=INK)
    a.set_title(f"field stars: |mu - mu_field| < {args.rwin:g}, |mu| > {args.rcl:g} mas/yr, chi2 > {args.chi2:g}; PM err < {args.emax}",
                fontsize=9.5, color=INK, loc="left")
    a.axvline(340, color=MUTED, lw=0.8, ls=":"); a.text(343, a.get_ylim()[0]*1.3 if a.get_ylim()[0] > 0 else 1, "no F625W/F814W beyond ~380''",
                                                        fontsize=7.5, color=MUTED, rotation=90, va="bottom")
    a.grid(True, color=GRID, lw=0.6); a.legend(fontsize=8, frameon=False)
    # coverage
    a = ax[1]
    ann = np.pi*(edges[1:]**2-edges[:-1]**2)/3600.
    a.plot(mid, area/ann, "o-", color=MUTED); a.set_ylim(0, 1.05)
    a.set_xlabel("R [arcsec]", color=INK); a.set_ylabel("footprint coverage of annulus", color=INK)
    a.set_title(f"coverage ({c:g}'' cell occupancy)", fontsize=10, color=INK, loc="left"); a.grid(True, color=GRID, lw=0.6)
    # vector-point diagram
    a = ax[2]
    sub = parent & (np.random.default_rng(1).random(len(R)) < 1.0)
    H, xe, ye = np.histogram2d(ma[sub], md[sub], bins=(np.linspace(-10, 10, 201), np.linspace(-10, 10, 201)))
    a.pcolormesh(xe, ye, H.T, cmap=RAMP, norm=LogNorm(vmin=1, vmax=H.max()), rasterized=True)
    fx, fy = np.median(ma[field]), np.median(md[field])
    th = np.linspace(0, 2*np.pi, 200)
    a.plot(FC[0]+args.rwin*np.cos(th), FC[1]+args.rwin*np.sin(th), "-", color="#2a78d6", lw=1.5, label="field window")
    a.plot(-FC[0]+args.rwin*np.cos(th), -FC[1]+args.rwin*np.sin(th), "--", color="#6b6a64", lw=1.5, label="control window")
    a.plot(args.rcl*np.cos(th), args.rcl*np.sin(th), ":", color="#eb6834", lw=1.2, label=f"|mu| = {args.rcl:g}")
    a.set_xlabel(r"$\mu_{\alpha*}$ rel. to cluster [mas/yr]", color=INK); a.set_ylabel(r"$\mu_\delta$ rel. [mas/yr]", color=INK)
    a.set_title("vector-point diagram, parent sample (log counts)", fontsize=10, color=INK, loc="left"); a.legend(fontsize=8, frameon=False)
    out["field_median_pm"] = [float(fx), float(fy)]
    for a in ax:
        a.tick_params(colors=MUTED, labelsize=9)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.suptitle("oMEGACat: radial surface density of proper-motion-selected field stars, with a mirrored-window control", fontsize=11.5, color=INK)
    fig.tight_layout()
    fig.savefig(ROOT/"plots/omegacat_field_counts.png", dpi=130)
    (ROOT/"results/plot_data/omegacat_field_counts.json").write_text(json.dumps(out, indent=1))
    print("parent", out["n_parent"], "field", out["n_field"], "field median PM", out["field_median_pm"])
    for lab, pr in out["profiles"].items():
        print(lab, np.round(pr["density_per_arcmin2"], 1).tolist())


if __name__ == "__main__":
    main()
