#!/usr/bin/env python3
"""Inner 10 arcsec of oMEGACat: field-window and control-window stars, same selection as bin/plot_field_counts.py.

Selection (identical to the full-radius plot): PM errors < 0.5 mas/yr; field window = within 3 mas/yr of the field
centroid (-2.58, +5.37) mas/yr (relative PM), > 4 mas/yr from the cluster and chi2 > 25 against the cluster
(intrinsic sigma 0.81 mas/yr, the 0-30'' value); control = the same window mirrored through the cluster.
With ~19 field stars per arcmin^2 (field minus control at R < 370'') only ~1.7 field stars are expected inside
10'' (0.087 arcmin^2), so the inner region is shown as counts per 2'' annulus with the Poisson expectation,
a vector-point diagram and a sky map. Stars in either window are listed (they may also be fast cluster members).

Usage: python bin/plot_field_counts_inner.py
Output: plots/omegacat_field_counts_inner10.png, results/plot_data/omegacat_field_counts_inner10.json
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
import numpy as np
from astropy.io import fits

INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"
FIELD_DENSITY = 19.0          # arcmin^-2, field minus control at R < 370'' (results/plot_data/omegacat_field_counts.json)
FC = np.array([-2.58, 5.37])
RMAX, STEP = 10.0, 2.0


def main():
    t = fits.open(ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits", memmap=True)[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.); y = 0.04*(np.asarray(t["y"], float)-15000.)
    R = np.hypot(x, y)
    ma, md = np.asarray(t["pmra_corrected"], float), np.asarray(t["pmdec_corrected"], float)
    ea, ed = np.asarray(t["pmra_corrected_err"], float), np.asarray(t["pmdec_corrected_err"], float)
    f6, f8 = np.asarray(t["f625w"], float), np.asarray(t["f814w"], float)
    inner = R < RMAX
    par = inner & np.isfinite(ma) & (ea < 0.5) & (ed < 0.5)
    chi2 = ma**2/(0.81**2+ea**2) + md**2/(0.81**2+ed**2)
    away = par & (chi2 > 25) & (np.hypot(ma, md) > 4)
    fld = away & (np.hypot(ma-FC[0], md-FC[1]) < 3)
    ctl = away & (np.hypot(ma+FC[0], md+FC[1]) < 3)
    edges = np.arange(0, RMAX+1e-9, STEP)
    mid = 0.5*(edges[1:]+edges[:-1])
    area = np.pi*(edges[1:]**2-edges[:-1]**2)/3600.
    nf = np.histogram(R[fld], edges)[0]; nc = np.histogram(R[ctl], edges)[0]; nall = np.histogram(R[par], edges)[0]
    expect = FIELD_DENSITY*area
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))
    a = ax[0]
    w = 0.35*STEP
    a.bar(mid-w/2, nf, w, color="#2a78d6", label=f"field window (N = {nf.sum()})")
    a.bar(mid+w/2, nc, w, color="#6b6a64", label=f"control window (N = {nc.sum()})")
    a.step(edges, np.r_[expect, expect[-1]], where="post", color="#eb6834", lw=1.8,
           label=f"uniform field, {FIELD_DENSITY:g} arcmin$^{{-2}}$ (total {expect.sum():.1f})")
    a.set_xlabel("R [arcsec]", color=INK); a.set_ylabel("stars per 2'' annulus", color=INK)
    a.set_title("inner 10'': counts vs uniform-field expectation", fontsize=10, color=INK, loc="left")
    a.legend(fontsize=8, frameon=False); a.grid(True, color=GRID, lw=0.6, axis="y")
    a.set_ylim(0, max(3.5, nf.max()+1, nc.max()+1))
    a2 = a.twinx()
    a2.plot(mid, nall, "o:", color=MUTED, ms=4)
    a2.set_ylabel("all stars with PM err < 0.5 mas/yr (dotted)", color=MUTED); a2.tick_params(colors=MUTED, labelsize=9)
    a = ax[1]
    a.plot(ma[par], md[par], ".", color="#9ec5f4", ms=3, alpha=0.6, rasterized=True, label=f"all inner stars (N = {par.sum():,})")
    a.plot(ma[fld], md[fld], "o", color="#2a78d6", ms=8, mec="white", mew=0.8, label="field window")
    a.plot(ma[ctl], md[ctl], "s", color="#6b6a64", ms=8, mec="white", mew=0.8, label="control window")
    th = np.linspace(0, 2*np.pi, 200)
    a.plot(FC[0]+3*np.cos(th), FC[1]+3*np.sin(th), "-", color="#2a78d6", lw=1.2)
    a.plot(-FC[0]+3*np.cos(th), -FC[1]+3*np.sin(th), "--", color="#6b6a64", lw=1.2)
    a.plot(4*np.cos(th), 4*np.sin(th), ":", color="#eb6834", lw=1.2)
    a.set_xlim(-10, 10); a.set_ylim(-10, 10); a.set_aspect("equal")
    a.set_xlabel(r"$\mu_{\alpha*}$ rel. to cluster [mas/yr]", color=INK); a.set_ylabel(r"$\mu_\delta$ rel. [mas/yr]", color=INK)
    a.set_title("vector-point diagram, R < 10''", fontsize=10, color=INK, loc="left"); a.legend(fontsize=8, frameon=False, loc="lower left")
    a = ax[2]
    a.plot(x[par], y[par], ".", color="#9ec5f4", ms=2.5, alpha=0.6, rasterized=True)
    a.plot(x[fld], y[fld], "o", color="#2a78d6", ms=9, mec="white", mew=0.8)
    a.plot(x[ctl], y[ctl], "s", color="#6b6a64", ms=9, mec="white", mew=0.8)
    for r_ in edges[1:]:
        a.add_patch(plt.Circle((0, 0), r_, fill=False, ec=GRID, lw=0.8))
    a.set_xlim(RMAX, -RMAX); a.set_ylim(-RMAX, RMAX); a.set_aspect("equal")
    a.set_xlabel("x (east) [arcsec]", color=INK); a.set_ylabel("y (north) [arcsec]", color=INK)
    a.set_title("positions (circles every 2'')", fontsize=10, color=INK, loc="left")
    for a in (ax[0], ax[1], ax[2]):
        a.tick_params(colors=MUTED, labelsize=9)
        for sp in ("top",):
            a.spines[sp].set_visible(False)
    fig.suptitle("oMEGACat inner 10'': proper-motion field window vs mirrored control (same selection as the full-radius plot)",
                 fontsize=11.5, color=INK)
    fig.tight_layout()
    fig.savefig(ROOT/"plots/omegacat_field_counts_inner10.png", dpi=130)
    rows = []
    for tag, sel in (("field", fld), ("control", ctl)):
        for i in np.where(sel)[0]:
            rows.append(dict(window=tag, id=int(t["ID"][i]), R_arcsec=float(R[i]), pmra=float(ma[i]), pmdec=float(md[i]),
                             pm_total=float(np.hypot(ma[i], md[i])), e_pmra=float(ea[i]), e_pmdec=float(ed[i]),
                             f625w=float(f6[i]), f814w=float(f8[i])))
    out = dict(edges_arcsec=edges.tolist(), n_field=nf.tolist(), n_control=nc.tolist(), n_all=nall.tolist(),
               expected_uniform_field=expect.tolist(), field_density_assumed=FIELD_DENSITY, stars=rows)
    (ROOT/"results/plot_data/omegacat_field_counts_inner10.json").write_text(json.dumps(out, indent=1))
    for r_ in rows:
        print(r_)


if __name__ == "__main__":
    main()
