#!/usr/bin/env python3
"""HST colour-magnitude diagrams of omega Cen members in radial annuli (oMEGACat catalogue).

Data: data/raw/omegacat_vi_kinematics/catalog_and_selections.fits (oMEGACat; Haeberle+2025,
arXiv:2503.04903; Zenodo 10.5281/zenodo.14978551). Radius from the catalogue pixel frame
(centre 15000, 15000; 0.04 arcsec/px), as in bin/build_count_profiles.py.
Membership by proper motion ONLY. Parent sample: `selection_hq_astrometry` (oMEGACat's quality
selection: good astrometry and photometry in both filters; ends near 340 arcsec). In each annulus
the intrinsic PM dispersion sigma (per component; PMs are relative to the cluster) is measured
robustly (1.4826 MAD of stars within 3 mas/yr, error-deconvolved with the median error), and a star
is a member if mu_a^2/(sigma^2+e_a^2) + mu_d^2/(sigma^2+e_d^2) < 11.83 (2 dof, 99.73%).
NB: the published flag `selection_hq_astrometry_and_membership` is NOT PM-only: 97% of the
hq-astrometry stars it rejects have |mu| < 2 mas/yr, and they form sharp-edged sequences either
side of the main ridge (blue MS, red populations/binaries), i.e. it includes a CMD selection.
Rows: (1) PM members (Hess diagram); (2) zoom on the white-dwarf region (F625W - F814W < 0.6,
F625W 19.5-26) as points. Row 2 uses a different parent, because the hq-astrometry selection removes
almost all faint blue stars (1096 of 26492 with PMs at colour < 0.35, F625W > 21): the photometric
quality flags only (selection_hq_f625w & selection_hq_f814w), same PM membership cut (sigma of the
annulus from the hq stars). A looser parent (any photometry, PM errors < 0.5 mas/yr) gives a diffuse
cloud of poorly measured stars and no sequence (tested 2026-10-08). Estimated start of the WD
cooling sequence: (m-M)_0 = 13.67 (5.43 kpc), A_F625W ~ 0.32 (E(B-V) = 0.12) -> F625W ~ 23.5-24 for
M ~ 9.5-10, i.e. within ~1 mag of the depth of this catalogue (photometric-quality sample 99.9th
percentile F625W = 25.2), so at most the top of the sequence is accessible.
Absolute magnitudes and dereddened colour: M_F625W = F625W - (m-M)_0 - A_F625W, colour
(F625W - F814W)_0, with D = 5.43 kpc and E(B-V) = 0.12 (Harris), A_V = 3.1 E(B-V), filter extinctions
from the MIST bolometric-correction tables (A_band = BC(A_V=0) - BC(A_V), median over the main
sequence), exactly as in bin/build_mass_function.py: A_F625W = 0.316, A_F814W = 0.219. Overlaid:
MIST v1.2 main sequence (ACS F625W/F814W; age 12.5 Gyr, [Fe/H] = -1.53, the mass-function ridge fit).

Usage: python bin/plot_cmd_radial.py [--edges 0 30 60 120 200 340] [--distance 5.43] [--ebv 0.12]
Output: plots/cmd_radial_omegacat.png, results/plot_data/cmd_radial_omegacat.json (counts per annulus).
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
from astropy.io import fits

RAMP = LinearSegmentedColormap.from_list("blue", ["#e8f0fb", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
INK, MUTED = "#1f1f1e", "#6b6a64"
CBINS = np.arange(-0.9, 1.9, 0.01)
MBINS = np.arange(-2.0, 12.0, 0.04)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--edges", nargs="+", type=float, default=[0, 30, 60, 120, 200, 340])
    p.add_argument("--distance", type=float, default=5.43, help="kpc")
    p.add_argument("--ebv", type=float, default=0.12)
    p.add_argument("--feh", type=float, default=-1.53, help="isochrone [Fe/H] (overlay only)")
    p.add_argument("--age", type=float, default=12.5, help="isochrone age [Gyr] (overlay only)")
    args = p.parse_args()
    sys.path.insert(0, str(ROOT/"bin"))
    from build_mass_function import BANDS, Isochrone
    iso = Isochrone(args.age, args.feh, args.distance, args.ebv)
    dm0, a6, a8 = iso.dm0, iso.a625, iso.a814
    iso_M = iso.m625 - dm0 - a6
    iso_c = (iso.m625 - iso.m814) - (a6 - a8)
    t = fits.open(ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits", memmap=True)[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.)
    y = 0.04*(np.asarray(t["y"], float)-15000.)
    R = np.hypot(x, y)
    f6 = np.asarray(t["f625w"], float) - dm0 - a6           # absolute, dereddened
    f8 = np.asarray(t["f814w"], float) - dm0 - a8
    hq = np.asarray(t["selection_hq_astrometry"]) == 1
    flag = np.asarray(t["selection_hq_astrometry_and_membership"]) == 1
    ma, md = np.asarray(t["pmra_corrected"], float), np.asarray(t["pmdec_corrected"], float)
    ea, ed = np.asarray(t["pmra_corrected_err"], float), np.asarray(t["pmdec_corrected_err"], float)
    good = np.isfinite(f6) & np.isfinite(f8)
    edges = args.edges
    n = len(edges)-1
    fig, axes = plt.subplots(2, n, figsize=(3.4*n, 10.5))
    loose = good & np.isfinite(ma) & np.isfinite(md) & (np.asarray(t["selection_hq_f625w"]) == 1) & (np.asarray(t["selection_hq_f814w"]) == 1)
    out = dict(source="oMEGACat catalog_and_selections.fits", edges_arcsec=edges, distance_kpc=args.distance, ebv=args.ebv,
               dm0=dm0, A_F625W=a6, A_F814W=a8, isochrone=dict(age_gyr=args.age, feh=args.feh, bands=BANDS), annuli=[])
    for i in range(n):
        ann = (R >= edges[i]) & (R < edges[i+1]) & good & hq & np.isfinite(ma) & np.isfinite(md)
        core = ann & (np.hypot(ma, md) < 3)
        mad = 1.4826*np.median(np.abs(np.concatenate((ma[core]-np.median(ma[core]), md[core]-np.median(md[core])))))
        e2 = np.median(np.concatenate((ea[core], ed[core])))**2
        sig = float(np.sqrt(max(mad**2-e2, 1e-4)))
        chi2 = ma**2/(sig**2+ea**2) + md**2/(sig**2+ed**2)
        sel_m = ann & (chi2 < 11.83)
        ring = (R >= edges[i]) & (R < edges[i+1])
        sel_w = ring & loose & (chi2 < 11.83)
        wd = sel_w & (f6-f8 < 0.5) & (f6 > 5.5)
        out["annuli"].append(dict(r_in=edges[i], r_out=edges[i+1], sigma_pm_masyr=sig, n_hq=int(ann.sum()),
                                  n_pm_members=int(sel_m.sum()), n_pm_nonmembers=int((ann & ~(chi2 < 11.83)).sum()),
                                  n_loose_pm_members=int(sel_w.sum()), n_loose_members_in_wd_box=int(wd.sum())))
        ax = axes[0, i]
        H, _, _ = np.histogram2d(f6[sel_m]-f8[sel_m], f6[sel_m], bins=(CBINS, MBINS))
        ax.pcolormesh(CBINS, MBINS, H.T, cmap=RAMP, norm=LogNorm(vmin=1, vmax=max(H.max(), 2)), rasterized=True)
        ax.set_xlim(CBINS[0], CBINS[-1]); ax.set_ylim(MBINS[-1], MBINS[0])
        ax.set_title(f"{edges[i]:.0f}-{edges[i+1]:.0f} arcsec: N = {sel_m.sum():,}\nPM members, hq astrometry\n"
                     f"(sigma_PM = {sig:.2f} mas/yr)", fontsize=9.5, color=INK, loc="left")
        ax.add_patch(plt.Rectangle((-0.9, 5.5), 1.4, 6.5, fill=False, ec=MUTED, lw=0.8, ls="--"))
        ax.plot(iso_c, iso_M, "-", color="#eb6834", lw=1.2, label="MIST MS")
        if i == 0:
            ax.legend(fontsize=8, frameon=False, loc="lower left")
        ax = axes[1, i]
        ax.plot(f6[wd]-f8[wd], f6[wd], ".", color="#2a78d6", ms=2.2, alpha=0.6, mew=0, rasterized=True)
        ax.plot(iso_c, iso_M, "-", color="#eb6834", lw=1.2)
        ax.set_xlim(-0.9, 0.5); ax.set_ylim(12, 5.5)
        ax.axhspan(9.5, 10.0, color="#efeee8", zorder=0)
        if i == 0:
            ax.text(-0.85, 9.45, "rough top of WD sequence (estimate)", fontsize=8, color=MUTED, va="bottom")
        ax.set_title(f"{edges[i]:.0f}-{edges[i+1]:.0f} arcsec: N = {wd.sum():,} in box\nPM members, photometric-quality\n"
                     "flags only (no astrometric hq cut)", fontsize=9.5, color=INK, loc="left")
        ax.set_xlabel("(F625W - F814W)$_0$", color=INK)
        for row in (0, 1):
            a = axes[row, i]
            a.tick_params(colors=MUTED, labelsize=9)
            for sp in ("top", "right"):
                a.spines[sp].set_visible(False)
            if i == 0:
                a.set_ylabel("$M_{F625W}$", color=INK)
    fig.suptitle(f"omega Cen HST CMDs by radius (oMEGACat), PM-selected members; D = {args.distance} kpc, E(B-V) = {args.ebv} "
                 f"(A_F625W = {a6:.3f}, A_F814W = {a8:.3f})\nTop: hq-astrometry sample, log counts (dashed box = zoom). "
                 "Bottom: white-dwarf region, photometric-quality sample", fontsize=10.5, color=INK)
    fig.tight_layout()
    fig.savefig(ROOT/"plots/cmd_radial_omegacat.png", dpi=130)
    (ROOT/"results/plot_data/cmd_radial_omegacat.json").write_text(json.dumps(out, indent=1))
    for a in out["annuli"]:
        print(a)


if __name__ == "__main__":
    main()
