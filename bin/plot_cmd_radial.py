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
Rows: (1) PM members; (2) PM members that the published flag removes (its CMD cut);
(3) PM non-members (field).
Colour F625W - F814W, magnitude F625W (Vega, not extinction-corrected; E(B-V) ~ 0.12 for omega Cen).

Usage: python bin/plot_cmd_radial.py [--edges 0 30 60 120 200 340]
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
CBINS = np.arange(-0.2, 2.0, 0.01)
MBINS = np.arange(12.0, 25.5, 0.04)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--edges", nargs="+", type=float, default=[0, 30, 60, 120, 200, 340])
    args = p.parse_args()
    t = fits.open(ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits", memmap=True)[1].data
    x = -0.04*(np.asarray(t["x"], float)-15000.)
    y = 0.04*(np.asarray(t["y"], float)-15000.)
    R = np.hypot(x, y)
    f6, f8 = np.asarray(t["f625w"], float), np.asarray(t["f814w"], float)
    hq = np.asarray(t["selection_hq_astrometry"]) == 1
    flag = np.asarray(t["selection_hq_astrometry_and_membership"]) == 1
    ma, md = np.asarray(t["pmra_corrected"], float), np.asarray(t["pmdec_corrected"], float)
    ea, ed = np.asarray(t["pmra_corrected_err"], float), np.asarray(t["pmdec_corrected_err"], float)
    good = np.isfinite(f6) & np.isfinite(f8)
    edges = args.edges
    n = len(edges)-1
    fig, axes = plt.subplots(3, n, figsize=(3.4*n, 14.5), sharex=True, sharey=True)
    out = dict(source="oMEGACat catalog_and_selections.fits", edges_arcsec=edges, annuli=[])
    for i in range(n):
        ann = (R >= edges[i]) & (R < edges[i+1]) & good & hq & np.isfinite(ma) & np.isfinite(md)
        core = ann & (np.hypot(ma, md) < 3)
        mad = 1.4826*np.median(np.abs(np.concatenate((ma[core]-np.median(ma[core]), md[core]-np.median(md[core])))))
        e2 = np.median(np.concatenate((ea[core], ed[core])))**2
        sig = float(np.sqrt(max(mad**2-e2, 1e-4)))
        chi2 = ma**2/(sig**2+ea**2) + md**2/(sig**2+ed**2)
        sel_m = ann & (chi2 < 11.83)
        sel_c = sel_m & ~flag
        sel_f = ann & ~(chi2 < 11.83)
        out["annuli"].append(dict(r_in=edges[i], r_out=edges[i+1], sigma_pm_masyr=sig, n_hq=int(ann.sum()),
                                  n_pm_members=int(sel_m.sum()), n_pm_members_removed_by_published_flag=int(sel_c.sum()),
                                  n_pm_nonmembers=int(sel_f.sum())))
        for row, (sel, lab) in enumerate(((sel_m, f"PM members\n(sigma_PM = {sig:.2f} mas/yr)"),
                                          (sel_c, "PM members cut by\npublished flag"), (sel_f, "PM non-members\n(field)"))):
            ax = axes[row, i]
            H, _, _ = np.histogram2d(f6[sel]-f8[sel], f6[sel], bins=(CBINS, MBINS))
            if H.sum() > 0:
                ax.pcolormesh(CBINS, MBINS, H.T, cmap=RAMP, norm=LogNorm(vmin=1, vmax=max(H.max(), 2)), rasterized=True)
            ax.set_title(f"{edges[i]:.0f}-{edges[i+1]:.0f} arcsec: N = {sel.sum():,}\n{lab}", fontsize=9.5, color=INK, loc="left")
            ax.tick_params(colors=MUTED, labelsize=9)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if row == 2:
                ax.set_xlabel("F625W - F814W", color=INK)
            if i == 0:
                ax.set_ylabel("F625W", color=INK)
    axes[0, 0].set_xlim(CBINS[0], CBINS[-1]); axes[0, 0].set_ylim(MBINS[-1], MBINS[0])
    fig.suptitle("omega Cen HST CMDs by radius (oMEGACat, hq astrometry): PM-only members (top), PM members that the published "
                 "membership flag removes (middle), PM non-members (bottom)\ncolour = log counts per 0.01 x 0.04 mag cell",
                 fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(ROOT/"plots/cmd_radial_omegacat.png", dpi=130)
    (ROOT/"results/plot_data/cmd_radial_omegacat.json").write_text(json.dumps(out, indent=1))
    for a in out["annuli"]:
        print(a)


if __name__ == "__main__":
    main()
