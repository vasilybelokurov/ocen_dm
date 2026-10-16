#!/usr/bin/env python3
"""How sensitive is the central HST count deficit to the adopted centre?

Recounts the 8 innermost HST F625W<19 bins (2-12.9 arcsec, the edges of ocen_counts_hst_f625w19) about
centres offset by -2..+2 arcsec in RA and Dec from the adopted Anderson & van der Marel (2010) centre
(oMEGACat pixel 15000, 15000), plus the oMEGACat VI kinematic centre (Haeberle et al. 2025: dRA, dDec =
-1.31, +1.44 arcsec). Coverage is complete inside 30 arcsec, so only the star assignment changes.
Compares with the fixed model expectations saved by bin/check_hst_count_completeness.py (no refit):
Poisson deviance over the 8 bins and the pooled lower-tail probability of the two innermost bins.

Outputs: plots/hst_count_centre_20261015.png, results/plot_data/hst_count_centre_20261015.json.
Usage: python bin/check_count_centre.py
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "ocen-df-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.stats import poisson

ROOT = Path(__file__).resolve().parents[1]
CATALOGUE = ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits"
KINEMATIC = (-1.31, 1.44)        # oMEGACat VI kinematic centre offset (RA, Dec), arcsec


def deviance(n, mu):
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(n > 0, n*np.log(n/mu), 0.)
    return float(2*np.sum(mu-n+term))


def main():
    models = json.loads((ROOT/"results/plot_data/hst_counts_completeness_20261015.json").read_text())["models"]
    with fits.open(CATALOGUE, memmap=True) as h:
        t = h[1].data
        x = -0.04*(np.asarray(t["x"], float)-15000.)       # arcsec, +RA (as in build_count_profiles.py)
        y = 0.04*(np.asarray(t["y"], float)-15000.)
        mag = np.asarray(t["f625w"], float)
    keep = np.isfinite(mag) & (mag < 19.) & (np.hypot(x, y) < 20.)
    x, y = x[keep], y[keep]
    edges = np.geomspace(2., 250., 21)[:9]
    offsets = [(dx, dy) for dx in (-2., -1., 0., 1., 2.) for dy in (-2., -1., 0., 1., 2.)]+[KINEMATIC]
    rows = []
    for dx, dy in offsets:
        n, _ = np.histogram(np.hypot(x-dx, y-dy), edges)
        row = dict(dRA=dx, dDec=dy, counts=n.tolist())
        for label, m in models.items():
            mu = np.asarray(m["mu"], float)
            row[label] = dict(deviance_8=deviance(n, mu), P_inner2=float(poisson.cdf(n[:2].sum(), mu[:2].sum())))
        rows.append(row)
    labels = list(models)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    grid = rows[:25]
    for k, (label, title) in enumerate([(labels[0], "no BH"), (labels[-1], "BH 40,000 Msun")]):
        img = np.array([r[label]["deviance_8"] for r in grid]).reshape(5, 5).T
        im = ax[k].imshow(img, origin="lower", extent=(-2.5, 2.5, -2.5, 2.5), cmap="viridis")
        for r in grid:
            ax[k].text(r["dRA"], r["dDec"], f"{r[label]['deviance_8']:.1f}\nP={r[label]['P_inner2']:.2g}",
                       ha="center", va="center", fontsize=7, color="w")
        ax[k].plot(*KINEMATIC, "r*", ms=12, label=f"kinematic centre: dev {rows[-1][label]['deviance_8']:.1f}, "
                                                   f"P={rows[-1][label]['P_inner2']:.2g}")
        ax[k].set(xlabel="centre offset dRA [arcsec]", ylabel="dDec [arcsec]",
                  title=f"{title}: deviance of 8 inner bins (P: inner 2 pooled)")
        ax[k].legend(fontsize=7, loc="lower left")
        fig.colorbar(im, ax=ax[k], shrink=.8)
    rc = np.sqrt(edges[:-1]*edges[1:])
    area = np.pi*np.diff(edges**2)
    for r in grid:
        ax[2].plot(rc, np.array(r["counts"])/area, color="0.75", lw=.8)
    ax[2].plot(rc, np.array(rows[12]["counts"])/area, "ko-", label="adopted centre")
    ax[2].plot(rc, np.array(rows[-1]["counts"])/area, "r*-", label="kinematic centre")
    for label, c in zip(labels, ("C0", "C1", "C2")):
        ax[2].plot(rc, np.asarray(models[label]["mu"])/area, color=c, ls="--", label=f"model: {label}")
    ax[2].set(xscale="log", xlabel="R [arcsec]", ylabel="stars per arcsec$^2$",
              title="inner counts for all 25 offset centres (grey)")
    ax[2].legend(fontsize=7)
    fig.savefig(ROOT/"plots/hst_count_centre_20261015.png", dpi=140)
    (ROOT/"results/plot_data/hst_count_centre_20261015.json").write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(f"dRA {r['dRA']:+.2f} dDec {r['dDec']:+.2f} inner counts {r['counts'][:3]} | "
              + " | ".join(f"{l}: dev {r[l]['deviance_8']:.1f} P {r[l]['P_inner2']:.2g}" for l in labels))


if __name__ == "__main__":
    main()
