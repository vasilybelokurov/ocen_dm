#!/usr/bin/env python3
"""Omega Cen rotation vs projected radius: every rotation curve available to the Jeans fits, in km/s.

1. MUSE line of sight, PUBLISHED: oMEGACat VI (Haeberle et al. 2025, arXiv:2503.04903) los_profile.fits, Table 3: v_rot per annulus (omegacat_vi_los_rotation.ecsv)
   (sigma_los fitted jointly with vbar = v_rot sin(phi - theta_0)). Used in the fits for the MUSE dataset.
2. Proper-motion rotation, Vasiliev & Baumgardt (2021) EDR3: mean tangential PM per annulus, median and 16-84%
   (vasiliev2021_ocen_pm_profiles.ecsv). Used in the fits for HST tangential and (default) Gaia tangential.
3. Proper-motion rotation, ours: |mean tangential PM| per annulus from our Gaia EDR3 mixture fit
   (outer_gaia.load_edr3_profile, column mean_pmt). Available as --gaia-rotation ours; not used in current runs.
PMs converted with 4.74047 * D[kpc] km/s per mas/yr at D = 5.43 kpc. Line-of-sight and plane-of-sky rotation are
different projections of the same rotation (inclination ~50 deg), so curve 1 need not equal curves 2-3.
Usage: PYTHONPATH=src python bin/plot_rotation_curves.py   ->  plots/rotation_curves_20261113.png
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

from ocen_dm.kinematics.outer_gaia import load_edr3_profile
from ocen_dm.paths import processed_dir

KMS = 4.74047*5.43          # km/s per mas/yr


def main():
    los = Table.read(processed_dir()/"kinematics"/"omegacat_vi_los_rotation.ecsv")
    vb = Table.read(processed_dir()/"kinematics"/"vasiliev2021_ocen_pm_profiles.ecsv")
    g = load_edr3_profile()
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.errorbar(los["r_median"], los["v_rot"], [los["v_rot_err_lo"], los["v_rot_err_hi"]], fmt="o", ms=3.5,
                color="C3", lw=0.8, label="line of sight, oMEGACat VI (Häberle+2025, MUSE) — used for MUSE")
    ax.fill_between(vb["r"], vb["vrot_pm_p16"]*KMS, vb["vrot_pm_p84"]*KMS, color="C0", alpha=0.25, lw=0)
    ax.plot(vb["r"], vb["vrot_pm"]*KMS, color="C0", lw=1.8,
            label="proper motion, Vasiliev & Baumgardt 2021 — used for HST and Gaia tangential")
    ax.plot(g["r_median"], np.abs(g["mean_pmt"])*KMS, "s", color="C2", ms=5,
            label="proper motion, ours (Gaia mixture fit) — not used")
    ax.axvspan(0, 100, color="0.9", zorder=0)
    ax.text(8, 0.4, "PM curve not measured here\n(Gaia crowding): model extrapolation", fontsize=7, color="0.35")
    ax.set_xscale("log"); ax.set_xlim(4, 3000); ax.set_ylim(0, 13)
    ax.set_xlabel("projected radius R [arcsec]   (100\" = 2.6 pc)"); ax.set_ylabel("rotation velocity [km/s]")
    sec = ax.secondary_xaxis("top", functions=(lambda x: x*5430/206264.8, lambda p: p*206264.8/5430))
    sec.set_xlabel("R [pc]")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    out = ROOT/"plots/rotation_curves_20261113.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
