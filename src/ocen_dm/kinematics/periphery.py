"""Beyond the Gaia catalogue's edge: does the cluster stop, and where?

The Vasiliev & Baumgardt member catalogue is truncated at **0.67 deg (2413 arcsec, 63 pc)**,
which is their retrieval radius for this cluster, not a physical boundary. The Jacobi radius
is larger: computed in the McMillan (2017) Milky Way potential for a cluster mass of
3.55e6 Msun, it is **90 pc (0.95 deg) at pericentre** and 185 pc at the present galactocentric
radius, so the Gaia sample stops at roughly 0.7 of the pericentric Jacobi radius. Everything
about "where the cluster ends" therefore has to come from other data.

Two datasets in this project reach further:

* ``kuzma2025_periphery`` -- Pristine CaHK plus Gaia, 157481 stars to 5.1 deg, with
  membership probabilities. **Its membership uses the proper motions**, so measuring a PM
  dispersion from it is circular: the answer tracks the selection window. This module
  therefore always reports the profile at two fixed windows, and the spread between them is
  the size of the systematic.
* ``kuzma2026_spectroscopy`` -- 592 line-of-sight velocities to 3.15 deg, 157 flagged
  members. These carry **no proper-motion selection**, so they are the honest probe.

Measured 2026-09-19. The proper-motion dispersion at a fixed 0.8 mas/yr window runs
6.7, 6.6, 5.7, 6.5, 8.9 km/s from 44 to 377 pc: flat, not falling. Widening the window to
1.2 mas/yr turns that into a steep rise to 12.6 km/s, which is the selection and not the
cluster. The spectroscopy gives 6.12 +- 0.48 km/s at 67 pc and 8.36 +- 1.07 at 83 pc.

So the cluster does not stop at the edge of the Gaia data. The dispersion **stops falling**
near 50-100 pc and settles at 6-8 km/s. Whether that plateau is unbound debris and potential
escapers, or a bound halo, is the question this project exists to answer, and it cannot be
settled without the Kuzma selection function (WP6).
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from ..paths import processed_dir
from .outer_profile import OCEN_DEC, OCEN_RA, dispersion_ml

__all__ = ["KMS_PER_MASYR_KPC", "GAIA_EDGE_DEG", "R_JACOBI_PERI_PC", "R_JACOBI_NOW_PC",
           "radius_deg", "periphery_pm_profile", "periphery_los_profile"]

KMS_PER_MASYR_KPC = 4.740470446
#: retrieval radius of the Vasiliev & Baumgardt catalogue, not a physical edge
GAIA_EDGE_DEG = 0.67
#: McMillan (2017) potential, M_cluster = 3.55e6 Msun (Baumgardt & Hilker 2018)
R_JACOBI_PERI_PC = 89.9
R_JACOBI_NOW_PC = 185.1


def radius_deg(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    return np.hypot((np.asarray(ra, float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC)),
                    np.asarray(dec, float) - OCEN_DEC)


def periphery_pm_profile(window_masyr: float = 0.8, prob_min: float = 0.5,
                         edges_deg=(0.20, 0.40, 0.67, 1.00, 1.60, 2.60, 5.20),
                         distance_kpc: float = 5.43, mu_sys=(-3.257, -6.730)) -> Table:
    """PM dispersion of the Pristine periphery sample inside a **fixed** PM window.

    The window must be fixed and stated, because the catalogue's own membership probability
    is proper-motion based: letting it set the sample makes the measured dispersion a
    readout of the selection. Compare two windows and treat the difference as systematic.
    """
    t = Table.read(processed_dir() / "tails" / "kuzma2025_periphery.ecsv")
    r = radius_deg(t["ra"], t["dec"])
    pa, pd = np.asarray(t["pmra"], float), np.asarray(t["pmdec"], float)
    ea, ed = np.asarray(t["pmra_error"], float), np.asarray(t["pmdec_error"], float)
    good = (np.asarray(t["membership_prob"], float) > prob_min) \
        & (np.hypot(pa - mu_sys[0], pd - mu_sys[1]) < window_masyr)
    rows = []
    for lo, hi in zip(edges_deg[:-1], edges_deg[1:]):
        m = good & (r >= lo) & (r < hi)
        if m.sum() < 5:
            continue
        a = dispersion_ml(pa[m] - mu_sys[0], ea[m]); b = dispersion_ml(pd[m] - mu_sys[1], ed[m])
        s = float(np.sqrt(0.5 * (a[0] ** 2 + b[0] ** 2)))
        e = float(0.5 * np.hypot(a[1], b[1]))
        rpc = float(np.median(r[m]) * np.pi / 180 * distance_kpc * 1e3)
        k = KMS_PER_MASYR_KPC * distance_kpc
        rows.append((lo, hi, rpc, int(m.sum()), s, e, s * k, e * k))
    return Table(rows=rows, names=("r_lower_deg", "r_upper_deg", "r_pc", "n_stars",
                                   "sigma_pm", "sigma_pm_err", "sigma_kms", "sigma_kms_err"),
                 meta={"window_masyr": window_masyr, "prob_min": prob_min,
                       "warning": "membership is PM-based; the window is part of the measurement"})


def periphery_los_profile(edges_deg=(0.0, 0.8, 1.2, 3.2), distance_kpc: float = 5.43) -> Table:
    """Line-of-sight dispersion of the flagged spectroscopic members: no PM selection."""
    t = Table.read(processed_dir() / "tails" / "kuzma2026_spectroscopy.ecsv")
    r = radius_deg(t["ra"], t["dec"])
    v, ev = np.asarray(t["vlos"], float), np.asarray(t["vlos_error"], float)
    ok = (np.asarray(t["member_flag"]) == "True") & np.isfinite(v) & np.isfinite(ev)
    v0 = float(np.median(v[ok]))
    rows = []
    for lo, hi in zip(edges_deg[:-1], edges_deg[1:]):
        m = ok & (r >= lo) & (r < hi)
        if m.sum() < 5:
            continue
        a = dispersion_ml(v[m] - v0, ev[m])
        rows.append((lo, hi, float(np.median(r[m]) * np.pi / 180 * distance_kpc * 1e3),
                     int(m.sum()), float(a[0]), float(a[1])))
    return Table(rows=rows, names=("r_lower_deg", "r_upper_deg", "r_pc", "n_stars",
                                   "sigma_kms", "sigma_kms_err"),
                 meta={"v_sys_kms": v0, "selection": "member_flag, no proper-motion cut"})
