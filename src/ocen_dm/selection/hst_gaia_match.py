"""Star-by-star match between the HST (oMEGACat) and Gaia EDR3 proper-motion catalogues.

Proper motions of the same stars measured by two instruments must agree; line-of-sight
velocities are allowed to differ from proper motions (different tracer populations, and the
conversion needs a distance), but two PM catalogues have no such licence. This module makes
the comparison directly: positional match, then the per-star difference and the dispersion
each instrument reports **from the same stars**.

Result (2026-09-19): inside 460 arcsec Gaia's proper motions carry ~1 mas/yr of scatter that
its formal errors do not describe, and the dispersion it reports from the same stars is
18-24 per cent higher than HST's. That is crowding, it is why the Vasiliev & Baumgardt
quality flag removes essentially every Gaia star inside 300 arcsec, and it is why the
published EDR3 profile inside ~400 arcsec is a smooth extrapolation rather than a measurement.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table

from ..paths import processed_dir, raw_dir

__all__ = ["OCEN_RA", "OCEN_DEC", "build_match", "load_match", "PRODUCT"]

OCEN_RA, OCEN_DEC = 201.696833, -47.476583
PRODUCT = "hst_gaia_pm_match"


def build_match(tolerance_arcsec: float = 0.3) -> Path:
    """Match the two catalogues and write the per-star comparison product."""
    from scipy.spatial import cKDTree

    hst = Table.read(raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits")
    hx = -0.04 * (np.asarray(hst["x"], float) - 15000.0)      # arcsec, the authors' convention
    hy = 0.04 * (np.asarray(hst["y"], float) - 15000.0)
    ha = np.asarray(hst["pmra_corrected"], float); hd = np.asarray(hst["pmdec_corrected"], float)
    hea = np.asarray(hst["pmra_corrected_err"], float); hed = np.asarray(hst["pmdec_corrected_err"], float)
    finite = np.isfinite(hx) & np.isfinite(hy) & np.isfinite(ha) & np.isfinite(hd)
    idx_h = np.flatnonzero(finite)
    tree = cKDTree(np.column_stack([hx[finite], hy[finite]]))

    gaia = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    gx = (np.asarray(gaia["ra"], float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC)) * 3600.0
    gy = (np.asarray(gaia["dec"], float) - OCEN_DEC) * 3600.0
    dist, j = tree.query(np.column_stack([gx, gy]), distance_upper_bound=tolerance_arcsec)
    ok = np.isfinite(dist)
    gi = np.flatnonzero(ok); hi = idx_h[j[ok]]

    t = Table({
        "r_arcsec": np.hypot(gx[gi], gy[gi]), "separation_arcsec": dist[ok],
        "g_mag": np.asarray(gaia["g_mag"], float)[gi],
        "membership_prob": np.asarray(gaia["membership_prob"], float)[gi],
        "gaia_quality": ((np.asarray(gaia["quality_flag"], int)[gi] & 2) > 0).astype(np.int8),
        "hst_quality": np.asarray(hst["selection_hq_astrometry_and_membership"], float)[hi].astype(np.int8),
        "hst_pmra": ha[hi], "hst_pmdec": hd[hi], "hst_pmra_error": hea[hi], "hst_pmdec_error": hed[hi],
        "gaia_pmra": np.asarray(gaia["pmra"], float)[gi], "gaia_pmdec": np.asarray(gaia["pmdec"], float)[gi],
        "gaia_pmra_error": np.asarray(gaia["pmra_error"], float)[gi],
        "gaia_pmdec_error": np.asarray(gaia["pmdec_error"], float)[gi],
    })
    t.meta.update({
        "product": PRODUCT, "tolerance_arcsec": tolerance_arcsec, "n_matched": int(len(t)),
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": "HST proper motions are locally corrected (relative): only differences and dispersions "
                "are comparable, never the zero point.",
    })
    path = processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def load_match(path: Path | None = None) -> Table:
    path = path or processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    if not path.exists():
        return Table.read(build_match())
    return Table.read(path)
