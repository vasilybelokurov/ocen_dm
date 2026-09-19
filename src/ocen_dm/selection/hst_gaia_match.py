"""Star-by-star match between the HST (oMEGACat) and Gaia EDR3 proper-motion catalogues.

Proper motions of the same stars measured by two instruments must agree; line-of-sight
velocities are allowed to differ from proper motions (different tracer populations, and the
conversion needs a distance), but two PM catalogues have no such licence. This module makes
the comparison directly: positional match, then the per-star difference and the dispersion
each instrument reports **from the same stars**.

Result (2026-09-19, after the match was corrected twice -- see JOURNAL): on identical stars
the **unflagged** Gaia catalogue carries 0.44-0.90 mas/yr of proper-motion scatter that its
quoted errors do not describe, growing inwards, and reports a dispersion 25-43 per cent
higher than HST. Applying the Vasiliev & Baumgardt astrometric quality flag removes those
stars: in the one annulus where both catalogues are usable (300-380 arcsec) flagged Gaia
gives 0.482 +- 0.034 against HST's 0.524 +- 0.002 mas/yr, a ratio of 0.92 +- 0.07. The two
instruments agree; the flag is what makes them agree, and it passes no star inside 200
arcsec, so the published EDR3 profile there is extrapolation rather than measurement.

"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table

from ..paths import processed_dir, raw_dir

__all__ = ["OCEN_RA", "OCEN_DEC", "build_match", "load_match", "PRODUCT"]

from ..cluster import OCEN_DEC, OCEN_RA  # single definition, derived from the data
PRODUCT = "hst_gaia_pm_match"


def build_match(tolerance_arcsec: float = 0.06) -> Path:
    """Match the two catalogues and write the per-star comparison product."""
    from scipy.spatial import cKDTree

    hst = Table.read(raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits")
    # match on the catalogue's own sky coordinates. The pixel-based x, y are good enough for
    # plots but not for a 0.05 arcsec match: using them the separations peaked at 0.18 arcsec
    # and the match rate fell from 93 to 1 per cent between the centre and 460 arcsec, so most
    # "matches" beyond the core were neighbouring stars (2026-09-19).
    hx = (np.asarray(hst["RA"], float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC)) * 3600.0
    hy = (np.asarray(hst["DEC"], float) - OCEN_DEC) * 3600.0
    ha = np.asarray(hst["pmra_corrected"], float); hd = np.asarray(hst["pmdec_corrected"], float)
    hea = np.asarray(hst["pmra_corrected_err"], float); hed = np.asarray(hst["pmdec_corrected_err"], float)
    finite = np.isfinite(hx) & np.isfinite(hy) & np.isfinite(ha) & np.isfinite(hd)
    idx_h = np.flatnonzero(finite)
    tree = cKDTree(np.column_stack([hx[finite], hy[finite]]))

    gaia = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    gx = (np.asarray(gaia["ra"], float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC)) * 3600.0
    gy = (np.asarray(gaia["dec"], float) - OCEN_DEC) * 3600.0
    # the two catalogues are at different epochs, so the cluster's systemic proper motion
    # displaces them by ~0.10 arcsec (7.5 mas/yr over ~13 years). Measure that offset from a
    # generous first pass and remove it, then match tightly.
    d0, j0 = tree.query(np.column_stack([gx, gy]), distance_upper_bound=0.4)
    ok0 = np.isfinite(d0)
    off_x = float(np.median(gx[ok0] - hx[finite][j0[ok0]]))
    off_y = float(np.median(gy[ok0] - hy[finite][j0[ok0]]))
    dist, j = tree.query(np.column_stack([gx - off_x, gy - off_y]), distance_upper_bound=tolerance_arcsec)
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
        "epoch_offset_arcsec": [off_x, off_y],
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
