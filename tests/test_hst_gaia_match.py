"""HST against Gaia EDR3 proper motions.

The match is the fragile part: it was wrong twice (pixel coordinates, then the ~0.10 arcsec
epoch offset from the cluster's systemic motion), and both times the wrong pairs inflated
the apparent disagreement. These tests pin the match quality first, then the physics:
unflagged Gaia carries undeclared scatter in the crowded region, and the astrometric
quality flag removes it.
"""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir, raw_dir

_HAS = ((raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits").exists()
        and (processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv").exists())
pytestmark = pytest.mark.skipif(not _HAS, reason="catalogues not present")


def test_match_is_tight_and_epoch_offset_removed():
    """A good match has separations far below the tolerance and a physical epoch offset."""
    from ocen_dm.selection.hst_gaia_match import load_match
    t = load_match()
    assert len(t) > 5000
    sep = np.asarray(t["separation_arcsec"])
    assert np.median(sep) < 0.02, "separations near the tolerance mean neighbouring stars"
    off = np.hypot(*t.meta["epoch_offset_arcsec"])
    # 7.5 mas/yr systemic PM over the ~13 yr between the HST and Gaia epochs
    assert 0.05 < off < 0.20


def test_unflagged_gaia_carries_undeclared_scatter():
    from ocen_dm.plotting.constraints import hst_gaia_excess_table
    ex = hst_gaia_excess_table()
    assert len(ex) >= 3
    # the excess is real, larger than the quoted errors, and falls outwards (crowding)
    assert np.all(ex["excess_noise"] > 1.5 * ex["quoted"])
    assert ex["excess_noise"][0] > ex["excess_noise"][-1]
    # and Gaia therefore reports a larger dispersion than HST from the identical stars
    assert np.all(ex["ratio"] > 1.1)


def test_quality_flag_reconciles_the_two_instruments():
    """In the only annulus where both catalogues are usable, they agree."""
    from ocen_dm.plotting.constraints import hst_gaia_overlap_table
    ov = hst_gaia_overlap_table()
    row = ov[ov["sample"] == "Gaia, quality flag"]
    assert len(row) == 1
    ratio, err = float(row["ratio"][0]), float(row["ratio_err"][0])
    assert abs(ratio - 1.0) < 2.0 * err, "flagged Gaia should agree with HST"
    unflagged = float(ov[ov["sample"] == "Gaia, no quality cut"]["ratio"][0])
    assert unflagged > 1.2, "and the unflagged sample should not"


def test_gaia_has_no_quality_stars_in_the_hst_region():
    """Which is why the published EDR3 profile inside ~400 arcsec is extrapolation."""
    from ocen_dm.kinematics.outer_profile import load_members
    s = load_members()
    q = s.select((s.quality_flag & 2) > 0)
    assert (q.r_arcsec < 200).sum() == 0
    assert ((q.r_arcsec >= 200) & (q.r_arcsec < 300)).sum() < 20
    assert ((q.r_arcsec >= 460) & (q.r_arcsec < 700)).sum() > 5000
