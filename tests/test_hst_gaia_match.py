"""HST and Gaia proper motions of the same stars: they must agree, and inside 460 arcsec they do not."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir, raw_dir

_HAS = ((raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits").exists()
        and (processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv").exists())


@pytest.mark.skipif(not _HAS, reason="catalogues not present")
def test_match_product_and_disagreement():
    from ocen_dm.kinematics.outer_profile import dispersion_ml
    from ocen_dm.selection.hst_gaia_match import load_match
    t = load_match()
    assert len(t) > 5000 and t.meta["tolerance_arcsec"] == 0.3
    r = np.asarray(t["r_arcsec"]); g = np.asarray(t["g_mag"])
    sel = (np.asarray(t["membership_prob"]) > 0.9) & (np.asarray(t["hst_quality"]) > 0) & (g < 18)
    assert sel.sum() > 500
    # per-star differences far exceed what the quoted errors allow
    da = np.asarray(t["gaia_pmra"])[sel] - np.asarray(t["hst_pmra"])[sel]
    da = da - np.median(da)
    quoted = np.hypot(np.median(np.asarray(t["hst_pmra_error"])[sel]), np.median(np.asarray(t["gaia_pmra_error"])[sel]))
    assert np.std(da) > 3 * quoted
    # and Gaia reports a larger dispersion than HST from the identical stars
    m = sel & (r > 150) & (r < 300)
    sh = dispersion_ml(np.asarray(t["hst_pmra"])[m], np.asarray(t["hst_pmra_error"])[m])[0]
    sg = dispersion_ml(np.asarray(t["gaia_pmra"])[m], np.asarray(t["gaia_pmra_error"])[m])[0]
    assert sg > 1.10 * sh


@pytest.mark.skipif(not _HAS, reason="catalogues not present")
def test_gaia_has_no_quality_stars_in_the_hst_region():
    """Which is why the published EDR3 profile inside ~400 arcsec is extrapolation."""
    from ocen_dm.kinematics.outer_profile import load_members
    s = load_members()
    q = s.select((s.quality_flag & 2) > 0)
    assert ((q.r_arcsec < 200)).sum() == 0
    assert ((q.r_arcsec >= 200) & (q.r_arcsec < 300)).sum() < 20
    assert ((q.r_arcsec >= 460) & (q.r_arcsec < 700)).sum() > 5000
