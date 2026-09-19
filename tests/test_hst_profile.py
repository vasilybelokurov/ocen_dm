"""HST measured in Gaia's annuli: the overlap the published profile hid."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import raw_dir

_HAS = (raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits").exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="oMEGACat catalogue not present")


def test_catalogue_reaches_far_beyond_its_published_profile():
    from ocen_dm.kinematics.likelihood import load_profile
    from ocen_dm.kinematics.hst_profile import HST_R_MAX_ARCSEC, load_hst_sample
    s = load_hst_sample()
    assert s.r_arcsec.max() > 460, "the stars go well past 300 arcsec"
    assert load_profile("hst_pm_combined").r.max() < 320, "but the published profile stops"
    assert HST_R_MAX_ARCSEC >= 460


def test_quality_flag_ends_long_before_the_data_do():
    from ocen_dm.kinematics.hst_profile import load_hst_sample
    s = load_hst_sample()
    flagged = s.quality_flag > 0
    assert ((s.r_arcsec >= 300) & (s.r_arcsec < 340) & flagged).sum() > 10000
    assert ((s.r_arcsec >= 380) & flagged).sum() == 0
    # yet thousands of unflagged stars remain out there
    assert ((s.r_arcsec >= 380) & (s.r_arcsec < 466)).sum() > 5000


def test_mixture_tames_the_unflagged_stars():
    """Raw, the unflagged stars give a dispersion that rises outwards; modelled, it falls."""
    from ocen_dm.kinematics.hst_profile import hst_profile
    from ocen_dm.kinematics.outer_profile import dispersion_ml, load_members  # noqa: F401
    from ocen_dm.kinematics.hst_profile import load_hst_sample
    s = load_hst_sample()
    raw = []
    for lo, hi in ((340, 380), (380, 420), (420, 466)):
        m = (s.r_arcsec >= lo) & (s.r_arcsec < hi) & (s.quality_flag == 0)
        x = dispersion_ml(s.mu_r[m], s.err_r[m])[0]
        y = dispersion_ml(s.mu_t[m], s.err_t[m])[0]
        raw.append(np.sqrt(0.5 * (x ** 2 + y ** 2)))
    assert raw[-1] > raw[0], "raw, it rises outwards, which no cluster does"
    t = hst_profile(edges_arcsec=(340., 380., 420., 466.))
    assert len(t) == 3
    assert t["sigma_pm"][0] > t["sigma_pm"][1], "modelled, it falls"
    assert np.all(np.asarray(t["f_field"]) < 0.06), "the field is a few per cent, not the signal"


def test_overlap_exists_and_the_instruments_agree():
    from ocen_dm.plotting.constraints import hst_gaia_overlap_profile
    ov = hst_gaia_overlap_profile()
    assert len(ov) == 2, "300-380 and 380-460 arcsec"
    assert np.all(np.asarray(ov["n_hst"]) > 5000) and np.all(np.asarray(ov["n_gaia"]) > 40)
    # the two samples sit at different radii inside the bin, hence the correction
    assert np.all(np.asarray(ov["r_gaia"]) > np.asarray(ov["r_hst"]))
    assert np.all(np.abs(np.asarray(ov["ratio"]) - np.asarray(ov["ratio_raw"])) > 0.005)
    r, e = np.asarray(ov["ratio"]), np.asarray(ov["ratio_err"])
    w = np.sum(r / e ** 2) / np.sum(1 / e ** 2)
    we = 1 / np.sqrt(np.sum(1 / e ** 2))
    assert abs(w - 1.0) < 2 * we, "the instruments agree once HST's unflagged stars are fixed"
    assert 0.92 < w < 1.08


def test_plot_renders(tmp_path):
    from ocen_dm.plotting.constraints import plot_hst_gaia_overlap
    out = plot_hst_gaia_overlap(tmp_path / "o.png")
    assert out.exists() and out.stat().st_size > 60_000


def test_unflagged_stars_are_biased_and_the_correction_is_applied():
    """The bias the published oMEGACat profile exposed: unflagged HST errors are too small."""
    from ocen_dm.kinematics.hst_profile import UNFLAGGED_BIAS, hst_profile, unflagged_bias
    ratio, err = unflagged_bias()
    assert 1.03 < ratio < 1.13 and err < 0.03
    assert abs(ratio - UNFLAGGED_BIAS[0]) < 4 * UNFLAGGED_BIAS[1], "stored value must track"
    edges = (300., 380., 460.)
    on = hst_profile(edges_arcsec=edges, correct_unflagged=True)
    off = hst_profile(edges_arcsec=edges, correct_unflagged=False)
    assert np.all(np.asarray(on["sigma_pm"]) < np.asarray(off["sigma_pm"]))
    # outside 340 arcsec every star is unflagged, so the full correction applies there
    assert on["f_unflagged"][-1] > 0.99
    assert np.isclose(off["sigma_pm"][-1] / on["sigma_pm"][-1], UNFLAGGED_BIAS[0], rtol=0.01)
    assert np.all(np.asarray(on["sigma_sys"]) > 0)


def test_flagged_measurement_reproduces_the_published_profile():
    """Validation of the whole chain against oMEGACat's own numbers."""
    from astropy.table import Table
    from ocen_dm.kinematics.hst_profile import hst_profile
    from ocen_dm.paths import processed_dir
    pub = Table.read(processed_dir() / "kinematics" / "omegacat_vi_pm_combined.ecsv")
    rp = np.asarray(pub["r_median"], float); sp = np.asarray(pub["sigma_pmc"], float)
    t = hst_profile(edges_arcsec=(150., 200., 250., 300.), require_flag=True,
                    correct_unflagged=False)
    ratio = np.asarray(t["sigma_pm"]) / np.interp(np.asarray(t["r_median"]), rp, sp)
    assert np.all(np.abs(ratio - 1.0) < 0.035), "our method must reproduce theirs to ~3 per cent"
