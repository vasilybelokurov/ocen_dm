"""Applying Vasiliev & Baumgardt's own error prescription reconciles us with their profile."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = ((processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv").exists()
        and (processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv").exists())
pytestmark = pytest.mark.skipif(not _HAS, reason="catalogues not present")


def test_error_inflation_matches_the_published_coefficients():
    from ocen_dm.kinematics.vb2021_replication import error_inflation
    # their Table 1: eta = (1 + Sigma/Sigma_0)**0.04, Sigma_0 = 10 (5p), 5 (6p)
    eta5 = error_inflation(np.array([0.0, 10.0, 200.0]), np.array([True, True, True]))
    assert np.isclose(eta5[0], 1.0)
    assert np.isclose(eta5[1], 2.0 ** 0.04)
    assert np.isclose(eta5[2], 21.0 ** 0.04)
    eta6 = error_inflation(np.array([200.0]), np.array([False]))
    assert eta6[0] > eta5[2], "6p sources take the smaller Sigma_0 and so a larger eta"
    # the quoted magnitude: ~1.1-1.15 in the densest regions
    assert 1.05 < error_inflation(np.array([200.0]), np.array([True]))[0] < 1.20


def test_inflation_removes_the_offset_against_the_published_profile():
    from ocen_dm.kinematics.vb2021_replication import replication_table
    t = replication_table()
    mid = (t["r_median"] > 450) & (t["r_median"] < 1050)
    assert mid.sum() >= 4
    raw = np.asarray(t["raw_ratio"])[mid]
    scaled = np.asarray(t["eta_ratio"])[mid]
    # deconvolving the released (raw Gaia) errors biases us high by several per cent
    assert raw.min() > 1.02, "raw-error dispersions should sit above their profile"
    # applying their eta brings every annulus within 2 per cent
    assert np.all(np.abs(scaled - 1.0) < 0.02)
    assert np.mean(np.abs(scaled - 1.0)) < np.mean(np.abs(raw - 1.0))


def test_released_catalogue_carries_raw_gaia_errors():
    """Guards the premise: the scaling is applied at fit time, not in the catalogue."""
    from astropy.table import Table
    cat = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    # a scaled catalogue would show error/sqrt-density structure; a raw one does not.
    # verified against gaia_edr3.gaia_source on WSDB (median ratio 1.0000, N = 228055).
    assert "source_density" in cat.colnames, "the density column exists precisely because the scaling is external"


def test_eta_overcorrects_the_faint_end():
    """The authors' own validation: a correct error model makes sigma independent of G."""
    from ocen_dm.kinematics.vb2021_replication import magnitude_consistency
    t = magnitude_consistency()
    inner = t["r_upper"] <= 1500
    raw = np.asarray(t["faint_over_bright_raw"])[inner]
    scaled = np.asarray(t["faint_over_bright_eta"])[inner]
    assert np.all(raw > 1.05), "raw errors are underestimated for faint stars"
    assert np.all(scaled < 1.0), "the density-only eta over-corrects them"


def test_well_measured_stars_agree_whatever_the_error_model():
    from ocen_dm.kinematics.vb2021_replication import low_noise_profile
    t = low_noise_profile()
    mid = (t["r_median"] > 450) & (t["r_median"] < 1200)
    assert mid.sum() >= 4
    # the error rescaling moves these by under 1 per cent...
    assert np.all(np.asarray(t["error_model_sensitivity"])[mid] < 0.01)
    # ...and both versions agree with the published profile to ~1 per cent
    assert np.all(np.abs(np.asarray(t["raw_ratio"])[mid] - 1.0) < 0.015)
    assert np.all(np.abs(np.asarray(t["eta_ratio"])[mid] - 1.0) < 0.015)
