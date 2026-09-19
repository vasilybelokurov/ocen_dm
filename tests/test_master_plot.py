"""The master figure must draw exactly the datasets the likelihood is given."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = (processed_dir() / "kinematics" / "ocen_pm_dispersion_edr3_ours.ecsv").exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="products not built")


def test_every_default_dataset_has_a_drawing_style():
    from ocen_dm.cli import DEFAULT_DATASETS
    from ocen_dm.plotting.style import DATASET_KEY_MAP, dataset_style
    for key in DEFAULT_DATASETS.split(","):
        assert key in DATASET_KEY_MAP, f"{key} would be drawn with a fallback style"
        dataset_style(key)


def test_point_count_matches_the_likelihood():
    from ocen_dm.cli import DEFAULT_DATASETS
    from ocen_dm.kinematics.likelihood import KinematicData, load_profile
    keys = DEFAULT_DATASETS.split(",")
    data = KinematicData.load(keys)
    assert data.n_points == sum(len(load_profile(k).r) for k in keys)
    assert data.n_points == 127


def test_master_plot_renders(tmp_path):
    from ocen_dm.plotting.constraints import plot_master_datasets
    out = plot_master_datasets(tmp_path / "m.png")
    assert out.exists() and out.stat().st_size > 80_000


def test_master_plot_follows_a_custom_dataset_list(tmp_path):
    """Passing the old published EDR3 key must still work, for comparison figures."""
    from ocen_dm.plotting.constraints import plot_master_datasets
    out = plot_master_datasets(tmp_path / "m2.png",
                               datasets="hst_pm_combined,muse_los_dispersion,gaia_edr3_pm")
    assert out.exists()


def test_master_plot_shows_the_components_the_fit_does_not_use():
    """Our Gaia measurement is fitted as a combination but measured as R and T."""
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    from ocen_dm.kinematics.likelihood import load_profile
    g = load_edr3_profile()
    assert {"sigma_pmr", "sigma_pmt"} <= set(g.colnames)
    combined = np.sqrt(0.5 * (np.asarray(g["sigma_pmr"]) ** 2 + np.asarray(g["sigma_pmt"]) ** 2))
    assert np.allclose(combined, np.asarray(g["sigma_pm"]), rtol=1e-6)
    assert load_profile("gaia_edr3_ours").kind == "pmc", "still fitted as the combination"
    # the anisotropy is a real signal, not noise: radial inside, tangential outside
    ratio = np.asarray(g["sigma_pmt"]) / np.asarray(g["sigma_pmr"])
    assert ratio[:4].mean() < 0.95 and ratio[-3:].mean() > 1.05


def test_usable_edr3_exists_between_300_and_380_arcsec():
    """The inner edge is set by where stars exist, and the innermost bin uses them."""
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    g = load_edr3_profile()
    first = g[0]
    assert 300 <= first["r_lower"] < 380
    assert first["n_stars"] >= 40, "these are real stars, not an extrapolation"
    assert first["sigma_pm_err"] / first["sigma_pm"] < 0.12
