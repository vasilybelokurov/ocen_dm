"""The master figure must draw exactly the datasets the likelihood is given."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = (processed_dir() / "kinematics" / "ocen_pm_dispersion_edr3_ours.ecsv").exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="products not built")


def test_every_default_dataset_has_a_drawing_style():
    from ocen_dm.cli import DEFAULT_DATASETS
    from ocen_dm.plotting.constraints import _MASTER_STYLE
    for key in DEFAULT_DATASETS.split(","):
        assert key in _MASTER_STYLE, f"{key} would be drawn with a fallback style"


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
