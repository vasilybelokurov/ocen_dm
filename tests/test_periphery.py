"""Beyond the Gaia edge: the PM answer depends on the selection window, the LOS answer does not."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = ((processed_dir() / "tails" / "kuzma2025_periphery.ecsv").exists()
        and (processed_dir() / "tails" / "kuzma2026_spectroscopy.ecsv").exists())
pytestmark = pytest.mark.skipif(not _HAS, reason="periphery catalogues not present")


def test_gaia_catalogue_is_truncated_not_bounded():
    """The member catalogue stops well inside the Jacobi radius."""
    from astropy.table import Table
    from ocen_dm.kinematics.periphery import (GAIA_EDGE_DEG, R_JACOBI_PERI_PC, radius_deg)
    t = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    r = radius_deg(t["ra"], t["dec"])
    assert r.max() < 0.68 and abs(r.max() - GAIA_EDGE_DEG) < 0.01
    edge_pc = GAIA_EDGE_DEG * 3600 * 5.43e3 / 206264.806
    assert edge_pc < R_JACOBI_PERI_PC, "the data stop before the cluster can"
    assert edge_pc / R_JACOBI_PERI_PC < 0.8


def test_periphery_pm_dispersion_tracks_the_selection_window():
    """Which is why it must never be quoted as a single number."""
    from ocen_dm.kinematics.periphery import periphery_pm_profile
    narrow, wide = periphery_pm_profile(0.8), periphery_pm_profile(1.2)
    outer = np.asarray(narrow["r_pc"]) > 100
    assert outer.sum() >= 2
    ratio = np.asarray(wide["sigma_kms"])[outer] / np.asarray(narrow["sigma_kms"])[outer]
    assert np.all(ratio > 1.3), "widening the window inflates the apparent dispersion"
    assert narrow.meta["window_masyr"] == 0.8


def test_line_of_sight_dispersion_does_not_keep_falling():
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    from ocen_dm.kinematics.periphery import KMS_PER_MASYR_KPC, periphery_los_profile
    los = periphery_los_profile()
    inner = los[los["r_pc"] < 80][0]
    assert inner["n_stars"] > 50
    # the outermost Gaia point, in km/s, is the value the profile would keep falling from
    g = load_edr3_profile()
    last = float(g["sigma_pm"][-1]) * KMS_PER_MASYR_KPC * 5.43
    assert 4.0 < last < 8.0
    assert abs(float(inner["sigma_kms"]) - last) < 3 * float(inner["sigma_kms_err"]) + 1.5


def test_plot_renders(tmp_path):
    from ocen_dm.plotting.constraints import plot_periphery
    out = plot_periphery(tmp_path / "p.png")
    assert out.exists() and out.stat().st_size > 40_000
