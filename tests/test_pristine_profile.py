"""Pristine stars through our own cluster-plus-field mixture, with no proper-motion cut."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = ((processed_dir() / "tails" / "kuzma2025_periphery.ecsv").exists()
        and (processed_dir() / "tails" / "kuzma2025_periphery_gaia_covariance.ecsv").exists())
pytestmark = pytest.mark.skipif(not _HAS, reason="periphery catalogues not present")


def test_sample_joins_the_covariance_by_source_id():
    """The two products are NOT in the same row order; a positional join would be wrong."""
    from astropy.table import Table
    from ocen_dm.kinematics.pristine_profile import load_pristine_sample
    a = Table.read(processed_dir() / "tails" / "kuzma2025_periphery.ecsv")
    b = Table.read(processed_dir() / "tails" / "kuzma2025_periphery_gaia_covariance.ecsv")
    assert not np.array_equal(np.asarray(a["source_id"], np.int64),
                              np.asarray(b["source_id"], np.int64))
    s = load_pristine_sample()
    assert s.err_corr is not None and np.all(np.abs(s.err_corr) <= 1.0)
    assert np.median(np.abs(s.err_corr)) > 0.01, "correlations must actually arrive"


def test_agrees_with_gaia_on_identical_annuli():
    """Two catalogues, two selections, two field models: the answer must not depend on them."""
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    from ocen_dm.kinematics.pristine_profile import pristine_profile
    g = load_edr3_profile()
    edges = np.array([300., 380., 460., 582., 737., 934., 1182., 1497., 1895., 2413.])
    p = pristine_profile(edges_deg=tuple(edges / 3600.0), min_stars=25)
    assert len(p) >= 8
    gr = np.asarray(g["r_median"]) * 5.43e3 / 206264.806
    ratio, err = [], []
    for row in p:
        j = int(np.argmin(np.abs(gr - row["r_pc"])))
        gs, ge = float(g["sigma_pm"][j]), float(g["sigma_pm_err"][j])
        ps, pe = float(row["sigma_pm"]), float(row["sigma_pm_err"])
        ratio.append(ps / gs); err.append((ps / gs) * np.hypot(pe / ps, ge / gs))
    ratio, err = np.array(ratio), np.array(err)
    w = np.sum(ratio / err ** 2) / np.sum(1 / err ** 2)
    assert abs(w - 1.0) < 0.05, "the two catalogues must agree to a few per cent"
    assert np.all(np.abs(ratio - 1.0) < 0.12)


def test_mixture_beats_a_fixed_window_outside_the_cluster():
    """The window method inflates where the cluster fraction collapses; the mixture does not."""
    from ocen_dm.kinematics.periphery import periphery_pm_profile
    from ocen_dm.kinematics.pristine_profile import pristine_profile
    wide = periphery_pm_profile(1.2)
    mix = pristine_profile()
    outer_mix = mix[(mix["r_pc"] > 150) & (mix["r_pc"] < 260)]
    outer_win = wide[(wide["r_pc"] > 150) & (wide["r_pc"] < 260)]
    assert len(outer_mix) and len(outer_win)
    assert float(outer_mix["sigma_kms"][0]) < float(outer_win["sigma_kms"][0])


def test_unreliable_bins_are_flagged_not_hidden():
    from ocen_dm.kinematics.pristine_profile import pristine_profile
    t = pristine_profile()
    assert "reliable" in t.colnames
    good = np.asarray(t["reliable"], bool)
    assert good.sum() >= 3 and (~good).sum() >= 1
    assert np.all(np.asarray(t["n_cluster"])[good] >= 10)
    # the cluster fraction collapses beyond ~100 pc, which is why those bins are flagged
    assert np.all(np.asarray(t["f_field"])[~good] > 0.9)


def test_plot_renders(tmp_path):
    from ocen_dm.plotting.constraints import plot_extended_profile
    out = plot_extended_profile(tmp_path / "e.png")
    assert out.exists() and out.stat().st_size > 60_000
