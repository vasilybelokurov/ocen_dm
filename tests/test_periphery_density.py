"""Are the periphery stars really there? A background-subtracted star-count test."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = (processed_dir() / "tails" / "kuzma2025_periphery.ecsv").exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="periphery catalogue not present")


def test_footprint_is_uniform_so_counts_need_no_correction():
    from ocen_dm.selection.periphery_density import load_periphery
    _, r, _ = load_periphery()
    dens = []
    for lo, hi in ((0.67, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.1)):
        dens.append(((r >= lo) & (r < hi)).sum() / (np.pi * (hi ** 2 - lo ** 2)))
    dens = np.array(dens)
    assert dens.std() / dens.mean() < 0.08, "coverage must be uniform for raw counts to work"


def test_overdensity_is_detected_where_the_spectroscopy_sits():
    from ocen_dm.selection.periphery_density import density_profile
    t = density_profile()
    inner = t[t["r_pc"] < 80]
    assert np.all(np.asarray(inner["significance"]) > 3.5), "57 and 76 pc must be solid"
    at100 = t[(t["r_pc"] > 90) & (t["r_pc"] < 110)][0]
    assert at100["significance"] > 2.0
    # and it truncates: beyond 150 pc the circular average is consistent with the field
    outer = t[t["r_pc"] > 150]
    assert np.all(np.abs(np.asarray(outer["significance"])) < 2.0)


def test_the_profile_truncates_steeply():
    from ocen_dm.selection.periphery_density import density_profile
    t = density_profile()
    r = np.asarray(t["r_pc"]); e = np.asarray(t["excess"])
    m = (r < 110) & (e > 0)
    slope = np.polyfit(np.log(r[m]), np.log(e[m]), 1)[0]
    assert slope < -4.0, "a truncating cluster, not an extended halo of stars"


def test_control_windows_are_a_real_background():
    from ocen_dm.selection.periphery_density import candidate_masks, load_periphery
    t, r, _ = load_periphery()
    sig, ctrl = candidate_masks(t)
    assert len(ctrl) == 8
    far = r > 3.6
    # far from the cluster the signal window should look like the controls
    n_sig = (sig & far).sum()
    n_ctrl = np.array([(c & far).sum() for c in ctrl], float)
    assert abs(n_sig - n_ctrl.mean()) < 3 * max(n_ctrl.std(ddof=1), 1.0)


def test_plot_renders(tmp_path):
    from ocen_dm.plotting.constraints import plot_periphery_density
    out = plot_periphery_density(tmp_path / "d.png")
    assert out.exists() and out.stat().st_size > 60_000
