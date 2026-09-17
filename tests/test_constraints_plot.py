"""The constraint-map / tracer-audit figures and the measurement behind them."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir, results_dir

_HAS_RUNS = (results_dir() / "fits" / "K1_noDM_composite" / "summary.json").exists()
_HAS_MEMBERS = (processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv").exists()


@pytest.mark.skipif(not _HAS_MEMBERS, reason="member catalogue not present")
def test_our_outer_profile_reproduces_the_published_one():
    """With the authors' quality flag our independent measurement matches their profile."""
    from astropy.table import Table
    from ocen_dm.plotting.constraints import our_outer_profile
    ours = our_outer_profile()
    pub = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    published = np.interp(ours["r_median"], pub["r"], pub["sigma_pm"])
    outer = np.asarray(ours["r_median"]) > 500.0
    frac = np.abs(np.asarray(ours["sigma_pm"])[outer] / published[outer] - 1)
    assert np.all(frac < 0.10) and np.median(frac) < 0.05
    # and the quality flag matters: without it the crowded inner annuli are inflated
    raw = our_outer_profile(quality_mask=None)
    inner = np.asarray(raw["r_median"]) < 600.0
    assert np.all(np.asarray(raw["sigma_pm"])[inner] > 1.15 * np.asarray(ours["sigma_pm"])[inner])


@pytest.mark.skipif(not _HAS_MEMBERS, reason="member catalogue not present")
def test_bright_subsample_resolves_the_outer_dispersion():
    """G < 18.5 members have per-star errors well below the dispersion they measure."""
    from ocen_dm.plotting.constraints import our_outer_profile
    bright = our_outer_profile(g_range=(0.0, 18.5))
    outer = np.asarray(bright["r_median"]) > 1000.0
    assert np.all(np.asarray(bright["median_err"])[outer] < 0.6 * np.asarray(bright["sigma_pm"])[outer])
    # and they agree with the full sample, whose errors exceed the signal
    full = our_outer_profile()
    assert np.all(np.abs(np.asarray(bright["sigma_pm"])[outer] / np.asarray(full["sigma_pm"])[outer] - 1) < 0.10)


@pytest.mark.skipif(not (_HAS_RUNS and _HAS_MEMBERS), reason="runs or members missing")
def test_figures_render(tmp_path):
    from ocen_dm.plotting.constraints import plot_constraint_map, plot_outer_tracer_audit
    a = plot_constraint_map(tmp_path / "map.png")
    b = plot_outer_tracer_audit(tmp_path / "audit.png")
    assert a.exists() and a.stat().st_size > 50_000
    assert b.exists() and b.stat().st_size > 50_000
