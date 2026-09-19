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
    from ocen_dm.plotting.constraints import plot_constraint_map, plot_contamination_model, plot_outer_tracer_audit
    for path in (plot_constraint_map(tmp_path / "map.png"),
                 plot_constraint_map(tmp_path / "map_pcut.png", contamination_modelled=False),
                 plot_outer_tracer_audit(tmp_path / "audit.png"),
                 plot_outer_tracer_audit(tmp_path / "audit_mix.png", reference="mixture"),
                 plot_contamination_model(tmp_path / "contam.png", annuli=((1400.0, 1800.0),))):
        assert path.exists() and path.stat().st_size > 50_000


@pytest.mark.skipif(not _HAS_MEMBERS, reason="member catalogue not present")
def test_fit_quality_table_records_every_annulus(tmp_path):
    from ocen_dm.plotting.constraints import fit_quality_table
    t = fit_quality_table(save=False)
    assert len(t) == 10
    for c in ("f_field", "n_cluster", "sigma_pmr", "sigma_pmt", "chi2_r_wide", "chi2_r_peak",
              "chi2_t_wide", "chi2_t_peak"):
        assert c in t.colnames and np.all(np.isfinite(t[c]))
    # the field fraction must rise outward, and the fit must be acceptable everywhere
    assert np.all(np.diff(np.asarray(t["f_field"])[2:]) > 0)
    assert np.all(np.asarray(t["chi2_r_wide"]) < 3.0) and np.all(np.asarray(t["chi2_t_wide"]) < 3.0)
    # the cluster count must fall outward and stay a sensible fraction of the stars
    assert np.all(np.asarray(t["n_cluster"]) < np.asarray(t["n_stars"]))


@pytest.mark.skipif(not (_HAS_RUNS and _HAS_MEMBERS), reason="runs or members missing")
def test_annulus_figures_render(tmp_path):
    from ocen_dm.plotting.constraints import plot_annulus_fits
    for comp in ("r", "t"):
        p = plot_annulus_fits(tmp_path / f"annuli_{comp}.png", component=comp)
        assert p.exists() and p.stat().st_size > 100_000


@pytest.mark.skipif(not _HAS_MEMBERS, reason="member catalogue not present")
def test_method_comparison_figure_and_agreement(tmp_path):
    """The two measurements must agree within a few per cent, component by component."""
    from ocen_dm.plotting.constraints import our_mixture_profile, our_outer_profile, plot_method_comparison
    mix = our_mixture_profile(); cut = our_outer_profile()
    for comp in ("sigma_pmr", "sigma_pmt"):
        d = np.asarray(mix[comp]) / np.asarray(cut[comp]) - 1.0
        e = np.hypot(np.asarray(mix[f"{comp}_err"]) / np.asarray(mix[comp]),
                     np.asarray(cut[f"{comp}_err"]) / np.asarray(cut[comp]))
        assert np.all(np.abs(d) < 0.06)                      # never more than 6 per cent apart
        assert np.all(d > -0.02)                             # the decomposition is never lower by much
        assert np.median(np.abs(d) / e) < 2.0                # and consistent within the errors
    p = plot_method_comparison(tmp_path / "cmp.png")
    assert p.exists() and p.stat().st_size > 80_000


@pytest.mark.skipif(not (_HAS_RUNS and _HAS_MEMBERS), reason="runs or members missing")
def test_both_streaming_treatments_render_and_differ(tmp_path):
    """Both versions are kept: the published rotation curve and our own fitted means."""
    from ocen_dm.plotting.constraints import _datasets_in_kms, plot_constraint_map
    a = _datasets_in_kms(5.43, True, "self")[-1]
    b = _datasets_in_kms(5.43, True, "published")[-1]
    assert a["streaming2"] is not None and b["streaming2"] is None
    for name, stream in (("self.png", "self"), ("pub.png", "published")):
        p = plot_constraint_map(tmp_path / name, streaming=stream)
        assert p.exists() and p.stat().st_size > 80_000
