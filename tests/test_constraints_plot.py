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


@pytest.mark.skipif(not (_HAS_RUNS and _HAS_MEMBERS), reason="runs or members missing")
def test_residual_significance_figure_and_flatness(tmp_path):
    """Between 460 and 1500 arcsec the residual is a constant offset, and inside 460 arcsec
    Gaia keeps too few stars for the measurement to mean anything."""
    import json
    from scipy.stats import chi2 as chi2_dist
    from ocen_dm.kinematics.report import _family_for
    from ocen_dm.plotting.constraints import (KMS_PER_MASYR_KPC, _sigma_1d_kms, our_mixture_profile,
                                              plot_residual_significance)
    from ocen_dm.kinematics.outer_profile import load_members
    from ocen_dm.paths import results_dir

    edges = np.array([300., 380, 460, 540, 630, 730, 850, 1000, 1200, 1500, 1900, 2400.])
    mix = our_mixture_profile(edges)
    summary = json.loads((results_dir() / "fits" / "K1_noDM_composite" / "summary.json").read_text())
    fam = _family_for(summary)
    x = np.array([summary["parameters"][n]["ml"] for n in fam.names])
    jeans, D, scales = fam.build(fam.to_dict(x))
    r = np.asarray(mix["r_median"])
    total = np.sqrt(0.5 * (np.asarray(mix["sigma_pmr"]) ** 2 + np.asarray(mix["mean_pmr"]) ** 2
                           + np.asarray(mix["sigma_pmt"]) ** 2 + np.asarray(mix["mean_pmt"]) ** 2))
    model = _sigma_1d_kms(jeans, D, r) / (KMS_PER_MASYR_KPC * D) * scales.get("GaiaEDR3", 1.0)
    res = 100 * (total / model - 1); err = 100 * np.asarray(mix["sigma_pm_err"]) / model
    good = (r > 460) & (r < 1500)
    w = 1 / err[good] ** 2
    c = np.sum(w * res[good]) / np.sum(w)
    chi2 = np.sum(w * (res[good] - c) ** 2)
    assert 1 - chi2_dist.cdf(chi2, good.sum() - 1) > 0.05      # flat: no structure at 13-36 pc
    assert 2.0 < c < 6.0                                        # a few per cent offset
    assert np.all(np.abs(res[r < 460] / err[r < 460]) < 2.5)    # the "wiggle" points are < 2.5 sigma
    assert res[-1] / err[-1] > 3.0                              # the outermost rise is real
    # and the reason the inner points are useless: almost nothing passes the quality flag
    s = load_members()
    q = s.select((s.quality_flag & 2) > 0)
    inner_all = ((s.r_arcsec >= 300) & (s.r_arcsec < 460)).sum()
    inner_q = ((q.r_arcsec >= 300) & (q.r_arcsec < 460)).sum()
    assert inner_q / inner_all < 0.05
    p = plot_residual_significance(tmp_path / "sig.png")
    assert p.exists() and p.stat().st_size > 60_000


@pytest.mark.skipif(not (_HAS_RUNS and _HAS_MEMBERS), reason="runs or members missing")
def test_hst_and_gaia_do_not_join_smoothly(tmp_path):
    """HST declines to -5 per cent of the model at its edge; Gaia sits at +10 to +12 per cent
    where it becomes usable. The step is the feature, not a bump inside the Gaia data."""
    import json
    from astropy.table import Table
    from ocen_dm.kinematics.report import _family_for
    from ocen_dm.paths import processed_dir, results_dir
    from ocen_dm.plotting.constraints import (KMS_PER_MASYR_KPC, _sigma_1d_kms, our_mixture_profile,
                                              plot_dataset_step)

    summary = json.loads((results_dir() / "fits" / "K1_noDM_composite" / "summary.json").read_text())
    fam = _family_for(summary)
    x = np.array([summary["parameters"][n]["ml"] for n in fam.names])
    jeans, D, _ = fam.build(fam.to_dict(x))
    hr = Table.read(processed_dir() / "kinematics" / "omegacat_vi_pm_radial.ecsv")
    ht = Table.read(processed_dir() / "kinematics" / "omegacat_vi_pm_tangential.ecsv")
    rh = np.asarray(hr["r_median"])
    hst = np.sqrt(0.5 * (np.asarray(hr["sigma_pmr"]) ** 2 + np.asarray(ht["sigma_pmt"]) ** 2))
    res_h = hst / (_sigma_1d_kms(jeans, D, rh) / (KMS_PER_MASYR_KPC * D)) - 1
    assert res_h[-1] < -0.03 and res_h[rh > 150][0] > -0.02      # HST declines from ~0 to below -3 %
    edges = np.array([460., 540, 630, 730, 850, 1000])
    mix = our_mixture_profile(edges)
    rg = np.asarray(mix["r_median"])
    tot = np.sqrt(0.5 * (np.asarray(mix["sigma_pmr"]) ** 2 + np.asarray(mix["mean_pmr"]) ** 2
                         + np.asarray(mix["sigma_pmt"]) ** 2 + np.asarray(mix["mean_pmt"]) ** 2))
    res_g = tot / (_sigma_1d_kms(jeans, D, rg) / (KMS_PER_MASYR_KPC * D)) - 1
    assert np.all(res_g > 0.05)                                   # Gaia sits well above the same model
    assert res_g.min() - res_h[-1] > 0.10                         # a step of more than 10 per cent
    # and the two barely overlap: HST ends at 346", Gaia is under 15 per cent usable until ~500"
    assert hr["r_upper"][-1] < 350
    p = plot_dataset_step(tmp_path / "step.png")
    assert p.exists() and p.stat().st_size > 60_000


@pytest.mark.skipif(not (_HAS_RUNS and _HAS_MEMBERS), reason="runs or members missing")
def test_model_anisotropy_disagrees_with_the_data_in_the_outskirts(tmp_path):
    """The K1 model is radial at every radius; the measurement turns tangential beyond ~1 kpc
    of arcsec. That mismatch is where the outer excess lives."""
    import json
    from ocen_dm.kinematics.report import _family_for
    from ocen_dm.paths import results_dir
    from ocen_dm.plotting.constraints import OUTER_EDGES, our_mixture_profile, plot_offset_explained

    summary = json.loads((results_dir() / "fits" / "K1_noDM_composite" / "summary.json").read_text())
    fam = _family_for(summary)
    x = np.array([summary["parameters"][n]["ml"] for n in fam.names])
    jeans, D, _ = fam.build(fam.to_dict(x))
    pc = D * 1e3 / 206264.806
    mix = our_mixture_profile(OUTER_EDGES)
    r = np.asarray(mix["r_median"])
    model = np.array([jeans.dispersions_kms(np.array([rr * pc]))["pmt"][0]
                      / jeans.dispersions_kms(np.array([rr * pc]))["pmr"][0] for rr in r])
    meas = np.asarray(mix["sigma_pmt"]) / np.asarray(mix["sigma_pmr"])
    inner = (r > 400) & (r < 800); outer = r > 1300
    assert np.all(np.abs(model[inner] - meas[inner]) < 0.06)     # they agree where HST anchors the fit
    assert np.all(model[outer] < 0.92) and np.all(meas[outer] > 1.0)
    assert (meas[outer] - model[outer]).min() > 0.10             # and disagree by > 0.1 in the outskirts
    p = plot_offset_explained(tmp_path / "offset.png")
    assert p.exists() and p.stat().st_size > 60_000
