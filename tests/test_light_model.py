"""Tests for the MGE light model: recovery, half-light radius, deprojection identity."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.light_model import (
    MGEFit,
    SurfaceBrightnessProfile,
    arcsec_to_pc,
    build_stellar_mge,
    fit_mge_projected,
    projected_half_light_radius,
)


def _synthetic_profile(sigmas, fractions, n=60, noise_mag=0.0, rng=np.random.default_rng(3),
                       depth_mag=12.0):
    """A noiseless MGE profile, cut where it has faded by ``depth_mag`` -- the
    dynamic range of the real Trager data. Extending a Gaussian to 30 mag below
    its centre is unphysical and only tests floating-point conditioning."""
    fit = MGEFit(np.asarray(sigmas, float), np.asarray(fractions, float), 0.0, 0.0, 0)
    r = np.geomspace(5.0, 5000.0, 400)
    with np.errstate(divide="ignore"):
        mu = 20.0 - 2.5 * np.log10(fit.surface_intensity(r) / fit.surface_intensity(r[0])[0])
    keep = mu - mu.min() <= depth_mag
    r, mu = r[keep][:: max(1, keep.sum() // n)], mu[keep][:: max(1, keep.sum() // n)]
    mu = mu + rng.normal(0, noise_mag, len(mu)) if noise_mag else mu
    return SurfaceBrightnessProfile(r, mu, np.ones(len(mu)), "synthetic")


def test_arcsec_to_pc():
    assert arcsec_to_pc(206264.806, 1.0) == pytest.approx(1000.0, rel=1e-9)   # 1 rad at 1 kpc
    assert arcsec_to_pc(1.0, 5.43) == pytest.approx(0.02633, rel=1e-3)


def test_single_gaussian_is_recovered_to_the_centi_magnitude():
    """Fixed-grid NNLS: recovery is limited by grid spacing (see the docstring)."""
    prof = _synthetic_profile([100.0], [1.0])
    fit = fit_mge_projected(prof)                 # default grid
    assert fit.rms_mag < 0.03
    fit128 = fit_mge_projected(prof, n_grid=128)
    assert fit128.rms_mag < 0.01                  # and improves with density
    # the light must sit in components close to the true width
    weighted_sigma = np.exp(np.sum(fit.fractions * np.log(fit.sigmas_arcsec)))
    assert weighted_sigma == pytest.approx(100.0, rel=0.15)


def test_two_gaussian_mixture_reprojects_to_the_data():
    prof = _synthetic_profile([40.0, 400.0], [0.6, 0.4])
    fit = fit_mge_projected(prof)
    assert fit.rms_mag < 0.02
    assert fit.max_abs_resid_mag < 0.1


def test_fractions_sum_to_one_and_are_non_negative():
    fit = fit_mge_projected(_synthetic_profile([40.0, 400.0], [0.6, 0.4]))
    assert fit.fractions.sum() == pytest.approx(1.0)
    assert np.all(fit.fractions >= 0)


def test_half_light_radius_of_one_gaussian_is_1p1774_sigma():
    fit = MGEFit(np.array([100.0]), np.array([1.0]), 0.0, 0.0, 1)
    assert projected_half_light_radius(fit) == pytest.approx(100.0 * np.sqrt(2 * np.log(2)), rel=1e-6)


def test_deprojection_keeps_widths_and_total_mass():
    fit = MGEFit(np.array([50.0, 500.0]), np.array([0.3, 0.7]), 0.0, 0.0, 2)
    mge = build_stellar_mge(fit, distance_kpc=5.43, total_mass=2.0e6)
    assert mge.total_mass == pytest.approx(2.0e6)
    assert np.allclose(mge._sigmas, arcsec_to_pc(fit.sigmas_arcsec, 5.43))
    assert mge.verify().passed


def test_projected_half_light_exceeds_nothing_unphysical():
    """3D half-mass radius must exceed the projected half-light radius (spherical)."""
    fit = MGEFit(np.array([50.0, 500.0]), np.array([0.3, 0.7]), 0.0, 0.0, 2)
    D = 5.43
    mge = build_stellar_mge(fit, D, total_mass=1.0)
    r_half_3d = mge.radius_enclosing(0.5)
    assert r_half_3d > arcsec_to_pc(projected_half_light_radius(fit), D)


def test_weights_downweight_bad_points():
    """A low-weight outlier must not pull the fit off the well-weighted points.

    Weights act as inverse variances, so a +3 mag point at weight 0.01 still
    carries chi^2 = 0.09 and may bend a noiseless fit by ~0.1 mag; the property
    that matters is that the fit through the unit-weight points stays good."""
    prof = _synthetic_profile([100.0], [1.0])
    clean = fit_mge_projected(prof)
    mu = prof.mu.copy(); mu[-1] += 3.0            # one wildly wrong outer point
    w = prof.weight.copy(); w[-1] = 0.01
    bad = SurfaceBrightnessProfile(prof.r_arcsec, mu, w, "synthetic")
    fit = fit_mge_projected(bad)
    heavy = w >= 0.5
    model = -2.5 * np.log10(fit.surface_intensity(prof.r_arcsec))
    zp = np.average(mu[heavy] - model[heavy])
    resid_heavy = mu[heavy] - (zp + model[heavy])
    assert np.sqrt(np.mean(resid_heavy**2)) < 0.15
    assert fit.rms_mag < 0.2
    # and with the outlier weight at the floor the fit is essentially the clean one
    w[-1] = 1e-3
    floor = fit_mge_projected(SurfaceBrightnessProfile(prof.r_arcsec, mu, w, "synthetic"))
    assert abs(projected_half_light_radius(floor) - projected_half_light_radius(clean)) < 2.0


def test_real_trager_profile_fits_if_present():
    from ocen_dm.light_model import load_trager_profile
    from ocen_dm.paths import processed_dir

    path = processed_dir() / "literature" / "trager1995_ocen_sbp.ecsv"
    if not path.is_file():
        pytest.skip("real Trager product not built")
    fit = fit_mge_projected(load_trager_profile(path))
    assert fit.rms_mag < 0.3
    assert 100.0 < projected_half_light_radius(fit) < 600.0     # arcsec; Harris r_h = 300"


# ------------------------------------------------------- star counts / composite ---
def _fake_profile(r, mu, w=None, source="x"):
    from ocen_dm.light_model import SurfaceBrightnessProfile
    r = np.asarray(r, float); mu = np.asarray(mu, float)
    return SurfaceBrightnessProfile(r, mu, np.ones_like(r) if w is None else np.asarray(w, float), source)


def test_composite_profile_matches_zero_point_and_splices():
    from ocen_dm.light_model import composite_profile
    r_out = np.geomspace(10, 1000, 30); mu_out = 15 + 2.5 * np.log10(1 + (r_out / 100) ** 2)   # Plummer-like
    r_in = np.geomspace(1, 80, 25); mu_in = 2.5 * np.log10(1 + (r_in / 100) ** 2) + 7.0          # same shape, shifted 8 mag
    comp = composite_profile(_fake_profile(r_out, mu_out), _fake_profile(r_in, mu_in), r_switch_arcsec=25.0,
                             anchor_arcsec=(30.0, 80.0))
    assert np.all(np.diff(comp.r_arcsec) > 0)
    assert comp.r_arcsec.min() == pytest.approx(1.0) and comp.r_arcsec.max() == pytest.approx(1000.0)
    # after the zero-point shift the inner points fall on the outer curve
    inner = comp.r_arcsec < 25
    expected = 15 + 2.5 * np.log10(1 + (comp.r_arcsec[inner] / 100) ** 2)
    np.testing.assert_allclose(comp.mu[inner], expected, atol=3e-3)   # limited by log-linear interpolation of the anchor
    assert not np.any((comp.r_arcsec >= 25) & np.isin(comp.r_arcsec, r_in))    # no inner points past the switch
    with pytest.raises(ValueError, match="anchor"):
        composite_profile(_fake_profile(r_out, mu_out), _fake_profile(r_in, mu_in), anchor_arcsec=(2000.0, 3000.0))
    # low-weight outliers in the outer profile must not move the zero-point
    mu_bad = mu_out.copy(); w_out = np.ones_like(r_out)
    bad = (r_out > 30) & (r_out < 80) & (np.arange(len(r_out)) % 2 == 0)
    mu_bad[bad] -= 0.8; w_out[bad] = 0.03
    comp2 = composite_profile(_fake_profile(r_out, mu_bad, w_out), _fake_profile(r_in, mu_in), r_switch_arcsec=25.0,
                              anchor_arcsec=(30.0, 80.0))
    np.testing.assert_allclose(comp2.mu[comp2.r_arcsec < 25], expected, atol=3e-3)


from ocen_dm.paths import raw_dir

_HAS_CAT = (raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits").exists()


@pytest.mark.skipif(not _HAS_CAT, reason="oMEGACat catalogue not present")
def test_star_counts_below_trager_light_in_the_core():
    """The documented finding: number density is 0.2-0.3 mag below the V-band light inside 20 arcsec."""
    from ocen_dm.light_model import load_tracer_profile, star_count_profile, fit_mge_projected
    counts = star_count_profile(mag_cut=19.0)
    assert counts.r_arcsec.min() < 3.0 and counts.r_arcsec.max() > 100.0
    comp = load_tracer_profile("composite")
    trager = fit_mge_projected(load_tracer_profile("trager"), sigma_range_arcsec=(10.5, 3000.0))
    # zero-point: weighted mean offset between the Trager MGE and the composite over 30-100"
    anchor = (comp.r_arcsec >= 30.0) & (comp.r_arcsec <= 100.0)
    mu_ref = np.average(comp.mu[anchor] + 2.5 * np.log10(trager.surface_intensity(comp.r_arcsec[anchor])),
                        weights=comp.weight[anchor])
    inner = (comp.r_arcsec > 3.0) & (comp.r_arcsec < 15.0)          # where N >= 15 per bin and the deficit lives
    deficit = comp.mu[inner] - trager.surface_brightness(comp.r_arcsec[inner], mu_ref)
    assert 0.1 < np.average(deficit, weights=comp.weight[inner]) < 0.5   # counts fainter (positive) than the light


def test_load_tracer_profile_rejects_unknown_kind():
    from ocen_dm.light_model import load_tracer_profile
    with pytest.raises(ValueError):
        load_tracer_profile("noyola")
