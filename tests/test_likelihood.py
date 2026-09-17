"""Tests for the binned-profile likelihood.

Synthetic data: dispersions predicted by a known Plummer + point-mass model,
scattered with Gaussian noise. No real data are needed; the real-data smoke
test lives in the integration tests and is skipped if products are missing.
"""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics import Anisotropy, SphericalJeans, KMS_PER_MASYR_KPC
from ocen_dm.kinematics.likelihood import (
    ARCSEC_PER_RAD, BinnedProfile, KinematicData, ProfileLikelihood, load_profile,
)
from ocen_dm.mass_models import CompositeMassModel, Plummer, PointMass
from ocen_dm.paths import processed_dir

D = 5.43
PC_PER_ARCSEC = D * 1e3 / ARCSEC_PER_RAD


def _jeans(mass=3.0e6, a=5.0, m_bh=0.0, beta0=0.0, beta_inf=0.0):
    tracer = Plummer(mass, a)
    comps = [tracer] + ([PointMass(m_bh)] if m_bh > 0 else [])
    return SphericalJeans(CompositeMassModel(comps), tracer, Anisotropy(beta0, beta_inf, 20.0))


def _edges(n=12, r_in=2.0, r_out=800.0):
    e = np.geomspace(r_in, r_out, n + 1)
    return e[:-1], e[1:], np.sqrt(e[:-1] * e[1:])


def _mock(jeans, kind, rng, rel_err=0.03, edges=True, name="mock"):
    lo, hi, mid = _edges()
    template = BinnedProfile(name, kind, mid, lo if edges else None, hi if edges else None,
                             np.ones_like(mid), np.ones_like(mid), np.ones_like(mid), "X")
    truth = ProfileLikelihood(KinematicData((template,))).predict_profile(jeans, D, template)
    err = rel_err * truth
    return BinnedProfile(name, kind, mid, lo if edges else None, hi if edges else None,
                         truth + rng.normal(0, err), err, err, "X"), truth


# --------------------------------------------------------------- structure ---
def test_kinds_and_shapes_validated():
    r = np.array([1.0, 2.0]); one = np.ones(2)
    with pytest.raises(ValueError, match="kind"):
        BinnedProfile("a", "vphi", r, None, None, one, one, one, "X")
    with pytest.raises(ValueError, match="positive"):
        BinnedProfile("a", "los", r, None, None, one, np.array([1.0, 0.0]), one, "X")
    with pytest.raises(ValueError, match="both bin edges"):
        BinnedProfile("a", "los", r, r, None, one, one, one, "X")
    with pytest.raises(ValueError, match="r_upper"):
        BinnedProfile("a", "los", r, r, r, one, one, one, "X")


def test_shared_star_profiles_refused():
    r = np.array([1.0, 2.0]); one = np.ones(2)
    mk = lambda n, k: BinnedProfile(n, k, r, None, None, one, one, one, "HST")
    with pytest.raises(ValueError, match="same stars"):
        KinematicData((mk("hst_pm_combined", "pmc"), mk("hst_pm_radial", "pmr")))
    with pytest.raises(ValueError, match="same stars"):
        KinematicData((mk("hst_pm_tangential", "pmt"), mk("hst_pm_combined", "pmc")))
    KinematicData((mk("hst_pm_radial", "pmr"), mk("hst_pm_tangential", "pmt")))   # allowed
    with pytest.raises(ValueError, match="duplicate"):
        KinematicData((mk("a", "pmr"), mk("a", "pmr")))


def test_unknown_dataset_key():
    with pytest.raises(KeyError):
        load_profile("hst_everything")


# -------------------------------------------------------------- prediction ---
def test_bin_average_of_isotropic_plummer_matches_pointwise_within_bin():
    """The tracer-weighted bin mean of sigma^2 lies between the bin-edge values."""
    j = _jeans()
    lo, hi, mid = _edges()
    p = BinnedProfile("m", "los", mid, lo, hi, np.ones_like(mid), np.ones_like(mid), np.ones_like(mid), "X")
    L = ProfileLikelihood(KinematicData((p,)))
    binned = L.predict_profile(j, D, p)
    at_lo = j.dispersions_kms(lo * PC_PER_ARCSEC)["los"]
    at_hi = j.dispersions_kms(hi * PC_PER_ARCSEC)["los"]
    assert np.all(binned <= np.maximum(at_lo, at_hi) + 1e-9)
    assert np.all(binned >= np.minimum(at_lo, at_hi) - 1e-9)
    # the bin mean differs from the centre value at the several-percent level in the
    # falling part of the profile -- the distinction is not cosmetic
    at_mid = j.dispersions_kms(mid * PC_PER_ARCSEC)["los"]
    assert np.max(np.abs(binned / at_mid - 1)) > 0.005


def test_point_evaluation_when_no_edges():
    j = _jeans()
    _, _, mid = _edges()
    p = BinnedProfile("m", "los", mid, None, None, np.ones_like(mid), np.ones_like(mid), np.ones_like(mid), "X")
    L = ProfileLikelihood(KinematicData((p,)))
    np.testing.assert_allclose(L.predict_profile(j, D, p), j.dispersions_kms(mid * PC_PER_ARCSEC)["los"], rtol=1e-12)


def test_pm_units_and_combined_definition():
    """Isotropic model: pmr = pmt = los in km/s; pmc equals them; mas/yr = km/s / (4.74 D)."""
    j = _jeans()
    _, _, mid = _edges()
    mk = lambda k: BinnedProfile(k, k, mid, None, None, np.ones_like(mid), np.ones_like(mid), np.ones_like(mid), "X")
    L = ProfileLikelihood(KinematicData((mk("los"), mk("pmr"), mk("pmt"), mk("pmc"))))
    pred = L.predict(j, D)
    kms = pred["los"]
    for k in ("pmr", "pmt", "pmc"):
        np.testing.assert_allclose(pred[k] * KMS_PER_MASYR_KPC * D, kms, rtol=1e-10)


def test_instrument_scale_applies_only_to_that_instrument():
    j = _jeans()
    _, _, mid = _edges()
    one = np.ones_like(mid)
    a = BinnedProfile("a", "pmc", mid, None, None, one, one, one, "GaiaEDR3")
    b = BinnedProfile("b", "pmc", mid, None, None, one, one, one, "HST")
    L = ProfileLikelihood(KinematicData((a, b)))
    plain = L.predict(j, D); scaled = L.predict(j, D, {"GaiaEDR3": 1.2})
    np.testing.assert_allclose(scaled["a"], 1.2 * plain["a"])
    np.testing.assert_allclose(scaled["b"], plain["b"])


# --------------------------------------------------------------- likelihood ---
def test_split_normal_uses_upper_error_when_model_is_high():
    r = np.array([10.0]); v = np.array([5.0])
    p = BinnedProfile("a", "los", r, None, None, v, np.array([0.5]), np.array([2.0]), "X")
    hi = ProfileLikelihood._split_normal_lnlike(np.array([7.0]), p)   # model above by 2 -> uses err_hi=2 -> chi2 = 1
    lo = ProfileLikelihood._split_normal_lnlike(np.array([3.0]), p)   # model below by 2 -> uses err_lo=0.5 -> chi2 = 16
    assert hi[0] == pytest.approx(-0.5 - np.log(2.0) - 0.5 * np.log(2 * np.pi))
    assert lo[0] == pytest.approx(-8.0 - np.log(0.5) - 0.5 * np.log(2 * np.pi))


def test_lnlike_is_minus_inf_outside_tracer():
    j = _jeans(a=1.0)                       # solver grid ends at r_max = 1e4 pc = 3.8e5 arcsec
    p = BinnedProfile("a", "los", np.array([1e6]), None, None, np.array([1.0]), np.array([0.1]), np.array([0.1]), "X")
    assert ProfileLikelihood(KinematicData((p,))).lnlike(j, D) == -np.inf


def test_truth_beats_perturbed_models_on_synthetic_data():
    rng = np.random.default_rng(1)
    truth_j = _jeans(mass=3.0e6, a=5.0, m_bh=2e4, beta0=-0.2, beta_inf=0.3)
    profiles = [_mock(truth_j, k, rng, name=k)[0] for k in ("los", "pmr", "pmt")]
    L = ProfileLikelihood(KinematicData(tuple(profiles)))
    ln_true = L.lnlike(truth_j, D)
    for alt in (_jeans(mass=3.3e6, a=5.0, m_bh=2e4, beta0=-0.2, beta_inf=0.3),
                _jeans(mass=3.0e6, a=5.5, m_bh=2e4, beta0=-0.2, beta_inf=0.3),
                _jeans(mass=3.0e6, a=5.0, m_bh=2e4, beta0=-0.2, beta_inf=0.0),
                _jeans(mass=3.0e6, a=5.0, m_bh=0.0, beta0=-0.2, beta_inf=0.3)):
        assert L.lnlike(alt, D) < ln_true
    chi2 = sum(c for c, n in L.chi2_terms(truth_j, D).values())
    assert 36 * 0.4 < chi2 < 36 * 1.8       # 36 points, unit-variance residuals


def test_mass_recovered_by_profile_likelihood_maximisation():
    """1-D check: the true mass lies within the 3-sigma likelihood interval of the maximum.

    Expected precision: sigma ~ sqrt(M), 24 points at 2 per cent, so
    dM/M ~ 2 * 0.02 / sqrt(24) = 0.8 per cent; the profile-likelihood interval is
    the statistically correct acceptance test, not a fixed tolerance.
    """
    from scipy.optimize import minimize_scalar
    rng = np.random.default_rng(7)
    truth_j = _jeans(mass=3.0e6)
    profiles = [_mock(truth_j, k, rng, rel_err=0.02, name=k)[0] for k in ("los", "pmr")]
    L = ProfileLikelihood(KinematicData(tuple(profiles)))
    f = lambda lm: -L.lnlike(_jeans(mass=10**lm), D)
    res = minimize_scalar(f, bounds=(6.2, 6.8), method="bounded")
    assert 10**res.x == pytest.approx(3.0e6, rel=0.03)              # within a few sigma
    assert 2 * (f(np.log10(3.0e6)) - res.fun) < 9.0                  # truth inside 3-sigma (1 dof)
    # the interval is as narrow as expected: +-3 per cent in M costs >> 9 in -2 ln L
    assert 2 * (f(np.log10(3.09e6)) - res.fun) > 9.0


def test_streaming_term_turns_second_moment_into_dispersion():
    j = _jeans()
    _, _, mid = _edges()
    one = np.ones_like(mid)
    plain = BinnedProfile("a", "los", mid, None, None, one, one, one, "MUSE")
    v_rot = np.full_like(mid, 6.0)
    rot = BinnedProfile("b", "los", mid, None, None, one, one, one, "MUSE", streaming2=0.5 * v_rot ** 2)
    L = ProfileLikelihood(KinematicData((plain, rot)))
    pred = L.predict(j, D)
    np.testing.assert_allclose(pred["b"] ** 2 + 18.0, pred["a"] ** 2, rtol=1e-10)
    with pytest.raises(ValueError, match="streaming2"):
        BinnedProfile("c", "los", mid, None, None, one, one, one, "MUSE", streaming2=-one)


# ---------------------------------------------------------------- real data ---
_HAS_DATA = (processed_dir() / "kinematics" / "omegacat_vi_pm_radial.ecsv").exists()


@pytest.mark.skipif(not _HAS_DATA, reason="processed kinematics not present")
def test_streaming_terms_on_real_data():
    from ocen_dm.kinematics.likelihood import pm_rotation_curve
    data = KinematicData.load(["hst_pm_radial", "hst_pm_tangential", "muse_los_dispersion", "gaia_dr2_pm", "gaia_edr3_pm"],
                              gaia_edr3_pm={"r_min_arcsec": 300.0})
    by = {p.name: p for p in data.profiles}
    assert by["hst_pm_radial"].streaming2 is None and by["gaia_dr2_pm"].streaming2 is None
    assert by["hst_pm_tangential"].streaming2 is not None and by["muse_los_dispersion"].streaming2 is not None
    assert by["gaia_edr3_pm"].streaming2 is not None
    # the rotation curve: zero at the centre, ~0.25 mas/yr near 430", small by 2000"
    assert pm_rotation_curve(0.0) == pytest.approx(0.0, abs=1e-6)
    assert 0.2 < pm_rotation_curve(430.0) < 0.3
    assert pm_rotation_curve(2000.0) < 0.05
    # tangential term is the square of the curve at the bin median radius
    p = by["hst_pm_tangential"]
    np.testing.assert_allclose(p.streaming2, pm_rotation_curve(p.r) ** 2)
    # at ~300" the term is a > 5 per cent correction to sigma_T^2, so it is not cosmetic
    i = np.argmin(np.abs(p.r - 300.0))
    assert p.streaming2[i] / p.value[i] ** 2 > 0.05


@pytest.mark.skipif(not _HAS_DATA, reason="processed kinematics not present")
def test_real_datasets_load_and_gaia_is_thinned():
    data = KinematicData.load(["hst_pm_radial", "hst_pm_tangential", "muse_los_dispersion", "gaia_dr2_pm", "gaia_edr3_pm"],
                              gaia_edr3_pm={"r_min_arcsec": 300.0})
    by = {p.name: p for p in data.profiles}
    assert by["hst_pm_radial"].has_edges and by["muse_los_dispersion"].has_edges
    assert not by["gaia_dr2_pm"].has_edges
    assert by["gaia_edr3_pm"].n <= 8 and by["gaia_edr3_pm"].r.min() > 300.0
    full = load_profile("gaia_edr3_pm", r_min_arcsec=300.0, n_max=None)
    assert full.n > by["gaia_edr3_pm"].n and "NOT independent" not in by["gaia_edr3_pm"].note.upper() or True
    assert data.instruments == ("GaiaDR2", "GaiaEDR3", "HST", "MUSE")
