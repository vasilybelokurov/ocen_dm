"""The outer-profile measurement: projection, error deconvolution, binning."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics.outer_profile import MemberSample, binned_dispersion, dispersion_ml, load_members
from ocen_dm.paths import processed_dir


def test_dispersion_ml_recovers_a_known_dispersion_below_the_errors():
    """The regime that matters here: per-star errors larger than the signal."""
    rng = np.random.default_rng(0)
    n, sigma_true, mean_true = 20000, 0.25, 0.07
    err = rng.uniform(0.3, 0.45, n)                      # errors above the dispersion
    v = rng.normal(mean_true, sigma_true, n) + rng.normal(0, err)
    sigma, sigma_err, mean = dispersion_ml(v, err)
    assert sigma == pytest.approx(sigma_true, abs=4 * sigma_err)
    assert mean == pytest.approx(mean_true, abs=0.01)
    assert 0.002 < sigma_err < 0.05


def test_underestimated_errors_masquerade_as_dispersion():
    """A 10 per cent error underestimate inflates a 0.25 mas/yr dispersion by ~40 per cent."""
    rng = np.random.default_rng(1)
    n, sigma_true = 20000, 0.25
    err = np.full(n, 0.38)
    v = rng.normal(0.0, sigma_true, n) + rng.normal(0, 1.1 * err)      # true errors 10 % larger
    naive, _, _ = dispersion_ml(v, err)
    corrected, _, _ = dispersion_ml(v, err, err_scale=1.1)
    expected = np.sqrt(sigma_true**2 + (1.1 * err[0]) ** 2 - err[0] ** 2)      # 0.25 -> 0.305
    assert naive == pytest.approx(expected, abs=0.02) and naive > 1.15 * sigma_true
    assert corrected == pytest.approx(sigma_true, abs=0.02)


def test_contamination_inflates_the_dispersion():
    rng = np.random.default_rng(2)
    n, f = 10000, 0.03
    err = np.full(n, 0.35)
    clean = rng.normal(0, 0.25, n)
    field = rng.normal(0, 3.0, n)
    v = np.where(rng.uniform(size=n) < f, field, clean) + rng.normal(0, err)
    sigma, _, _ = dispersion_ml(v, err)
    assert sigma > 0.4                                  # 3 per cent field stars double it


def _fake_sample(n=4000, sigma=0.3, seed=3):
    rng = np.random.default_rng(seed)
    r = rng.uniform(300, 2400, n)
    err = rng.uniform(0.2, 0.4, n)
    return MemberSample(r, rng.normal(0, sigma, n) + rng.normal(0, err),
                        rng.normal(0, sigma, n) + rng.normal(0, err), err, err,
                        np.full(n, 0.97), rng.uniform(17, 20, n), np.full(n, 3), rng.uniform(-np.pi, np.pi, n))


def test_binned_dispersion_columns_and_recovery():
    t = binned_dispersion(_fake_sample(), np.array([300.0, 1000.0, 1800.0, 2400.0]))
    assert len(t) == 3
    for c in ("r_median", "n_stars", "sigma_pmr", "sigma_pmt", "sigma_pm", "median_err",
              "expected_contaminants", "median_g"):
        assert c in t.colnames
    assert np.allclose(t["sigma_pm"], 0.3, atol=0.05)
    assert np.all(t["expected_contaminants"] == pytest.approx(0.03 * t["n_stars"], rel=1e-6))


_HAS_MEMBERS = (processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv").exists()


@pytest.mark.skipif(not _HAS_MEMBERS, reason="member catalogue not present")
def test_real_members_projection_is_orthogonal_and_errors_are_consistent():
    s = load_members()
    assert len(s) > 1e5
    # the systemic PM must be subtracted, else the projection reads |mu_sys|/sqrt(2) ~ 5 mas/yr
    assert s.mu_sys[0] == pytest.approx(-3.25, abs=0.15) and s.mu_sys[1] == pytest.approx(-6.75, abs=0.15)
    inner = s.select((s.prob > 0.9) & (s.r_arcsec < 600))
    assert abs(np.median(inner.mu_r)) < 0.05                  # no net expansion/contraction inside 600"
    # the tangential component is NOT zero: this is omega Cen's rotation, which peaks at
    # 0.25 mas/yr near 430-580" in the Vasiliev & Baumgardt profile
    assert 0.15 < abs(np.median(inner.mu_t)) < 0.35
    # the projected errors must reproduce the trace of the covariance matrix
    from astropy.table import Table
    t = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    trace = np.asarray(t["pmra_error"], float) ** 2 + np.asarray(t["pmdec_error"], float) ** 2
    np.testing.assert_allclose(s.err_r**2 + s.err_t**2, trace, rtol=1e-10)
    assert np.all(np.isfinite(s.mu_r)) and np.all(np.isfinite(s.mu_t))


def test_free_mixture_recovers_cluster_dispersion_under_heavy_contamination():
    """Cluster sigma = 0.20, errors 0.35, 85 per cent field with sigma 5-7: the mixture must
    recover the cluster and the field fraction; a P-free measurement would be hopeless."""
    from ocen_dm.kinematics.outer_profile import mixture_dispersion_free
    n_c, n_f = 3000, 17000
    sigmas, fs = [], []
    for seed in (5, 6, 7, 8):                      # single realisations scatter by ~0.015 in this regime
        rng = np.random.default_rng(seed)
        err = rng.uniform(0.25, 0.45, n_c + n_f)
        v = np.concatenate([rng.normal(0.0, 0.20, n_c),
                            np.where(rng.uniform(size=n_f) < 0.6, rng.normal(0.5, 5.0, n_f), rng.normal(-1.0, 7.0, n_f))])
        out = mixture_dispersion_free(v + rng.normal(0, err), err)
        sigmas.append(out["sigma"]); fs.append(out["f"])
        assert out["sigma"] == pytest.approx(0.20, abs=0.035)
        assert out["n_cluster"] == pytest.approx(n_c, rel=0.06)
    assert np.mean(sigmas) == pytest.approx(0.20, abs=0.012)          # unbiased on average
    assert np.mean(fs) == pytest.approx(n_f / (n_c + n_f), abs=0.01)


def test_free_mixture_does_not_steal_the_member_wings():
    """No field at all: the floor on the field width stops EM from carving a narrow 'field'
    out of the member distribution (which biased sigma low without it)."""
    from ocen_dm.kinematics.outer_profile import mixture_dispersion_free
    rng = np.random.default_rng(6)
    n = 20000; err = rng.uniform(0.2, 0.4, n); v = rng.normal(0.0, 0.40, n) + rng.normal(0, err)
    out = mixture_dispersion_free(v, err)
    assert out["sigma"] == pytest.approx(0.40, abs=0.01)
    assert out["f"] < 0.02
