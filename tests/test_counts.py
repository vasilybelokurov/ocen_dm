import numpy as np
import pytest

from ocen_dm.kinematics.counts import CountProfile, deviance_residuals, poisson_profile_fit


def _profile(rng, a=2.0, b=0.3, fit_field=True):
    edges = np.geomspace(10., 1000., 13)
    lo, hi = edges[:-1], edges[1:]
    area = np.pi*(hi**2-lo**2)*0.9
    nodes = np.sqrt(lo*hi)[:, None]*np.array([.7, .9, 1.0, 1.1, 1.3])[None, :]
    shape = 1./(1.+(nodes/100.)**2)          # a toy Sigma(R)
    mu = area*(a*shape.mean(axis=1)+b)
    counts = rng.poisson(mu)
    return CountProfile("toy", lo, hi, nodes, counts, area, "synthetic", "toy", fit_field=fit_field), shape, mu


def test_poisson_fit_recovers_amplitude_and_field():
    rng = np.random.default_rng(1)
    p, shape, mu_true = _profile(rng, a=2e-3, b=3e-4)
    out = p.compare(shape)
    assert out["amplitude"] == pytest.approx(2e-3, rel=.05)
    assert out["field"] == pytest.approx(3e-4, rel=.15)
    assert out["deviance"] == pytest.approx(out["residual"] @ out["residual"])
    assert out["deviance"] < 3*p.n  # a correct model has deviance ~ n


def test_field_is_pinned_at_zero_when_not_fitted_and_never_negative():
    rng = np.random.default_rng(2)
    p, shape, _ = _profile(rng, a=2e-3, b=0., fit_field=False)
    assert p.compare(shape)["field"] == 0.
    p2, shape2, _ = _profile(rng, a=2e-3, b=0., fit_field=True)
    assert p2.compare(shape2)["field"] >= 0.


def test_deviance_residuals_and_perfect_model():
    n = np.array([0, 3, 10, 50]); mu = np.array([1., 3., 8., 60.])
    r = deviance_residuals(n, mu)
    assert r[1] == pytest.approx(0.) and r[2] > 0 and r[3] < 0 and r[0] == pytest.approx(-np.sqrt(2.))
    np.testing.assert_allclose(deviance_residuals(mu, mu), 0., atol=1e-12)


def test_wrong_shape_gives_larger_deviance_and_roundtrip():
    rng = np.random.default_rng(3)
    p, shape, _ = _profile(rng, a=2e-3, b=3e-4)
    good = p.compare(shape)["deviance"]
    bad = p.compare(shape**2)["deviance"]
    assert bad > good+10
    q = CountProfile.from_dict(p.to_dict())
    np.testing.assert_array_equal(q.counts, p.counts)
    assert q.fit_field == p.fit_field and q.compare(shape)["deviance"] == pytest.approx(good)


def test_validation_rejects_bad_inputs():
    lo = np.array([1., 2.]); hi = np.array([2., 3.]); nodes = np.array([[1.5], [2.5]])
    with pytest.raises(ValueError):
        CountProfile("x", lo, hi, nodes, np.array([1, -1]), np.ones(2), "s", "sel")
    with pytest.raises(ValueError):
        CountProfile("x", lo, hi, nodes, np.array([1, 2]), np.array([1., 0.]), "s", "sel")
    with pytest.raises(ValueError):
        poisson_profile_fit([1, 2], [1., 1.], [-1., 1.])
