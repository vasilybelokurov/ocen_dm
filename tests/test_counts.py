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


def test_fixed_amplitude_matches_profiled_at_optimum_and_is_worse_elsewhere():
    rng = np.random.default_rng(4)
    p, shape, _ = _profile(rng, a=2e-3, b=3e-4)
    free = p.compare(shape)
    pinned = p.compare(shape, amplitude=free["amplitude"])
    assert pinned["amplitude_fixed"] and not free["amplitude_fixed"]
    assert pinned["deviance"] == pytest.approx(free["deviance"], rel=1e-6)
    assert pinned["field"] == pytest.approx(free["field"], rel=1e-4)
    worse = p.compare(shape, amplitude=1.5*free["amplitude"])
    assert worse["amplitude"] == pytest.approx(1.5*free["amplitude"]) and worse["deviance"] > free["deviance"]+10
    assert worse["field"] >= 0.       # the field is re-profiled but never negative
    assert p.compare(shape, amplitude=10*free["amplitude"])["field"] == 0.   # far too many model stars
    q, shape_q, _ = _profile(rng, a=2e-3, b=0., fit_field=False)
    assert q.compare(shape_q, amplitude=1e-3)["field"] == 0.
    with pytest.raises(ValueError):
        p.compare(shape, amplitude=-1.)


def _quadrature():
    import importlib.util
    from pathlib import Path
    spec = importlib.util.spec_from_file_location("bcp", Path(__file__).resolve().parents[1]/"bin"/"build_count_profiles.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.quadrature


@pytest.mark.parametrize("slope", [0., -1., -2., -3.])
def test_area_weighted_bin_average_matches_annulus_integral(slope):
    quadrature = _quadrature()
    edges = np.geomspace(2., 250., 21)          # the HST product's bins
    nodes, weights, areas = zip(*(quadrature(lo, hi) for lo, hi in zip(edges[:-1], edges[1:])))
    p = CountProfile("pl", edges[:-1], edges[1:], np.array(nodes), np.ones(20, int), np.array(areas), "s", "sel",
                     fit_field=False, node_weights=np.array(weights))
    got = p.bin_average(np.array(nodes)**slope)
    lo, hi = edges[:-1], edges[1:]
    exact = (2*np.log(hi/lo) if slope == -2 else 2*(hi**(slope+2)-lo**(slope+2))/(slope+2))/(hi**2-lo**2)
    np.testing.assert_allclose(got, exact, rtol=1e-3)   # vs a 1-6% star-weighted bias
    np.testing.assert_allclose(areas, np.pi*(hi**2-lo**2), rtol=1e-12)


def test_partial_coverage_weights_follow_covered_area():
    quadrature = _quadrature()
    cut = 100.*2**.5                      # sub-ring edge 4 of 8 in [100, 200)
    nodes, w, area = quadrature(100., 200., lambda r1, r2: np.pi*(r2**2-r1**2)*(r1 < cut-1e-9))
    assert np.all(w[4:] == 0.) and np.all(w[:4] > 0.) and w.sum() == pytest.approx(1.)
    assert area == pytest.approx(np.pi*(cut**2-100.**2), rel=1e-12)


def test_weighted_profile_roundtrip_and_legacy_replay():
    lo, hi = np.array([1., 2.]), np.array([2., 4.])
    nodes = np.array([[1.2, 1.8], [2.5, 3.5]])
    w = np.array([[.3, .7], [.5, .5]])
    p = CountProfile("x", lo, hi, nodes, np.array([10, 5]), np.ones(2), "s", "sel", fit_field=False, node_weights=w)
    q = CountProfile.from_dict(p.to_dict())
    assert p.to_dict()["schema_version"] == 2
    np.testing.assert_array_equal(q.node_weights, w)
    np.testing.assert_allclose(q.bin_average(nodes), [.3*1.2+.7*1.8, 3.])
    legacy = CountProfile("x", lo, hi, nodes, np.array([10, 5]), np.ones(2), "s", "sel", fit_field=False)
    assert legacy.to_dict()["schema_version"] == 1 and CountProfile.from_dict(legacy.to_dict()).node_weights is None
    np.testing.assert_allclose(legacy.bin_average(nodes), [1.5, 3.])
    with pytest.raises(ValueError):
        CountProfile("x", lo, hi, nodes, np.array([10, 5]), np.ones(2), "s", "sel", node_weights=np.array([[.3, .6], [.5, .5]]))


def _pair(rng, a=0.2, ratio=12., b=2e-5, c_num=None, c_den=None):
    """Two products of one population: 'inner' (no field) and 'outer' (field); complete-sample ratio ``ratio``."""
    def make(name, edges, amp, field, comp):
        lo, hi = edges[:-1], edges[1:]
        area = np.pi*(hi**2-lo**2)
        nodes = np.sqrt(lo*hi)[:, None]*np.array([.8, 1., 1.2])[None, :]
        shape = 1./(1.+(nodes/100.)**2)
        c = np.ones(len(lo)) if comp is None else comp
        mu = area*(amp*c*shape.mean(axis=1)+field)
        p = CountProfile(name, lo, hi, nodes, rng.poisson(mu), area, "synthetic",
                         "toy", fit_field=field > 0, completeness=comp)
        return p, shape
    num, s_num = make("inner", np.geomspace(5., 250., 11), a, 0., c_num)
    den, s_den = make("outer", np.geomspace(300., 2400., 9), a/ratio, b, c_den)
    return num, s_num, den, s_den


def test_completeness_scales_the_model_and_roundtrips():
    rng = np.random.default_rng(7)
    c = np.linspace(.8, 1., 8)
    _, _, den, s = _pair(rng, c_den=c)
    out = den.compare(s)
    assert out["amplitude"] == pytest.approx(0.2/12., rel=.05)          # complete-sample amplitude recovered
    q = CountProfile.from_dict(den.to_dict())
    np.testing.assert_array_equal(q.completeness, c)
    assert q.compare(s)["deviance"] == pytest.approx(out["deviance"])
    bare = CountProfile.from_dict({k: v for k, v in den.to_dict().items() if k != "completeness"})
    assert bare.completeness is None and bare.compare(s)["deviance"] > out["deviance"]
    with pytest.raises(ValueError):
        CountProfile.from_dict(dict(den.to_dict(), completeness=[1.2]*den.n))


def test_tied_compare_recovers_ratio_and_reduces_to_free_fits_with_a_wide_prior():
    from ocen_dm.kinematics.counts import CountTie, tied_compare
    rng = np.random.default_rng(11)
    num, s_num, den, s_den = _pair(rng, ratio=12., c_den=np.full(8, .95))
    tn, td, tt = tied_compare(num, den, s_num, s_den, CountTie("inner", "outer", 12., .2))
    assert tt["ratio"] == pytest.approx(12., abs=.5) and tn["amplitude"]/td["amplitude"] == pytest.approx(tt["ratio"])
    assert tt["chi2"] == pytest.approx(((tt["ratio"]-12.)/.2)**2)
    # a very wide prior reproduces the separately profiled amplitudes
    fn, fd = num.compare(s_num), den.compare(s_den)
    wn, wd, wt = tied_compare(num, den, s_num, s_den, CountTie("inner", "outer", fn["amplitude"]/fd["amplitude"], 1e3))
    assert wn["amplitude"] == pytest.approx(fn["amplitude"], rel=1e-5) and wd["amplitude"] == pytest.approx(fd["amplitude"], rel=1e-4)
    assert wn["deviance"]+wd["deviance"] == pytest.approx(fn["deviance"]+fd["deviance"], abs=1e-4)
    # a prior far from the free ratio costs deviance + prior
    _, _, bad = tied_compare(num, den, s_num, s_den, CountTie("inner", "outer", 8., .1))
    assert bad["chi2"] > 1. and bad["ratio"] < 12.
    with pytest.raises(ValueError):
        tied_compare(den, num, s_den, s_num, CountTie("outer", "inner", 1/12., .01))   # numerator with a field
