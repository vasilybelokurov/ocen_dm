"""Analytical and independent numerical controls for the compact stellar DF."""
from dataclasses import replace
import json

import numpy as np
import pytest
from scipy.special import gamma

from ocen_dm.kinematics.df_fit import FitCoordinate, config_at, model_config_from_dict
from ocen_dm.kinematics.regularized_df import (
    ContourNumerics, ExponentialParameters, PrescribedMatter, RegularizedDFConfig,
    RegularizedExponentialDF, SphericalFrequencyRatio,
)
from ocen_dm.kinematics.positive_df import agama_pc, spherical_velocity_moments


@pytest.fixture(scope="module")
def agama():
    pytest.importorskip("agama")
    return agama_pc()


@pytest.mark.parametrize("alpha", [.7, 1., 2.])
def test_harmonic_endpoint_and_normalization(alpha):
    p = ExponentialParameters(J0=80., alpha=alpha, b_out=0)
    df = RegularizedExponentialDF(2.3e6, p, frequency_ratio=lambda jr, l: 2+0*(jr+l))
    jr = np.geomspace(1e-8, 1e7, 71)
    l = jr[::-1]*.17
    np.testing.assert_allclose(df.contour_label(jr, l), 2*jr+l, rtol=3e-6)
    expected = p.J0**3*gamma(3/alpha)/(2*alpha)
    assert df.integral == pytest.approx(expected, rel=2e-6)
    assert df.normalization_integral(320, extend=2) == pytest.approx(expected, rel=3e-6)
    assert df.totalMass() == pytest.approx(2.3e6, rel=2e-6)
    assert df([0, 0, 0]) == df.amplitude


@pytest.mark.parametrize("kw", [dict(alpha=0), dict(alpha=.5), dict(J_a=0), dict(b_out=np.nan),
                              dict(b_outer=-.4), dict(b_outer=-.4, J_outer=20)])
def test_invalid_stellar_parameters(kw):
    with pytest.raises(ValueError):
        ExponentialParameters(**kw)


def test_configuration_and_optimizer_roundtrip():
    cfg = RegularizedDFConfig(matter=PrescribedMatter(rho20=2., gamma=1.))
    assert model_config_from_dict(json.loads(json.dumps(cfg.to_dict()))) == cfg
    coordinates = [FitCoordinate("stellar.b_out", -1, 1),
                   FitCoordinate("matter.rho20", .1, 10, log=True)]
    new = config_at(cfg, coordinates, [-.5, np.log(3)])
    assert new.stellar.b_out == -.5
    assert new.matter.rho20 == pytest.approx(3.)
    assert cfg.matter.rho20 == 2.
    with pytest.raises(ValueError):
        FitCoordinate("mixture_log_ratio.1", -3, 3).get(cfg)
    with pytest.raises(ValueError):
        model_config_from_dict(dict(family="unrecognized"))


def test_density_at_20pc_and_no_dm():
    for slope in (0., 1.):
        matter = PrescribedMatter(rho20=3.7, r_s=17, r_t=43, gamma=slope)
        assert matter.halo_density(20) == pytest.approx(3.7, rel=1e-14)
        assert np.all(matter.halo_density([.1, 1, 10, 100]) >= 0)
    np.testing.assert_array_equal(PrescribedMatter(gamma=1).halo_density([0, 20]), 0)
    with pytest.raises(ValueError):
        PrescribedMatter(gamma=2.)


def test_exact_isochrone_frequencies_and_ergodic_contours(agama):
    p = agama.Potential(type="Isochrone", mass=3e6, scaleRadius=5.)
    pars = ExponentialParameters(b_out=0)
    df = RegularizedExponentialDF(3e6, pars, potential=p)
    rng = np.random.default_rng(32)
    jr, l = 10**rng.uniform(-3, 3, (2, 80))
    k = agama.G*3e6*5
    np.testing.assert_allclose(df.frequency_ratio(jr, l),
                               2/(1+l/np.sqrt(l*l+4*k)), rtol=2e-10)
    h = jr+.5*(l+np.sqrt(l*l+4*k))
    expected = h-k/h
    np.testing.assert_allclose(df.contour_label(jr, l), expected, rtol=2e-5)
    # At fixed radius and speed, orientation must not change an ergodic f(E).
    theta = np.linspace(0, np.pi/2, 41)
    xv = np.zeros((len(theta), 6))
    xv[:, 0] = 8.
    xv[:, 3], xv[:, 4] = 12*np.cos(theta), 12*np.sin(theta)
    f = df(agama.ActionFinder(p)(xv))
    assert np.ptp(f)/np.mean(f) < 5e-5
    rho, pr, pt = spherical_velocity_moments(p, df, agama.ActionFinder(p), [.1, 1, 10, 60], 64).T
    np.testing.assert_allclose(pt/pr, 1., atol=1e-4)
    with pytest.raises(ValueError, match="point mass"):
        SphericalFrequencyRatio(agama.Potential(type="Plummer", mass=1e4, scaleRadius=0.))


@pytest.mark.parametrize("bias", [-.8, .8])
def test_contour_consistency_regularity_and_sign(agama, bias):
    pot = agama.Potential(type="Isochrone", mass=3e6, scaleRadius=5.)
    df = RegularizedExponentialDF(3e6, ExponentialParameters(b_out=bias), potential=pot)
    jr = np.geomspace(.1, 3000, 35)
    l = .23*jr
    np.testing.assert_allclose(df.contour_label(jr, l), df.direct_contour_label(jr, l), rtol=4e-5)
    eps = jr*1e-6
    dr = (df.contour_label(jr+eps, 0)-df.contour_label(jr-eps, 0))/(2*eps)
    dl = (df.contour_label(jr, eps)-df.contour_label(jr, 0))/eps
    np.testing.assert_allclose(dr/dl, 2., rtol=3e-5)
    # Spherical symmetry includes the sign of Jphi, not a positive-only wedge.
    a = np.column_stack((jr, l*.7, l*.3))
    b = np.column_stack((jr, l, l*0))
    np.testing.assert_allclose(df(a), df(b), rtol=1e-12)
    a[:, 2] *= -1
    np.testing.assert_allclose(df(a), df(b), rtol=1e-12)
    _, pr, pt = spherical_velocity_moments(pot, df, agama.ActionFinder(pot), [15., 40.], 64).T
    assert np.all((1-pt/pr)*bias > 0)


def test_two_transitions_have_correct_limits():
    p = ExponentialParameters(b_out=.8, J_a=100., b_outer=-.6, J_outer=3000.)
    assert abs(p.bias(1e-8)) < 1e-15
    assert p.bias(500) > .7
    assert p.bias(1e9) == pytest.approx(-.6, abs=1e-8)


def test_two_transition_velocity_anisotropy(agama):
    pot = agama.Potential(type="Isochrone", mass=3e6, scaleRadius=5.)
    p = ExponentialParameters(b_out=.8, J_a=50., b_outer=-.6, J_outer=500.)
    df = RegularizedExponentialDF(3e6, p, potential=pot)
    _, pr, pt = spherical_velocity_moments(pot, df, agama.ActionFinder(pot), [.1, 10., 80.], 64).T
    beta = 1-pt/pr
    assert abs(beta[0]) < .001 and beta[1] > .05 and beta[2] < -.2


def test_joint_origin_velocity_derivative():
    # At r=0 in a harmonic core, Jr=v^2/(4 Omega). alpha>1/2 makes
    # [f(v)-f(0)]/v tend to zero. This is separate from fixed-Jr, L->0.
    df = RegularizedExponentialDF(1., ExponentialParameters(alpha=.7, J0=1, b_out=.8),
                                  frequency_ratio=lambda jr, l: 2+0*(jr+l))
    v = np.array([1e-2, 1e-3, 1e-4])
    actions = np.column_stack((v*v/4, v*0, v*0))
    slopes = abs((df(actions)/df.amplitude-1)/v)
    assert slopes[2] < slopes[1] < slopes[0]
    np.testing.assert_allclose(slopes[1:]/slopes[:-1], 10**(-.4), rtol=.002)


def test_frequency_interpolation_failure_regression(agama):
    # A native interpolated frequency derivative is negative at this orbit.
    # The direct integral remains physical, including the circular limit.
    pot = agama.Potential(type="Plummer", mass=3e6, scaleRadius=5.)
    ratio = SphericalFrequencyRatio(pot)
    values = ratio(np.array([575567.514, 2.9e-12, 0]), np.array([874.366, 1.1e-5, 100.]))
    assert np.all(np.isfinite(values) & (values >= 1) & (values <= 2.00001))
