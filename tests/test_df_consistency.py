import numpy as np
import pytest

from ocen_dm.kinematics.anisotropy import Anisotropy, TurnoverAnisotropy
from ocen_dm.kinematics.df_consistency import (
    mge_log_slopes, beta_log_derivative, necessary_profiles,
    summarize_profiles, sampled_negative_intervals,
)
from ocen_dm.kinematics.jeans import SphericalJeans
from ocen_dm.mass_models import MGE, PointMass, CompositeMassModel, Plummer
from ocen_dm.mass_models.base import G


def test_single_gaussian_slopes_including_density_underflow():
    tracer = MGE([1.], [2.])
    r = np.array([.01, 1., 3., 1000.])
    gamma, dot = mge_log_slopes(tracer, r)
    np.testing.assert_allclose(gamma, r*r/4)
    np.testing.assert_allclose(dot, r*r/2)


@pytest.mark.parametrize("a", [Anisotropy(-.7, .4, 3.),
                              TurnoverAnisotropy(-.7, .4, -.9, .2, 20.)])
def test_analytic_derivatives_against_independent_finite_differences(a):
    r = np.geomspace(.03, 50., 50)
    h = 1e-4
    expected = (a.beta(r*np.exp(h))-a.beta(r*np.exp(-h)))/(2*h)
    np.testing.assert_allclose(beta_log_derivative(a, r), expected, rtol=2e-7, atol=1e-9)
    tracer = MGE([.2, .8], [.3, 4.])
    model = SphericalJeans(CompositeMassModel([tracer, PointMass(3.), Plummer(10., 2.)]), tracer, a)
    p = lambda x: tracer.density(x)*np.exp(a.log_integrating_factor(x))
    plus, minus, centre = p(r*np.exp(h)), p(r*np.exp(-h)), p(r)
    # Differentiation with respect to ln r, then d/dPsi = -r/(GM) d/dln r.
    plus2, minus2 = p(r*np.exp(2*h)), p(r*np.exp(-2*h))
    dp = (-plus2+8*plus-8*minus+minus2)/(12*h)
    ddp = (-plus2+16*plus-30*centre+16*minus-minus2)/(12*h*h)
    m = model.mass.enclosed_mass(r)
    q = 4*np.pi*r**3*model.mass.density(r)/m
    result = necessary_profiles(model, r)
    np.testing.assert_allclose(result['B1'], -dp/centre, rtol=3e-5, atol=1e-6)
    np.testing.assert_allclose(result['B2'], (ddp+(1-q)*dp)/centre, rtol=1e-4, atol=1e-5)


def test_central_bound_saturation_is_not_sufficient():
    # beta=-1/2 Gaussian tracer in a Kepler potential: P=nu/r,
    # P''/positive_factor = z(z-1). Negative for 0<r<sigma.
    r = np.array([.1, .5, 1., 2.])
    tracer = MGE([1.], [1.])
    a = Anisotropy(-.5, -.5, 1.)
    model = SphericalJeans(PointMass(1/G), tracer, a)
    out = necessary_profiles(model, r)
    np.testing.assert_allclose(out['B2'], r*r*(r*r-1), atol=1e-14)
    assert summarize_profiles(out, -.5)['separable_df_ruled_out']


def test_more_tangential_gaussian_passes_second_derivative_screen():
    r = np.geomspace(.01, 10., 50)
    tracer = MGE([1.], [1.])
    model = SphericalJeans(PointMass(1/G), tracer, Anisotropy(-1., -1., 1.))
    out = necessary_profiles(model, r)
    np.testing.assert_allclose(out['B2'], r**4+r*r+2, atol=1e-12)
    summary = summarize_profiles(out, -1.)
    assert not summary['separable_df_ruled_out']
    assert not summary['df_nonnegative_certified']


def test_second_derivative_not_misapplied_to_isotropic_case():
    r = np.array([.1, .5])
    tracer = MGE([1.], [1.])
    model = SphericalJeans(PointMass(1/G), tracer, Anisotropy())
    out = necessary_profiles(model, r)
    assert np.all(out['B2'] < 0)
    summary = summarize_profiles(out, 0.)
    assert not summary['B2']['necessary_for_separable_df']
    assert not summary['separable_df_ruled_out']
    assert not summary['df_nonnegative_certified']


def test_disjoint_negative_intervals_and_nan_rejection():
    assert sampled_negative_intervals([1,2,3,4,5], [-1,0,-2,-3,0]) == [[1.,1.],[3.,4.]]
    with pytest.raises(ValueError):
        sampled_negative_intervals([1], [np.nan])
