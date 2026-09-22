"""Independent King-limit, Jeans and embedded-potential controls."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.special import gammainc

from ocen_dm.kinematics.df_mock import LoweredIsothermalMock, MockPopulation, king_integrals
from ocen_dm.kinematics.df_fit import model_config_from_dict
from ocen_dm.kinematics.lowered_isothermal import (
    LoweredIsothermalConfig, LoweredIsothermalModel, lowered_exponential, lowered_moments,
)
from ocen_dm.kinematics.regularized_df import PrescribedMatter
from ocen_dm.kinematics.positive_df import agama_pc


@pytest.fixture(scope="module")
def agama():
    pytest.importorskip("agama")
    return agama_pc()


def test_lowering_function_limits():
    x = np.geomspace(1e-10, 10, 21)
    np.testing.assert_allclose(lowered_exponential(0, x), np.exp(x), rtol=1e-13)
    np.testing.assert_allclose(lowered_exponential(1, x), np.expm1(x), rtol=1e-13)
    for g in (.5, 1.5, 2.):
        np.testing.assert_allclose(lowered_exponential(g, x), np.exp(x)*gammainc(g, x), rtol=3e-14)
    np.testing.assert_array_equal(lowered_exponential(0, [-1, 0]), 0)


def test_king_density_and_pressure_analytic_integrals():
    W = np.array([1e-5, .1, 1, 7.])
    s = 15.
    numerical = lowered_moments(W, 0, s, 40, 1)
    rho, p = king_integrals(W)
    np.testing.assert_allclose(numerical[:, 0], (2*np.pi*s*s)**1.5*rho, rtol=1e-12)
    np.testing.assert_allclose(numerical[:, 1], (2*np.pi*s*s)**1.5*s*s*p, rtol=1e-12)
    np.testing.assert_allclose(numerical[:, 1], numerical[:, 2], rtol=1e-12)


def test_isotropic_independent_poisson_solution(agama):
    mock = LoweredIsothermalMock([MockPopulation("stars", 1., 1.)])
    config = LoweredIsothermalConfig(W0=mock.W0, rho0=mock.rho0,
                                     velocity_scale=mock.velocity_kms, anisotropy_radius=1e10)
    model = LoweredIsothermalModel(config)
    assert model_config_from_dict(config.to_dict()) == config
    assert model.M_star == pytest.approx(mock.total_mass, rel=2e-7)
    assert model.r_t == pytest.approx(mock.r_t, rel=2e-7)
    r = np.geomspace(.01, .95*mock.r_t, 31)
    actual = model.intrinsic_moments(r)
    expected = mock.intrinsic_components(r)[:, 0, :]
    for i, k in enumerate(("rho", "radial_pressure", "tangential_pressure")):
        np.testing.assert_allclose(actual[k], expected[:, i], rtol=1e-6)
    R = np.geomspace(.03, .8*mock.r_t, 12)
    for key, expected in mock.projected_moments(R).items():
        np.testing.assert_allclose(model.projected_moments(R)[key], expected, rtol=1e-6)


@pytest.mark.parametrize("slope", [0., 1.])
def test_embedded_control_jeans_and_boundary(agama, slope):
    cfg = LoweredIsothermalConfig(matter=PrescribedMatter(rho20=2., gamma=slope, M_rem=3e5))
    model = LoweredIsothermalModel(cfg)
    r = np.geomspace(.05, .8*model.r_t, 30)
    m = model.intrinsic_moments(r)
    eps = 1e-4
    gradient = (model.intrinsic_moments(r*(1+eps))["radial_pressure"]-
                model.intrinsic_moments(r*(1-eps))["radial_pressure"])/(2*eps*r)
    lhs = gradient+2*m["beta"]*m["radial_pressure"]/r
    rhs = m["rho"]*model.potential.force(np.column_stack((r, r*0, r*0)))[:, 0]
    np.testing.assert_allclose(lhs, rhs, rtol=4e-4)
    beta = model.intrinsic_moments([.001, 30., model.r_t*.9999])["beta"]
    assert abs(beta[0]) < 1e-7 and beta[1] > .05 and abs(beta[2]) < .001
    xv = [[1., 0, 0, 2, 3, 4], [1., 0, 0, -2, -3, -4], [2*model.r_t, 0, 0, 0, 0, 0]]
    f = model.phase_space_density(xv)
    assert f[0] > 0 and f[0] == pytest.approx(f[1]) and f[2] == 0
    # Refined velocity quadrature checks the anisotropic angular integral.
    high = model.intrinsic_moments(r, 128)
    np.testing.assert_allclose(m["radial_pressure"], high["radial_pressure"], rtol=1e-8)


@pytest.mark.parametrize("kw", [dict(W0=-1), dict(truncation_index=-1), dict(rho0=0),
                              dict(anisotropy_radius=np.nan)])
def test_invalid_config(kw):
    with pytest.raises(ValueError):
        LoweredIsothermalConfig(**kw)


def test_matched_mock_through_shared_joint_likelihood(agama):
    from ocen_dm.kinematics.df_fit import DFJointProblem, PhotometricData
    from ocen_dm.kinematics.likelihood import ARCSEC_PER_RAD, BinnedProfile, KinematicData
    cfg = LoweredIsothermalConfig(matter=PrescribedMatter(rho20=2., gamma=1.))
    model = LoweredIsothermalModel(cfg)
    angular = np.array([10., 100., 600.])
    projected = model.projected_moments(angular*cfg.distance_kpc*1000/ARCSEC_PER_RAD)
    profile = BinnedProfile("synthetic_los", "los", angular, None, None,
                            np.sqrt(projected["los"]/projected["Sigma"]),
                            np.ones(3)*.3, np.ones(3)*.3, "synthetic")
    photo = PhotometricData(angular, -2.5*np.log10(projected["Sigma"])+23,
                            np.full(3, .03), "matched mock", "known independent errors")
    # No supplied model: exercise configuration dispatch and rebuild as in a fit.
    evaluation = DFJointProblem(KinematicData((profile,)), photo).evaluate(cfg)
    assert evaluation["objective"] < 1e-18
