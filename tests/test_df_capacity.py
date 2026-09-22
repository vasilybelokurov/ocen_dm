"""Independent physics, derivatives, binning and recovery controls."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import simpson
from scipy.optimize import least_squares

from ocen_dm.kinematics.df_mock import king_integrals, challenge_mock
from ocen_dm.kinematics.df_capacity import (
    dpl_shape, FixedPotentialProjector, CapacityObservations, capacity_bounds,
)
from ocen_dm.kinematics.df_fit import PhotometricData
from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData, ProfileLikelihood
from ocen_dm.kinematics.positive_df import ActionComponent, agama_pc


@pytest.fixture(scope="module")
def truth():
    pytest.importorskip("agama")
    return challenge_mock("mixed_with_dm")


@pytest.fixture(scope="module")
def observations(truth):
    r = np.geomspace(4., 1500., 8)
    profiles = []
    for i, kind in enumerate(("los", "pmr", "pmt")):
        profiles.append(BinnedProfile(kind, kind, r, r*.9, r*1.1, r*0+1, r*0+.03,
                                     r*0+.04, "mock", streaming2=r*0+10,
                                     r_nodes=np.array([r*.96, r*1.04]).T if i == 1 else None))
    photo = PhotometricData(r, r*0, r*0+.1, "test", "adopted")
    return CapacityObservations(KinematicData(tuple(profiles)), photo, truth)


@pytest.fixture(scope="module")
def kernel(truth, observations):
    return FixedPotentialProjector(truth.agama_potential(), observations.radii,
                                   radial_nodes=200, velocity_nodes=48, projection_nodes=100)


def parameters(count):
    shapes = [[np.log(150*(i+1)), .1-i*.2, 7+i, np.log(2.), 1.3, 1.1+i*.3]
              for i in range(count)]
    return np.r_[np.ravel(shapes), np.arange(1, count)*-.4]


def test_king_integrals_direct_velocity_quadrature():
    x, w = np.polynomial.legendre.leggauss(160)
    for W in (1e-8, .001, .1, 1., 7., 20.):
        v = np.sqrt(2*W)*(x+1)/2
        f = np.expm1(W-v*v/2)
        expected = np.array([4*np.pi*np.sum(w*f*v*v)*np.sqrt(2*W)/2,
                             4*np.pi/3*np.sum(w*f*v**4)*np.sqrt(2*W)/2])/(2*np.pi)**1.5
        np.testing.assert_allclose(king_integrals(W), expected, rtol=2e-13)


def test_mock_df_moments_and_anisotropy(truth):
    x, w = np.polynomial.legendre.leggauss(180)
    u, w = (x+1)/2, w/2
    for r in (.1, 10., 50.):
        vmax = np.sqrt(2*truth._state(np.array(r))[0])*truth.velocity_kms
        numerical = []
        mu = u[None, :]
        for i, ra in enumerate(truth.ra):
            # Integrate to the actual Q=0 ellipsoid; no unresolved internal cut.
            limit = vmax/np.sqrt(mu*mu+(1+(r/ra)**2)*(1-mu*mu))
            v = limit*u[:, None]
            f = truth.phase_space_components(r, v*mu, v*np.sqrt(1-mu*mu))[..., i]
            assert np.all(f >= 0)
            weight = 4*np.pi*v*v*limit*w[:, None]*w[None, :]
            numerical.append([np.sum(f*weight*k) for k in (1., (v*mu)**2,
                                                           v*v*(1-mu*mu)/2)])
        np.testing.assert_allclose(numerical, truth.intrinsic_components(np.array(r)), rtol=3e-6)
    assert np.all(truth.phase_space_components(truth.r_t*1.01, 0., 0.) == 0)


def test_independent_mass_integral_and_potential(truth):
    r = np.geomspace(1e-5, truth.r_t, 4001)
    rho = truth.intrinsic_components(r)[..., 0]
    integrated = simpson(4*np.pi*r[:, None]**3*rho, x=np.log(r), axis=0)
    np.testing.assert_allclose(integrated, truth.masses, rtol=1e-7)
    p = truth.agama_potential()
    r = np.geomspace(.01, truth.r_t*2, 120)
    xyz = np.column_stack((r, r*0, r*0))
    np.testing.assert_allclose(-p.force(xyz)[:, 0]*r*r/truth.G, truth.enclosed_mass(r), rtol=2e-6)
    np.testing.assert_allclose(p.potential(xyz), truth.potential_value(r), rtol=2e-7)


def test_shape_matches_native_agama_and_derivatives():
    agama = agama_pc()
    rng = np.random.default_rng(74)
    actions = 10**rng.uniform(-2, 5, (120, 3))
    theta = parameters(1)
    f, grad = dpl_shape(actions, theta)
    df = agama.DistributionFunction(type="DoublePowerLaw", norm=(2*np.pi)**3,
         J0=np.exp(theta[0]), slopeIn=theta[1], slopeOut=theta[2], steepness=np.exp(theta[3]),
         coefJrIn=theta[4], coefJzIn=(3-theta[4])/2, coefJrOut=theta[5], coefJzOut=(3-theta[5])/2)
    np.testing.assert_allclose(f, df(actions), rtol=2e-13)
    for i in range(6):
        delta = np.eye(6)[i]*1e-5
        finite = (np.log(dpl_shape(actions, theta+delta, False))-
                  np.log(dpl_shape(actions, theta-delta, False)))/(2e-5)
        np.testing.assert_allclose(grad[:, i], finite, rtol=2e-6, atol=2e-8)


@pytest.mark.parametrize("count", [1, 2, 3])
def test_full_observable_jacobian(kernel, observations, count):
    x = parameters(count)
    value, jac = kernel.project(x, count)
    residual, deriv = observations.residual(value, jac)
    for i in range(len(x)):
        delta = np.eye(len(x))[i]*1e-5
        finite = (observations.residual(kernel.project(x+delta, count, False))-
                  observations.residual(kernel.project(x-delta, count, False)))/(2e-5)
        np.testing.assert_allclose(deriv[:, i], finite, rtol=3e-5, atol=2e-6)


def test_saved_binning_and_unit_conversion(observations, truth):
    expected = observations.likelihood.predict(truth, observations.distance_kpc)
    np.testing.assert_allclose(observations.truth[:observations.nkin],
                               np.concatenate(list(expected.values())), rtol=1e-13)
    # Changing one common photometric zero point has no effect.
    projected = observations.truth_projected.copy()
    projected *= 4.
    np.testing.assert_allclose(observations.residual(projected), 0., atol=1e-12)
    assert all(p.streaming2 is None for p in observations.data.profiles)


def test_direct_agama_projection(kernel):
    x = parameters(2)
    df = kernel.native_df(x, 2)
    gm = agama_pc().GalaxyModel(kernel.potential, df)
    r = kernel.R[::8]
    rho, tensor = gm.moments(np.column_stack((r, r*0)))
    direct = np.column_stack((rho, rho*tensor[:, 2], rho*tensor[:, 0], rho*tensor[:, 1]))
    np.testing.assert_allclose(kernel.project(x, 2, False)[::8], direct, rtol=.003)


def test_local_same_family_recovery(kernel, observations):
    # Only an optimizer sanity control; the production challenge uses independent DFs.
    target = parameters(1)
    projected = kernel.project(target, 1, False)
    target_residual = observations.residual(projected)
    def objective(x):
        m, dm = kernel.project(x, 1)
        r, j = observations.residual(m, dm)
        return r-target_residual, j
    x0 = target+np.array([.02, .03, -.03, .02, -.03, .03])
    result = least_squares(lambda x: objective(x)[0], x0, jac=lambda x: objective(x)[1],
                           bounds=capacity_bounds(1), max_nfev=100, gtol=1e-8, xtol=1e-9, ftol=1e-9)
    assert result.success
    assert np.max(abs(result.fun)) < 1e-5
