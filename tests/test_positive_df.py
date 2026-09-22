"""Physical and numerical controls for the new positive-DF branch."""
from dataclasses import replace
import json

import numpy as np
import pytest
from scipy.integrate import cumulative_simpson
from scipy.interpolate import CubicSpline

from ocen_dm.kinematics.positive_df import (
    ActionComponent, DFConvergenceError, DFModelConfig, DFNumerics,
    PositiveDFModel, agama_pc, spherical_velocity_moments,
)
from ocen_dm.kinematics.df_fit import (
    DFJointProblem, FitCoordinate, PhotometricData, config_at,
)
from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData, ProfileLikelihood


@pytest.mark.parametrize("kw", [dict(J0=0), dict(fraction=-1), dict(h_r=3), dict(g_r=0),
                              dict(slope_in=3), dict(slope_out=3), dict(J0=np.nan)])
def test_invalid_action_df_rejected(kw):
    with pytest.raises(ValueError):
        ActionComponent(**kw)


def test_configuration_roundtrip_and_normalized_mixture():
    cfg = DFModelConfig(components=(ActionComponent(fraction=.3),
                                   ActionComponent(fraction=.7, J0=1000, g_r=2)))
    assert DFModelConfig.from_dict(json.loads(json.dumps(cfg.to_dict()))) == cfg
    with pytest.raises(ValueError, match="sum to one"):
        replace(cfg, components=(ActionComponent(fraction=.8),))
    with pytest.raises(ValueError):
        DFNumerics(velocity_nodes=4)


@pytest.fixture(scope="module")
def agama():
    pytest.importorskip("agama")
    return agama_pc()


@pytest.fixture(scope="module")
def model(agama):
    return PositiveDFModel(DFModelConfig())


def test_positive_df_spherical_invariance_and_mixture_mass(agama):
    a, b = ActionComponent(fraction=.3, slope_in=-1), ActionComponent(fraction=.7, g_r=2)
    da, db = a.build(3e6), b.build(3e6)
    composite = agama.DistributionFunction(da, db)
    rng = np.random.default_rng(6819)
    actions = 10**rng.uniform(-8, 8, (400, 3))
    actions[:, 2] *= rng.choice([-1, 1], len(actions))
    f = composite(actions)
    assert np.all(np.isfinite(f) & (f >= 0))
    np.testing.assert_allclose(f, da(actions)+db(actions), rtol=1e-13)
    equivalent = actions.copy()
    equivalent[:, 1] += abs(equivalent[:, 2])
    equivalent[:, 2] = 0
    np.testing.assert_allclose(f, composite(equivalent), rtol=1e-13)
    # Native composite totalMass performs a fresh numerical action integral.
    assert composite.totalMass() == pytest.approx(3e6, rel=1e-7)


def test_point_mass_is_not_softened(model, agama):
    r = np.geomspace(1e-7, 100, 12)
    xyz = np.column_stack((r, r*0, r*0))
    np.testing.assert_allclose(model.static_potentials[0].force(xyz)[:, 0],
                               -agama.G*model.config.M_bh/r**2, rtol=1e-12)


def test_kepler_isotropic_control(agama):
    # In a Kepler potential Jr+L = GM/sqrt(-2E), so f(Jr+L) is isotropic.
    pot = agama.Potential(type="Plummer", mass=1e6, scaleRadius=0.)
    df = ActionComponent(J0=100, slope_in=1, h_r=1, g_r=1).build(1e4)
    r = np.array([.1, 1., 10.])
    rho, pr, pt = spherical_velocity_moments(pot, df, agama.ActionFinder(pot), r, 96).T
    # AGAMA's spherical action finder interpolates Jr even in this analytic potential.
    np.testing.assert_allclose(pr, pt, rtol=3e-5)
    x, w = np.polynomial.legendre.leggauss(160)
    escape = np.sqrt(2*agama.G*1e6/r)
    v = escape[:, None]*(x[None, :]+1)/2
    action_sum = agama.G*1e6/np.sqrt(escape[:, None]**2-v*v)
    actions = np.column_stack((action_sum.ravel(), np.zeros((action_sum.size, 2))))
    f = df(actions).reshape(v.shape)
    expected_rho = 4*np.pi*np.sum(f*v*v*w[None, :]*escape[:, None]/2, axis=1)
    expected_pressure = 4*np.pi/3*np.sum(f*v**4*w[None, :]*escape[:, None]/2, axis=1)
    np.testing.assert_allclose(rho, expected_rho, rtol=5e-5)
    np.testing.assert_allclose(pr, expected_pressure, rtol=5e-5)
    # Independent native AGAMA adaptive integration, including velocity moments.
    density, tensor = agama.GalaxyModel(pot, df).moments(np.column_stack((r, r*0, r*0)))
    np.testing.assert_allclose(rho, density, rtol=.002)
    np.testing.assert_allclose(pr/rho, tensor[:, 0], rtol=.002)


def test_self_consistency_and_projection(model):
    assert model.diagnostics["density_closure"] < .001
    assert model.diagnostics["mass_error"] < .001
    r = np.geomspace(.02, 80, 14)
    a, b = model.projected_moments(r), model.direct_projected_moments(r)
    for key in a:
        np.testing.assert_allclose(a[key], b[key], rtol=.003)


def test_independent_jeans_equation_and_poisson_mass(model, agama):
    # DF velocity integration must satisfy Jeans without using it to construct f.
    t = model.moment_table
    lr = np.log(t["r"])
    r = np.geomspace(.03, 100, 100)
    m = model.intrinsic_moments(r)
    pressure_slope = CubicSpline(lr, np.log(t["radial_pressure"]))(np.log(r), 1)
    pressure_gradient = m["radial_pressure"]*(pressure_slope+2*m["beta"])/r
    force = model.potential.force(np.column_stack((r, r*0, r*0)))[:, 0]
    np.testing.assert_allclose(pressure_gradient, m["rho"]*force, rtol=.012)
    # Independent radial integration of DF density supplies the stellar gravity.
    enclosed = cumulative_simpson(4*np.pi*t["r"]**3*t["rho"], x=lr, initial=0)
    interp_mass = np.interp(np.log(r), lr, enclosed)+model.config.M_bh
    np.testing.assert_allclose(-force*r*r/agama.G, interp_mass, rtol=.006)


@pytest.mark.parametrize("gamma", [0., 1.])
def test_halo_mass_normalization_and_force(gamma, agama):
    cfg = DFModelConfig(M_dm_100=2e6, gamma=gamma, M_rem=4e5)
    m = PositiveDFModel(cfg)
    assert m.halo.enclosed_mass([100])[0] == pytest.approx(2e6, rel=1e-12)
    r = np.geomspace(.02, 500, 30)
    force = m.static_potentials[-1].force(np.column_stack((r, r*0, r*0)))[:, 0]
    np.testing.assert_allclose(-force*r*r/agama.G, m.halo.enclosed_mass(r), rtol=.003)
    assert np.all(np.isfinite(m.projected_moments([.1, 10, 50])["los"]))


def test_nonconvergence_and_out_of_domain_are_errors(model):
    with pytest.raises(DFConvergenceError, match="did not converge"):
        PositiveDFModel(replace(model.config, numerics=DFNumerics(max_iterations=2)))
    for radius in (0, -1, np.nan, 1e5):
        with pytest.raises(ValueError):
            model.projected_moments([radius])


def test_changed_global_units_are_rejected(model, agama, monkeypatch):
    monkeypatch.setattr(agama, "getUnits", lambda: dict(length=1, mass=1, velocity=1))
    with pytest.raises(RuntimeError, match="units"):
        model.projected_moments([1])


def test_phase_space_density_bound_unbound_and_reflection(model):
    a = np.array([[1., 0, 0, 1., 2., 3.], [1., 0, 0, -1., -2., -3.],
                  [1., 0, 0, 1e5, 0, 0]])
    values = model.phase_space_density(a)
    assert values[0] > 0 and values[2] == 0
    assert values[0] == pytest.approx(values[1], rel=1e-12)
    with pytest.raises(ValueError):
        model.phase_space_density([[0., 0, 0, 0, 0, 0]])


def photo_data():
    return PhotometricData(np.array([1., 10., 100.]), np.array([10., 11., 13.]),
                           np.array([.1, .2, .3]), "synthetic", "known errors")


def test_photometric_zero_point_and_snapshot():
    p = photo_data()
    assert PhotometricData.from_dict(json.loads(json.dumps(p.to_dict()))).to_dict() == p.to_dict()
    intensity = 10**(-.4*p.mu)
    assert p.compare(intensity)["chi2"] < 1e-20
    assert p.compare(100*intensity)["chi2"] < 1e-20
    with pytest.raises(ValueError):
        p.compare([1, 0, 1])


def test_existing_binning_likelihood_and_infeasible_streaming(model):
    p = BinnedProfile("test", "pmt", np.array([10., 50.]), None, None,
                      np.array([.6, .5]), np.array([.01, .02]), np.array([.02, .03]),
                      "HST", streaming2=np.array([.01, .02]),
                      r_nodes=np.array([[8., 10., 12.], [40., 50., 60.]]))
    data = KinematicData((p,))
    result = DFJointProblem(data, photo_data()).evaluate(model.config, model)
    expected = ProfileLikelihood(data).predict(model, model.config.distance_kpc)
    np.testing.assert_allclose(result["predictions"][p.name], expected[p.name], rtol=1e-14)
    impossible = KinematicData((replace(p, streaming2=np.array([1e6, 1e6])),))
    with pytest.raises(ValueError, match="streaming exceeds"):
        DFJointProblem(impossible, photo_data()).evaluate(model.config, model)


def test_fit_coordinate_roundtrip_and_input_preservation():
    cfg = DFModelConfig()
    coordinates = [FitCoordinate("M_star", 1e6, 1e7, log=True),
                   FitCoordinate("components.0.g_r", .3, 2.8)]
    result = config_at(cfg, coordinates, [np.log(4e6), 2.])
    assert result.M_star == pytest.approx(4e6)
    assert result.components[0].g_r == 2.
    assert cfg.M_star == 3e6 and cfg.components[0].g_r == 1.
    with pytest.raises(ValueError):
        config_at(cfg, coordinates, [np.log(4e6), 3.])


def test_fitted_mixture_weights_stay_positive_and_normalized():
    cfg = DFModelConfig(components=(ActionComponent(fraction=.6),
                                   ActionComponent(fraction=.3, J0=500),
                                   ActionComponent(fraction=.1, J0=100)))
    coords = [FitCoordinate("mixture_log_ratio.1", -8., 8.),
              FitCoordinate("mixture_log_ratio.2", -8., 8.)]
    reconstructed = config_at(cfg, coords, [c.get(cfg) for c in coords])
    np.testing.assert_allclose([c.fraction for c in reconstructed.components], [.6, .3, .1])
    changed = config_at(cfg, coords, [3., -5.])
    weights = np.array([c.fraction for c in changed.components])
    assert np.all(weights > 0) and weights.sum() == pytest.approx(1.)
    assert weights[1]/weights[0] == pytest.approx(np.exp(3))
    assert weights[2]/weights[0] == pytest.approx(np.exp(-5))
    assert cfg.components[0].fraction == .6


def test_invalid_mixture_coordinates_rejected():
    for path in ("mixture_log_ratio.0", "mixture_log_ratio.-1", "components.0.fraction"):
        with pytest.raises(ValueError):
            FitCoordinate(path, -1., 1.)
    c = FitCoordinate("mixture_log_ratio.1", -8., 8.)
    with pytest.raises(ValueError, match="index"):
        c.get(DFModelConfig())
    with pytest.raises(ValueError, match="index"):
        config_at(DFModelConfig(), [c], [0.])
