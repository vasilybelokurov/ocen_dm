"""Analytical and independent numerical controls for the compact stellar DF."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.special import gamma

from ocen_dm.kinematics.df_fit import FitCoordinate, config_at, model_config_from_dict
from ocen_dm.kinematics.regularized_df import (
    ExactSphericalActions,
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
    # A point mass is supported in the exact mode only; the epicycle reference assumes a finite centre.
    with pytest.raises(ValueError, match="finite centre"):
        SphericalFrequencyRatio(agama.Potential(type="Plummer", mass=1e4, scaleRadius=0.), mode="epicycle")


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


def test_frequency_representative_preserves_angular_momentum(agama, monkeypatch):
    # Original polar mapping reached R=0 with nonfinite velocities. Exporting
    # the potential softens this to a 4e-4 angular-momentum error, still wrong.
    fixture = Path(__file__).parent/"fixtures/compact_frequency"
    record = json.loads((fixture/"orbit.json").read_text())
    jr, angular = record["action_angles"][:2]
    pot = agama.Potential(str(fixture/"potential.ini"))
    original = agama.actions
    recovered = []

    def capture(potential, xv, **kwargs):
        assert np.all(np.isfinite(xv))
        result = original(potential, xv, **kwargs)
        recovered.append(result[0])
        return result

    monkeypatch.setattr(agama, "actions", capture)
    ratio = SphericalFrequencyRatio(pot)(jr, angular)
    action = np.concatenate(recovered)
    np.testing.assert_allclose(action[:, 0], jr, rtol=2e-6)
    np.testing.assert_allclose(action[:, 1]+abs(action[:, 2]), angular, rtol=1e-6)
    assert ratio == pytest.approx(2., rel=2e-5)


# ----------------------------------------------------------------------------- central point mass
def test_point_mass_kepler_frequencies_contours_and_isotropy(agama):
    """Pure Kepler: Omega_r = Omega_t, H depends on Jr+L only, so ergodic contours are L_c = Jr+L."""
    p = agama.Potential(type="Plummer", mass=3e6, scaleRadius=0.)
    ratio = SphericalFrequencyRatio(p)
    assert ratio.point_mass and ratio.radial_limit == 1.
    rng = np.random.default_rng(7)
    jr, l = 10**rng.uniform(-2, 3, (2, 60))
    np.testing.assert_allclose(ratio(jr, l), 1., rtol=3e-8)   # ~1e-8 rounding floor on near-circular orbits
    np.testing.assert_allclose(ratio(jr, 0*l), 1.)
    df = RegularizedExponentialDF(3e6, ExponentialParameters(b_out=0), potential=p)
    np.testing.assert_allclose(df.contour_label(jr, l), jr+l, rtol=2e-5)
    np.testing.assert_allclose(df.direct_contour_label(jr, l), jr+l, rtol=1e-6)
    rho, pr, pt = spherical_velocity_moments(p, df, agama.ActionFinder(p), [.01, .1, 1, 10], 64).T
    np.testing.assert_allclose(pt/pr, 1., atol=1e-4)


def test_point_mass_sets_the_radial_orbit_limit(agama):
    """With a finite core plus a point mass, near-radial orbits swing 2 pi at pericentre: ratio -> 1, not 2."""
    core = agama.Potential(type="Isochrone", mass=3e6, scaleRadius=5.)
    both = agama.Potential(core, agama.Potential(type="Plummer", mass=1e4, scaleRadius=0.))
    assert SphericalFrequencyRatio(core).radial_limit == 2.
    ratio = SphericalFrequencyRatio(both)
    assert ratio.radial_limit == 1.
    jr = np.array([1., 10., 100., 1000.])
    small = ratio(jr, 1e-6*jr)
    np.testing.assert_allclose(small, 1., atol=2e-2)
    np.testing.assert_allclose(ratio(jr, 0*jr), 1.)


def test_uniform_df_moments_around_a_point_mass(agama):
    """f ~ const near a point mass (q << J0): rho = 4 pi f v_esc^3/3, sigma_r^2 = v_esc^2/5 (escape-shell accuracy)."""
    gm = agama.G*1e4
    p = agama.Potential(type="Plummer", mass=1e4, scaleRadius=0.)
    df = RegularizedExponentialDF(1., ExponentialParameters(J0=1e9, alpha=1., b_out=0), potential=p)
    r = np.array([.01, .1, 1.])
    rho, pr, pt = spherical_velocity_moments(p, df, agama.ActionFinder(p), r, 32).T
    vesc = np.sqrt(2*gm/r)
    f0 = df.amplitude
    np.testing.assert_allclose(rho, 4*np.pi*f0*vesc**3/3, rtol=2e-3)
    np.testing.assert_allclose(pr/rho, vesc**2/5, rtol=2e-3)
    np.testing.assert_allclose(pt/pr, 1., atol=2e-3)


@pytest.mark.parametrize("b_out,b_outer", [(2., None), (-2., None), (.8, -.6), (-1.5, 1.5)])
def test_contour_label_lower_bound_used_for_underflow(agama, b_out, b_outer):
    """Lc >= q exp(1-e^B), |bias| <= B: the bound that keeps near-escape point-mass orbits out of the mapper."""
    core = agama.Potential(type="Isochrone", mass=3e6, scaleRadius=5.)
    both = agama.Potential(core, agama.Potential(type="Plummer", mass=4e4, scaleRadius=0.))
    p = ExponentialParameters(b_out=b_out, J_a=50., b_outer=b_outer, J_outer=None if b_outer is None else 500.)
    bound = abs(b_out)+(abs(b_outer-b_out) if b_outer is not None else 0.)
    for pot in (core, both):
        df = RegularizedExponentialDF(3e6, p, potential=pot)
        q = np.geomspace(1e-2, 1e5, 15)[:, None]
        c = np.array([0., 1e-3, .1, .5, .9])[None, :]
        lc = df.contour_label(q*(1-c), q*c+0*q)
        assert np.all(lc >= q*np.exp(1-np.exp(bound)))
    huge = np.array([[2.6e12, 0., 7.9e9], [3e20, 0., 1e-35]])
    np.testing.assert_array_equal(df(huge), 0.)


def test_prescribed_point_mass_config_and_potential(agama):
    from ocen_dm.kinematics.positive_df import DFNumerics
    n = DFNumerics()
    assert len(PrescribedMatter().potentials(n)) == 0
    pots = PrescribedMatter(M_bh=1e4).potentials(n)
    assert len(pots) == 1
    force = -pots[0].force([2., 0., 0.])[0]
    assert force == pytest.approx(agama.G*1e4/4, rel=1e-12)
    with pytest.raises(ValueError):
        PrescribedMatter(M_bh=-1.)
    cfg = RegularizedDFConfig(matter=PrescribedMatter(M_rem=2e5, a_rem=2., M_bh=8200.))
    back = RegularizedDFConfig.from_dict(json.loads(json.dumps(cfg.to_dict())))
    assert back.matter.M_bh == 8200.
    legacy = cfg.to_dict()
    del legacy["matter"]["M_bh"]           # configurations written before the point mass existed
    assert RegularizedDFConfig.from_dict(legacy).matter.M_bh == 0.


def test_mass_profiles_include_the_point_mass(agama):
    from types import SimpleNamespace
    from ocen_dm.kinematics.compact_recovery import mass_profiles
    star = agama.Potential(type="Plummer", mass=3e6, scaleRadius=5.)
    model = SimpleNamespace(stellar_potential=star, config=SimpleNamespace(matter=PrescribedMatter(M_bh=8200.)),
                            diagnostics=dict(stellar_mass=3e6))
    out = mass_profiles(model, [.1, 1., 10.])
    np.testing.assert_allclose(out["bh"], 8200.)
    np.testing.assert_allclose(np.array(out["total"]), np.array(out["stars"])+8200., rtol=1e-12)


def test_apsidal_quadrature_matches_the_exact_isochrone_ratio(agama):
    """The point-mass frequency path (direct apsidal angle) reproduces the analytic isochrone ratio."""
    p = agama.Potential(type="Isochrone", mass=3e6, scaleRadius=5.)
    ratio = SphericalFrequencyRatio(p)
    rng = np.random.default_rng(11)
    jr, l = 10**rng.uniform(-2, 3, (2, 128))
    xv = ratio.mapper(np.column_stack((jr, 0*jr, l, np.ones((len(jr), 3)))))
    r0 = np.linalg.norm(xv[:, :3], axis=1)
    energy = ratio._phi_radial(r0)+.5*np.sum(xv[:, 3:]**2, axis=1)
    k = agama.G*3e6*5
    np.testing.assert_allclose(ratio._apsidal_ratio(energy, l, r0), 2/(1+l/np.sqrt(l*l+4*k)), rtol=1e-7)


def test_radial_action_inversion_matches_the_exact_isochrone_energy(agama):
    """Jr(E, L) inversion (fallback where the mapper fails) vs H = -(GM)^2/(2 [Jr + (L + sqrt(L^2+4GMb))/2]^2)."""
    gm, b = agama.G*3e6, 5.
    p = agama.Potential(type="Isochrone", mass=3e6, scaleRadius=b)
    ratio = SphericalFrequencyRatio(p)
    rng = np.random.default_rng(12)
    jr, l = 10**rng.uniform(-2, 9, 96), 10**rng.uniform(-3, 7, 96)
    exact = -.5*(gm/(jr+.5*(l+np.sqrt(l*l+4*gm*b))))**2
    r_circ = np.atleast_1d(p.Rcirc(L=l))
    np.testing.assert_allclose(ratio._energy_from_actions(jr, l, r_circ), exact, rtol=1e-9)


def test_radial_action_inversion_is_exact_in_kepler(agama):
    """Jr + L = GM/sqrt(-2E) for all eccentricities, including near-radial orbits."""
    gm = agama.G*8200.
    ratio = SphericalFrequencyRatio(agama.Potential(type="Plummer", mass=8200., scaleRadius=0.))
    q = np.geomspace(1e-2, 1e8, 21)
    for c in (1e-4, 1e-3, .3, .95):
        jr, l = q*(1-c), q*c
        np.testing.assert_allclose(ratio._energy_from_actions(jr, l, l**2/gm), -.5*(gm/q)**2, rtol=1e-11)


def test_point_mass_frequencies_survive_mapper_failure(agama, monkeypatch):
    """Around a point mass AGAMA's mapper can return NaN at large Jr: those rows use the Jr(E, L) inversion,
    which AGAMA's (independent) action finder confirms, and the ratio matches the mapper path where both work."""
    multi = agama.Potential(type="Multipole", density=agama.Density(type="Plummer", mass=3e6, scaleRadius=3.),
                            lmax=0, gridSizeR=80, rmin=1e-4, rmax=1e5)
    both = agama.Potential(multi, agama.Potential(type="Plummer", mass=2.7e5, scaleRadius=2.2),
                           agama.Potential(type="Plummer", mass=8200., scaleRadius=0.))
    ratio = SphericalFrequencyRatio(both)
    q = np.geomspace(1e-2, 1e8, 33)
    for c in (1e-3, .3, .95):
        jr, l = q*(1-c), q*c
        r_circ = np.atleast_1d(both.Rcirc(L=l))
        energy = ratio._energy_from_actions(jr, l, r_circ)
        vr = np.sqrt(np.maximum(2*(energy-ratio._phi_radial(r_circ))-(l/r_circ)**2, 0.))
        xv = np.column_stack((r_circ, 0*q, 0*q, vr, l/r_circ, 0*q))
        # AGAMA's action finder is itself off at huge actions and (4.5e-4 in exact Kepler) at c = 1e-3
        if c > .1:
            small = q < 1e5
            np.testing.assert_allclose(agama.actions(both, xv[small])[:, 0], jr[small], rtol=1e-6)
        reference = ratio(jr, l)
        mapper = ratio.mapper
        monkeypatch.setattr(ratio, "mapper", lambda aa: np.where(np.arange(len(aa))[:, None] % 2 == 0, np.nan, mapper(aa)))
        mixed = ratio(jr, l)
        monkeypatch.setattr(ratio, "mapper", mapper)
        assert np.all(np.isfinite(mixed)) and np.all((mixed > 1-1e-9) & (mixed < 2))
        np.testing.assert_allclose(mixed, reference, atol=1e-4)


@pytest.mark.parametrize("kind", ["isochrone", "kepler"])
def test_exact_spherical_actions_match_analytic_actions(agama, kind):
    """Jr = GM/sqrt(-2E) - (L + sqrt(L^2 + 4 GM b))/2 (b = 0: Kepler), incl. near-radial and near-escape points."""
    gm, b = agama.G*3e6, (3. if kind == "isochrone" else 0.)
    p = (agama.Potential(type="Isochrone", mass=3e6, scaleRadius=b) if b else
         agama.Potential(type="Plummer", mass=3e6, scaleRadius=0.))
    rng = np.random.default_rng(5)
    n = 400
    r = 10**rng.uniform(-3, 4, n)
    vesc = np.sqrt(-2*p.potential(np.column_stack((r, 0*r, 0*r))))
    v = vesc*np.concatenate([rng.uniform(0, 1, n-50), 1-10**rng.uniform(-6, -2, 50)])
    mu = np.concatenate([rng.uniform(-1, 1, n-100), 1-10**rng.uniform(-8, -3, 100)])
    pos = r[:, None]*np.array([[0.6, 0., 0.8]])
    vel = v[:, None]*np.column_stack((0.6*mu, np.sqrt(1-mu**2), 0.8*mu))
    actions = ExactSphericalActions(p)(np.column_stack((pos, vel)))
    energy = p.potential(pos)+.5*v**2
    L = r*v*np.sqrt(1-mu**2)
    jr = gm/np.sqrt(-2*energy)-.5*(L+np.sqrt(L*L+4*gm*b))
    np.testing.assert_allclose(actions[:, 0], jr, rtol=1e-8, atol=1e-12*gm/np.sqrt(-2*energy).min())
    np.testing.assert_allclose(actions[:, 1]+np.abs(actions[:, 2]), L, rtol=1e-12)
