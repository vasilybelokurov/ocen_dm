"""Rung-2 profiles: independent Jeans quadrature, limits and replay."""
import json

import numpy as np
import pytest
from scipy.integrate import quad

from ocen_dm.kinematics.anisotropy import Anisotropy, TurnoverAnisotropy
from ocen_dm.kinematics.fit import NoDarkMatterModel, DarkMatterModel
from ocen_dm.kinematics.jeans import SphericalJeans
from ocen_dm.kinematics.run_io import family_config, family_from_config
from ocen_dm.light_model import MGEFit
from ocen_dm.mass_models import Plummer
from ocen_dm.mass_models.base import G

MGE = MGEFit(np.array([50., 250.]), np.array([.7, .3]), 0., 0., 2)


def test_turnover_bounds_and_nonmonotonic_freedom():
    r = np.geomspace(1e-5, 1e6, 1000)
    for levels in [(-.75, .4, -.2), (-.6, -.9, .5), (-.8, -.2, .4)]:
        a = TurnoverAnisotropy(*levels, 1., 29.)
        b = a.beta(r)
        assert b.min() >= min(levels) and b.max() <= max(levels)
        assert a.beta(0.) == levels[0]
        assert b[-1] == pytest.approx(levels[-1], abs=1e-8)
    peak = TurnoverAnisotropy(-.75, .4, -.2, 1., 29.).beta(r)
    assert peak.max() > .3 and np.any(np.diff(peak) > 0) and np.any(np.diff(peak) < 0)


@pytest.mark.parametrize("levels", [(-.75, .4, -.2), (-.6, -.9, .5)])
def test_integrating_factor_matches_direct_quadrature(levels):
    a = TurnoverAnisotropy(*levels, .1, 29.9)
    for radius in [.002, .04, 1., 20., 300.]:
        expected = quad(lambda u: 2*float(a.beta(np.exp(u))), 0, np.log(radius), epsabs=1e-10)[0]
        assert a.log_integrating_factor(radius)-a.log_integrating_factor(1.) == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("levels", [(-.5, -.5, -.5), (-.8, .2, .2), (-.8, -.8, .2)])
def test_constant_and_single_transition_limits_reproduce_projections(levels):
    a = TurnoverAnisotropy(*levels, 1., 19.)
    radius = 20. if levels[0] == levels[1] else 1.
    legacy = Anisotropy(levels[0], levels[-1], radius)
    r = np.geomspace(.01, 100., 30)
    np.testing.assert_allclose(a.beta(r), legacy.beta(r), atol=1e-15)
    p = Plummer(3e6, 5.)
    actual = SphericalJeans(p, p, a).dispersions_kms(r)
    expected = SphericalJeans(p, p, legacy).dispersions_kms(r)
    for key in actual:
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-12)


def test_turnover_jeans_solution_matches_independent_integral():
    p = Plummer(3e6, 5.)
    a = TurnoverAnisotropy(-.75, .4, -.3, 1., 25.)
    model = SphericalJeans(p, p, a)
    for r in [.1, 1., 5., 20., 70.]:
        # Numerically integrate beta itself, independently of the closed form.
        def integrand(u):
            s = np.exp(u)
            log_ratio = quad(lambda v: 2*float(a.beta(np.exp(v))), np.log(r), u, epsabs=1e-9)[0]
            return np.exp(log_ratio)*p.density(s).item()*G*p.enclosed_mass(s).item()/s
        expected = quad(integrand, np.log(r), np.log(1e4), epsabs=1e-8, epsrel=1e-8)[0]
        assert model.nu_sigma_r2(r)[0] == pytest.approx(expected, rel=2e-4)


@pytest.mark.parametrize("name", ["beta_0", "beta_mid", "beta_inf", "r_beta", "delta_r_beta"])
@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_nonfinite_parameters_rejected(name, value):
    with pytest.raises(ValueError):
        TurnoverAnisotropy(**{name: value})


@pytest.mark.parametrize("kw", [dict(beta_mid=1.), dict(r_beta=0.), dict(delta_r_beta=-1.)])
def test_nonphysical_levels_or_unordered_scales_rejected(kw):
    with pytest.raises(ValueError):
        TurnoverAnisotropy(**kw)


@pytest.mark.parametrize("cls", [NoDarkMatterModel, DarkMatterModel])
def test_saved_turnover_model_replays_exactly(cls):
    family = cls(mge_fit=MGE, instruments=(), anisotropy_profile="turnover")
    assert len(family.names) == (10 if cls is NoDarkMatterModel else 12)
    assert family.parameters[family.names.index("beta_0")].prior.hi == -.5
    config = json.loads(json.dumps(family_config(family)))
    restored = family_from_config(config)
    assert family_config(restored) == config
    x = family.transform(np.full(len(family.names), .4))
    before = family.build(family.to_dict(x))[0]
    after = restored.build(restored.to_dict(x))[0]
    for key, values in before.dispersions_kms(np.array([.1, 3., 20., 40.])).items():
        np.testing.assert_array_equal(values, after.dispersions_kms(np.array([.1, 3., 20., 40.]))[key])
    with pytest.raises(ValueError, match="jeans"):
        family_from_config(config, backend="jam")


def test_legacy_snapshot_without_profile_flag_remains_monotonic():
    family = NoDarkMatterModel(mge_fit=MGE, instruments=())
    config = family_config(family)
    config.pop("anisotropy_profile")
    assert family_from_config(config).anisotropy_profile == "monotonic"


@pytest.mark.parametrize("kw", [dict(constant_beta=True), dict(backend="jam"), dict(backend="agama")])
def test_incompatible_family_configuration_rejected(kw):
    with pytest.raises(ValueError, match="turnover"):
        NoDarkMatterModel(mge_fit=MGE, anisotropy_profile="turnover", **kw)


def test_cli_turnover_flag():
    from ocen_dm.cli import build_parser
    args = build_parser().parse_args(["fit", "--anisotropy-profile", "turnover", "--no-scales"])
    assert args.anisotropy_profile == "turnover" and not args.constant_beta
