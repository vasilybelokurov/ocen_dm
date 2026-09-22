"""The simple companion keeps the rung-1 problem and its constant-beta limit."""
from dataclasses import replace
from pathlib import Path
import runpy

import numpy as np
import pytest

from ocen_dm.kinematics.fit import DarkMatterModel, NoDarkMatterModel, Prior
from ocen_dm.kinematics.run_io import family_config, family_from_config
from ocen_dm.light_model import MGEFit

DRIVER = runpy.run_path(str(Path(__file__).resolve().parents[1] / "bin/run_rung2_fits.py"))
make_family = DRIVER["companion_family"]
MGE = MGEFit(np.array([50., 250.]), np.array([.7, .3]), 0., 0., 2)


def parent(gamma=None, **kw):
    cls = NoDarkMatterModel if gamma is None else DarkMatterModel
    return cls(mge_fit=MGE, instruments=(), constant_beta=True, tracer="composite",
               **({} if gamma is None else dict(gamma=gamma, r_t=800.)), **kw)


@pytest.mark.parametrize("gamma", [None, 0., 1.])
def test_simple_priors_and_saved_model_replay(gamma):
    old = parent(gamma)
    simple = make_family(old, "simple")
    assert len(simple.names) == (8 if gamma is None else 10)
    assert set(simple.names) - set(old.names) == {"beta_inf", "r_beta"}
    priors = {p.name: p.prior for p in simple.parameters}
    assert priors["beta_0"] == priors["beta_inf"] == Prior("uniform", -1., .5)
    assert priors["r_beta"] == Prior("loguniform", .5, 100.)
    for p in old.parameters:
        assert priors[p.name] == p.prior
    restored = family_from_config(family_config(simple))
    assert family_config(restored) == family_config(simple)
    x = simple.transform(np.linspace(.15, .85, len(simple.names)))
    np.testing.assert_array_equal(restored.transform(np.linspace(.15, .85, len(simple.names))), x)
    before, _, _ = simple.build(simple.to_dict(x))
    after, _, _ = restored.build(restored.to_dict(x))
    for component in ("los", "pmr", "pmt"):
        np.testing.assert_array_equal(before.projected_moments(np.array([.1, 1., 10., 40.]))[component],
                                      after.projected_moments(np.array([.1, 1., 10., 40.]))[component])


@pytest.mark.parametrize("gamma", [None, 0., 1.])
def test_constant_limit_preserves_all_three_projected_moments(gamma):
    old = parent(gamma)
    simple = make_family(old, "simple")
    theta = old.to_dict(old.transform(np.full(len(old.names), .5)))
    radii = np.array([.1, 1., 5., 20., 40.])
    for beta in (-.7, .1, .45):
        theta["beta_0"] = beta
        reference, _, _ = old.build(theta)
        expected = reference.projected_moments(radii)
        for scale in (.5, 10., 100.):
            trial, _, _ = simple.build({**theta, "beta_inf": beta, "r_beta": scale})
            actual = trial.projected_moments(radii)
            for component in ("los", "pmr", "pmt"):
                np.testing.assert_allclose(actual[component], expected[component], rtol=1e-12, atol=1e-12)


def test_custom_parent_mass_distance_and_fixed_values_are_preserved():
    old = parent(0., fixed={"a_rem": 2.}, distance_kpc=5.2,
                 distance_prior=Prior("normal", mu=5.2, sigma=.08))
    old.parameters = tuple(replace(p, prior=Prior("uniform", 1e4, 9e5))
                           if p.name == "M_rem" else p for p in old.parameters)
    simple = make_family(old, "simple")
    assert simple.fixed == old.fixed
    assert simple.fixed_distance_kpc == old.fixed_distance_kpc
    assert simple.distance_prior == old.distance_prior
    assert simple.gamma == old.gamma and simple.r_t == old.r_t
    for p in old.parameters:
        assert next(q for q in simple.parameters if q.name == p.name).prior == p.prior


def test_turnover_recipe_retains_its_original_priors():
    family = make_family(parent(), "turnover")
    reference = NoDarkMatterModel(mge_fit=MGE, tracer="composite", instruments=(),
                                  anisotropy_profile="turnover")
    assert family_config(family) == family_config(reference)


def test_launch_refuses_existing_batch_and_fit_before_starting_any_worker(tmp_path, monkeypatch):
    launch = DRIVER["launch"]
    monkeypatch.setitem(launch.__globals__, "ROOT", tmp_path)
    with pytest.raises(FileExistsError, match="Batch already exists"):
        launch(tmp_path, "simple")
    existing = tmp_path / "results/fits/rung2_K2_nfw_simple_beta"
    existing.mkdir(parents=True)
    marker = existing / "untouched.txt"
    marker.write_text("preserve")
    with pytest.raises(FileExistsError, match="Fit already exists"):
        launch(tmp_path / "new_batch", "simple")
    assert marker.read_text() == "preserve"
    assert not (tmp_path / "new_batch").exists()


def test_unexpected_parent_beta_prior_fails_explicitly():
    old = parent()
    old.parameters = tuple(replace(p, prior=Prior("uniform", -1., -.5))
                           if p.name == "beta_0" else p for p in old.parameters)
    with pytest.raises(ValueError, match="agreed rung-1 beta prior"):
        make_family(old, "simple")
