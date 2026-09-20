"""Tests for the model families, priors and fitting driver."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics import (
    Anisotropy, DarkMatterModel, FitProblem, KinematicData, NoDarkMatterModel, Prior,
    SphericalJeans, maximum_likelihood,
)
from ocen_dm.kinematics.likelihood import BinnedProfile
from ocen_dm.light_model import MGEFit
from ocen_dm.mass_models import CompositeMassModel, Plummer


# A tiny two-Gaussian "light model" so the tests need no data files.
_MGE = MGEFit(sigmas_arcsec=np.array([60.0, 250.0]), fractions=np.array([0.7, 0.3]),
              rms_mag=0.0, max_abs_resid_mag=0.0, n_grid=2)


# ----------------------------------------------------------------- priors ---
def test_prior_transforms_hit_bounds_and_are_monotonic():
    u = np.linspace(1e-6, 1 - 1e-6, 101)
    p = Prior("uniform", -1.0, 2.0)
    assert p.transform(0.0) == -1.0 and p.transform(1.0) == 2.0
    p = Prior("loguniform", 1e2, 1e6)
    assert p.transform(0.0) == pytest.approx(1e2) and p.transform(1.0) == pytest.approx(1e6)
    assert p.transform(0.5) == pytest.approx(1e4)
    p = Prior("normal", mu=5.43, sigma=0.05)
    assert p.transform(0.5) == pytest.approx(5.43)
    assert p.transform(0.8413) == pytest.approx(5.48, abs=1e-3)      # +1 sigma
    p = Prior("truncnormal", 0.0, 1.0, mu=0.0, sigma=0.5)
    x = p.transform(u)
    assert np.all(np.diff(x) > 0) and x[0] >= 0.0 and x[-1] <= 1.0


def test_prior_validation():
    with pytest.raises(ValueError):
        Prior("flat", 0, 1)
    with pytest.raises(ValueError):
        Prior("uniform", 1, 0)
    with pytest.raises(ValueError):
        Prior("loguniform", 0, 1)
    with pytest.raises(ValueError):
        Prior("normal", sigma=0)


# --------------------------------------------------------------- families ---
def test_k1_parameter_vector_and_build():
    fam = NoDarkMatterModel(mge_fit=_MGE)
    assert fam.names[:7] == ("M_star", "M_rem", "a_rem", "M_bh", "beta_0", "beta_inf", "r_beta")
    assert fam.names[7:] == ("s_MUSE", "s_GaiaDR2", "s_GaiaEDR3", "distance")   # distance is never fixed by default
    x = fam.transform(np.full(len(fam.names), 0.5))
    theta = fam.to_dict(x)
    # An & Evans with a central point mass: the prior is [-1, -1/2], so the midpoint of
    # the unit cube maps to -0.75 (narrowed from [-1, 0] on 2026-09-20)
    assert theta["beta_0"] == pytest.approx(-0.75)
    jeans, D, scales = fam.build(theta)
    assert D == pytest.approx(5.43) and scales == {"MUSE": 1.0, "GaiaDR2": 1.0, "GaiaEDR3": 1.0}
    fixed = NoDarkMatterModel(mge_fit=_MGE, fix_distance=True, instruments=("GaiaDR2", "GaiaEDR3"))
    assert fixed.names[7:] == ("s_GaiaDR2", "s_GaiaEDR3")
    assert isinstance(jeans, SphericalJeans)
    total = float(jeans.mass.enclosed_mass(1e4)[0])
    assert total == pytest.approx(theta["M_star"] + theta["M_rem"] + theta["M_bh"], rel=1e-2)


def test_k1_distance_becomes_a_parameter_when_given_a_prior():
    fam = NoDarkMatterModel(mge_fit=_MGE, distance_prior=Prior("normal", mu=5.43, sigma=0.05))
    assert fam.names[-1] == "distance"
    assert NoDarkMatterModel(mge_fit=_MGE, fix_distance=True).names[-1] != "distance"
    theta = fam.to_dict(fam.transform(np.full(len(fam.names), 0.5)))
    assert fam.build(theta)[1] == pytest.approx(5.43)


def test_k2_halo_mass_parametrisation_is_exact():
    for gamma in (0.0, 1.0):
        fam = DarkMatterModel(gamma=gamma, mge_fit=_MGE)
        assert "M_dm_100" in fam.names and "r_s" in fam.names
        theta = fam.to_dict(fam.transform(np.full(len(fam.names), 0.37)))
        halo = fam.halo(theta)
        assert float(np.asarray(halo.enclosed_mass(100.0)).reshape(-1)[0]) == pytest.approx(theta["M_dm_100"], rel=1e-6)
        assert halo.gamma == gamma
        jeans, _, _ = fam.build(theta)
        assert any(c.name.startswith("dm") for c in jeans.mass.components)
    assert DarkMatterModel(gamma=0.0, mge_fit=_MGE).label == "K2_cored"
    assert DarkMatterModel(gamma=1.0, mge_fit=_MGE).label == "K2_nfw"


# ----------------------------------------------------------------- problem ---
def _synthetic_data(fam, x_true, rng, rel_err=0.03):
    """Bins on 3 to 800 arcsec, three projections plus one no-edge Gaia-like profile."""
    e = np.geomspace(3.0, 800.0, 13)
    lo, hi, mid = e[:-1], e[1:], np.sqrt(e[:-1] * e[1:])
    one = np.ones_like(mid)
    template = KinematicData((
        BinnedProfile("hst_pm_radial", "pmr", mid, lo, hi, one, one, one, "HST"),
        BinnedProfile("hst_pm_tangential", "pmt", mid, lo, hi, one, one, one, "HST"),
        BinnedProfile("muse_los_dispersion", "los", mid, lo, hi, one, one, one, "MUSE"),
        BinnedProfile("gaia_edr3_pm", "pmc", mid[6:], None, None, one[6:], one[6:], one[6:], "GaiaEDR3"),
    ))
    P0 = FitProblem(fam, template)
    pred = P0.predict(x_true)
    profiles = []
    for p in template.profiles:
        err = rel_err * pred[p.name]
        from dataclasses import replace
        profiles.append(replace(p, value=pred[p.name] + rng.normal(0, err), err_lo=err, err_hi=err))
    return KinematicData(tuple(profiles))


_K1_TEST = dict(mge_fit=_MGE, fix_distance=True, instruments=("GaiaDR2", "GaiaEDR3"))


def test_problem_returns_finite_or_floor_never_raises():
    fam = NoDarkMatterModel(**_K1_TEST)
    rng = np.random.default_rng(0)
    data = _synthetic_data(fam, fam.transform(np.full(9, 0.5)), rng)
    P = FitProblem(fam, data)
    assert np.isfinite(P.loglike_unit(np.full(9, 0.5)))
    assert P.loglike_unit(np.full(9, 1e-9)) <= 0.0          # corner of the cube evaluates, no exception
    assert P.n_calls == 2


def test_mock_data_keeps_bins_and_errors():
    fam = NoDarkMatterModel(**_K1_TEST)
    rng = np.random.default_rng(3)
    x = fam.transform(np.full(9, 0.5))
    data = _synthetic_data(fam, x, rng)
    P = FitProblem(fam, data)
    mock = P.mock_data(x, np.random.default_rng(4))
    for a, b in zip(data.profiles, mock.profiles):
        assert a.name == b.name and np.array_equal(a.r, b.r) and np.array_equal(a.err_lo, b.err_lo)
        assert not np.array_equal(a.value, b.value)


def test_maximum_likelihood_recovers_injected_stellar_mass_and_scale():
    """Injection: with M_bh at the prior floor and no remnant weight, the multi-start
    optimiser must find a model whose -2 ln L is within a few units of the truth's
    (the likelihood surface is degenerate in M_rem/a_rem/M_bh, so only the total
    dynamical mass and the Gaia scale are checked as recovered quantities)."""
    fam = NoDarkMatterModel(**_K1_TEST)
    rng = np.random.default_rng(11)
    theta_true = {"M_star": 3.0e6, "M_rem": 1.5e5, "a_rem": 1.5, "M_bh": 1e2,
                  "beta_0": -0.7, "beta_inf": 0.2, "r_beta": 10.0, "s_GaiaDR2": 1.0, "s_GaiaEDR3": 1.15}
    x_true = np.array([theta_true[n] for n in fam.names])
    data = _synthetic_data(fam, x_true, rng, rel_err=0.02)
    P = FitProblem(fam, data)
    x_ml, ln_ml = maximum_likelihood(P, n_starts=4, seed=1, maxiter=600)
    ln_true = P.loglike_vector(x_true)
    assert ln_ml >= ln_true - 3.0                                 # optimiser at least as good as truth (allow noise)
    th = fam.to_dict(x_ml)
    m_true = theta_true["M_star"] + theta_true["M_rem"]
    assert (th["M_star"] + th["M_rem"] + th["M_bh"]) == pytest.approx(m_true, rel=0.05)
    assert th["s_GaiaEDR3"] == pytest.approx(1.15, abs=0.05)


def test_backend_option_changes_parametrisation_and_label():
    fam = NoDarkMatterModel(mge_fit=_MGE, backend="agama", fix_distance=True, instruments=())
    assert fam.names == ("M_star", "M_rem", "a_rem", "M_bh", "beta_0", "r_a") and fam.label.endswith("_agama")
    jam = NoDarkMatterModel(mge_fit=_MGE, backend="jam", fix_distance=True, instruments=())
    assert jam.names == NoDarkMatterModel(mge_fit=_MGE, fix_distance=True, instruments=()).names
    with pytest.raises(ValueError):
        NoDarkMatterModel(mge_fit=_MGE, backend="galpy")


def test_jam_and_jeans_backends_agree_on_a_k1_vector():
    pytest.importorskip("jampy")
    fam_j = NoDarkMatterModel(mge_fit=_MGE, fix_distance=True, instruments=())
    fam_m = NoDarkMatterModel(mge_fit=_MGE, fix_distance=True, instruments=(), backend="jam")
    x = fam_j.transform(np.full(len(fam_j.names), 0.5))
    R = np.geomspace(0.5, 20.0, 6)
    a = fam_j.build(fam_j.to_dict(x))[0].dispersions_kms(R)
    b = fam_m.build(fam_m.to_dict(x))[0].dispersions_kms(R)
    for k in ("los", "pmr", "pmt"):
        np.testing.assert_allclose(b[k], a[k], rtol=3e-3)
