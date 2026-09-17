"""Spherical Jeans solver: closed forms, identities, JamPy and AGAMA cross-checks."""

from __future__ import annotations

import io
import contextlib

import numpy as np
import pytest

from ocen_dm.kinematics import Anisotropy, SphericalJeans
from ocen_dm.kinematics.jeans import KMS_PER_MASYR_KPC
from ocen_dm.mass_models import MGE, CompositeMassModel, Plummer, PointMass, TruncatedGNFW
from ocen_dm.mass_models.base import G

R = np.geomspace(0.05, 60.0, 12)


# ------------------------------------------------------------- closed forms --
def test_isotropic_plummer_radial_dispersion():
    p = Plummer(3.0e6, 5.0)
    j = SphericalJeans(p, p, Anisotropy(0.0, 0.0, 10.0))
    expected = np.sqrt(G * 3.0e6 / (6.0 * np.sqrt(R**2 + 25.0)))
    assert np.allclose(j.sigma_r(R), expected, rtol=1e-4)


def test_isotropic_plummer_projected_dispersion():
    """sigma_p^2(R) = 3 pi G M / (64 sqrt(R^2 + a^2)) (Dejonghe 1987)."""
    p = Plummer(3.0e6, 5.0)
    j = SphericalJeans(p, p, Anisotropy(0.0, 0.0, 10.0))
    expected = np.sqrt(3 * np.pi * G * 3.0e6 / (64.0 * np.sqrt(R**2 + 25.0)))
    d = j.dispersions_kms(R)
    for k in ("los", "pmr", "pmt"):
        assert np.allclose(d[k], expected, rtol=1e-4), k


def test_isotropy_makes_the_three_projections_identical():
    p = Plummer(1.0e6, 3.0)
    d = SphericalJeans(p, p, Anisotropy()).dispersions_kms(R)
    assert np.allclose(d["pmr"], d["los"]) and np.allclose(d["pmt"], d["los"])


def test_surface_density_of_a_plummer_is_the_closed_form():
    p = Plummer(1.0e6, 3.0)
    j = SphericalJeans(p, p)
    expected = 1.0e6 * 9.0 / (np.pi * (R**2 + 9.0) ** 2)
    assert np.allclose(j.surface_density(R), expected, rtol=1e-4)


# --------------------------------------------------------------- anisotropy --
def test_radial_anisotropy_orders_the_projections():
    """beta > 0 outside: sigma_pmr > sigma_los > sigma_pmt at large R (Binney & Mamon)."""
    p = Plummer(3.0e6, 5.0)
    d = SphericalJeans(p, p, Anisotropy(0.0, 0.6, 5.0)).dispersions_kms(np.array([20.0, 40.0]))
    assert np.all(d["pmr"] > d["los"]) and np.all(d["los"] > d["pmt"])


def test_tangential_anisotropy_reverses_the_order():
    p = Plummer(3.0e6, 5.0)
    d = SphericalJeans(p, p, Anisotropy(0.0, -1.0, 5.0)).dispersions_kms(np.array([20.0, 40.0]))
    assert np.all(d["pmt"] > d["los"]) and np.all(d["los"] > d["pmr"])


def test_sigma_t_over_sigma_r_is_one_minus_beta():
    p = Plummer(3.0e6, 5.0)
    a = Anisotropy(-0.3, 0.4, 6.0)
    j = SphericalJeans(p, p, a)
    ratio2 = (j.sigma_t(R) / j.sigma_r(R)) ** 2
    assert np.allclose(ratio2, 1.0 - a.beta(R), rtol=1e-10)


def test_anisotropy_profile_and_integrating_factor():
    a = Anisotropy(-0.2, 0.5, 4.0)
    assert a.beta(0.0) == pytest.approx(-0.2)
    assert a.beta(1e6) == pytest.approx(0.5, abs=1e-9)
    assert a.beta(4.0) == pytest.approx(0.15)                        # halfway at r = r_beta
    r = np.array([1.0, 2.0, 5.0])
    numeric = np.array([2 * __import__("scipy").integrate.quad(lambda t: a.beta(t) / t, 1.0, x)[0] for x in r])
    analytic = a.log_integrating_factor(r) - a.log_integrating_factor(1.0)
    assert np.allclose(numeric, analytic, atol=1e-8)


def test_anisotropy_validation():
    with pytest.raises(ValueError):
        Anisotropy(beta_inf=1.0)
    with pytest.raises(ValueError):
        Anisotropy(r_beta=0.0)


def test_jampy_logistic_form_matches_the_spec_form():
    """The spec's beta(r) equals JamPy's logistic anisotropy with alpha = 2."""
    a = Anisotropy(0.1, 0.6, 7.0)
    r_a, b0, binf, alpha = a.jampy_logistic
    r = np.geomspace(0.1, 100, 20)
    jampy_beta = b0 + (binf - b0) / (1 + (r_a / r) ** alpha)
    assert np.allclose(a.beta(r), jampy_beta)


# -------------------------------------------------------------------- units --
def test_observed_units_convert_pm_to_mas_per_year():
    p = Plummer(3.0e6, 5.0)
    j = SphericalJeans(p, p)
    D = 5.43
    obs = j.dispersions_observed(np.array([100.0, 300.0]), D)
    kms = j.dispersions_kms(np.array([100.0, 300.0]) * D * 1e3 / 206264.806)
    assert np.allclose(obs["sigma_pmr"] * KMS_PER_MASYR_KPC * D, kms["pmr"])
    assert np.allclose(obs["sigma_los"], kms["los"])
    # 1 mas/yr at 1 kpc is 4.74 km/s
    assert KMS_PER_MASYR_KPC == pytest.approx(4.7405, rel=1e-4)


# ------------------------------------------------------------------ hygiene --
def test_tracer_edge_is_truncated_not_extrapolated():
    """Gaussian tracers underflow long before r_max; the table must stop first."""
    mge = MGE([1e6], [10.0])
    j = SphericalJeans(mge, mge, r_max=1e4)
    assert 50.0 < j.r_cut < 200.0
    inside = np.geomspace(0.1, 0.9 * j.r_cut, 40)
    d = j.dispersions_kms(inside)
    assert np.all(np.isfinite(d["los"])) and np.all(d["los"] < 100.0)
    with pytest.raises(ValueError, match="tracer's edge"):
        j.dispersions_kms(np.array([1.5 * j.r_cut]))


def test_a_model_costs_milliseconds():
    import time

    p = Plummer(3.0e6, 5.0)
    start = time.perf_counter()
    for i in range(10):
        SphericalJeans(p, p, Anisotropy(0.0, 0.3 + 1e-3 * i, 8.0)).dispersions_kms(np.geomspace(0.05, 60, 70))
    assert (time.perf_counter() - start) / 10 < 0.05


# ---------------------------------------------------------- external checks --
jampy = pytest.importorskip("jampy", reason="JamPy not installed")


def test_matches_jampy_spherical_with_black_hole_and_varying_beta():
    from jampy.sph.jam_sph_proj import jam_sph_proj

    D_kpc = 5.43
    ppa = D_kpc * 1e3 / 206264.806
    frac, sig, Mtot, mbh = np.array([0.3, 0.5, 0.2]), np.array([2.0, 6.0, 15.0]), 3.0e6, 1.0e4
    mge = MGE((frac * Mtot).tolist(), sig.tolist(), name="stars")
    model = CompositeMassModel([mge, PointMass(mbh)])
    ani = Anisotropy(0.05, 0.5, 8.0)
    ours = SphericalJeans(model, mge, ani, r_max=2e3).dispersions_kms(np.geomspace(2, 800, 15) * ppa)
    surf = frac * Mtot / (2 * np.pi * sig**2)
    for tensor in ("los", "pmr", "pmt"):
        with contextlib.redirect_stdout(io.StringIO()):
            out = jam_sph_proj(surf, sig / ppa, surf, sig / ppa, mbh, D_kpc / 1e3, np.geomspace(2, 800, 15),
                               beta=[ani.r_beta / ppa, ani.beta_0, ani.beta_inf, 2.0], logistic=True,
                               tensor=tensor, epsrel=1e-5, plot=False, quiet=True)
        assert np.allclose(ours[tensor], np.asarray(out.model), rtol=3e-3), tensor


agama = pytest.importorskip("agama", reason="AGAMA not installed")


def _agama_moments(pot, dens, beta0, r):
    agama.setUnits(mass=1, length=1e-3, velocity=1)
    df = agama.DistributionFunction(type="QuasiSpherical", potential=pot, density=dens, beta0=beta0)
    rho, v2 = agama.GalaxyModel(pot, df).moments(np.column_stack([r, 0 * r, 0 * r]), dens=True, vel=False, vel2=True)
    return rho, np.sqrt(v2[:, 0]), np.sqrt(v2[:, 1])


def test_matches_agama_distribution_function_isotropic_and_tangential():
    agama.setUnits(mass=1, length=1e-3, velocity=1)
    pot = agama.Potential(type="Plummer", mass=3.0e6, scaleRadius=6.0)
    dens = agama.Density(type="Plummer", mass=3.0e6, scaleRadius=6.0)
    stars = Plummer(3.0e6, 6.0)
    r = np.geomspace(0.3, 100.0, 8)
    for beta0 in (0.0, -0.5):
        _, sr, st = _agama_moments(pot, dens, beta0, r)
        j = SphericalJeans(stars, stars, Anisotropy(beta0, beta0, 10.0))
        assert np.allclose(j.sigma_r(r), sr, rtol=1e-3), beta0
        assert np.allclose(j.sigma_t(r), st, rtol=1e-3), beta0


def test_an_evans_no_positive_df_for_radial_anisotropy_in_a_core():
    """A constant beta > 0 needs a tracer cusp gamma >= 2 beta; a Plummer core has
    gamma = 0. AGAMA's realised density then departs from the input at the centre,
    which is the sign that the requested model has no positive distribution
    function -- a constraint the K1 prior must respect (beta_0 <= 0 for a cored MGE)."""
    agama.setUnits(mass=1, length=1e-3, velocity=1)
    pot = agama.Potential(type="Plummer", mass=3.0e6, scaleRadius=6.0)
    dens = agama.Density(type="Plummer", mass=3.0e6, scaleRadius=6.0)
    stars = Plummer(3.0e6, 6.0)
    rho, _, _ = _agama_moments(pot, dens, 0.3, np.array([0.3]))
    assert rho[0] / stars.density(0.3)[0] > 1.5
