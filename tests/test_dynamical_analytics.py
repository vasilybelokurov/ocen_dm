"""Independent analytic limits for the pilot comparison quadrature."""
import importlib.util
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

spec = importlib.util.spec_from_file_location(
    "dynamical_analytics", Path(__file__).resolve().parents[1]/"bin/compare_dynamical_analytics.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_constant_df_recovers_exact_kepler_energy_cut():
    # f(E)=constant gives rho proportional to psi**(3/2), hence a Kepler
    # r**(-3/2) cusp with exact retained fraction (1-r/r_t)**(3/2).
    psi = np.geomspace(.2, 100, 21)
    rho = comparison.density_integral(psi, 0., np.ones_like)
    assert_allclose(rho, 4*np.pi/3*(2*psi)**1.5, rtol=1e-12)
    cut = comparison.density_integral(psi, 1., np.ones_like)
    assert_allclose(cut/rho, np.maximum(1-1/psi, 0)**1.5, atol=1e-13)


def test_linear_df_recovers_nontrivial_kepler_formula():
    psi = np.geomspace(1.01, 100., 20)
    retained = comparison.density_integral(psi, 1., lambda e: e)/comparison.density_integral(psi, 0., lambda e: e)
    expected = comparison.kepler_truncation_factor(1/psi, gamma=2.5)
    assert_allclose(retained, expected, rtol=1e-12)


def test_shell_averaging_uses_mass_weights_and_matches_exact_integrals():
    radii = np.array([3., 20., 70.])
    # rho=r^-1 and an artificial retention 1-r: ratio involves <r> weighted
    # by r^2 rho dr; this catches pointwise or volume-only shell averaging.
    lo, hi = radii/1000*np.exp(-.1), radii/1000*np.exp(.1)
    expected = 1-(hi**3-lo**3)/3/((hi**2-lo**2)/2)
    retained = comparison.shell_retention(radii, lambda r: 1/r, lambda r: 1/r-1)
    assert_allclose(retained, expected, rtol=1e-12)
