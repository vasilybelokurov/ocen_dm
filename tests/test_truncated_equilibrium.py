import numpy as np
from ocen_dm.tails.truncated_equilibrium import kepler_truncation_factor, shock_heating


def test_truncation_factor_limits_and_closed_form():
    assert kepler_truncation_factor(0.0, 1.0) == 1.0 and kepler_truncation_factor(1.0, 1.0) < 1e-12
    x = np.linspace(0, 1, 50)
    assert np.all(np.diff(kepler_truncation_factor(x, 1.0)) <= 1e-12)
    assert np.allclose(kepler_truncation_factor(x, 1.5), (1 - x) ** 1.5, atol=1e-10)
    assert np.all(kepler_truncation_factor(x[1:-1], 1.5) > kepler_truncation_factor(x[1:-1], 1.0))   # steeper cusps keep more


def test_shock_heating_scalings():
    s1 = shock_heating(20.0, 1.36e10, 1570.0, 387.0); s2 = shock_heating(40.0, 1.36e10, 1570.0, 387.0)
    assert s2["dv"] / s1["dv"] == 2.0 and s2["omega_tau"] < s1["omega_tau"] and s2["dE_E"] > s1["dE_E"]
    assert 0 < s1["A"] < 1
