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


def test_truncation_formula_against_monte_carlo():
    """N1 of the Section-11 programme: sample p(v) ~ v^2 (-E)^(gamma-3/2) in the Kepler well by
    numerical inverse-CDF (substitution v = v_esc sin(theta) removes the endpoint singularity),
    truncate at Phi(r_J) and compare the retained fraction with 1 - I_x(gamma-1/2, 3/2).
    Uses no beta-function identity, so it is an independent check of the closed form."""
    rng = np.random.default_rng(7)
    th = np.linspace(0, np.pi / 2, 4001)
    for g in (1.0, 1.5):
        w = np.sin(th) ** 2 * np.cos(th) ** (2 * g - 2)
        c = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(th))]); c /= c[-1]
        for x in (0.1, 0.3, 0.6):                      # r in units of r_J; G = M = r_J = 1
            v = np.sqrt(2.0 / x) * np.sin(np.interp(rng.random(40000), c, th))
            f = np.mean(0.5 * v**2 - 1.0 / x < -1.0)
            e = np.sqrt(f * (1 - f) / v.size)
            assert abs(f - kepler_truncation_factor(x, g)) < 4 * e
