import numpy as np
import pytest
tt = pytest.importorskip("ocen_dm.tails.tidal_tracks")


def test_track_and_inversion():
    for x in (0.9, 0.3, 0.05, 1e-3):
        assert abs(tt.x_from_T(tt.track_T(x)) - x) / x < 1e-6
    assert abs(tt.track_v(1.0) - 1.0) < 1e-12 and abs(tt.track_M(1.0) - 1.0) < 1e-12
    x = np.logspace(-4, 0, 50); assert np.all(np.diff(tt.track_T(x)) > 0) and np.all(np.diff(tt.track_M(x)) > 0)


def test_nfw_numbers():
    h = tt.NFW(1e10, 5.5, z=2.0)
    assert 20 < h.r200 < 23 and 8.0 < h.rmx0 < 9.0 and 45 < h.Vmx0 < 51
    assert abs(h.mass(h.r200) / 1e10 - 1) < 1e-6
    # EN21 asymptotic criterion: dense host at pericentre leaves a tiny self-bound remnant
    r, m = tt.asymptotic_density_radius(h, 7e8)
    assert 0.002 < r < 0.01 and m < 1e4


def test_walk_monotonic_and_delayed_by_eccentricity():
    h = tt.NFW(1e10, 5.5)
    n = 30; tp = np.linspace(9, 0, n); rp = np.full(n, 10.0); Tp = np.full(n, 0.3); To = np.full(n, 0.3)
    w_circ = tt.walk(h, tp, rp, rp, Tp, To)
    w_ecc = tt.walk(h, tp, rp, rp * 10, Tp, To)
    assert np.all(np.diff(w_circ["mfrac"]) <= 1e-12)
    assert w_ecc["mfrac"][-1] > w_circ["mfrac"][-1]          # eccentric orbits are delayed (eq. 4)
