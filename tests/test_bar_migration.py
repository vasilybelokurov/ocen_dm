"""Class-3 bar-migration set-up (Dillamore+2026 via the oCen_bar pipeline)."""
import numpy as np
import pytest

bm = pytest.importorskip("ocen_dm.tails.bar_migration")
pytestmark = pytest.mark.skipif(not bm.POT_DIR.exists(), reason="oCen_bar clone with Hunter24 potentials not found")


def test_bar_history_matches_the_paper():
    h = bm.BarHistory(24.0)
    assert abs(h.omega[-1] - 24.0) < 1e-6                      # ends at the requested present-day speed
    assert abs(h.omega_1 - 45.0) < 1.0                          # paper: initial ~45 km/s/kpc for eta=0.003
    assert np.all(np.diff(h.omega) <= 1e-9)                     # monotonic deceleration
    assert h.frac[0] == 0 and h.frac[-1] == 1 and abs(bm.bar_fraction(0.5) - 0.5) < 1e-12
    assert abs(h.phi[-1]) < 1e-12                               # bar angle zero today
    eta_late = -np.gradient(h.omega, h.t)[-100:] / h.omega[-100:] ** 2
    assert np.allclose(eta_late, 0.003, atol=1e-4)


def test_present_day_state_is_retrograde_and_outside_gse():
    pa = bm.axisymmetric_potential()
    E, Lz = bm.energy_lz(bm.present_day_samples(1), pa)
    assert Lz[0] < -450 and -1.5e5 < E[0] < -1.4e5
    assert not bm.inside_gse(E, Lz, bm.gse_contours(pa))[0]
    s = bm.present_day_samples(200, seed=1)
    assert s.shape == (200, 6) and np.std(np.linalg.norm(s[:, :3], axis=1)) < 0.1


def test_gse_contours_nested():
    cont = bm.gse_contours(bm.axisymmetric_potential())
    assert len(cont) == 5 and cont[0].shape[0] > cont[1].shape[0]
    inner_centre = cont[3].mean(0)
    assert bm.inside_gse([inner_centre[1]], [inner_centre[0]], cont, level=0)[0]


def test_short_backward_integration_runs_and_migrates():
    run = bm.back_integrate(24.0, n_samples=20, n_times=41)
    assert run.traj.shape == (41, 20, 6) and abs(run.t[0]) < 1e-9 and abs(run.t[-1] - 8.0) < 1e-9
    # today (t = tf) retrograde; at t = 0 the samples have moved to higher E and less retrograde Lz
    assert np.median(run.Lz[-1]) < -450 and np.median(run.Lz[0]) > np.median(run.Lz[-1]) + 300
    assert np.median(run.E[0]) > np.median(run.E[-1])
    assert 0.0 <= run.frac_inside <= 1.0
