"""Massive-remnant migration test: the vectorised leapfrog in the time-dependent barred potential."""
import numpy as np
import pytest

bm = pytest.importorskip("ocen_dm.tails.bar_migration")
ft = pytest.importorskip("ocen_dm.tails.friction_test")
pytestmark = pytest.mark.skipif(not bm.POT_DIR.exists(), reason="oCen_bar clone not found")


def test_leapfrog_without_friction_matches_agama_orbit():
    ref = bm.back_integrate(24.0, n_samples=10, n_times=161)
    res = ft.integrate_with_friction(24.0, 0.0, n_samples=10, dt_myr=0.25)
    # the same samples 8 Gyr back: energies agree to the level set by chaotic divergence near the resonance
    assert np.median(np.abs(res["E"][0] - ref.E[0]) / np.abs(ref.E[0])) < 0.05
    assert abs(res["frac"] - ref.frac_inside) <= 0.2


def test_friction_changes_the_endpoint_monotonically_with_mass():
    e0 = [np.median(ft.integrate_with_friction(24.0, m, n_samples=10)["E"][0]) for m in (0.0, 1e7, 1e8)]
    assert e0[0] < e0[1] < e0[2]          # backward friction pumps energy; more mass, more pumping
