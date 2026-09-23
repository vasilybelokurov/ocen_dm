"""Driver helpers for compact-DF batches (bin/run_compact_df_recovery.py)."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

_PATH = Path(__file__).resolve().parents[1]/"bin"/"run_compact_df_recovery.py"
_spec = importlib.util.spec_from_file_location("run_compact_df_recovery", _PATH)
driver = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(driver)

COORDS = [dict(path="M_star", lower=1e6, upper=6e6, log=True),
          dict(path="stellar.alpha", lower=.6, upper=3.)]


def test_override_bounds_changes_only_named_paths_and_rejects_unknown():
    out = driver.override_bounds(COORDS, {"stellar.alpha": (.3, 4.)})
    assert out[1]["lower"] == .3 and out[1]["upper"] == 4.
    assert out[0] == COORDS[0] and COORDS[1]["lower"] == .6  # input not mutated
    with pytest.raises(ValueError):
        driver.override_bounds(COORDS, {"stellar.J0": (1., 2.)})


def test_latin_starts_are_reproducible_inside_central_box_and_stratified():
    a, b = driver.latin_starts(COORDS, 5), driver.latin_starts(COORDS, 5)
    assert a == b and len(a) == 5
    m = np.array([s[0] for s in a])
    u = (np.log(m)-np.log(1e6))/(np.log(6e6)-np.log(1e6))
    assert np.all((u > .2) & (u < .8))
    # one start per stratum of the central box
    assert sorted(np.floor((u-.2)/.12).astype(int).tolist()) == [0, 1, 2, 3, 4]
    assert driver.latin_starts(COORDS, 0) == []


def _profile(name, err):
    from ocen_dm.kinematics.likelihood import BinnedProfile
    return BinnedProfile(name, "pmr", np.array([10., 20.]), None, None, np.array([.5, .4]),
                         np.full(2, err), np.full(2, err), "x")


def test_observed_variant_subsets_and_floors_only_gaia():
    from ocen_dm.kinematics.likelihood import KinematicData
    data = KinematicData((_profile("hst_pm_radial_ours", .01), _profile("gaia_edr3_ours_radial", .01)))
    v = driver.observed_variant(data, gaia_error_floor=.02)
    np.testing.assert_allclose(v.profiles[0].err_lo, .01)
    np.testing.assert_allclose(v.profiles[1].err_hi, np.hypot(.01, .02))
    only = driver.observed_variant(data, datasets=["gaia_edr3_ours_radial"])
    assert [p.name for p in only.profiles] == ["gaia_edr3_ours_radial"]
    with pytest.raises(ValueError):
        driver.observed_variant(data, datasets=["muse_los_dispersion"])
