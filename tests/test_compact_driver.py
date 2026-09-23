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
