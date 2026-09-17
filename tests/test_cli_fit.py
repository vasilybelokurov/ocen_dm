"""End-to-end smoke test of ``ocen fit`` on mock data with a deliberately tiny sampler budget."""

from __future__ import annotations

import json

import numpy as np
import pytest

from ocen_dm.cli import main
from ocen_dm.paths import processed_dir, results_dir

_HAS_DATA = (processed_dir() / "kinematics" / "omegacat_vi_pm_radial.ecsv").exists()


@pytest.mark.skipif(not _HAS_DATA, reason="processed kinematics not present")
def test_fit_cli_writes_all_products(tmp_path, monkeypatch):
    from ocen_dm.kinematics import NoDarkMatterModel
    fam = NoDarkMatterModel()
    x = fam.transform(np.full(len(fam.names), 0.5))
    np.save(tmp_path / "x.npy", x)
    monkeypatch.setattr("ocen_dm.paths.results_dir", lambda: tmp_path / "results")   # results -> tmp, data untouched
    rc = main(["fit", "--family", "K1", "--datasets", "hst_pm_radial,muse_los_dispersion",
               "--mock-from", str(tmp_path / "x.npy"), "--n-live", "40", "--max-ncalls", "1200",
               "--label", "smoke", "--seed", "3"])
    assert rc == 0
    out = tmp_path / "results" / "fits" / "smoke"
    for f in ("posterior.ecsv", "summary.json", "profiles.npz", "run.yaml", "ml_x.npy"):
        assert (out / f).exists(), f
    summary = json.loads((out / "summary.json").read_text())
    assert summary["family"] == "K1_noDM_composite" and set(summary["parameters"]) == set(fam.names)   # default tracer
    assert summary["n_points"] == 69                       # 40 + 29 bins
    prof = np.load(out / "profiles.npz")
    assert prof["M_total"].shape[1] == 60 and np.all(np.isfinite(prof["v_circ"]))


def test_fit_cli_rejects_unknown_family():
    with pytest.raises(SystemExit):
        main(["fit", "--family", "K9"])
