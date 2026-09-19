"""Our rebuilt Gaia EDR3 profile: independent of the error model, and only where data exist."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = (processed_dir() / "kinematics" / "ocen_pm_dispersion_edr3_ours.ecsv").exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="product not built")


def test_starts_where_gaia_becomes_usable():
    from ocen_dm.kinematics.outer_gaia import R_MIN_ARCSEC, load_edr3_profile
    t = load_edr3_profile()
    assert t["r_lower"].min() >= R_MIN_ARCSEC
    assert R_MIN_ARCSEC >= 460.0, "inside this the quality flag leaves too few stars"
    assert len(t) >= 6 and np.all(t["n_stars"] > 500)


def test_error_model_barely_moves_it():
    """The point of the err < 0.4 sigma cut: raw and inflated errors must agree."""
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    t = load_edr3_profile()
    spread = np.abs(np.asarray(t["sigma_raw"]) / np.asarray(t["sigma_eta"]) - 1.0)
    assert np.all(spread < 0.02), "an error-model-independent measurement moves by < 2 per cent"
    # and that spread is carried, not discarded
    assert np.allclose(t["sigma_sys"], 0.5 * np.abs(t["sigma_raw"] - t["sigma_eta"]))
    assert np.all(t["sigma_pm_err"] >= t["sigma_stat"])


def test_profile_declines_and_carries_rotation():
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    t = load_edr3_profile()
    assert np.all(np.diff(np.asarray(t["sigma_pm"])) < 0), "dispersion must fall outwards"
    # streaming is the mean-square of the fitted means, to remove from the model second moment
    assert np.allclose(t["streaming2"], 0.5 * (np.asarray(t["mean_pmr"]) ** 2
                                               + np.asarray(t["mean_pmt"]) ** 2))
    assert np.all(t["streaming2"] >= 0)


def test_registered_and_excludes_the_published_spline():
    from ocen_dm.kinematics.likelihood import KinematicData, load_profile
    p = load_profile("gaia_edr3_ours")
    assert p.instrument == "GaiaEDR3" and p.kind == "pmc"
    assert p.r.min() > 460 and p.streaming2 is not None
    with pytest.raises(ValueError):
        KinematicData.load(["gaia_edr3_ours", "gaia_edr3_pm"])


def test_is_the_default_dataset():
    from ocen_dm.cli import DEFAULT_DATASETS
    assert "gaia_edr3_ours" in DEFAULT_DATASETS
    assert "gaia_edr3_pm" not in DEFAULT_DATASETS.replace("gaia_edr3_ours", "")


def test_plot_renders(tmp_path):
    from ocen_dm.plotting.constraints import plot_pm_datasets
    out = plot_pm_datasets(tmp_path / "pm.png")
    assert out.exists() and out.stat().st_size > 50_000
