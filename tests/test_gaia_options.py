"""The two modelling choices are explicit options, not averaged over."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.paths import processed_dir

_HAS = (processed_dir() / "kinematics" / "ocen_pm_dispersion_edr3_ours.ecsv").exists()
pytestmark = pytest.mark.skipif(not _HAS, reason="product not built")


def test_default_is_raw_errors_and_the_published_rotation():
    import inspect
    from ocen_dm.kinematics.likelihood import _edr3_ours_component
    sig = inspect.signature(_edr3_ours_component)
    assert sig.parameters["error_model"].default == "raw"
    assert sig.parameters["rotation"].default == "published"


def test_the_product_carries_both_error_models():
    from ocen_dm.kinematics.outer_gaia import load_edr3_profile
    t = load_edr3_profile()
    for c in ("sigma_pmr", "sigma_pmt"):
        assert c in t.colnames and c + "_eta" in t.colnames
    # the inflated fit is systematically the lower of the two: it removes more noise
    assert np.mean(np.asarray(t["sigma_pmr_eta"]) <= np.asarray(t["sigma_pmr"])) > 0.7
    # and the recorded systematic is no longer folded into the quoted error
    assert np.allclose(t["sigma_pm_err"], t["sigma_stat"])


def test_the_options_change_the_data_they_are_supposed_to():
    from ocen_dm.kinematics.likelihood import load_profile
    base = load_profile("gaia_edr3_ours_tangential")
    eta = load_profile("gaia_edr3_ours_tangential", error_model="eta")
    ours = load_profile("gaia_edr3_ours_tangential", rotation="ours")
    assert not np.allclose(base.value, eta.value), "the error model must move the values"
    assert np.allclose(base.value, ours.value), "the rotation must not move the values"
    assert not np.allclose(base.streaming2, ours.streaming2), "it must move the streaming"
    # the radial component carries no published rotation, but does carry ours
    r_pub = load_profile("gaia_edr3_ours_radial")
    r_our = load_profile("gaia_edr3_ours_radial", rotation="ours")
    assert np.allclose(r_pub.streaming2, 0.0)
    assert np.all(r_our.streaming2 >= 0) and np.any(r_our.streaming2 > 0)


def test_gaia_dr2_is_out_of_the_default_but_still_loadable():
    from ocen_dm.cli import DEFAULT_DATASETS
    from ocen_dm.kinematics.likelihood import KinematicData, load_profile
    assert "gaia_dr2_pm" not in DEFAULT_DATASETS
    assert KinematicData.load(DEFAULT_DATASETS.split(",")).n_points == 89
    assert len(load_profile("gaia_dr2_pm").r) == 9, "still available by key"


def test_the_cli_exposes_both_switches():
    from ocen_dm.cli import build_parser
    a = build_parser().parse_args(["fit", "--gaia-errors", "eta", "--gaia-rotation", "ours"])
    assert (a.gaia_errors, a.gaia_rotation) == ("eta", "ours")
    d = build_parser().parse_args(["fit"])
    assert (d.gaia_errors, d.gaia_rotation) == ("raw", "published")
