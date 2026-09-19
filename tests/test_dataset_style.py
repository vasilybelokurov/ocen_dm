"""One symbol and one colour per dataset; the component is carried by the fill only."""

from __future__ import annotations

import pytest

from ocen_dm.plotting.style import (DATASET_KEY_MAP, DATASET_STYLE, dataset_label,
                                    dataset_style)


def test_components_share_symbol_and_colour():
    for key in DATASET_STYLE:
        base = dataset_style(key, "combined")
        for comp in ("radial", "tangential"):
            st = dataset_style(key, comp)
            assert st["marker"] == base["marker"], f"{key}/{comp} changed symbol"
            assert st["color"] == base["color"], f"{key}/{comp} changed colour"
            assert st["fillstyle"] != base["fillstyle"], f"{key}/{comp} is not distinguishable"
    assert dataset_style("hst", "radial")["fillstyle"] == "left"
    assert dataset_style("hst", "tangential")["fillstyle"] == "right"


def test_radial_and_tangential_datasets_map_to_one_instrument():
    r, t = DATASET_KEY_MAP["hst_pm_radial"], DATASET_KEY_MAP["hst_pm_tangential"]
    assert r[0] == t[0] == "hst"
    assert (r[1], t[1]) == ("radial", "tangential")
    assert dataset_style("hst_pm_radial")["marker"] == dataset_style("hst_pm_tangential")["marker"]


def test_every_symbol_and_colour_is_unique_across_datasets():
    pairs = [(v["color"], v["marker"]) for v in DATASET_STYLE.values()]
    assert len(set(pairs)) == len(pairs), "two datasets would be indistinguishable"
    assert len({v["marker"] for v in DATASET_STYLE.values()}) == len(DATASET_STYLE)


def test_unfitted_data_keeps_its_identity():
    a = dataset_style("pristine", fitted=True)
    b = dataset_style("pristine", fitted=False)
    assert a["marker"] == b["marker"] and a["color"] == b["color"]
    assert b["alpha"] < a["alpha"]
    assert "not fitted" in dataset_label("pristine", fitted=False)


def test_every_likelihood_key_has_a_style():
    from ocen_dm.cli import DEFAULT_DATASETS
    from ocen_dm.kinematics.likelihood import DATASETS
    for key in list(DATASETS) + DEFAULT_DATASETS.split(","):
        assert key in DATASET_KEY_MAP, f"{key} has no agreed symbol"
        dataset_style(key)
