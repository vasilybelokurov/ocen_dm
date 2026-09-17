"""Every radius axis carries the other unit: arcsec plots get pc, pc plots get arcsec."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from ocen_dm.plotting.style import ARCSEC_PER_RAD, add_arcsec_axis, add_pc_axis

D = 5.43
PC_PER_ARCSEC = D * 1e3 / ARCSEC_PER_RAD          # 0.02633 pc per arcsec at 5.43 kpc


def _forward(sec):
    """The secondary axis's parent -> secondary conversion."""
    return sec._functions[0]


def test_pc_axis_converts_arcsec_and_arcmin():
    fig, ax = plt.subplots()
    sec = add_pc_axis(ax, D)
    assert _forward(sec)(1.0) == pytest.approx(PC_PER_ARCSEC, rel=1e-12)
    assert _forward(sec)(300.0) == pytest.approx(7.9, abs=0.05)        # 300" ~ 7.9 pc
    assert "pc" in sec.get_xlabel() and "5.43" in sec.get_xlabel()   # two decimals, not 5.43000001
    sec_min = add_pc_axis(ax, D, per_unit_arcsec=60.0)
    assert _forward(sec_min)(1.0) == pytest.approx(60 * PC_PER_ARCSEC, rel=1e-12)
    plt.close(fig)


def test_arcsec_axis_is_the_inverse_conversion():
    fig, ax = plt.subplots()
    sec = add_arcsec_axis(ax, D)
    assert _forward(sec)(PC_PER_ARCSEC) == pytest.approx(1.0, rel=1e-12)
    assert _forward(sec)(100.0) == pytest.approx(3798.0, abs=1.0)      # 100 pc ~ 63 arcmin
    sec_min = add_arcsec_axis(ax, D, unit="arcmin")
    assert _forward(sec_min)(100.0) == pytest.approx(63.3, abs=0.1)
    assert "arcmin" in sec_min.get_xlabel()
    plt.close(fig)


def test_round_trip_through_both_helpers():
    fig, ax = plt.subplots()
    f_pc = _forward(add_pc_axis(ax, D))
    f_as = _forward(add_arcsec_axis(ax, D))
    for a in (1.0, 37.0, 2400.0):
        assert f_as(f_pc(a)) == pytest.approx(a, rel=1e-12)
    plt.close(fig)


def test_distance_enters_the_conversion():
    fig, ax = plt.subplots()
    near, far = _forward(add_pc_axis(ax, 4.59)), _forward(add_pc_axis(ax, 5.43))
    assert far(100.0) / near(100.0) == pytest.approx(5.43 / 4.59, rel=1e-12)
    plt.close(fig)
