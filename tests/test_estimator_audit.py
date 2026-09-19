"""The before/after audit must reproduce both defects exactly, and only inside its block."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics import outer_profile as op
from ocen_dm.kinematics.estimator_audit import with_legacy


def _sample(**kw):
    import sys
    sys.path.insert(0, "tests")
    from test_mixture_recovery import _fit, _sample as mk
    return mk(**kw), _fit


def test_legacy_flags_are_off_by_default_and_restored():
    assert op.LEGACY == set()
    with with_legacy("signs"):
        assert "signs" in op.LEGACY
    assert op.LEGACY == set(), "the context manager must not leak"
    with pytest.raises(ZeroDivisionError):
        with with_legacy("interval"):
            raise ZeroDivisionError
    assert op.LEGACY == set(), "must be restored even on an exception"


def test_legacy_signs_reproduce_the_biased_streaming():
    (samp, ad), fit = _sample(mt=-0.20)
    good = fit(samp, ad)
    with with_legacy("signs"):
        bad = fit(samp, ad)
    assert abs(good["mean_t"] + 0.20) < 0.015
    assert abs(bad["mean_t"] + 0.20) > 0.05, "the defect must actually reappear"
    assert abs(bad["mean_t"]) < abs(good["mean_t"]), "it biased the streaming low"


def test_legacy_interval_reproduces_the_floor():
    (samp, ad), fit = _sample(n=60000)
    good = fit(samp, ad)
    with with_legacy("interval"):
        bad = fit(samp, ad)
    assert abs(bad["sigma_r_err"] / bad["sigma_r"] - 0.02) < 0.004, "2 per cent per component"
    assert good["sigma_r_err"] < 0.4 * bad["sigma_r_err"]


def test_audit_plot_renders(tmp_path):
    from ocen_dm.plotting.constraints import plot_estimator_audit
    out = plot_estimator_audit(tmp_path / "audit.png")
    assert out.exists() and out.stat().st_size > 80_000
    assert op.LEGACY == set()
