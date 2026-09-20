"""The split-normal likelihood must be a continuous, normalised density."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from ocen_dm.kinematics.likelihood import BinnedProfile, ProfileLikelihood


def _profile(lo, hi):
    return BinnedProfile(name="t", kind="pmc", r=np.array([1.0]), r_lower=None, r_upper=None,
                         value=np.array([1.0]), err_lo=np.array([lo]), err_hi=np.array([hi]),
                         instrument="X")


@pytest.mark.parametrize("lo,hi", [(0.5, 2.0), (2.0, 0.5), (1.0, 1.0), (0.1, 0.3)])
def test_continuous_across_the_datum(lo, hi):
    """It jumped by err_hi/err_lo until 2026-09-20 -- a factor of 4 for (0.5, 2.0)."""
    p = _profile(lo, hi)
    f = ProfileLikelihood._split_normal_lnlike
    below = f(np.array([1.0 - 1e-9]), p)[0]
    above = f(np.array([1.0 + 1e-9]), p)[0]
    assert abs(below - above) < 1e-6


@pytest.mark.parametrize("lo,hi", [(0.5, 2.0), (2.0, 0.5), (1.0, 1.0)])
def test_integrates_to_one(lo, hi):
    p = _profile(lo, hi)
    dens = lambda x: float(np.exp(ProfileLikelihood._split_normal_lnlike(np.array([x]), p)[0]))
    total, _ = quad(dens, 1.0 - 40 * lo, 1.0 + 40 * hi, limit=400)
    assert abs(total - 1.0) < 1e-4


def test_reduces_to_a_normal_when_the_errors_are_equal():
    p = _profile(0.7, 0.7)
    x = np.array([0.3, 1.0, 2.2])
    ref = -0.5 * ((x - 1.0) / 0.7) ** 2 - np.log(0.7) - 0.5 * np.log(2 * np.pi)
    assert np.allclose(ProfileLikelihood._split_normal_lnlike(x, p), ref, atol=1e-12)


def test_the_agama_prior_stays_inside_what_the_backend_accepts():
    from ocen_dm.kinematics.fit import AGAMA_BETA0_MIN, NoDarkMatterModel
    fam = NoDarkMatterModel(backend="agama")
    b0 = [q for q in fam.parameters if q.name == "beta_0"][0]
    assert b0.prior.lo >= AGAMA_BETA0_MIN, "AGAMA rejects beta_0 below -0.5"
    plain = NoDarkMatterModel()
    assert [q for q in plain.parameters if q.name == "beta_0"][0].prior.lo == -1.0


def test_run_records_every_switch_that_changes_the_inputs():
    import inspect
    from ocen_dm.kinematics.fit import run_nested
    assert "dataset_options" in inspect.signature(run_nested).parameters
    from ocen_dm.kinematics.report import comparison_table
    src = inspect.getsource(comparison_table)
    assert "dataset_options" in src, "the fingerprint must include the switches"
    assert '"seed"' in src or "'seed'" in src, "and the mock seed"
