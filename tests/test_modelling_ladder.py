"""The implemented ladder and legacy scaled diagnostic build with the intended priors.

The case numbered 3 below exercises the optional scale diagnostic; it is not an
adopted rung for evidence comparison."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics.fit import (BETA0_MAX_CORED_TRACER_WITH_BH, CONSTANT_BETA_RANGE,
                                    DarkMatterModel, FitProblem, NoDarkMatterModel)
from ocen_dm.kinematics.likelihood import KinematicData
from ocen_dm.paths import processed_dir

_HAS = (processed_dir() / "kinematics" / "ocen_pm_dispersion_edr3_ours.ecsv").exists()

RUNGS = {
    0: dict(constant_beta=True, fixed={"beta_0": 0.0}, instruments=()),
    1: dict(constant_beta=True, instruments=()),
    2: dict(instruments=(), anisotropy_profile="turnover"),
    3: dict(),
}
EXPECTED_PARAMS = {0: (5, 7), 1: (6, 8), 2: (10, 12), 3: (11, 13)}


@pytest.mark.parametrize("rung", sorted(RUNGS))
def test_parameter_counts_match_the_plan(rung):
    k1 = NoDarkMatterModel(**RUNGS[rung]); k2 = DarkMatterModel(gamma=0.0, **RUNGS[rung])
    assert (len(k1.names), len(k2.names)) == EXPECTED_PARAMS[rung]
    # K1 is nested in K2: same parameters plus the two halo ones
    assert set(k2.names) - set(k1.names) == {"M_dm_100", "r_s"}


def test_an_evans_bound_applies_only_where_there_is_a_distinct_centre():
    """Pin the adopted diagnostic prior policy, not an exemption from the physical DF bound."""
    varying = [p for p in NoDarkMatterModel(instruments=()).parameters if p.name == "beta_0"][0]
    assert varying.prior.hi == BETA0_MAX_CORED_TRACER_WITH_BH == -0.5
    constant = [p for p in NoDarkMatterModel(constant_beta=True, instruments=()).parameters
                if p.name == "beta_0"][0]
    assert (constant.prior.lo, constant.prior.hi) == CONSTANT_BETA_RANGE == (-1.0, 0.5)
    # the range the data imply (beta ~ -0.4 to +0.3) must be inside the constant-beta prior
    assert constant.prior.lo < -0.4 and constant.prior.hi > 0.3
    assert "beta_0" not in NoDarkMatterModel(**RUNGS[0]).names, "rung 0 has no anisotropy parameter"


@pytest.mark.skipif(not _HAS, reason="products not built")
@pytest.mark.parametrize("rung", sorted(RUNGS))
def test_every_rung_evaluates_a_finite_likelihood(rung):
    from ocen_dm.cli import DEFAULT_DATASETS
    data = KinematicData.load(DEFAULT_DATASETS.split(","))
    for fam in (NoDarkMatterModel(**RUNGS[rung]), DarkMatterModel(gamma=0.0, **RUNGS[rung])):
        prob = FitProblem(fam, data)
        for frac in (0.25, 0.5, 0.75):
            ln = prob.loglike_unit(np.full(len(fam.names), frac))
            assert np.isfinite(ln) and ln > -1e50


def test_cli_exposes_every_rung():
    from ocen_dm.cli import build_parser
    p = build_parser()
    a = p.parse_args(["fit", "--isotropic", "--no-scales"])
    assert a.isotropic and a.no_scales and not a.constant_beta
    b = p.parse_args(["fit", "--constant-beta", "--no-scales"])
    assert b.constant_beta and not b.isotropic
