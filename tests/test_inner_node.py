"""Fixed inner beta node for the spline anisotropy (2026-11-19, docs/PLAN_2026-11-19_limits.md step 1)."""
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import integrate

from ocen_dm.kinematics.anisotropy import SPLINE_NODES_ARCSEC, SplineAnisotropy
from ocen_dm.kinematics.fit import NoDarkMatterModel
from ocen_dm.kinematics.run_io import family_config, family_from_config

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT/"results/fits/recover_real_K1_spline_nw"


def test_inner_node_log_g_exact():
    a = SplineAnisotropy((0.05, 0.26, 2.6, 6.4), (-0.5, 0.0, 0.3, -0.4))
    r = np.geomspace(0.01, 30, 40)
    q = [2*integrate.quad(lambda x: float(a.beta(np.exp(x))), np.log(0.05), np.log(ri), limit=200)[0] for ri in r]
    assert np.allclose(a.log_integrating_factor(r), q, atol=1e-10)


def test_inner_node_priors_and_beta():
    fam = NoDarkMatterModel(anisotropy_profile="spline", inner_node_arcsec=2.0, instruments=())
    pri = {p.name: p.prior for p in fam.parameters}
    assert (pri["b10"].lo, pri["b10"].hi) == (-1.0, 0.0) and pri["b100"].hi == 0.95
    legacy = NoDarkMatterModel(anisotropy_profile="spline", instruments=(), mge_fit=fam.mge_fit)
    assert {p.name: p.prior for p in legacy.parameters}["b10"].hi == -0.5
    th = dict(M_star=2e6, M_rem=5e5, a_rem=4., M_bh=3e4, distance=5.43,
              **{f"b{int(r)}": 0.0 for r in SPLINE_NODES_ARCSEC})
    model, D, _ = fam.build(th)
    pc = D*1e3/206264.806
    assert model.anisotropy.beta(2*pc) == pytest.approx(-0.5) and model.anisotropy.beta(0.5*pc) == pytest.approx(-0.5)
    assert model.anisotropy.beta(10*pc) == pytest.approx(0.0)
    with pytest.raises(ValueError):
        NoDarkMatterModel(anisotropy_profile="turnover", inner_node_arcsec=2.0, instruments=())


def test_inner_node_config_round_trip():
    fam = NoDarkMatterModel(anisotropy_profile="spline", inner_node_arcsec=2.0, instruments=())
    cfg = family_config(fam)
    assert cfg["inner_node_arcsec"] == 2.0
    assert family_from_config(cfg).inner_node_arcsec == 2.0
    assert "inner_node_arcsec" not in family_config(NoDarkMatterModel(anisotropy_profile="spline", instruments=(),
                                                                      mge_fit=fam.mge_fit))


@pytest.mark.skipif(not RUN.exists(), reason="production spline run not present")
def test_legacy_spline_likelihood_unchanged():
    from ocen_dm.kinematics.report import load_run_metadata, problem_for
    P = problem_for(load_run_metadata(RUN))
    assert getattr(P.family, "inner_node_arcsec", None) is None
    s = json.loads((RUN/"summary.json").read_text())
    assert float(np.sum(P.loglike_vector(np.load(RUN/"ml_x.npy")))) == pytest.approx(s["lnL_max"], abs=1e-6)
