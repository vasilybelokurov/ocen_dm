"""Check fixed-density refits against analytic masses and the fitted likelihood."""
import importlib.util
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest

from ocen_dm.kinematics.fit import DarkMatterModel, FitProblem
from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData
from ocen_dm.light_model import MGEFit
from ocen_dm.mass_models import TruncatedGNFW


spec = importlib.util.spec_from_file_location(
    "density_profile", Path(__file__).resolve().parents[1]/"bin/profile_dm_density_limit.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def test_fixed_density_mass_matches_analytic_nfw():
    for rs in [3., 20., 1000.]:
        rho = 2.
        rho_s = rho*(20/rs)*(1+20/rs)**2
        expected = 4*np.pi*rho_s*rs**3*(np.log1p(100/rs)-(100/rs)/(1+100/rs))
        assert_allclose(analysis.mass100_for_density(rho, rs, 1, np.inf),
                        expected, rtol=1e-11)


@pytest.mark.parametrize("gamma", [0., 1.])
def test_fixed_density_mass_matches_production_halo_across_prior(gamma):
    for rs in np.geomspace(3, 1000, 9):
        unit = TruncatedGNFW(rho_s=1., r_s=rs, gamma=gamma, r_t=1000.)
        m100 = analysis.mass100_for_density(5., rs, gamma, 1000.)
        actual_rho = m100/unit.enclosed_mass(100.)*unit.density(20.)
        assert_allclose(actual_rho, 5., rtol=1e-7)


def test_parameter_elimination_and_objective_match_full_likelihood():
    mge = MGEFit(np.array([30.,300.,1500.]), np.array([.1,.8,.1]), 0., 0., 3)
    family = DarkMatterModel(gamma=0., mge_fit=mge, instruments=(),
                             anisotropy_profile="turnover")
    data = KinematicData(tuple(
        BinnedProfile(name=kind, kind=kind, r=np.array([30.,100.,500.]),
                      r_lower=None, r_upper=None, value=np.array([.3,.8,.2]),
                      err_lo=np.full(3,.02), err_hi=np.full(3,.05), instrument="mock")
        for kind in ("pmr", "pmt")))
    problem = FitProblem(family, data)
    profile = analysis.ProfileProblem(problem, 1e6)
    x = family.transform(np.full(len(family.names), .5))
    encoded = profile.encode(x, 2.)
    xx = profile.decode(encoded, 2.)
    assert_allclose(xx[profile.keep], x[profile.keep])
    assert_allclose(family.halo(family.to_dict(xx)).density(20.), 2., rtol=1e-7)
    residual = profile.residual(encoded, 2.)
    chi = sum(v[0] for v in problem.chi2(xx).values())
    assert_allclose(residual@residual, chi, rtol=1e-12)  # distance is at its mean
    x2 = xx.copy()
    x2[family.names.index("distance")] += .03
    e2 = profile.encode(x2, 2.)
    r2 = profile.residual(e2, 2.)
    expected_change = -2*(problem.loglike_vector(x2)-problem.loglike_vector(xx)) + (.03/.05)**2
    assert_allclose(r2@r2-residual@residual, expected_change, rtol=1e-10)
    # A density requiring a halo mass beyond the archived support is rejected.
    with pytest.raises(ValueError, match="DM-mass support"):
        profile.residual(encoded, 1e6)
