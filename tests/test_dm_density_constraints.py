"""Analytic checks on likelihood tempering and model-averaged sensitivity."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
from numpy.testing import assert_allclose
from numpy.polynomial.hermite import hermgauss

BIN = Path(__file__).resolve().parents[1]/"bin"
sys.path.insert(0,str(BIN))
spec = importlib.util.spec_from_file_location("dm_density_constraints",BIN/"analyse_dm_density_constraints.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def test_tempering_matches_analytic_gaussian_refit():
    # N(0,1) prior and y=2 with sigma=1 give N(1,1/2) posterior.
    nodes,base = hermgauss(100)
    x = 1+nodes
    base = base/np.sqrt(np.pi)
    ll = -.5*(x-2)**2-.5*np.log(2*np.pi)
    for fraction in (0.,.2,1.):
        w,_ = analysis.tempered_weights(base,ll,fraction)
        mean = 2*(1-fraction)/(2-fraction)
        var = 1/(2-fraction)
        assert_allclose(w@x,mean,atol=1e-12)
        assert_allclose(w@((x-mean)**2),var,atol=1e-12)


def toy_models():
    return [dict(rho=np.array([rho]),counts=np.array([3]),meta={"logz":0.},
                 loglike=np.array([[-i,-2*i]],float)) for i,rho in enumerate([0.,.5,2.])]


def test_model_evidence_updates_when_data_are_weakened():
    rows = toy_models()
    result = analysis.distribution(rows,np.array([True,False]),.2)
    expected = analysis.model_probabilities([0,.2,.4],analysis.PRIORS)
    assert_allclose(result["model_probabilities"],expected)
    assert_allclose(result["probability_rho_gt_1"],expected[2])
    # A data-only normalisation shift cancels between the models.
    for row in rows:
        row["loglike"] += 100
    shifted = analysis.distribution(rows,np.array([True,False]),.2)
    assert_allclose(shifted["model_probabilities"],expected)


def test_bin_derivative_matches_reweighting_and_adds_by_dataset():
    rows = toy_models()
    full,conditional = analysis.local_tail_sensitivity(rows)
    base = analysis.distribution(rows,np.array([False,False]),0.)
    eps = 1e-6
    for mask in (np.array([True,False]),np.array([False,True]),np.array([True,True])):
        perturbed = analysis.distribution(rows,mask,eps)
        fd = (perturbed["probability_rho_gt_1"]-base["probability_rho_gt_1"])/eps
        assert_allclose(fd,full[mask].sum(),atol=2e-6)
        p0 = base["probability_rho_gt_1"]/base["probability_dm"]
        p1 = perturbed["probability_rho_gt_1"]/perturbed["probability_dm"]
        assert_allclose((p1-p0)/eps,conditional[mask].sum(),atol=2e-6)
