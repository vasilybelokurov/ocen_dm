"""Analytic and discrete-mixture checks for the density posterior report."""
import importlib.util
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest

spec = importlib.util.spec_from_file_location("dm_density_posterior",
    Path(__file__).resolve().parents[1]/"bin/plot_dm_density_posterior.py")
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def test_density_matches_analytic_untruncated_nfw():
    rs = np.array([3.,20.,1000.])
    m100 = np.array([1e3,1e5,3e7])
    x = 100/rs
    norm = m100/(4*np.pi*rs**3*(np.log1p(x)-x/(1+x)))
    expected = norm/(20/rs*(1+20/rs)**2)
    assert_allclose(analysis.density_at_radius(m100,rs,1.,np.inf),expected,rtol=1e-11)


def test_mixture_keeps_zero_atom_and_is_independent_of_draw_counts():
    values = [np.array([0.]),np.array([1.,2.]),np.array([4.,4.,4.,4.])]
    probability = np.array([.5,.25,.25])
    sample,weight = analysis.mixture_samples(values,probability)
    assert_allclose(weight.sum(),1.)
    assert_allclose(weight[sample==0].sum(),.5)
    assert_allclose(analysis.weighted_quantiles(sample,weight,[.1,.5,.6,.75,.8,.95]),[0,0,1,2,4,4])
    doubled = [np.tile(v,2) for v in values]
    a,b = analysis.mixture_samples(doubled,probability)
    assert_allclose(analysis.weighted_quantiles(a,b),analysis.weighted_quantiles(sample,weight))


def test_evidence_weights_respect_priors_and_log_offset():
    prior = np.array([.5,.25,.25])
    assert_allclose(analysis.model_probabilities([0,0,0],prior),prior)
    assert_allclose(analysis.model_probabilities([1000,1000,1000],prior),prior)
    assert_allclose(analysis.model_probabilities([0,np.log(2),np.log(2)],prior),np.ones(3)/3)
    with pytest.raises(ValueError):
        analysis.model_probabilities([0,0,0],[1,1,1])
