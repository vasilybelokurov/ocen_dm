"""The mixture estimator must recover injected means and dispersions, including streaming.

Written after a Codex review found that the ``cad`` terms in the mean-update normal equations
carried the wrong sign (2026-09-19). The existing mocks all had zero streaming and isotropic
dispersions, where ``cad`` vanishes, so nothing caught it. These tests inject both.
"""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics.outer_profile import MemberSample, dispersion_2d


def _sample(n=5000, sr=0.50, st=0.25, mr=0.05, mt=-0.20, err=0.05, seed=7):
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, 2 * np.pi, n)
    c, s = np.cos(phi), np.sin(phi)
    vr = rng.normal(mr, sr, n); vt = rng.normal(mt, st, n)
    a = vr * c - vt * s + rng.normal(0, err, n)
    d = vr * s + vt * c + rng.normal(0, err, n)
    samp = MemberSample(np.full(n, 500.0), vr, vt, np.full(n, err), np.full(n, err),
                        np.ones(n), np.full(n, 18.0), np.full(n, 2, int), np.arctan2(c, s),
                        mu_a=a, mu_d=d, err_a=np.full(n, err), err_d=np.full(n, err),
                        err_corr=np.zeros(n), sys_a=np.zeros(n), sys_d=np.zeros(n),
                        mu_sys=(0.0, 0.0), exact=True)
    return samp, (a, d)


def _fit(samp, ad):
    return dispersion_2d(samp, np.ones(len(samp), bool),
                         lambda x, y: np.full(np.shape(x), 1e-8), depth_var=0.0, field_at=ad)


@pytest.mark.parametrize("mt", [-0.20, 0.0, 0.30])
def test_streaming_is_recovered_under_anisotropy(mt):
    """The case the sign error broke: sigma_R != sigma_T with a non-zero tangential mean."""
    samp, ad = _sample(mt=mt)
    o = _fit(samp, ad)
    assert abs(o["mean_t"] - mt) < 0.015, "tangential mean must be recovered"
    assert abs(o["sigma_r"] - 0.50) < 0.02 and abs(o["sigma_t"] - 0.25) < 0.02


def test_anisotropy_is_recovered_both_ways():
    for sr, st in ((0.50, 0.25), (0.25, 0.50)):
        samp, ad = _sample(sr=sr, st=st, mt=-0.20)
        o = _fit(samp, ad)
        assert abs(o["sigma_r"] - sr) < 0.02 and abs(o["sigma_t"] - st) < 0.02
        assert abs(o["sigma_t"] / o["sigma_r"] - st / sr) < 0.06


def test_matches_a_direct_numerical_maximisation():
    """Independent of the analytic update: brute-force the same Gaussian likelihood."""
    from scipy.optimize import minimize
    samp, (a, d) = _sample()
    c, s = np.sin(samp.phi), np.cos(samp.phi)
    sr2, st2, e2 = 0.50 ** 2, 0.25 ** 2, 0.05 ** 2
    caa = sr2 * c ** 2 + st2 * s ** 2 + e2
    cdd = sr2 * s ** 2 + st2 * c ** 2 + e2
    cad = (sr2 - st2) * c * s
    det = caa * cdd - cad ** 2

    def nll(p):
        ma, md = p[0] * c - p[1] * s, p[0] * s + p[1] * c
        da, dd = a - ma, d - md
        return 0.5 * np.sum((cdd * da ** 2 - 2 * cad * da * dd + caa * dd ** 2) / det)

    ref = minimize(nll, [0.0, 0.0], method="Nelder-Mead",
                   options={"xatol": 1e-8, "fatol": 1e-10}).x
    o = _fit(samp, (a, d))
    assert abs(o["mean_r"] - ref[0]) < 0.01
    assert abs(o["mean_t"] - ref[1]) < 0.01


@pytest.mark.parametrize("n,tol", [(2000, 0.25), (20000, 0.20)])
def test_uncertainty_scales_with_sample_size(n, tol):
    """Guards the 2 per cent step floor that made every large bin report the same error."""
    samp, ad = _sample(n=n)
    o = _fit(samp, ad)
    expected = 0.50 / np.sqrt(2 * n)
    assert abs(o["sigma_r_err"] / expected - 1.0) < tol
    # and the floor itself is gone: 2 per cent of sigma would be 0.01
    if n >= 20000:
        assert o["sigma_r_err"] < 0.006
