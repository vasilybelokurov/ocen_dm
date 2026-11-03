"""Tests for ocen_dm.centre_marg (docs/CENTRE_MARGINALISATION.md, T1, T3 and helpers).

Synthetic tests always run; tests marked ``needs_data`` need the oMEGACat catalogue and a frozen
DF run (results/df/dftwo_counts_area_20261013) and are skipped otherwise.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from ocen_dm import centre_marg as cm

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT/"results/df/dftwo_counts_area_20261013/observed"
CAT = ROOT/"data/raw/omegacat_vi_kinematics/catalog_and_selections.fits"
needs_data = pytest.mark.skipif(not (RUN.exists() and CAT.exists()), reason="frozen run or catalogue missing")


def test_grid_moments():
    for s in (1., 2., 3.):
        c, w = cm.centre_grid(s, 6)
        assert w.shape == (36,) and np.isclose(w.sum(), 1.)
        assert np.allclose(w@c, 0., atol=1e-12)
        assert np.allclose(w@(c**2), s**2)                 # exact for degree <= 11
    c, w = cm.centre_grid(0., 6)
    assert c.shape == (1, 2) and w[0] == 1.


def test_marginal_delta_and_constant():
    lnl = np.array([-3., -1., -2.])
    assert np.isclose(cm.marginal_loglike(lnl, [0., 1., 0.]), -1.)
    assert np.isclose(cm.marginal_loglike(np.full(4, -7.), np.full(4, .25)), -7.)


def test_ml_sigma_recovery():
    rng = np.random.default_rng(1)
    n, sig = 4000, 18.
    phi = rng.uniform(-np.pi, np.pi, n)
    ev = rng.uniform(1., 6., n)
    v = 3.+4.*np.sin(phi-.7)+rng.normal(0, sig, n)+rng.normal(0, ev)
    est = cm._ml_sigma(v, ev, phi)
    assert abs(est-sig) < 3*sig/np.sqrt(2*n)


def test_quadrature_matches_builder():
    spec = importlib.util.spec_from_file_location("bcp", ROOT/"bin/build_count_profiles.py")
    bcp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bcp)
    f = lambda r1, r2: .9*np.pi*(r2**2-r1**2)
    for args in ((2., 2.55, None), (58., 74., f)):
        a, b = bcp.quadrature(*args), cm._quadrature(*args)
        assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]) and a[2] == b[2]


def _frozen():
    from ocen_dm.kinematics.counts import CountProfile
    from ocen_dm.kinematics.run_io import read_data_snapshot
    rec = json.loads((RUN/"problem.json").read_text())
    counts = [CountProfile.from_dict(json.loads((RUN/f"{n}.json").read_text())) for n in rec["counts"]]
    return read_data_snapshot(RUN, rec["data_snapshot"]), counts


@needs_data
def test_T1_counts_reproduce_frozen():
    _, counts = _frozen()
    for c in counts:
        new = cm.hst_counts(c) if c.name.startswith(cm.HST_NAME) else cm.gaia_counts(c)
        assert np.array_equal(new.counts, c.counts), c.name
        assert np.allclose(new.area_arcsec2, c.area_arcsec2, rtol=1e-12, atol=0), c.name
        assert np.allclose(new.r_nodes, c.r_nodes, rtol=1e-12, atol=0), c.name
        assert np.allclose(new.node_weights, c.node_weights, rtol=1e-12, atol=1e-15), c.name


@needs_data
def test_T1_kinematics_reproduce_frozen():
    data, counts = _frozen()
    cd = cm.CentreData(data, counts)
    assert cd.n_hst == 6 and cd.n_muse == 3
    new = {p.name: p for p in cd.kinematics((0., 0.)).profiles}
    for p in data.profiles:
        q = new[p.name]
        for k in ("r", "value"):
            assert np.allclose(getattr(q, k), getattr(p, k), rtol=1e-12, atol=0), (p.name, k)
        for k in ("err_lo", "err_hi"):           # mixture-fit Hessian errors: ~5e-10 relative noise
            assert np.allclose(getattr(q, k), getattr(p, k), rtol=1e-8, atol=0), (p.name, k)
        if p.r_nodes is not None:
            assert np.allclose(q.r_nodes, p.r_nodes, rtol=1e-12, atol=0), p.name


@needs_data
def test_T5_T6_offset_counts_conserve_and_cover():
    _, counts = _frozen()
    hst = next(c for c in counts if c.name.startswith(cm.HST_NAME))
    a = cm.hst_counts(hst, (2., -2.))
    cov0 = hst.area_arcsec2/(np.pi*(hst.r_upper**2-hst.r_lower**2))
    cov1 = a.area_arcsec2/(np.pi*(a.r_upper**2-a.r_lower**2))
    assert np.max(np.abs(cov1-cov0)) < 5e-3
    # total within 2-250": only stars crossing the 2" and 250" circles can change it
    assert abs(int(a.counts.sum())-int(hst.counts.sum())) < 200
    assert not np.array_equal(a.counts, hst.counts)
