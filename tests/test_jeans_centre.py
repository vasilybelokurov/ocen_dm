"""Sampled model centre for the Jeans fits (FitProblem + centre_model.CentreGeometry).

Synthetic data and star positions, so no catalogue is needed. The real-data accuracy of the
speed options (cells, n_moments, cont_kwargs) is recorded in JOURNAL 2026-11-12.
"""
from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.centre_model import CentreGeometry, PointSet, projected_variance, ModelMoments
from ocen_dm.kinematics import DarkMatterModel, FitProblem, KinematicData, NoDarkMatterModel
from ocen_dm.kinematics.likelihood import BinnedProfile
from ocen_dm.kinematics.run_io import family_config, family_from_config
from ocen_dm.light_model import MGEFit

_MGE = MGEFit(sigmas_arcsec=np.array([60.0, 250.0]), fractions=np.array([0.7, 0.3]),
              rms_mag=0.0, max_abs_resid_mag=0.0, n_grid=2)
EDGES = np.array([1., 5., 15., 40., 80.])        # bins 0-1 start inside r_cont = 15": continuous


def _profile(name, kind, edges, value, inst, rng):
    lo, hi = edges[:-1], edges[1:]
    r = 0.5*(lo+hi)
    nodes = np.sort(rng.uniform(lo[:, None], hi[:, None], (len(r), 8)), axis=1)
    v = np.full(len(r), value)
    return BinnedProfile(name, kind, r, lo, hi, v, 0.05*v, 0.06*v, inst, r_nodes=nodes)


def _data():
    rng = np.random.default_rng(0)
    return KinematicData((
        _profile("hst_pm_radial_ours", "pmr", EDGES, 0.6, "HST", rng),
        _profile("hst_pm_tangential_ours", "pmt", EDGES, 0.6, "HST", rng),
        _profile("muse_los_dispersion", "los", EDGES, 15., "MUSE", rng),
        _profile("gaia_edr3_ours_radial", "pmr", np.array([300., 450., 600.]), 0.4, "GaiaEDR3", rng),
    ))


def _stars(n=60000, seed=1):
    """Uniform disc to 85" with a 60-degree wedge missing (an incomplete footprint)."""
    rng = np.random.default_rng(seed)
    r = 85*np.sqrt(rng.uniform(size=n))
    p = rng.uniform(0, 2*np.pi, n)
    keep = ~((p > 1.0) & (p < 1.0+np.pi/3))
    return r[keep]*np.sin(p[keep]), r[keep]*np.cos(p[keep])


def _geometry(data, **kw):
    s = _stars()
    return CentreGeometry(data, (), hst_stars=s, muse_stars=_stars(8000, 2), **kw)


def _problem(sigma=2.0, **opts):
    data = _data()
    fam = NoDarkMatterModel(mge_fit=_MGE, centre_sigma_arcsec=sigma, centre_options=opts or None)
    o = fam.centre_options
    return FitProblem(fam, data, geometry=_geometry(data, cells=o["cells"], cont_kwargs=o["cont_kwargs"]))


def _x(fam, u=0.5, c=(0., 0.)):
    x = fam.transform(np.full(len(fam.names), u))
    x[fam.names.index("centre_x")], x[fam.names.index("centre_y")] = c
    return x


def test_default_family_has_fixed_centre_and_legacy_config():
    fam = NoDarkMatterModel(mge_fit=_MGE)
    assert "centre_x" not in fam.names and fam.centre({}) is None
    cfg = family_config(fam)
    assert "centre_sigma_arcsec" not in cfg and "centre_options" not in cfg


def test_centre_parameters_priors_label_and_roundtrip():
    fam = NoDarkMatterModel(mge_fit=_MGE, centre_sigma_arcsec=2.0)
    assert fam.names[-2:] == ("centre_x", "centre_y") and fam.label.endswith("_centre2")
    for p in fam.parameters[-2:]:
        assert (p.prior.kind, p.prior.mu, p.prior.sigma) == ("normal", 0.0, 2.0)
    assert family_config(family_from_config(family_config(fam))) == family_config(fam)
    dm = DarkMatterModel(gamma=1.0, mge_fit=_MGE, centre_sigma_arcsec=3.0)
    assert dm.label == "K2_nfw_centre3" and family_from_config(family_config(dm)).names == dm.names
    with pytest.raises(ValueError):
        NoDarkMatterModel(mge_fit=_MGE, centre_sigma_arcsec=2.0, backend="jam")
    with pytest.raises(ValueError):
        NoDarkMatterModel(mge_fit=_MGE, centre_sigma_arcsec=2.0, centre_options={"bogus": 1})


def test_zero_centre_reproduces_fixed_likelihood_exactly():
    pb = _problem()
    fixed = FitProblem(NoDarkMatterModel(mge_fit=_MGE), pb.data)
    for u in (0.2, 0.5, 0.8):
        x = _x(pb.family, u)
        assert pb.loglike_vector(x) == fixed.loglike_vector(x[:-2])
        assert pb.loglike_vector(x) > -1e50


def test_centre_moves_hst_and_muse_not_gaia_and_is_continuous():
    pb = _problem()
    p0 = pb.predict(_x(pb.family))
    p1 = pb.predict(_x(pb.family, c=(3., -1.)))
    assert np.array_equal(p0["gaia_edr3_ours_radial"], p1["gaia_edr3_ours_radial"])
    for n in ("hst_pm_radial_ours", "hst_pm_tangential_ours", "muse_los_dispersion"):
        assert np.max(np.abs(p1[n]/p0[n]-1)) > 1e-4
    l0 = pb.loglike_vector(_x(pb.family))
    assert abs(pb.loglike_vector(_x(pb.family, c=(1e-7, 0.))) - l0) < 1e-5


def test_cell_compression_matches_every_star():
    data = _data()
    full, cells = _geometry(data), _geometry(data, cells=(8, 128))
    mom = ModelMoments.from_functions(Sigma=lambda r: (1+(r/30)**2)**-2, s2r=lambda r: 1+50/(r+1),
                                      s2t=lambda r: 0.8+40/(r+1), s2los=lambda r: 0.9+45/(r+1))
    for c in ((2., 1.), (-4., 3.)):
        ff, fc = full.factors(mom, c), cells.factors(mom, c)
        for n in ff:
            k = slice(full.cont[n].n_bins, None)               # the discrete bins
            change = np.abs(ff[n][k]-1)
            assert np.all(np.abs(fc[n][k]-ff[n][k]) < 0.01*change+1e-6), n


def test_default_geometry_reproduces_full_star_means():
    """The DF path (no cells): discrete bins are the plain means over every star, as before."""
    data = _data()
    g = _geometry(data)
    mom = ModelMoments.from_functions(Sigma=lambda r: 1/(1+r), s2r=lambda r: 1+10/(r+1),
                                      s2t=lambda r: 1+5/(r+1), s2los=lambda r: 1+7/(r+1))
    c = (1.5, -2.)
    for n, comp in (("hst_pm_radial_ours", "r"), ("muse_los_dispersion", "los")):
        ps = g.sets[n]
        with np.errstate(invalid="ignore", divide="ignore"):
            expect = ps.mean(projected_variance(mom, ps.x, ps.y, c, comp))
        m = g.cont[n].n_bins
        assert np.array_equal(g._means(mom, n, c)[m:], expect[m:])
