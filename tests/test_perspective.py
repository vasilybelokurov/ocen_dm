"""Perspective/projection effects: exact field, and a zero-dispersion mock that must read zero."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics.perspective import (KMS_PER_MASYR_KPC, depth_dispersion, systemic_pm_field,
                                            systemic_velocity_vector, unit_vectors)
from ocen_dm.kinematics.outer_profile import MemberSample, binned_dispersion, dispersion_ml
from ocen_dm.mass_models import Plummer

RA0, DEC0, D, VLOS = 201.696833, -47.476583, 5.43, 232.6
PM0 = (-3.248, -6.746)


def test_field_matches_astropy_exactly():
    astropy = pytest.importorskip("astropy")
    from astropy.coordinates import CartesianDifferential, CartesianRepresentation, ICRS, SkyCoord
    import astropy.units as u
    v = systemic_velocity_vector(RA0, DEC0, *PM0, VLOS, D)
    rng = np.random.default_rng(0)
    ra = RA0 + rng.uniform(-1, 1, 50); dec = DEC0 + rng.uniform(-0.7, 0.7, 50)
    n, _, _ = unit_vectors(ra, dec)
    pos = CartesianRepresentation(n.T * D * u.kpc, differentials=CartesianDifferential(np.tile(v, (50, 1)).T * u.km / u.s))
    sc = SkyCoord(ICRS(pos))
    pma, pmd = systemic_pm_field(ra, dec, v, D)
    np.testing.assert_allclose(pma, sc.pm_ra_cosdec.to_value(u.mas / u.yr), atol=1e-12)
    np.testing.assert_allclose(pmd, sc.pm_dec.to_value(u.mas / u.yr), atol=1e-12)


def test_field_reduces_to_the_perspective_term_at_small_angle():
    v = systemic_velocity_vector(RA0, DEC0, *PM0, VLOS, D)
    theta = 0.5                                                   # deg, due north
    pma, pmd = systemic_pm_field(np.array([RA0]), np.array([DEC0 + theta]), v, D)
    expected = -VLOS * np.radians(theta) / (KMS_PER_MASYR_KPC * D)   # apparent contraction (approaching? sign from v_los > 0)
    assert pmd[0] - PM0[1] == pytest.approx(expected, rel=0.02)


def _mock_rigid_cluster(n=60000, sigma_int=0.0, seed=1, a_pc=15.0):
    """Stars in a Plummer sphere moving rigidly with the systemic velocity (plus an optional
    isotropic internal dispersion), observed with exact geometry including their depth."""
    rng = np.random.default_rng(seed)
    # Plummer positions in pc, centre at distance D along the line of sight
    u3 = rng.uniform(size=n); r = a_pc / np.sqrt(u3 ** (-2 / 3) - 1)
    r = r[r < 300]; n = len(r)
    dirs = rng.normal(size=(n, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    n0, e0, nth0 = unit_vectors(RA0, DEC0)
    pos = D * 1e3 * n0 + r[:, None] * dirs                          # pc, ICRS cartesian
    dist = np.linalg.norm(pos, axis=1); nhat = pos / dist[:, None]
    ra = np.degrees(np.arctan2(nhat[:, 1], nhat[:, 0])) % 360; dec = np.degrees(np.arcsin(nhat[:, 2]))
    v = systemic_velocity_vector(RA0, DEC0, *PM0, VLOS, D)
    vel = np.tile(v, (n, 1)) + rng.normal(0, sigma_int * KMS_PER_MASYR_KPC * D, size=(n, 3))
    _, e, nth = unit_vectors(ra, dec)
    vt = KMS_PER_MASYR_KPC * dist / 1e3
    pmra = np.sum(vel * e, axis=1) / vt; pmdec = np.sum(vel * nth, axis=1) / vt
    err = np.full(n, 0.02); pmra += rng.normal(0, err); pmdec += rng.normal(0, err)
    return ra, dec, pmra, pmdec, err


def _project(ra, dec, pmra, pmdec, err, exact):
    cosd = np.cos(np.radians(DEC0))
    x = (ra - RA0) * cosd * 3600; y = (dec - DEC0) * 3600; rr = np.hypot(x, y)
    cp, sp = x / rr, y / rr
    if exact:
        v = systemic_velocity_vector(RA0, DEC0, *PM0, VLOS, D)
        ea, ed = systemic_pm_field(ra, dec, v, D)
    else:
        ea, ed = PM0
    pa, pd = pmra - ea, pmdec - ed
    mu_r = pa * cp + pd * sp; mu_t = -pa * sp + pd * cp
    return MemberSample(rr, mu_r, mu_t, err, err, np.ones_like(rr), np.full_like(rr, 18.0),
                        np.full(len(rr), 3), np.arctan2(x, y), mu_sys=PM0, exact=exact)


def test_rigid_cluster_reads_zero_after_exact_correction_and_depth_term():
    ra, dec, pmra, pmdec, err = _mock_rigid_cluster(sigma_int=0.0)
    edges = np.array([600.0, 1200.0, 2000.0, 3000.0])
    tracer = Plummer(1.0, 15.0)
    exact = binned_dispersion(_project(ra, dec, pmra, pmdec, err, True), edges, prob_min=0.0,
                              depth_tracer=tracer, distance_kpc=D)
    naive = binned_dispersion(_project(ra, dec, pmra, pmdec, err, False), edges, prob_min=0.0)
    # exact + depth term: consistent with zero internal dispersion (errors are 0.02 mas/yr)
    assert np.all(np.asarray(exact["sigma_pm"]) < 0.01), exact["sigma_pm"]
    # the depth term alone is real: without it the exact projection still shows |mu_sys| sigma_z / D
    no_depth = binned_dispersion(_project(ra, dec, pmra, pmdec, err, True), edges, prob_min=0.0)
    expected = depth_dispersion(tracer, np.asarray(no_depth["r_median"]) * D * 1e3 / 206264.806, np.hypot(*PM0), D)
    assert np.allclose(np.asarray(no_depth["sigma_pm"]), expected / np.sqrt(2), rtol=0.35)
    # the naive constant subtraction leaves a spurious signal that grows with radius
    assert naive["sigma_pm"][-1] > 3 * exact["sigma_pm"][-1]
    assert naive["sigma_pm"][-1] > 0.02


def test_internal_dispersion_recovered_with_the_full_correction():
    ra, dec, pmra, pmdec, err = _mock_rigid_cluster(sigma_int=0.25, seed=2)
    edges = np.array([600.0, 1200.0, 2000.0, 3000.0])
    t = binned_dispersion(_project(ra, dec, pmra, pmdec, err, True), edges, prob_min=0.0,
                          depth_tracer=Plummer(1.0, 15.0), distance_kpc=D)
    assert np.allclose(np.asarray(t["sigma_pm"]), 0.25, atol=0.02)
