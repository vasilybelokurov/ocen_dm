"""The JamPy and AGAMA back-ends must reproduce our Jeans solver where the physics coincides."""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.kinematics import Anisotropy, SphericalJeans
from ocen_dm.kinematics.backends import mge_from_component, project_density
from ocen_dm.mass_models import MGE, CompositeMassModel, PointMass, RemnantPlummer, TruncatedGNFW

D = 5.43
STARS = MGE([0.3e6, 1.5e6, 1.2e6], [1.0, 3.5, 12.0])
REM = RemnantPlummer(2e5, 1.0)
DM = TruncatedGNFW.cored(rho_s=0.05, r_s=30.0, r_t=1000.0)
MASS = CompositeMassModel([STARS, PointMass(1e4), REM, DM])
R = np.geomspace(0.3, 60, 8)


def test_projection_of_plummer_is_analytic():
    Rp = np.array([0.5, 1.0, 3.0, 10.0])
    exact = REM.mass * REM.scale**2 / (np.pi * (Rp**2 + REM.scale**2) ** 2)
    np.testing.assert_allclose(project_density(REM, Rp), exact, rtol=1e-6)


def test_projected_mge_of_non_gaussian_components_conserves_enclosed_mass():
    r = np.geomspace(0.1, 1000, 8)
    for comp in (REM, DM):
        mge = mge_from_component(comp, np.geomspace(0.05, 500.0, 40))
        ratio = np.asarray(mge.enclosed_mass(r)) / np.asarray(comp.enclosed_mass(r))
        assert np.all(np.abs(ratio - 1) < 0.01), (comp.name, ratio)


jampy = pytest.importorskip("jampy", reason="JamPy not installed")


def test_jam_backend_matches_our_solver_on_a_full_composite_model():
    from ocen_dm.kinematics.backends import JamBackend
    ani = Anisotropy(-0.2, 0.4, 5.0)
    ours = SphericalJeans(MASS, STARS, ani).dispersions_kms(R)
    jam = JamBackend(MASS, STARS, ani, D)
    theirs = jam.dispersions_kms(R)
    for k in ("los", "pmr", "pmt"):
        np.testing.assert_allclose(theirs[k], ours[k], rtol=3e-3)
    m = jam.projected_moments(R)
    np.testing.assert_allclose(m["Sigma"], SphericalJeans(MASS, STARS, ani).surface_density(R), rtol=2e-3)


agama = pytest.importorskip("agama", reason="AGAMA not installed")


def test_agama_df_backend_matches_our_solver_for_constant_beta():
    from ocen_dm.kinematics.backends import AgamaDFBackend
    for b0 in (0.0, -0.3):
        ours = SphericalJeans(MASS, STARS, Anisotropy(b0, b0, 1.0)).dispersions_kms(R)
        df = AgamaDFBackend(MASS, STARS, beta0=b0, r_a=np.inf, distance_kpc=D)
        theirs = df.dispersions_kms(R)
        for k in ("los", "pmr", "pmt"):
            np.testing.assert_allclose(theirs[k], ours[k], rtol=5e-3)
        # a positive DF exists: the realised density is the input density
        np.testing.assert_allclose(df.density_check(np.geomspace(0.3, 30, 6)), 1.0, rtol=2e-2)


def test_agama_density_check_flags_an_unphysical_jeans_model():
    """Radial anisotropy in a cored tracer has no positive DF (An & Evans 2006)."""
    from ocen_dm.kinematics.backends import AgamaDFBackend
    df = AgamaDFBackend(MASS, STARS, beta0=0.3, r_a=np.inf, distance_kpc=D)
    ratio = df.density_check(np.array([0.3, 0.5]))
    assert np.any(np.abs(ratio - 1) > 0.05)
