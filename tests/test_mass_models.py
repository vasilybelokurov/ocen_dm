"""Tests for the mass-profile library (specification Milestone 2, section 13).

Every analytic profile is checked against independent numerics, and the two
components with an external implementation are checked against AGAMA.
"""

from __future__ import annotations

import numpy as np
import pytest

from ocen_dm.mass_models import (
    MGE,
    Burkert,
    CompositeMassModel,
    Plummer,
    PointMass,
    RemnantPlummer,
    TruncatedGNFW,
)
from ocen_dm.mass_models.base import G, MassComponent, dphi_dr

RADII = np.geomspace(1e-3, 1e4, 40)


# --------------------------------------------------------------------- units --
def test_g_matches_astropy():
    """The module constant must not drift from CODATA via astropy."""
    from astropy import constants as const
    from astropy import units as u

    reference = (const.G).to(u.pc * u.km**2 / u.s**2 / u.Msun).value
    assert abs(G - reference) / reference < 1e-12


def test_circular_velocity_of_a_point_mass_is_keplerian():
    bh = PointMass(1e4)
    r = np.array([1.0, 10.0, 100.0])
    assert np.allclose(bh.circular_velocity(r), np.sqrt(G * 1e4 / r), rtol=1e-12)


def test_units_are_self_consistent_in_km_per_s():
    """A 1e6 Msun Plummer with a 5 pc scale must give tens of km/s, not m/s or kpc/Gyr."""
    v = Plummer(1e6, 5.0).circular_velocity(5.0)[0]
    assert 10.0 < v < 40.0


# ------------------------------------------------------------------- Plummer --
def test_plummer_analytic_matches_quadrature():
    report = Plummer(1e6, 5.0).verify()
    assert report.passed, report


def test_plummer_enclosed_mass_matches_the_closed_form():
    p = Plummer(1e6, 5.0)
    r = np.array([1.0, 5.0, 50.0])
    expected = 1e6 * r**3 / (r**2 + 25.0) ** 1.5
    assert np.allclose(p.enclosed_mass(r), expected, rtol=1e-14)


def test_plummer_half_mass_radius():
    p = Plummer(1e6, 5.0)
    assert np.isclose(p.enclosed_mass(p.half_mass_radius)[0], 0.5e6, rtol=1e-12)


def test_plummer_tends_to_its_total_mass():
    p = Plummer(1e6, 5.0)
    assert np.isclose(p.enclosed_mass(1e8)[0], 1e6, rtol=1e-8)
    assert p.total_mass == 1e6


def test_plummer_rejects_unphysical_parameters():
    with pytest.raises(ValueError):
        Plummer(-1.0, 5.0)
    with pytest.raises(ValueError):
        Plummer(1e6, 0.0)


# ----------------------------------------------------------------------- MGE --
def test_mge_single_gaussian_matches_the_closed_form():
    from scipy import special

    mge = MGE([1e6], [10.0])
    r = np.array([1.0, 10.0, 100.0])
    x = r / 10.0
    expected = 1e6 * (special.erf(x / np.sqrt(2)) - np.sqrt(2 / np.pi) * x * np.exp(-0.5 * x**2))
    assert np.allclose(mge.enclosed_mass(r), expected, rtol=1e-12)


def test_mge_is_additive():
    a, b = MGE([6e5], [2.0]), MGE([4e5], [20.0])
    both = MGE([6e5, 4e5], [2.0, 20.0])
    assert np.allclose(both.enclosed_mass(RADII), a.enclosed_mass(RADII) + b.enclosed_mass(RADII))
    assert np.allclose(both.potential(RADII), a.potential(RADII) + b.potential(RADII))


def test_mge_mass_to_light_scales_the_whole_profile():
    base = MGE([1e6], [5.0])
    scaled = MGE([1e6], [5.0], mass_to_light=2.0)
    assert np.allclose(scaled.enclosed_mass(RADII), 2.0 * base.enclosed_mass(RADII))
    assert scaled.total_mass == 2e6


def test_mge_verifies():
    assert MGE([6e5, 4e5], [2.0, 20.0]).verify().passed


def test_mge_potential_is_finite_at_the_centre():
    assert np.isfinite(MGE([1e6], [5.0]).potential(0.0)[0])


def test_mge_rejects_mismatched_inputs():
    with pytest.raises(ValueError):
        MGE([1e6, 1e5], [5.0])
    with pytest.raises(ValueError):
        MGE([], [])
    with pytest.raises(ValueError):
        MGE([1e6], [-5.0])


# ---------------------------------------------------------------- point mass --
def test_point_mass_encloses_its_mass_everywhere_outside_the_origin():
    bh = PointMass(1e4)
    assert np.allclose(bh.enclosed_mass([1e-6, 1.0, 1e6]), 1e4)
    assert bh.enclosed_mass(0.0)[0] == 0.0


def test_point_mass_potential_and_force():
    bh = PointMass(1e4)
    r = np.array([2.0, 20.0])
    assert np.allclose(bh.potential(r), -G * 1e4 / r, rtol=1e-14)
    assert np.allclose(bh.force(r), -G * 1e4 / r**2, rtol=1e-14)


def test_zero_imbh_is_allowed_and_contributes_nothing():
    bh = PointMass(0.0)
    assert np.all(bh.enclosed_mass(RADII) == 0.0)
    assert np.all(bh.potential(RADII) == 0.0)


def test_negative_imbh_is_rejected():
    with pytest.raises(ValueError):
        PointMass(-1.0)


# ------------------------------------------------------------- dark matter ----
def test_nfw_closed_form_matches_quadrature():
    assert TruncatedGNFW.nfw(0.05, 50.0).verify().passed


def test_gnfw_gamma_one_equals_the_analytic_nfw_mass():
    nfw = TruncatedGNFW.nfw(0.05, 50.0)
    r = np.array([5.0, 50.0, 500.0])
    x = r / 50.0
    expected = 4 * np.pi * 0.05 * 50.0**3 * (np.log1p(x) - x / (1 + x))
    assert np.allclose(nfw.enclosed_mass(r), expected, rtol=1e-12)


def test_cored_and_cuspy_differ_in_the_centre_but_not_the_outskirts():
    cusp = TruncatedGNFW.nfw(0.05, 50.0, r_t=500.0)
    core = TruncatedGNFW.cored(0.05, 50.0, r_t=500.0)
    assert cusp.density(0.1)[0] > 10 * core.density(0.1)[0]
    assert cusp.density(400.0)[0] == pytest.approx(core.density(400.0)[0], rel=0.3)


def test_truncation_removes_mass_and_makes_the_total_finite():
    free = TruncatedGNFW.nfw(0.05, 50.0)
    cut = TruncatedGNFW.nfw(0.05, 50.0, r_t=200.0)
    assert cut.enclosed_mass(1000.0)[0] < free.enclosed_mass(1000.0)[0]
    assert np.isinf(free.total_mass)
    assert np.isfinite(cut.total_mass)


def test_truncated_profiles_verify():
    for component in (
        TruncatedGNFW.nfw(0.05, 50.0, r_t=300.0),
        TruncatedGNFW.cored(0.05, 50.0, r_t=300.0),
        Burkert(0.1, 30.0, r_t=300.0),
    ):
        assert component.verify().passed, component.verify()


def test_burkert_closed_form_matches_quadrature():
    burkert = Burkert(0.1, 30.0)
    r = np.array([3.0, 30.0, 300.0])
    numeric = MassComponent.enclosed_mass(burkert, r)
    assert np.allclose(burkert.enclosed_mass(r), numeric, rtol=1e-6)


def test_burkert_has_a_flat_core():
    burkert = Burkert(0.1, 30.0)
    assert burkert.density(1e-3)[0] == pytest.approx(0.1, rel=1e-3)


def test_dm_parameters_are_validated():
    with pytest.raises(ValueError):
        TruncatedGNFW(0.05, 50.0, gamma=3.0)
    with pytest.raises(ValueError):
        TruncatedGNFW(0.05, -1.0)
    with pytest.raises(ValueError):
        Burkert(0.1, 0.0)


# ------------------------------------------------------------------ composite --
def _model(with_dm: bool = True, with_bh: bool = True) -> CompositeMassModel:
    parts = [MGE([2.5e6, 1.0e6], [3.0, 25.0]), RemnantPlummer(2e5, 2.0)]
    if with_bh:
        parts.append(PointMass(5e3))
    if with_dm:
        parts.append(TruncatedGNFW.nfw(0.02, 60.0, r_t=300.0))
    return CompositeMassModel(parts)


def test_composite_mass_and_potential_are_additive():
    model = _model()
    total_mass = sum(c.enclosed_mass(RADII) for c in model.components)
    total_phi = sum(c.potential(RADII) for c in model.components)
    assert np.allclose(model.enclosed_mass(RADII), total_mass)
    assert np.allclose(model.potential(RADII), total_phi)


def test_composite_verifies_including_its_force():
    reports = _model().verify_all()
    assert all(r.passed for r in reports), [str(r) for r in reports]


def test_mass_breakdown_sums_to_the_total():
    breakdown = _model().mass_breakdown(RADII)
    parts = sum(v for k, v in breakdown.items() if k != "total")
    assert np.allclose(parts, breakdown["total"])


def test_dark_matter_mass_selects_only_dm_components():
    model = _model()
    assert np.allclose(model.dark_matter_mass(RADII), model["dm_nfw"].enclosed_mass(RADII))
    assert np.all(model.dark_matter_fraction(RADII) <= 1.0)
    assert np.all(model.dark_matter_fraction(RADII) >= 0.0)


def test_a_model_without_dark_matter_reports_zero():
    model = _model(with_dm=False)
    assert np.all(model.dark_matter_mass(RADII) == 0.0)
    assert np.all(model.dark_matter_fraction(RADII) == 0.0)


def test_omitting_the_imbh_equals_including_a_zero_mass_one():
    without = CompositeMassModel([MGE([1e6], [5.0])])
    with_zero = CompositeMassModel([MGE([1e6], [5.0]), PointMass(0.0)])
    assert np.allclose(without.enclosed_mass(RADII), with_zero.enclosed_mass(RADII))
    assert np.allclose(without.potential(RADII), with_zero.potential(RADII))


def test_duplicate_component_names_are_rejected():
    with pytest.raises(ValueError, match="unique"):
        CompositeMassModel([Plummer(1e6, 5.0), Plummer(1e5, 1.0)])


def test_empty_composite_is_rejected():
    with pytest.raises(ValueError):
        CompositeMassModel([])


def test_component_lookup():
    model = _model()
    assert model["remnants"].total_mass == 2e5
    assert model.get("nonexistent") is None
    with pytest.raises(KeyError):
        model["nonexistent"]


# ------------------------------------------------------- derived quantities ---
def test_escape_speed_is_positive_and_decreasing():
    v = _model().escape_speed(RADII)
    assert np.all(v > 0)
    assert np.all(np.diff(v) < 0)


def test_escape_speed_to_a_finite_radius_is_smaller_than_to_infinity():
    model = _model()
    r = np.array([1.0, 10.0])
    assert np.all(model.escape_speed(r, r_esc=100.0) < model.escape_speed(r))


def test_escape_speed_vanishes_at_the_reference_radius():
    model = _model()
    assert model.escape_speed(100.0, r_esc=100.0)[0] == pytest.approx(0.0, abs=1e-9)


def test_jacobi_radius_solves_its_defining_equation():
    model = _model()
    R, M_gal = 6000.0, 5.0e10          # 6 kpc, a plausible enclosed Galactic mass
    r_j = model.jacobi_radius(R, M_gal)
    predicted = R * (model.enclosed_mass(r_j)[0] / (3 * M_gal)) ** (1 / 3)
    assert r_j == pytest.approx(predicted, rel=1e-8)


def test_jacobi_radius_grows_with_galactocentric_distance():
    model = _model()
    near = model.jacobi_radius(4000.0, 4.0e10)
    far = model.jacobi_radius(8000.0, 8.0e10)
    assert far > near


def test_jacobi_radius_rejects_bad_input():
    model = _model()
    with pytest.raises(ValueError):
        model.jacobi_radius(-1.0, 1e10)
    with pytest.raises(ValueError):
        model.jacobi_radius(6000.0, 1e10, bracket=(1e-6, 1e-5))


def test_radial_profile_carries_the_quantities_the_spec_asks_for():
    profile = _model().radial_profile(RADII)
    for key in ("r", "M_total", "M_dm", "f_dm", "rho_total", "v_circ", "v_esc"):
        assert key in profile
        assert np.all(np.isfinite(profile[key]))
    assert "M_dm_nfw" in profile and "M_remnants" in profile and "M_imbh" in profile


def test_radius_enclosing_inverts_the_mass_profile():
    model = _model()
    target = 1.0e6
    r = model.radius_enclosing(target)
    assert model.enclosed_mass(r)[0] == pytest.approx(target, rel=1e-8)


def test_radius_enclosing_reports_an_unreachable_mass():
    with pytest.raises(ValueError, match="encloses only"):
        Plummer(1e6, 5.0).radius_enclosing(1e9)


# ------------------------------------------------------------------- hygiene --
def test_negative_radii_are_rejected():
    with pytest.raises(ValueError):
        Plummer(1e6, 5.0).density(-1.0)


def test_all_components_are_finite_over_six_decades():
    for component in (
        MGE([1e6], [5.0]), RemnantPlummer(1e5, 2.0), PointMass(1e4),
        TruncatedGNFW.nfw(0.02, 60.0, r_t=300.0), Burkert(0.1, 30.0, r_t=300.0),
    ):
        mass = component.enclosed_mass(RADII)
        assert np.all(np.isfinite(mass)) and np.all(mass >= 0)


# ---------------------------------------------------------- external check ----
agama = pytest.importorskip("agama", reason="AGAMA not installed")


def test_plummer_matches_agama():
    """Independent implementation, no shared code."""
    agama.setUnits(mass=1, length=1e-3, velocity=1)   # Msun, pc, km/s
    r = np.geomspace(0.1, 500.0, 15)
    xyz = np.column_stack([r, np.zeros_like(r), np.zeros_like(r)])
    ours = Plummer(1e6, 5.0)
    theirs = agama.Potential(type="Plummer", mass=1e6, scaleRadius=5.0)
    # AGAMA carries a slightly different value of G, which sets the floor here.
    assert np.allclose(ours.enclosed_mass(r), theirs.enclosedMass(r), rtol=1e-10)
    assert np.allclose(ours.potential(r), theirs.potential(xyz), rtol=1e-8)


def test_nfw_matches_agama():
    agama.setUnits(mass=1, length=1e-3, velocity=1)
    r = np.geomspace(0.1, 500.0, 15)
    xyz = np.column_stack([r, np.zeros_like(r), np.zeros_like(r)])
    ours = TruncatedGNFW.nfw(0.05, 50.0)
    theirs = agama.Potential(
        type="Spheroid", densityNorm=0.05, scaleRadius=50.0, gamma=1, beta=3, alpha=1
    )
    assert np.allclose(ours.density(r), theirs.density(xyz), rtol=1e-5)
    assert np.allclose(ours.enclosed_mass(r), theirs.enclosedMass(r), rtol=1e-5)
    assert np.allclose(ours.potential(r), theirs.potential(xyz), rtol=1e-6)


# ------------------------------------------------------------ fast path ------
class TestProfileTables:
    """The spline tables must match quadrature and be cheap enough for inference."""

    def test_tables_match_quadrature_for_every_numeric_profile(self):
        for component in (
            TruncatedGNFW.nfw(0.05, 50.0, r_t=300.0),
            TruncatedGNFW.cored(0.05, 50.0, r_t=300.0),
            TruncatedGNFW(0.05, 50.0, gamma=0.5, r_t=300.0),
            TruncatedGNFW(0.05, 50.0, gamma=1.5, r_t=1000.0),
            Burkert(0.1, 30.0, r_t=300.0),
            Burkert(0.1, 30.0),
        ):
            r = np.geomspace(1e-3, 1e4, 30)
            fast = component.enclosed_mass(r)
            slow = component.quad_enclosed_mass(r)
            assert np.allclose(fast, slow, rtol=1e-6), component
            assert np.allclose(component.potential(r), component.quad_potential(r), rtol=1e-6)

    def test_tables_are_built_once_per_instance(self):
        component = TruncatedGNFW.nfw(0.05, 50.0, r_t=300.0)
        assert component.tables is component.tables

    def test_fresh_parameters_get_fresh_tables(self):
        a = TruncatedGNFW.nfw(0.05, 50.0, r_t=300.0)
        b = TruncatedGNFW.nfw(0.06, 50.0, r_t=300.0)
        assert a.enclosed_mass(100.0)[0] < b.enclosed_mass(100.0)[0]

    def test_total_mass_of_a_truncated_halo_is_finite_and_consistent(self):
        component = TruncatedGNFW.nfw(0.05, 50.0, r_t=300.0)
        total = component.total_mass
        assert np.isfinite(total)
        assert component.enclosed_mass(1e5)[0] == pytest.approx(total, rel=1e-4)
        assert component.enclosed_mass(np.inf)[0] == total

    def test_radii_below_the_grid_follow_the_inner_power_law(self):
        core = TruncatedGNFW.cored(0.05, 50.0, r_t=300.0)   # rho -> const, M ~ r^3
        m = core.enclosed_mass([1e-7, 2e-7])
        assert m[1] / m[0] == pytest.approx(8.0, rel=1e-3)

    def test_radii_above_the_grid_fall_back_to_quadrature(self):
        component = Burkert(0.1, 30.0, r_t=300.0)
        r = np.array([5e6])
        assert np.isfinite(component.enclosed_mass(r)[0])
        assert component.enclosed_mass(r)[0] == pytest.approx(component.quad_enclosed_mass(r)[0], rel=1e-6)

    def test_a_non_convergent_potential_is_refused(self):
        from ocen_dm.mass_models.base import ProfileTables

        with pytest.raises(ValueError, match="does not converge"):
            ProfileTables(lambda r: 1.0 / (1.0 + r))      # rho ~ r^-1 at large r

    def test_composite_evaluation_is_fast_enough_for_nested_sampling(self):
        """Fresh DM parameters every call, as a sampler would do: under 5 ms each."""
        import time

        r = np.geomspace(0.1, 300.0, 40)
        start = time.perf_counter()
        n = 10
        for i in range(n):
            model = CompositeMassModel([
                MGE([2.2e6, 1.3e6], [3.0, 12.0]),
                RemnantPlummer(1.5e5, 2.0),
                TruncatedGNFW.nfw(0.02 * (1 + 1e-3 * i), 60.0, r_t=300.0),
            ])
            model.enclosed_mass(r)
            model.potential(r)
        per_call = (time.perf_counter() - start) / n
        assert per_call < 5e-3, f"{per_call * 1e3:.1f} ms per evaluation"

    def test_analytic_components_bypass_the_tables(self):
        """Plummer and untruncated NFW have closed forms and must not pay for a table."""
        assert "tables" not in Plummer(1e6, 5.0).__dict__
        Plummer(1e6, 5.0).enclosed_mass(RADII)
        assert "tables" not in Plummer(1e6, 5.0).__dict__
