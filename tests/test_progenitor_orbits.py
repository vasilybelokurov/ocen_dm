"""Progenitor-orbit set: frame, integrators and the chosen orbits (all McMillan17, seconds)."""
import numpy as np
import pytest

po = pytest.importorskip("ocen_dm.tails.progenitor_orbits")


def test_present_day_frame_matches_astropy():
    """galpy's solarmotion is the peculiar motion only; the orbit must be retrograde."""
    import astropy.units as u, astropy.coordinates as ac
    c = ac.SkyCoord(ra=po.PRESENT["ra"] * u.deg, dec=po.PRESENT["dec"] * u.deg,
                    distance=po.PRESENT["distance_kpc"] * u.kpc,
                    pm_ra_cosdec=po.PRESENT["pmra"] * u.mas / u.yr, pm_dec=po.PRESENT["pmdec"] * u.mas / u.yr,
                    radial_velocity=po.PRESENT["vlos"] * u.km / u.s)
    vsun = ac.CartesianDifferential([-po.VSUN[0], po.VSUN[1] + po.V0_KMS, po.VSUN[2]] * u.km / u.s)
    g = c.transform_to(ac.Galactocentric(galcen_distance=po.R0_KPC * u.kpc, z_sun=po.ZSUN_KPC * u.kpc,
                                         galcen_v_sun=vsun))
    o = po.present_day_orbit()
    # galpy is left-handed in x
    assert abs(-o.x() - g.x.value) < 0.01 and abs(o.y() - g.y.value) < 0.01 and abs(o.z() - g.z.value) < 0.01
    assert abs(-o.vx() - g.v_x.value) < 0.5 and abs(o.vy() - g.v_y.value) < 0.5 and abs(o.vz() - g.v_z.value) < 0.5
    assert o.R() * o.vT() < -450          # retrograde, L_z ~ -529 kpc km/s


def test_class1_mcmillan_matches_agama_jacobi_run():
    o = po.class1_current_orbit("McMillan17", t_gyr=3.0, n=6001)
    assert abs(o.r_peri - 1.59) < 0.1 and abs(o.r_apo - 6.99) < 0.15   # 2026-09-19 AGAMA values


def test_fast_integrator_without_friction_reproduces_class1():
    f = po.class2_friction_backwards_fast(0.0, 1.0, t_gyr=3.0)
    c = po.class1_current_orbit("McMillan17", t_gyr=3.0, n=6001)
    assert abs(f.r_peri - c.r_peri) < 0.1 and abs(f.r_apo - c.r_apo) < 0.1
    assert not f.notes["friction"]


def test_backward_friction_pumps_the_orbit_outward_with_mass():
    apos = [po.class2_friction_backwards_fast(m, 1.0, t_gyr=3.0).r_apo for m in (0.0, 1e9, 1e10)]
    assert apos[0] < apos[1] < apos[2]


def test_exponential_stripping_floors_at_the_nucleus():
    m = po.exponential_stripping(1e10, 1.0, 10.0)
    assert m(11.0) == po.M_NUCLEUS_MSUN and m(10.0) == 1e10
    assert m(0.0) == po.M_NUCLEUS_MSUN and 1e9 < m(8.0) < 1e10


def test_class3_orbit_starts_at_apocentre_with_requested_lz():
    o = po.class3_gse_debris_orbit(15.5, -300.0, 60.0, t_gyr=1.0)
    assert abs(o.r_apo - 15.5) < 0.2 and o.eccentricity > 0.85   # later apocentres drift slightly (L not conserved)
    assert abs(o.R[0] * o.vT[0] + 300.0) < 1.0


def test_class3_bar_migration_orbits_are_gse_like_before_migration():
    bm = pytest.importorskip("ocen_dm.tails.bar_migration")
    if not bm.POT_DIR.exists():
        pytest.skip("oCen_bar clone not found")
    orbs = po.class3_bar_migration_orbits(n_samples=200, n_times=401)
    assert len(orbs) == 3
    for o in orbs:
        assert abs(o.t_gyr[0]) < 1e-9 and abs(o.t_gyr[-1] + 8.0) < 1e-9      # look-back axis, 0 today first
        assert abs(o.r[0] - 6.46) < 0.2                                      # today's radius
        early = o.notes["early_0_1_gyr"]
        assert early["ecc"] > 0.8 and 9 < early["r_apo"] < 14                # GSE-debris-like orbit
        assert o.notes["Lz0_kpc_kms"] > -100                                # no longer strongly retrograde
