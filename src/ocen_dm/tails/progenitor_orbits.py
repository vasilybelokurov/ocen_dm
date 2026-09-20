"""Candidate orbit histories for omega Cen's host dwarf, in three classes.

The plan (JOURNAL 2026-09-20) is to constrain the cluster's dark matter from the other end:
simulate the tidal disruption of the host and see what the nucleus is left with. That needs
the host's orbit over the last ~10 Gyr, which is not observable. Three classes bracket it.

* **Class 1 -- today's orbit, static axisymmetric potential.** The dwarf is assumed to have
  always been on the present orbit. This is the destructive limit, not a realistic history:
  a pericentre of 1.6 kpc from the start strips a dwarf to its nucleus in a few passages.
  Its value is as a bound. Variants: different Milky Way potentials.
* **Class 2 -- today's orbit integrated backwards with dynamical friction.** The dwarf sank
  from a wider orbit. Backward integration with friction is well posed *given* the satellite
  mass history; the mass history is not known, so it is scanned. Variants: satellite mass.
* **Class 3 -- deposited by GSE, then transported inward by the bar.** Dillamore, Zhang &
  Belokurov (2026, arXiv:2606.12516) show this works only for a present-day bar pattern speed
  below ~26 km/s/kpc. The early orbit is then GSE's own infall orbit (Naidu et al. 2021,
  arXiv:2103.03251): from the virial radius at z ~ 2, circularity 0.5, inclination 15 deg,
  retrograde, with M_* = 5e8 and M_DM = 2e11 Msun. This class is built from the literature
  rather than integrated here, because it needs a time-dependent barred potential.

Units: kpc, km/s, Msun, Gyr. Present-day phase space from :mod:`ocen_dm.cluster`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..cluster import MU_SYS, OCEN_DEC, OCEN_RA, OCEN_VSYS_KMS


def _import_galpy_safely() -> None:
    """galpy imports astroquery for Simbad name lookups, which we never use, and in this
    environment astroquery crashes on import because astropy's package metadata reports no
    version (``importlib.metadata.version('astropy')`` is None while ``astropy.__version__``
    is 8.0.1). Blocking the optional dependency lets galpy fall through cleanly. Done here,
    once, rather than by editing the environment."""
    import sys
    try:
        import astroquery.simbad  # noqa: F401  -- the sub-module is what fails
    except Exception:
        sys.modules["astroquery"] = None
        sys.modules["astroquery.simbad"] = None
    # galpy <= 1.9.2 (checkout of 2024-03) imports scipy's private ``vectorize1`` in
    # galpy/util/quadpack.py and then shadows it with its own definition on the next lines;
    # scipy >= 1.15 removed the name. Provide a placeholder so the import succeeds
    # (upstream fix: galpy commit 4fc5db9e).
    import scipy.integrate._quadrature as _q
    if not hasattr(_q, "vectorize1"):
        _q.vectorize1 = None


_import_galpy_safely()

__all__ = ["PRESENT", "OrbitSummary", "present_day_orbit", "class1_current_orbit",
           "class2_friction_backwards", "class2_friction_backwards_fast",
           "class3_gse_initial_conditions", "class3_gse_debris_orbit", "plausible_orbit_set",
           "exponential_stripping", "CLASS2_SET", "CLASS3_SET", "M_NUCLEUS_MSUN",
           "class3_bar_migration_orbits", "CLASS3_OMEGA_B"]

#: present-day observables (D and v_los as used throughout the project)
PRESENT = dict(ra=OCEN_RA, dec=OCEN_DEC, distance_kpc=5.43,
               pmra=MU_SYS[0], pmdec=MU_SYS[1], vlos=OCEN_VSYS_KMS)
#: solar parameters used for every conversion (McMillan 2017 / galpy defaults)
R0_KPC, V0_KMS, ZSUN_KPC = 8.21, 233.1, 0.0208
#: solar peculiar motion in galpy's ``[-U, V, W]`` convention (Schoenrich et al. 2010);
#: galpy adds v_circ itself, so V0 must NOT be included here (checked against astropy:
#: identical galactocentric x, v to 0.1 km/s, L_z = -529 kpc km/s, retrograde).
VSUN = [-11.1, 12.24, 7.25]


@dataclass
class OrbitSummary:
    label: str
    potential: str
    t_gyr: np.ndarray
    R: np.ndarray            # cylindrical radius, kpc
    z: np.ndarray            # height, kpc
    r: np.ndarray            # spherical radius, kpc
    vR: np.ndarray
    vT: np.ndarray
    vz: np.ndarray
    notes: dict = field(default_factory=dict)

    @property
    def r_peri(self) -> float:
        return float(self.r.min())

    @property
    def r_apo(self) -> float:
        return float(self.r.max())

    @property
    def eccentricity(self) -> float:
        return (self.r_apo - self.r_peri) / (self.r_apo + self.r_peri)

    def at(self, t_lookback_gyr: float) -> dict:
        """State at a given look-back time (positive = past)."""
        i = int(np.argmin(np.abs(self.t_gyr + t_lookback_gyr)))
        return dict(t=float(self.t_gyr[i]), R=float(self.R[i]), z=float(self.z[i]),
                    r=float(self.r[i]), vR=float(self.vR[i]), vT=float(self.vT[i]),
                    vz=float(self.vz[i]))


def _potential(name: str):
    from galpy.potential import MWPotential2014, mwpotentials
    if name == "McMillan17":
        return mwpotentials.McMillan17
    if name == "MWPotential2014":
        return MWPotential2014
    if name == "Irrgang13I":
        return mwpotentials.Irrgang13I
    raise ValueError(name)


def present_day_orbit(ro: float = R0_KPC, vo: float = V0_KMS):
    """galpy Orbit at the present day from the project's observables.

    ``ro, vo`` must match the potential the orbit is integrated in; galpy's bundled
    potentials carry their own (Irrgang13I: 8.4 kpc, 242 km/s).
    """
    from galpy.orbit import Orbit
    return Orbit([PRESENT["ra"], PRESENT["dec"], PRESENT["distance_kpc"],
                  PRESENT["pmra"], PRESENT["pmdec"], PRESENT["vlos"]],
                 radec=True, ro=ro, vo=vo, zo=ZSUN_KPC, solarmotion=VSUN)


def _orbit_in(pot):
    """Present-day orbit with the physical scaling of ``pot``."""
    from galpy.util.conversion import get_physical
    ph = get_physical(pot)
    return present_day_orbit(ro=ph.get("ro", R0_KPC), vo=ph.get("vo", V0_KMS))


def _summarise(o, ts, label, potname, notes=None) -> OrbitSummary:
    return OrbitSummary(label, potname, np.asarray(ts.to_value("Gyr")) if hasattr(ts, "to_value") else ts,
                        o.R(ts), o.z(ts), o.r(ts), o.vR(ts), o.vT(ts), o.vz(ts), notes or {})


def class1_current_orbit(potential: str = "McMillan17", t_gyr: float = 10.0,
                         n: int = 20001) -> OrbitSummary:
    """The present orbit integrated for ``t_gyr`` in a static potential, no friction."""
    import astropy.units as u
    pot = _potential(potential)
    o = _orbit_in(pot)
    ts = np.linspace(0.0, -t_gyr, n) * u.Gyr
    o.integrate(ts, pot, method="dop853_c")
    return _summarise(o, ts, "class1 %s" % potential, potential,
                      {"friction": False, "t_gyr": t_gyr})


def class2_friction_backwards(m_sat_msun: float, r_half_kpc: float,
                              potential: str = "McMillan17", t_gyr: float = 10.0,
                              n: int = 20001, mass_history: str = "constant") -> OrbitSummary:
    """The present orbit integrated **backwards** with Chandrasekhar dynamical friction.

    Time-reversed friction pumps the orbit outward, recovering where a satellite of mass
    ``m_sat_msun`` would have had to be to sink to the present orbit. ``mass_history``
    ``"constant"`` keeps the mass fixed (an upper bound on the friction late in the orbit,
    when the real satellite was lighter). The satellite's half-mass radius sets the
    softening of the friction force.
    """
    import astropy.units as u
    from galpy.potential import ChandrasekharDynamicalFrictionForce
    pot = _potential(potential)
    from galpy.util.conversion import get_physical
    ph = get_physical(pot)
    cdf = ChandrasekharDynamicalFrictionForce(GMs=m_sat_msun * u.Msun, rhm=r_half_kpc * u.kpc,
                                              dens=pot, ro=ph.get("ro", R0_KPC), vo=ph.get("vo", V0_KMS))
    o = _orbit_in(pot)
    ts = np.linspace(0.0, -t_gyr, n) * u.Gyr
    o.integrate(ts, pot + cdf, method="dop853_c")
    return _summarise(o, ts, "class2 M=%.1e" % m_sat_msun, potential,
                      {"friction": True, "m_sat_msun": m_sat_msun, "r_half_kpc": r_half_kpc,
                       "mass_history": mass_history, "t_gyr": t_gyr})


def class3_gse_initial_conditions() -> dict:
    """GSE's infall orbit and structure at z ~ 2, from Naidu et al. (2021), Tables 1-3 and
    Section 4.2, with the bar-migration requirement from Dillamore et al. (2026)."""
    return {
        "source": ["Naidu et al. 2021, arXiv:2103.03251 (fiducial model M1)",
                   "Dillamore, Zhang & Belokurov 2026, arXiv:2606.12516"],
        "lookback_gyr": 10.0, "redshift": 2.0,
        "start": "at the Milky Way's virial radius, E = E_circ(R_vir)",
        "circularity": 0.5, "inclination_deg": 15.0, "sense": "retrograde",
        "M_star_msun": 5e8, "M_dm_msun": 2e11, "c200": 4.0,
        "disk_r50_kpc": 1.7, "dm_profile": "NFW mapped to Hernquist (GalIC)",
        "mw_at_z2": {"M200_msun": 5e11, "c200": 3.8, "disk_msun": 6e9, "bulge_msun": 1.4e10},
        "bar_requirement": "present-day pattern speed <~ 26 km/s/kpc (decelerating bar)",
        "caveat": "GSE membership of omega Cen is contested; oMEGACat X (arXiv:2603.23589) "
                  "argues for Sequoia/Thamnos as the earliest-stripped debris and leaves the "
                  "GSE connection unsure.",
    }


# ---------------------------------------------------------------------------------------
# Fast friction integrator. galpy's ChandrasekharDynamicalFrictionForce is a Python force,
# so galpy falls back to its Python ODE solver (tens of minutes per orbit). This leapfrog
# evaluates the host with AGAMA (McMillan 2017 as shipped in agama/data) and applies the same
# Chandrasekhar formula galpy uses:
#   a_df = -4 pi G^2 M rho(r) ln(Lambda) [erf(X) - 2X/sqrt(pi) exp(-X^2)] v / v^3,
#   X = v / (sqrt(2) sigma(r)),  sigma = v_circ(r) / sqrt(2)  (isothermal-sphere estimate),
#   ln Lambda = ln( r / max(r_hm, G M / v^2) )   (galpy's default, gamma = 1).
# Integrating backwards flips the sign of the friction term: the orbit is pumped outward.
# Units: kpc, km/s, Msun. The leapfrog runs in AGAMA's natural time unit kpc/(km/s) =
# 0.977792 Gyr, in which x' = v and v' = a with a in (km/s)^2/kpc; times are converted to Gyr
# only for output. (An earlier version divided the acceleration by 1.02271 instead of
# multiplying: kicks 4.4% too weak, pericentre 4% too large; caught by Codex on 2026-09-20.)
# ---------------------------------------------------------------------------------------
_G_KPC_KMS2_MSUN = 4.30091e-6  # kpc (km/s)^2 / Msun
TIME_UNIT_GYR = 0.977792       # 1 kpc/(km/s) in Gyr (AGAMA natural unit for kpc, km/s)
_KPC_PER_KMS_GYR = 1.0 / TIME_UNIT_GYR   # = 1.02271, kept for external callers


def _agama_potential(name: str = "McMillan17"):
    import agama, os
    agama.setUnits(mass=1, length=1, velocity=1)   # Msun, kpc, km/s
    if name != "McMillan17":
        raise ValueError("AGAMA host available for McMillan17 only")
    return agama.Potential(os.path.join(os.path.dirname(agama.__file__), "data", "McMillan17.ini"))


def _friction_accel(pot, x, v, m_sat, r_hm):
    """Chandrasekhar deceleration [(km/s)^2/kpc = km/s per time unit] in the *forward* sense of time."""
    from scipy.special import erf
    r = np.linalg.norm(x); vmag = np.linalg.norm(v)
    rho = float(pot.density(x))                                   # Msun/kpc^3
    R = np.hypot(x[0], x[1])
    vc = np.sqrt(-R * float(pot.force([R, 0.0, 0.0])[0]))        # km/s
    sigma = vc / np.sqrt(2.0)
    X = vmag / (np.sqrt(2.0) * sigma)
    bmin = max(r_hm, _G_KPC_KMS2_MSUN * m_sat / vmag**2)
    lnL = max(np.log(r / bmin), 0.0)
    coeff = 4.0 * np.pi * _G_KPC_KMS2_MSUN**2 * m_sat * rho * lnL * (
        erf(X) - 2.0 * X / np.sqrt(np.pi) * np.exp(-X**2))       # (km/s)^2 / kpc
    return -coeff * v / vmag**3


def class2_friction_backwards_fast(m_sat_msun: float, r_half_kpc: float,
                                   potential: str = "McMillan17", t_gyr: float = 10.0,
                                   dt_myr: float = 0.25, mass_history=None) -> OrbitSummary:
    """Leapfrog version of :func:`class2_friction_backwards` (seconds instead of tens of minutes).

    ``mass_history(t_lookback_gyr) -> M_sat`` may replace the constant mass, e.g. to mimic
    tidal mass loss (a lighter satellite in the past experiences less friction).
    """
    pot = _agama_potential(potential)
    o = present_day_orbit()
    # galpy is left-handed (x towards the Sun); AGAMA/astropy right-handed: flip x and vx
    x = np.array([-o.x(), o.y(), o.z()], float)
    v = np.array([-o.vx(), o.vy(), o.vz()], float)
    dt = -dt_myr * 1e-3 / TIME_UNIT_GYR       # negative: backwards; natural time units
    n = int(round(t_gyr / (dt_myr * 1e-3)))
    m_of = (lambda t: m_sat_msun) if mass_history is None else mass_history

    def accel(x, v, t_nat):
        a = np.asarray(pot.force(x), float)                        # (km/s)^2/kpc = km/s per time unit
        if m_sat_msun > 0:
            a = a + _friction_accel(pot, x, v, m_of(-t_nat * TIME_UNIT_GYR), r_half_kpc)
        return a

    ts = np.empty(n + 1); X = np.empty((n + 1, 3)); V = np.empty((n + 1, 3))
    t = 0.0; ts[0] = t; X[0] = x; V[0] = v
    a = accel(x, v, t)
    for i in range(1, n + 1):
        v_half = v + 0.5 * dt * a
        x = x + dt * v_half
        t += dt
        a = accel(x, v_half, t)                # velocity-dependent force: use half-step v
        v = v_half + 0.5 * dt * a
        ts[i] = t; X[i] = x; V[i] = v
    R = np.hypot(X[:, 0], X[:, 1]); phi = np.arctan2(X[:, 1], X[:, 0])
    vR = V[:, 0] * np.cos(phi) + V[:, 1] * np.sin(phi)
    # right-handed (astropy) frame has disc rotation with v_phi < 0; store prograde-positive like galpy
    vT = -(-V[:, 0] * np.sin(phi) + V[:, 1] * np.cos(phi))
    return OrbitSummary("class2fast M=%.1e" % m_sat_msun, potential, ts * TIME_UNIT_GYR, R, X[:, 2],
                        np.linalg.norm(X, axis=1), vR, vT, V[:, 2],
                        {"friction": m_sat_msun > 0, "m_sat_msun": m_sat_msun,
                         "r_half_kpc": r_half_kpc,
                         "mass_history": "constant" if mass_history is None else "custom",
                         "t_gyr": t_gyr, "integrator": "leapfrog dt=%.2f Myr, AGAMA host, natural units" % dt_myr})


# ---------------------------------------------------------------------------------------
# The plausible set: 2-3 orbits per class (decided 2026-09-20, see docs/PROGENITOR_ORBITS.md)
# ---------------------------------------------------------------------------------------
M_NUCLEUS_MSUN = 3.55e6      # Baumgardt & Hilker 2018; friction is negligible at this mass
#: class 2: (label, M_infall [Msun], stripping e-folding time [Gyr], t_infall [Gyr ago]);
#: chosen so that the orbit is at 50-150 kpc (the virial radius of the young Milky Way) at t_infall
CLASS2_SET = (("Sequoia-like  1e10, tau 2.0", 1e10, 2.0, 10.0),
              ("intermediate  3e10, tau 1.25", 3e10, 1.25, 10.0),
              ("GSE-like      1e11, tau 0.75", 1e11, 0.75, 10.0))
#: class 3 (superseded static guess, kept for tests): (label, apocentre, L_z, inclination). The
#: proper class-3 orbits come from :mod:`bar_migration` (Dillamore+2026 set-up); see
#: :func:`plausible_orbit_set`.
CLASS3_SET = (("GSE debris apo 11.5, Lz -300", 11.5, -300.0, 60.0),
              ("GSE debris apo 15.5, Lz -300", 15.5, -300.0, 60.0),
              ("GSE debris apo 21,   Lz -300", 21.0, -300.0, 60.0))
#: present-day bar pattern speed for the class-3 back-integration (paper's fiducial; migration
#: into the GSE debris works for Omega_b,0 <~ 26 km/s/kpc)
CLASS3_OMEGA_B = 24.0


def exponential_stripping(m_inf: float, tau_gyr: float, t_inf_gyr: float):
    """Bound mass vs look-back time: M_inf exp(-(t_inf - t_lb)/tau), floored at the nucleus."""
    def m_of(t_lb):
        if t_lb > t_inf_gyr:
            return M_NUCLEUS_MSUN
        return max(M_NUCLEUS_MSUN, m_inf * np.exp(-(t_inf_gyr - t_lb) / tau_gyr))
    return m_of


def class3_gse_debris_orbit(r_apo_kpc: float, lz_kpc_kms: float, incl_deg: float = 60.0,
                            t_gyr: float = 4.0, dt_myr: float = 0.25) -> OrbitSummary:
    """A static-potential (McMillan 2017, AGAMA) orbit started at apocentre on the x-axis with
    the requested L_z and an out-of-plane velocity giving the requested orbital inclination.
    Represents the debris orbit omega Cen would have had **before** bar migration; the bar
    itself is not modelled here (that needs the Hunter et al. 2024 barred potential)."""
    pot = _agama_potential("McMillan17")
    vphi = lz_kpc_kms / r_apo_kpc
    vz = abs(vphi) * np.tan(np.radians(incl_deg))
    x = np.array([r_apo_kpc, 0.0, 0.0]); v = np.array([0.0, -vphi, vz])   # right-handed: disk L_z < 0
    dt = dt_myr * 1e-3 / TIME_UNIT_GYR; n = int(round(t_gyr / (dt_myr * 1e-3)))   # natural units
    ts = np.empty(n + 1); X = np.empty((n + 1, 3)); V = np.empty((n + 1, 3))
    ts[0] = 0; X[0] = x; V[0] = v
    a = np.asarray(pot.force(x))
    for i in range(1, n + 1):
        v_half = v + 0.5 * dt * a
        x = x + dt * v_half
        a = np.asarray(pot.force(x))
        v = v_half + 0.5 * dt * a
        ts[i] = i * dt * TIME_UNIT_GYR; X[i] = x; V[i] = v
    R = np.hypot(X[:, 0], X[:, 1]); phi = np.arctan2(X[:, 1], X[:, 0])
    E = 0.5 * (V[0]**2).sum() + float(pot.potential(X[0]))
    return OrbitSummary("class3 apo=%.1f Lz=%.0f" % (r_apo_kpc, lz_kpc_kms), "McMillan17", ts, R, X[:, 2],
                        np.linalg.norm(X, axis=1),
                        V[:, 0] * np.cos(phi) + V[:, 1] * np.sin(phi),
                        -(-V[:, 0] * np.sin(phi) + V[:, 1] * np.cos(phi)),   # prograde positive
                        V[:, 2], {"E_km2s2": E, "Lz_kpc_kms": lz_kpc_kms, "incl_deg": incl_deg,
                                  "bar_migration": "not modelled; Dillamore+2026 need Omega_b <~ 26"})


def class3_bar_migration_orbits(omega_b: float = CLASS3_OMEGA_B, n_samples: int = 1000,
                                seed: int = 42, n_times: int = 161) -> list:
    """Class 3 done properly: omega Cen's phase-space samples integrated back 8 Gyr through the
    growing, decelerating bar of Dillamore et al. (2026); the three samples at the 16/50/84th
    percentiles of E(t=0) among those ending inside the GSE debris contours. Time axis of the
    returned OrbitSummary is look-back time (negative = past), like classes 1 and 2."""
    from . import bar_migration as bm
    run = bm.back_integrate(omega_b, n_samples=n_samples, seed=seed, n_times=n_times)
    out = []
    for p in bm.pick_class3(run):                      # fine (0.25 Myr) re-integration of each pick
        i = p["sample"]; X = p["traj"][:, :3]; V = p["traj"][:, 3:]
        lb = p["t"] - run.hist.tf                      # 0 today, -8 at bar formation
        R = np.hypot(X[:, 0], X[:, 1]); phi = np.arctan2(X[:, 1], X[:, 0])
        vR = V[:, 0] * np.cos(phi) + V[:, 1] * np.sin(phi)
        vT = -V[:, 0] * np.sin(phi) + V[:, 1] * np.cos(phi)     # oCen_bar frame: disc L_z > 0
        o = OrbitSummary("class3 bar-migrated, E0 q%.2f (sample %d)" % (p["q"], i), "Hunter24+bar",
                         lb[::-1], R[::-1], X[::-1, 2], np.linalg.norm(X, axis=1)[::-1],
                         vR[::-1], vT[::-1], V[::-1, 2],
                         {"omega_b_present": omega_b, "sample": i, "E0_km2s2": p["E0"], "Lz0_kpc_kms": p["Lz0"],
                          "early_0_1_gyr": p["early"], "late_last_gyr": p["late"],
                          "frac_samples_inside_gse_t0": run.frac_inside,
                          "set_up": "Dillamore, Zhang & Belokurov 2026; Hunter et al. 2024 potential"})
        out.append(o)
    return out


def plausible_orbit_set(potentials=("McMillan17", "MWPotential2014", "Irrgang13I"),
                        class3: str = "bar") -> dict:
    """The 3 + 3 + 3 orbits. Keys ``class1``, ``class2``, ``class3`` -> list of OrbitSummary.
    ``class3="bar"`` (default) uses the Dillamore+2026 back-integration; ``"static"`` the
    superseded static GSE-debris guesses."""
    out = {"class1": [class1_current_orbit(p) for p in potentials], "class2": [], "class3": []}
    for label, m_inf, tau, t_inf in CLASS2_SET:
        rh = 1.0 * (m_inf / 1e10) ** (1.0 / 3.0)
        o = class2_friction_backwards_fast(m_inf, rh, t_gyr=t_inf,
                                           mass_history=exponential_stripping(m_inf, tau, t_inf))
        o.label = label; o.notes.update(m_infall=m_inf, tau_gyr=tau, t_infall_gyr=t_inf,
                                        r_at_infall_kpc=o.at(t_inf)["r"],
                                        apo_last_gyr_kpc=float(o.r[o.t_gyr < -(t_inf - 1)].max()))
        out["class2"].append(o)
    if class3 == "bar":
        out["class3"] = class3_bar_migration_orbits()
    else:
        for label, apo, lz, incl in CLASS3_SET:
            o = class3_gse_debris_orbit(apo, lz, incl); o.label = label
            out["class3"].append(o)
    return out
