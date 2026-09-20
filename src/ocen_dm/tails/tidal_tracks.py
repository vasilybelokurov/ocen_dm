"""Tidal tracks (Penarrubia+2010; Errani & Navarro 2021, EN21) walked along our nine orbits.

Model (EN21, arXiv:2011.07077):
  * track (eq. 5):  V_mx/V_mx0 = 2^a (r_mx/r_mx0)^b [1 + (r_mx/r_mx0)^2]^-a,  a = 0.4, b = 0.65;
    with T_mx = 2 pi r_mx / V_mx = (3 pi / G rho_mx)^1/2 this maps the crossing time to M_mx/M_mx0.
  * time evolution in the heavy-mass-loss regime (T_mx0/T_peri > 2/3), eqs 10-14:
      Y(t) = (T_mx - T_asy)/T_peri = Y0 [1 + (t/tau)^eta]^(-1/eta),  T_asy = 0.22 T_peri,
      tau = 0.65 T_orb / Y0,  eta = 1 - exp(-2.5 Y0);  eccentric orbits are delayed by
      f_ecc = [2x/(x+1)]^3.2, x = r_apo/r_peri (eq. 4).
    modest regime (T_mx0/T_peri < 2/3), eqs 15-16: T_asy' = T_mx0 / (1 + T_mx0/T_peri)^2.2,
      eta' = 0.67, tau' = 1.2 (T_mx0/T_peri)^-1/2 T_orb.
  * remnant profile (eqs 7-9): rho = rho_NFW exp(-r/r_cut) / (1 + r_s/r_cut)^0.3,
      r_cut/r_mx0 = 0.44 (M_mx/M_mx0)^0.44 [1 - (M_mx/M_mx0)^0.3]^-1.1.
  T_peri is the circular time at pericentre, here generalised to sqrt(3 pi / G rho_host(<r_peri))
  (EN21 state the remnant density is set by the host mean density at pericentre and that the
  model applies to other hosts "after proper scaling", with caution).
Validity: EN21 simulated 0.2 < T_mx0/T_peri < 2 and r_apo/r_peri <= 20. Our orbits are far outside
the upper end (see ``walk``), so the results are extrapolations and are flagged as such.

Orbit dependence enters only through the per-passage (T_peri, T_orb, f_ecc): the state (Y, hence
T_mx, M_mx) is advanced passage by passage with the local parameters, using the "equivalent time"
of the current state under the local law (Markov approximation).
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

G = 4.30091e-6            # kpc (km/s)^2 / Msun
TU = 0.977792             # kpc/(km/s) in Gyr
A_TR, B_TR = 0.4, 0.65
RHO_CRIT_Z0 = 277.5 * 0.7**2 / 1e9 * 1e9   # Msun/kpc^3 for h = 0.7: 2.775e11 h^2 Msun/Mpc^3 = 136 Msun/kpc^3
RHO_CRIT_Z0 = 2.775e11 * 0.7**2 / 1e9      # = 136.0 Msun/kpc^3


def rho_crit(z, om=0.3):
    return RHO_CRIT_Z0 * (om * (1 + z)**3 + 1 - om)


@dataclass
class NFW:
    m200: float
    c: float
    z: float = 2.0

    def __post_init__(self):
        self.r200 = (3 * self.m200 / (800 * np.pi * rho_crit(self.z))) ** (1 / 3)
        self.rs = self.r200 / self.c
        f = np.log(1 + self.c) - self.c / (1 + self.c)
        self.rho_s = self.m200 / (4 * np.pi * self.rs**3 * f)
        self.rmx0 = 2.163 * self.rs
        self.Mmx0 = self.mass(self.rmx0)
        self.Vmx0 = np.sqrt(G * self.Mmx0 / self.rmx0)
        self.Tmx0 = 2 * np.pi * self.rmx0 / self.Vmx0      # kpc/(km/s)

    def mass(self, r):
        x = np.asarray(r) / self.rs
        return 4 * np.pi * self.rho_s * self.rs**3 * (np.log(1 + x) - x / (1 + x))

    def rho(self, r):
        x = np.asarray(r) / self.rs
        return self.rho_s / (x * (1 + x)**2)


# ---------------------------------------------------------------- track: T_mx <-> M_mx
def track_v(x):           # x = r_mx / r_mx0
    return 2**A_TR * x**B_TR * (1 + x**2) ** (-A_TR)


def track_T(x):           # T_mx / T_mx0
    return x / track_v(x)


def track_M(x):           # M_mx / M_mx0 = (r/r0) (V/V0)^2
    return x * track_v(x) ** 2


def x_from_T(T_ratio):
    """invert T_mx/T_mx0 -> r_mx/r_mx0 (monotonic for x < 1)."""
    if T_ratio >= 1.0:
        return 1.0
    return brentq(lambda x: track_T(x) - T_ratio, 1e-12, 1.0)


# ---------------------------------------------------------------- remnant profile (eqs 7-9)
def r_cut(halo: NFW, mfrac):
    mfrac = min(max(mfrac, 1e-12), 1 - 1e-12)
    return halo.rmx0 * 0.44 * mfrac**0.44 * (1 - mfrac**0.3) ** (-1.1)


def remnant_rho(halo: NFW, mfrac, r):
    rc = r_cut(halo, mfrac)
    return halo.rho(r) * np.exp(-np.asarray(r) / rc) / (1 + halo.rs / rc) ** 0.3


def remnant_mass_within(halo: NFW, mfrac, r_out):
    from scipy.integrate import quad
    return quad(lambda r: 4 * np.pi * r**2 * remnant_rho(halo, mfrac, r), 0, r_out, limit=200)[0]


# ---------------------------------------------------------------- EN21 time evolution
def _law(Tmx0, Tperi, Torb):
    """Return (T_asy, tau, eta) for the local pericentre; Tmx0 is the *initial* crossing time."""
    ratio = Tmx0 / Tperi
    if ratio > 2 / 3:
        Tasy = 0.22 * Tperi
        Y0 = (Tmx0 - Tasy) / Tperi
        return Tasy, 0.65 * Torb / Y0, 1 - np.exp(-2.5 * Y0)
    Tasy = Tmx0 / (1 + ratio) ** 2.2
    return Tasy, 1.2 * ratio ** (-0.5) * Torb, 0.67


def _Y(t, Y0, tau, eta):
    return Y0 * (1 + (t / tau) ** eta) ** (-1 / eta)


def _t_of_Y(Y, Y0, tau, eta):
    if Y >= Y0:
        return 0.0
    return tau * ((Y0 / Y) ** eta - 1) ** (1 / eta)


def walk(halo: NFW, t_peri_gyr, r_peri, r_apo, Tperi_gyr, Torb_gyr):
    """Advance the remnant passage by passage. Inputs per passage: look-back time of the passage,
    pericentre, the following apocentre (for f_ecc), local T_peri and the radial period, all in
    Gyr/kpc. Returns dict with arrays of T_mx, M_mx/M_mx0, r_mx over the passages."""
    Tmx0 = halo.Tmx0 * TU                                    # Gyr
    Tmx = Tmx0
    out = dict(t=[], Tmx=[], mfrac=[], x=[], ratio0=[])
    for tp, rp, ra, Tp, To in zip(t_peri_gyr, r_peri, r_apo, Tperi_gyr, Torb_gyr):
        Tasy, tau, eta = _law(Tmx0, Tp, To)
        Y0 = (Tmx0 - Tasy) / Tp
        Ynow = max((Tmx - Tasy) / Tp, 1e-9)
        if Ynow > Y0:                    # gentler pericentre than any before: state cannot exceed Y0
            Ynow = Y0
        fecc = (2 * (ra / rp) / (ra / rp + 1)) ** 3.2
        t_eq = _t_of_Y(Ynow, Y0, tau, eta) + To / fecc
        Ynew = _Y(t_eq, Y0, tau, eta)
        Tmx = min(Tmx, Tasy + Ynew * Tp)
        x = x_from_T(Tmx / Tmx0)
        out["t"].append(tp); out["Tmx"].append(Tmx); out["mfrac"].append(track_M(x)); out["x"].append(x); out["ratio0"].append(Tmx0 / Tp)
    return {k: np.array(v) for k, v in out.items()}


def passages(o, pot):
    """Pericentre passages of an OrbitSummary (present first, t_gyr <= 0), chronological:
    look-back time, r_peri, the apocentre following each passage, T_peri from the host mean
    density inside r_peri, and the radial period."""
    r = o.r; t = -o.t_gyr
    ip = np.flatnonzero((r[1:-1] < r[:-2]) & (r[1:-1] < r[2:])) + 1        # indices, present first
    ia = np.flatnonzero((r[1:-1] > r[:-2]) & (r[1:-1] > r[2:])) + 1
    ip = ip[::-1]                                                          # chronological
    rp, tp = r[ip], t[ip]
    ra = np.array([r[ia[ia < i]].max() if np.any(ia < i) else r[ia].max() for i in ip])  # next apocentre (smaller index = later)
    Menc = np.array([pot.enclosedMass(x) for x in rp])
    rho = Menc / (4 / 3 * np.pi * rp**3)
    Tperi = np.sqrt(3 * np.pi / (G * rho)) * TU
    Torb = np.abs(np.gradient(tp)) if len(tp) > 1 else np.array([0.1])
    return tp, rp, ra, Tperi, Torb


# ---------------------------------------------------------------- the robust part: densities
def asymptotic_density_radius(halo: NFW, rho_host_peri, factor=20.0):
    """Radius inside which the *initial* NFW mean density exceeds ``factor`` x the host mean density
    inside the pericentre (EN21: T_asy = 0.22 T_peri <=> rho_asy ~ 1/0.22^2 ~ 20 rho_host), and the
    NFW mass inside it. This is the scale of what can survive as a self-bound DM remnant."""
    target = factor * rho_host_peri
    f = lambda r: halo.mass(r) / (4 / 3 * np.pi * r**3) - target
    if f(1e-6) < 0:
        return 0.0, 0.0
    r = brentq(f, 1e-6, halo.r200)
    return r, halo.mass(r)


def jacobi_radius(pot, r_gal, m):
    """Jacobi radius r (m / (3 M_host(<r)))^(1/3), the simple form used throughout the project."""
    return r_gal * (m / (3 * pot.enclosedMass(r_gal))) ** (1 / 3)


def budget_table(orbits_by_class, pots, haloes, m_nuc=3.55e6):
    """For each orbit and initial halo: host mean density at the median and minimum pericentre,
    the self-bound-remnant radius/mass (asymptotic-density criterion), and the NFW mass inside the
    nucleus's Jacobi radius at the minimum pericentre (what the nucleus's own potential can hold)."""
    from astropy.table import Table
    rows = []
    for cls, orbs in orbits_by_class.items():
        for o in orbs:
            pot = pots[o.potential]
            tp, rp, ra, Tp, To = passages(o, pot)
            rmin, rmed = rp.min(), np.median(rp)
            rho_min = pot.enclosedMass(rmin) / (4 / 3 * np.pi * rmin**3)
            rj = jacobi_radius(pot, rmin, m_nuc)
            for h in haloes:
                r_asy, m_asy = asymptotic_density_radius(h, rho_min)
                tmx_ratio = h.Tmx0 * TU / (np.sqrt(3 * np.pi / (G * rho_min)) * TU)
                rows.append(dict(cls=cls, orbit=o.label, potential=o.potential, n_peri=len(rp), r_peri_min=rmin, r_peri_med=rmed,
                                 rho_host_peri_min=rho_min, m200=h.m200, c=h.c, Tmx0_over_Tperi=tmx_ratio,
                                 r_selfbound_pc=1e3 * r_asy, m_selfbound=m_asy, rj_nucleus_pc=1e3 * rj,
                                 m_dm_in_rj_initial=h.mass(rj), m_dm_in_35pc=h.mass(0.035), m_dm_in_70pc=h.mass(0.07)))
    return Table(rows=rows)


def budget_figure(tab, haloes, path="plots/tidal_tracks_budget.png"):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    r = np.logspace(-3, 1, 300)
    for h in haloes:
        ax[0].plot(r * 1e3, h.mass(r), label="NFW M200 = %.0e, c = %.1f" % (h.m200, h.c))
    ax[0].axhline(3.55e6, color="k", ls=":", lw=0.8); ax[0].text(1.2, 4e6, "nucleus mass", fontsize=8)
    c3 = tab[tab["cls"] == "class3"]; c12 = tab[tab["cls"] != "class3"]
    ax[0].axvspan(c3["rj_nucleus_pc"].min(), c3["rj_nucleus_pc"].max(), color="C6", alpha=0.25, label="r_J(nucleus) at class-3 pericentre")
    ax[0].axvspan(c12["rj_nucleus_pc"].min(), c12["rj_nucleus_pc"].max(), color="C0", alpha=0.15, label="r_J(nucleus) at class-1/2 pericentre")
    ax[0].axvspan(tab["r_selfbound_pc"].min(), tab["r_selfbound_pc"].max(), color="grey", alpha=0.3, label="self-bound remnant (EN21 density criterion)")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlim(1, 1e4); ax[0].set_ylim(1e2, 1e11)
    ax[0].set_xlabel("r [pc]"); ax[0].set_ylabel(r"initial NFW mass inside r [M$_\odot$]"); ax[0].legend(fontsize=7, loc="lower right"); ax[0].grid(alpha=0.3)
    labels = ["%.0e, c%.0f" % (h.m200, h.c) for h in haloes]; xs = np.arange(len(haloes))
    for k, (cls, col) in enumerate((("class1", "C0"), ("class2", "C3"), ("class3", "C6"))):
        lo = [tab[(tab["cls"] == cls) & (tab["m200"] == h.m200) & (tab["c"] == h.c)]["m_dm_in_rj_initial"].min() for h in haloes]
        hi = [tab[(tab["cls"] == cls) & (tab["m200"] == h.m200) & (tab["c"] == h.c)]["m_dm_in_rj_initial"].max() for h in haloes]
        ax[1].bar(xs + (k - 1) * 0.25, hi, width=0.25, color=col, alpha=0.8, label=cls)
        ax[1].errorbar(xs + (k - 1) * 0.25, hi, yerr=[np.array(hi) - np.array(lo), np.zeros(len(hi))], fmt="none", ecolor="k", lw=0.8)
    ax[1].axhline(3.55e6, color="k", ls=":", lw=0.8); ax[1].set_yscale("log"); ax[1].set_xticks(xs); ax[1].set_xticklabels(labels, fontsize=8)
    ax[1].set_ylabel(r"initial DM inside the nucleus's Jacobi radius at pericentre [M$_\odot$]"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, axis="y")
    ax[1].set_title("upper envelope of DM the nucleus can hold (no contraction, no shocks)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
