"""Class 3 done properly: omega Cen back-integrated through a growing, decelerating bar.

Reproduces the set-up of Dillamore, Zhang & Belokurov (2026, arXiv:2606.12516), following the
authors' pipeline in the ``oCen_bar`` repository (github.com/adllmr/oCen_bar, local clone at
``~/Work/Code/oCen_bar``): Hunter et al. (2024) Milky Way potential split into an axisymmetric
part and a bar part; the bar amplitude grows with the Dehnen (2000) switch-on between t0 and t1,
the pattern speed decelerates smoothly between t1 and t2 and then with constant
eta = -dOmega/dt / Omega^2 = 0.003 until t_f = 8 Gyr; the bar length scales as
S = Omega_b,f / Omega_b(t) while the amplitude is held fixed. 10^3 samples of omega Cen's present
phase-space point are integrated backwards from t_f to 0 and compared with the GSE debris
contours of Belokurov et al. (2023) in the axisymmetric potential.

The Hunter et al. potential files and the GSE contour file are *referenced* from the oCen_bar
clone (``OCEN_BAR_DIR``), not copied.
"""
from __future__ import annotations
import os
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

from ..cluster import OCEN_RA, OCEN_DEC, OCEN_VSYS_KMS, MU_SYS
from .progenitor_orbits import PRESENT

OCEN_BAR_DIR = Path(os.environ.get("OCEN_BAR_DIR", "~/Work/Code/oCen_bar")).expanduser()
POT_DIR = OCEN_BAR_DIR / "agama_potentials"
GSE_CONTOURS = OCEN_BAR_DIR / "artifacts" / "gaiadr3_gse_elz.fits"

#: 1-sigma uncertainties on the present-day observables (Vasiliev & Baumgardt 2021 PMs,
#: Baumgardt & Vasiliev 2021 distance and v_los), as in oCen_bar/data/GC_catalogue.txt
PRESENT_ERR = dict(distance_kpc=0.05, pmra=0.025, pmdec=0.025, vlos=0.21)
BAR_ANGLE_DEG = 28.0          # present-day bar angle to the Sun-GC line (oCen_bar convention)
ETA = 0.003
T0, T1, T2, TF = 0.0, 1.0, 2.0, 8.0   # Gyr; paper text: growth to t1 ~ 1, smooth deceleration to t2 ~ 2


def _agama():
    import agama
    agama.setUnits(length=1, velocity=1, mass=1)
    return agama


# ----------------------------------------------------------------------------- bar history
def bar_fraction(t, t0=T0, t1=T1):
    """Dehnen (2000) eq. 4 smooth switch-on of the bar amplitude between t0 and t1."""
    xi = 2.0 * (np.asarray(t, float) - t0) / (t1 - t0) - 1.0
    return np.clip(3 / 16 * xi**5 - 5 / 8 * xi**3 + 15 / 16 * xi + 0.5, 0.0, 1.0)


def omega_b(t, omega_1, eta=ETA, t1=T1, t2=T2):
    """Pattern speed: constant omega_1 to t1, smoothly starting deceleration to t2, then
    constant eta (Omega = Omega_2 / (1 + eta Omega_2 (t - t2)))."""
    omega_2 = omega_1 / (1 + 0.5 * eta * omega_1 * (t2 - t1))
    t = np.asarray(t, float)
    out = np.where(t <= t1, omega_1,
                   np.where(t <= t2, omega_1 / (1 + 0.5 * eta * omega_1 * (t - t1) ** 2 / (t2 - t1)),
                            omega_2 / (1 + eta * omega_2 * (t - t2))))
    return out


def omega_1_for_final(omega_f, eta=ETA, t1=T1, t2=T2, tf=TF):
    """Initial pattern speed that ends at omega_f at tf (inverse of :func:`omega_b`)."""
    omega_2 = omega_f / (1 - eta * omega_f * (tf - t2))
    return omega_2 / (1 - 0.5 * eta * omega_2 * (t2 - t1))


@dataclass
class BarHistory:
    omega_f: float
    eta: float = ETA
    t0: float = T0
    t1: float = T1
    t2: float = T2
    tf: float = TF
    n: int = 1001
    t: np.ndarray = field(init=False)
    omega: np.ndarray = field(init=False)
    phi: np.ndarray = field(init=False)
    frac: np.ndarray = field(init=False)
    omega_1: float = field(init=False)

    def __post_init__(self):
        self.omega_1 = float(omega_1_for_final(self.omega_f, self.eta, self.t1, self.t2, self.tf))
        self.t = np.linspace(0.0, self.tf, self.n)
        self.omega = omega_b(self.t, self.omega_1, self.eta, self.t1, self.t2)
        # bar angle = integral of Omega, zero at present (t = tf)
        from scipy.integrate import cumulative_trapezoid
        phi = cumulative_trapezoid(self.omega, self.t, initial=0.0)
        self.phi = phi - phi[-1]
        self.frac = bar_fraction(self.t, self.t0, self.t1)


# ----------------------------------------------------------------------------- potentials
def axisymmetric_potential():
    agama = _agama()
    return agama.Potential(file=str(POT_DIR / "MWPotentialHunter24_axi.ini"))


def slowing_bar_potential(hist: BarHistory):
    """Axisymmetric Hunter24 + (scaled, rotating) baryonic bar part - (scaled) axisymmetrised
    baryonic part, exactly as ``oCen_bar.ocen_bar_common.make_potential``."""
    agama = _agama()
    S = hist.omega_f / hist.omega            # length scale ~ corotation radius
    A = hist.frac * S                        # amplitude switch-on times length scale
    pot_axi = axisymmetric_potential()
    pot_bar = agama.Potential(file=str(POT_DIR / "MWPotentialHunter24_baryon_full.ini"),
                              scale=np.column_stack([hist.t, A, S]),
                              rotation=np.column_stack([hist.t, hist.phi]))
    pot_bar_axi_neg = agama.Potential(file=str(POT_DIR / "MWPotentialHunter24_baryon_axi.ini"),
                                      scale=np.column_stack([hist.t, -A, S]))
    return agama.Potential(pot_axi, pot_bar, pot_bar_axi_neg), pot_axi


# ----------------------------------------------------------------------------- omega Cen samples
def present_day_samples(n: int = 1000, seed: int = 42, bar_angle_deg: float = BAR_ANGLE_DEG):
    """Galactocentric Cartesian samples in the bar frame (bar along x at t = tf), astropy default
    frame (R0 = 8.122 kpc, v_sun = (12.9, 245.6, 7.78), z_sun = 20.8 pc), as in the paper;
    left-handed x-flip as in oCen_bar so that the disc rotates with L_z > 0."""
    import astropy.units as u, astropy.coordinates as ac
    rng = np.random.default_rng(seed)
    d = rng.normal(PRESENT["distance_kpc"], PRESENT_ERR["distance_kpc"], n)
    pmra = rng.normal(PRESENT["pmra"], PRESENT_ERR["pmra"], n)
    pmdec = rng.normal(PRESENT["pmdec"], PRESENT_ERR["pmdec"], n)
    v = rng.normal(PRESENT["vlos"], PRESENT_ERR["vlos"], n)
    if n == 1:   # nominal
        d[:], pmra[:], pmdec[:], v[:] = PRESENT["distance_kpc"], PRESENT["pmra"], PRESENT["pmdec"], PRESENT["vlos"]
    c = ac.SkyCoord(ra=np.full(n, PRESENT["ra"]) * u.deg, dec=np.full(n, PRESENT["dec"]) * u.deg,
                    distance=d * u.kpc, pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr,
                    radial_velocity=v * u.km / u.s).transform_to(ac.Galactocentric())
    ic = np.column_stack([c.x.value, c.y.value, c.z.value, c.v_x.value, c.v_y.value, c.v_z.value])
    ic *= np.array([-1, 1, 1, -1, 1, 1])
    a = np.radians(bar_angle_deg)
    rot = np.array([[np.cos(a), np.sin(a), 0], [-np.sin(a), np.cos(a), 0], [0, 0, 1]])
    return np.hstack([ic[:, :3] @ rot.T, ic[:, 3:] @ rot.T])


def energy_lz(xv, pot_axi):
    xv = np.atleast_2d(xv)
    E = 0.5 * (xv[:, 3:] ** 2).sum(1) + pot_axi.potential(xv[:, :3])
    Lz = xv[:, 0] * xv[:, 4] - xv[:, 1] * xv[:, 3]
    return E, Lz


# ----------------------------------------------------------------------------- GSE contours
def gse_contours(pot_axi):
    """Belokurov et al. (2023) GSE contours in (L_z [kpc km/s], E [km^2/s^2]) shifted to the
    Hunter24 energy zero point (their potential: MWPotential2014_mod, matched at R = 8.2 kpc).
    Returns a list of (Lz, E) arrays, outermost first."""
    from astropy.io import fits
    agama = _agama()
    pot_b = agama.Potential(file=str(POT_DIR / "MWPotential2014_mod.ini"))
    offset = float(pot_axi.potential([[8.2, 0, 0]])[0] - pot_b.potential([[8.2, 0, 0]])[0])
    with fits.open(GSE_CONTOURS) as h:
        meta, pts = h[1].data, np.asarray(h[2].data, float)
    out = []
    for row in meta:
        seg = pts[row["OFFSET"]:row["OFFSET"] + row["N"]]
        out.append(np.column_stack([seg[:, 0], seg[:, 1] * 1e5 + offset]))
    return out


def inside_gse(E, Lz, contours, level: int = 0):
    from matplotlib.path import Path as MplPath
    poly = MplPath(np.column_stack([contours[level][:, 0], contours[level][:, 1] / 1e5]))
    return poly.contains_points(np.column_stack([Lz, np.asarray(E) / 1e5]))


# ----------------------------------------------------------------------------- integration
@dataclass
class MigrationRun:
    omega_f: float
    hist: BarHistory
    t: np.ndarray                 # time since bar formation, 0..tf (Gyr)
    traj: np.ndarray              # (n_times, n_samples, 6) in the bar frame at t = tf
    E: np.ndarray                 # (n_times, n_samples) in the axisymmetric potential
    Lz: np.ndarray
    inside_t0: np.ndarray         # bool per sample: inside outer GSE contour at t = 0
    ic: np.ndarray = None         # (n_samples, 6) present-day initial conditions (bar frame)

    @property
    def frac_inside(self) -> float:
        return float(self.inside_t0.mean())


def back_integrate(omega_f: float, n_samples: int = 1000, seed: int = 42, n_times: int = 161,
                   eta: float = ETA, t1: float = T1, t2: float = T2, tf: float = TF) -> MigrationRun:
    agama = _agama()
    hist = BarHistory(omega_f, eta=eta, t1=t1, t2=t2, tf=tf)
    pot, pot_axi = slowing_bar_potential(hist)
    ic = present_day_samples(n_samples, seed)
    orbits = agama.orbit(ic=ic, potential=pot, timestart=tf, time=-tf, trajsize=n_times)
    traj = np.stack([o for o in orbits[:, 1]]).transpose(1, 0, 2)   # (n_times, n, 6), t from tf down to 0
    t = np.asarray(orbits[0, 0])
    order = np.argsort(t); t, traj = t[order], traj[order]
    E = np.empty(traj.shape[:2]); Lz = np.empty_like(E)
    for j in range(traj.shape[0]):
        E[j], Lz[j] = energy_lz(traj[j], pot_axi)
    cont = gse_contours(pot_axi)
    inside = inside_gse(E[0], Lz[0], cont, level=0)
    return MigrationRun(omega_f, hist, t, traj, E, Lz, inside, ic)


def refine_trajectory(run: MigrationRun, i: int, n_times: int = 32001):
    """Re-integrate sample ``i`` alone with fine output (0.25 Myr for tf = 8 Gyr), so that
    pericentre passages (~1 Myr long at 0.5 kpc) are resolved. Returns (t ascending, traj)."""
    agama = _agama()
    pot, _ = slowing_bar_potential(run.hist)
    t, tr = agama.orbit(ic=run.ic[i], potential=pot, timestart=run.hist.tf, time=-run.hist.tf, trajsize=n_times)
    order = np.argsort(t)
    return np.asarray(t)[order], np.asarray(tr)[order]


# ----------------------------------------------------------------------------- grid, picks, products
OMEGA_GRID = tuple(np.arange(20.0, 31.0, 1.0)) + (32.0, 35.0, 40.0)
OMEGA_FIDUCIAL = 24.0          # the paper's forward-GSE case; migration succeeds for <~ 26


def run_grid(omegas=OMEGA_GRID, n_samples: int = 1000, seed: int = 42, n_times: int = 161):
    """Backward integration for each present-day pattern speed. Returns list of MigrationRun."""
    return [back_integrate(w, n_samples=n_samples, seed=seed, n_times=n_times) for w in omegas]


def orbit_elements(traj_1, t, t_lo, t_hi):
    """peri/apo/ecc/zmax of one sample's trajectory (n_times, 6) within t_lo <= t <= t_hi."""
    m = (t >= t_lo) & (t <= t_hi)
    r = np.linalg.norm(traj_1[m, :3], axis=1)
    return dict(r_peri=float(r.min()), r_apo=float(r.max()),
                ecc=float((r.max() - r.min()) / (r.max() + r.min())),
                z_max=float(np.abs(traj_1[m, 2]).max()))


def pick_class3(run: MigrationRun, quantiles=(0.16, 0.5, 0.84), refine: bool = True):
    """Among samples inside the GSE contour at t = 0, the ones at the given quantiles of E(t=0).
    Returns list of dicts (sample index, fine trajectory, elements at 0-1 Gyr and in the last
    Gyr, E0, Lz0). Elements come from a 0.25-Myr re-integration when ``refine`` is set."""
    idx = np.flatnonzero(run.inside_t0)
    if idx.size == 0:
        return []
    E0 = run.E[0, idx]
    picks = []
    for q in quantiles:
        target = np.quantile(E0, q)
        i = idx[np.argmin(np.abs(E0 - target))]
        t, tr = refine_trajectory(run, i) if refine else (run.t, run.traj[:, i])
        picks.append(dict(sample=int(i), q=q, E0=float(run.E[0, i]), Lz0=float(run.Lz[0, i]), t=t, traj=tr,
                          early=orbit_elements(tr, t, 0.0, 1.0),
                          late=orbit_elements(tr, t, run.hist.tf - 1.0, run.hist.tf)))
    return picks


def build_products(out_dir: Path = Path("results/tails"), plot_dir: Path = Path("plots"),
                   omegas=OMEGA_GRID, n_samples: int = 1000, seed: int = 42):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.table import Table
    out_dir.mkdir(parents=True, exist_ok=True); plot_dir.mkdir(parents=True, exist_ok=True)
    runs = run_grid(omegas, n_samples=n_samples, seed=seed, n_times=161)
    pot_axi = axisymmetric_potential()
    cont = gse_contours(pot_axi)
    fid = back_integrate(OMEGA_FIDUCIAL, n_samples=n_samples, seed=seed, n_times=161)
    picks = pick_class3(fid)
    E_now, Lz_now = energy_lz(present_day_samples(1), pot_axi)

    tab = Table(rows=[dict(omega_f=r.omega_f, omega_1=r.hist.omega_1, frac_inside_gse=r.frac_inside,
                           E0_median=float(np.median(r.E[0])), Lz0_median=float(np.median(r.Lz[0])),
                           E0_p16=float(np.quantile(r.E[0], 0.16)), E0_p84=float(np.quantile(r.E[0], 0.84)))
                      for r in runs])
    tab.meta.update(eta=ETA, t1=T1, t2=T2, tf=TF, n_samples=n_samples, seed=seed,
                    E_now=float(E_now[0]), Lz_now=float(Lz_now[0]),
                    note="Dillamore+2026 set-up; Hunter24 potential; GSE contours Belokurov+2023 (outer level)")
    tab.write(out_dir / "bar_migration_grid.ecsv", overwrite=True)
    np.savez_compressed(out_dir / "bar_migration_fiducial.npz", t=fid.t, traj=fid.traj, E=fid.E, Lz=fid.Lz,
                        inside_t0=fid.inside_t0, omega_f=fid.omega_f, picks=np.array([p["sample"] for p in picks]))
    ptab = Table(rows=[dict(q=p["q"], sample=p["sample"], E0=p["E0"], Lz0=p["Lz0"],
                            **{"early_" + k: v for k, v in p["early"].items()},
                            **{"late_" + k: v for k, v in p["late"].items()}) for p in picks])
    ptab.meta.update(omega_f=OMEGA_FIDUCIAL, early="0-1 Gyr after bar formation (pre-migration)",
                     late="last Gyr before today")
    ptab.write(out_dir / "bar_migration_class3_picks.ecsv", overwrite=True)

    # ---- figure 1: E-Lz at t=0 for several pattern speeds (paper Fig. 4 style) + fraction curve
    show = [r for r in runs if r.omega_f in (22.0, 24.0, 26.0, 30.0, 35.0)]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    axes = axes.ravel()
    for ax, r in zip(axes[:5], show):
        for c in cont:
            ax.plot(c[:, 0], c[:, 1] / 1e5, color="k", lw=0.6)
        ax.scatter(r.Lz[0], r.E[0] / 1e5, s=3, alpha=0.4, c=np.where(r.inside_t0, "C2", "C3"))
        ax.plot(Lz_now, E_now / 1e5, marker="*", ms=12, color="C1", mec="k", ls="none", label="omega Cen today")
        ax.set_title(r"$\Omega_{b,0}$ = %.0f km/s/kpc: %.0f%% inside GSE at t=0" % (r.omega_f, 100 * r.frac_inside), fontsize=10)
        ax.set_xlabel(r"$L_z$ [kpc km/s]"); ax.set_ylabel(r"$E$ [$10^5$ km$^2$ s$^{-2}$]")
        ax.set_xlim(-1500, 1500); ax.set_ylim(-1.7, -0.7)
    axes[0].legend(fontsize=8, loc="lower right")
    ax = axes[5]
    ax.plot(tab["omega_f"], tab["frac_inside_gse"], "o-")
    ax.axvline(26, color="grey", ls="--", lw=0.8); ax.text(26.3, 0.9, "paper: <~ 26", fontsize=8, color="grey")
    ax.set_xlabel(r"present-day $\Omega_{b,0}$ [km/s/kpc]"); ax.set_ylabel("fraction of samples inside GSE contour at t=0")
    ax.set_ylim(0, 1); ax.grid(alpha=0.3)
    fig.suptitle("omega Cen back-integrated 8 Gyr through a growing, decelerating bar (Dillamore+2026 set-up, Hunter+2024 potential)", fontsize=11)
    fig.tight_layout(); fig.savefig(plot_dir / "bar_migration_elz.png", dpi=150); plt.close(fig)

    # ---- figure 2: the three picked histories
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    lb = fid.hist.tf - fid.t     # look-back time
    for p in picks:
        i = p["sample"]; r = np.linalg.norm(p["traj"][:, :3], axis=1); lbp = fid.hist.tf - p["t"]
        lab = "sample %d (E0 quantile %.2f): early peri/apo %.1f/%.1f -> late %.1f/%.1f" % (
            i, p["q"], p["early"]["r_peri"], p["early"]["r_apo"], p["late"]["r_peri"], p["late"]["r_apo"])
        axes[0].plot(lbp, r, lw=0.4, label=lab)
        axes[1].plot(lb, fid.E[:, i] / 1e5, lw=0.8)
        axes[2].plot(lb, fid.Lz[:, i], lw=0.8)
    axes[0].set_yscale("log"); axes[0].set_ylabel("r [kpc]"); axes[0].legend(fontsize=7, loc="upper left")
    axes[1].set_ylabel(r"$E$ [$10^5$ km$^2$ s$^{-2}$]"); axes[2].set_ylabel(r"$L_z$ [kpc km/s]")
    axes[2].set_xlabel("look-back time [Gyr]  (bar forms at 8, grows to 7, decelerates from 7)")
    ax2 = axes[1].twinx(); ax2.plot(fid.hist.tf - fid.hist.t, fid.hist.omega, color="grey", ls=":", lw=1); ax2.set_ylabel(r"$\Omega_b$ [km/s/kpc]", color="grey")
    for ax in axes: ax.grid(alpha=0.3)
    axes[0].set_title(r"Class 3, $\Omega_{b,0}$ = %.0f: three back-integrated omega Cen samples ending inside GSE" % OMEGA_FIDUCIAL, fontsize=10)
    fig.tight_layout(); fig.savefig(plot_dir / "bar_migration_class3_picks.png", dpi=150); plt.close(fig)
    return tab, ptab


if __name__ == "__main__":
    t, p = build_products()
    t.pprint(max_width=200); p.pprint(max_width=250)
