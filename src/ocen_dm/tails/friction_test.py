"""The cheap decisive test: does a massive (nucleus + bound DM remnant) system migrate through
the decelerating bar like a test particle?

Class-3 back-integration (Dillamore+2026 set-up, :mod:`bar_migration`) repeated with
Chandrasekhar dynamical friction for a constant bound mass M, all samples integrated together
with a vectorised leapfrog in the time-dependent AGAMA potential. Friction is evaluated with the
total (bar-inclusive) density at the current time, sigma = v_c(R)/sqrt(2) of the axisymmetric
potential, ln Lambda = ln(r / max(r_h, GM/v^2)) (galpy's default), no background streaming.
Integrating backwards flips the friction sign. Constant mass is an upper bound on the friction
of a system that was in fact still losing mass.

Usage: PYTHONPATH=src python -m ocen_dm.tails.friction_test
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
from scipy.special import erf

from . import bar_migration as bm

_G = 4.30091e-6           # kpc (km/s)^2 / Msun
TU = bm.TIME_UNIT_GYR     # natural time unit, Gyr
MASSES = (0.0, 1e7, 1e8, 1e9)
OUT = Path("results/tails"); PLOTS = Path("plots")


def r_half(m):
    """Half-mass radius of the bound remnant: 0.1 kpc (M/1e8)^(1/3), well inside the Jacobi radius."""
    return 0.1 * (max(m, 1e6) / 1e8) ** (1.0 / 3.0)


def integrate_with_friction(omega_f: float, m_sat: float, n_samples: int = 100, seed: int = 42,
                            dt_myr: float = 0.25, n_out: int = 161):
    hist = bm.BarHistory(omega_f)
    pot, pot_axi = bm.slowing_bar_potential(hist)
    ic = bm.present_day_samples(n_samples, seed)
    x = ic[:, :3].copy(); v = ic[:, 3:].copy()
    tf = hist.tf; dt = -dt_myr * 1e-3 / TU; n = int(round(tf / (dt_myr * 1e-3 / TU)))   # natural units
    rh = r_half(m_sat)

    def accel(x, v, t):
        a = pot.force(x, t=t)
        if m_sat > 0:
            r = np.linalg.norm(x, axis=1); vm = np.linalg.norm(v, axis=1)
            rho = pot.density(x, t=t)
            R = np.hypot(x[:, 0], x[:, 1])
            xr = np.column_stack([R, np.zeros_like(R), np.zeros_like(R)])
            vc = np.sqrt(np.maximum(-R * pot_axi.force(xr)[:, 0], 1.0))
            X = vm / vc                                   # v / (sqrt2 sigma), sigma = vc/sqrt2
            bmin = np.maximum(rh, _G * m_sat / vm**2)
            lnL = np.maximum(np.log(r / bmin), 0.0)
            coeff = 4 * np.pi * _G**2 * m_sat * rho * lnL * (erf(X) - 2 * X / np.sqrt(np.pi) * np.exp(-X**2))
            a = a - (coeff / vm**3)[:, None] * v
        return a

    every = max(1, n // (n_out - 1))
    ts, X, V = [tf], [x.copy()], [v.copy()]
    t = tf; a = accel(x, v, t)
    for i in range(1, n + 1):
        vh = v + 0.5 * dt * a
        x = x + dt * vh
        t += dt
        a = accel(x, vh, t)
        v = vh + 0.5 * dt * a
        if i % every == 0 or i == n:
            ts.append(t); X.append(x.copy()); V.append(v.copy())
    ts = np.array(ts)[::-1]; traj = np.stack([np.hstack([xx, vv]) for xx, vv in zip(X, V)])[::-1]
    E = np.empty(traj.shape[:2]); Lz = np.empty_like(E)
    for j in range(traj.shape[0]):
        E[j], Lz[j] = bm.energy_lz(traj[j], pot_axi)
    cont = bm.gse_contours(pot_axi)
    inside = bm.inside_gse(E[0], Lz[0], cont)
    return dict(omega_f=omega_f, m_sat=m_sat, t=ts, traj=traj, E=E, Lz=Lz, inside=inside,
                frac=float(inside.mean()), hist=hist)


def summarise(res):
    r0 = np.linalg.norm(res["traj"][0, :, :3], axis=1)   # radius 8 Gyr ago
    tl = res["hist"].tf - res["t"]
    m = (tl >= 7) & (tl <= 8)                              # pre-migration window
    r = np.linalg.norm(res["traj"][m][:, :, :3], axis=2)
    return dict(frac=res["frac"], E0=float(np.median(res["E"][0])), Lz0=float(np.median(res["Lz"][0])),
                E0_16=float(np.quantile(res["E"][0], .16)), E0_84=float(np.quantile(res["E"][0], .84)),
                Lz0_16=float(np.quantile(res["Lz"][0], .16)), Lz0_84=float(np.quantile(res["Lz"][0], .84)),
                apo_pre=float(np.median(r.max(0))), peri_pre=float(np.median(r.min(0))))


def main(omegas=(24.0,), masses=MASSES, n_samples=100):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from astropy.table import Table
    OUT.mkdir(parents=True, exist_ok=True); PLOTS.mkdir(exist_ok=True)
    rows, keep = [], {}
    ref = bm.back_integrate(omegas[0], n_samples=n_samples, n_times=161)
    print("reference agama.orbit, Omega_f %.1f, M=0: frac %.2f  medE0 %.4e medLz0 %.0f" % (omegas[0], ref.frac_inside, np.median(ref.E[0]), np.median(ref.Lz[0])))
    for w in omegas:
        for m in masses:
            res = integrate_with_friction(w, m, n_samples=n_samples)
            s = summarise(res); rows.append(dict(omega_f=w, m_sat=m, r_half_kpc=r_half(m), **s)); keep[(w, m)] = res
            print("Omega_f %.1f  M %.0e : frac inside GSE %.2f | E0 med %.4e [%.4e, %.4e] | Lz0 med %5.0f [%5.0f, %5.0f] | pre-migration peri/apo %.2f/%.1f" % (
                w, m, s["frac"], s["E0"], s["E0_16"], s["E0_84"], s["Lz0"], s["Lz0_16"], s["Lz0_84"], s["peri_pre"], s["apo_pre"]), flush=True)
    tab = Table(rows=rows); tab.meta.update(n_samples=n_samples, seed=42, note="constant bound mass; Chandrasekhar friction in the time-dependent barred potential")
    tab.write(OUT / "friction_test.ecsv", overwrite=True)

    w = omegas[0]; pa = bm.axisymmetric_potential(); cont = bm.gse_contours(pa)
    E_now, Lz_now = bm.energy_lz(bm.present_day_samples(1), pa)
    fig, axes = plt.subplots(2, len(masses), figsize=(4.2 * len(masses), 8))
    for k, m in enumerate(masses):
        res = keep[(w, m)]; ax = axes[0, k]
        for c in cont: ax.plot(c[:, 0], c[:, 1] / 1e5, "k", lw=0.6)
        ax.scatter(res["Lz"][0], res["E"][0] / 1e5, s=4, c=np.where(res["inside"], "C2", "C3"), alpha=0.6)
        ax.plot(Lz_now, E_now / 1e5, "*", ms=12, color="C1", mec="k")
        ax.set_xlim(-1500, 1500); ax.set_ylim(-1.7, -0.7); ax.set_xlabel(r"$L_z$ [kpc km/s]"); ax.set_ylabel(r"$E$ [$10^5$ km$^2$ s$^{-2}$]")
        ax.set_title("M = %.0e M$_\\odot$: %.0f%% in GSE 8 Gyr ago" % (m, 100 * res["frac"]) if m > 0 else "test particle: %.0f%% in GSE" % (100 * res["frac"]), fontsize=10)
        ax = axes[1, k]; lb = res["hist"].tf - res["t"]
        for i in range(0, res["traj"].shape[1], 10):
            ax.plot(lb, res["Lz"][:, i], lw=0.4, c="C2" if res["inside"][i] else "C3", alpha=0.7)
        ax.set_xlabel("look-back time [Gyr]"); ax.set_ylabel(r"$L_z(t)$ [kpc km/s]"); ax.set_ylim(-1500, 800); ax.grid(alpha=0.3)
    fig.suptitle(r"Back-integration through the decelerating bar ($\Omega_{b,0}$ = %.0f) with Chandrasekhar friction for a constant bound mass" % w, fontsize=11)
    fig.tight_layout(); fig.savefig(PLOTS / "friction_test_elz.png", dpi=150); plt.close(fig)
    return tab


if __name__ == "__main__":
    main()
