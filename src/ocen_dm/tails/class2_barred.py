"""Class 2 in a barred potential: does the backward friction orbit change when the host has a
rotating (growing, decelerating) bar?

The class-2 picks (exponential stripping since infall 10 Gyr ago) are re-integrated in the
Hunter et al. (2024) potential in three versions: axisymmetrised (baseline, isolates the
potential change from McMillan 2017), with the Dillamore+2026 decelerating bar ending at the
mainstream Omega_b,0 = 37.5 km/s/kpc (Hunter's own value), and ending at the slow 24 km/s/kpc.
The bar is 8 Gyr old, so for look-back times > 8 Gyr the potential is axisymmetric (checked:
the AGAMA scale modifier holds A = 0 before t = 0). Friction uses the time-dependent total
density. 20 error samples per case give the spread.

Usage: PYTHONPATH=src python -m ocen_dm.tails.class2_barred
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
from scipy.special import erf

from . import bar_migration as bm
from .progenitor_orbits import CLASS2_SET, exponential_stripping

_G = 4.30091e-6; _KPC_PER_KMS_GYR = 1.02271
OUT = Path("results/tails"); PLOTS = Path("plots")
CASES = (("axisymmetric", None), ("bar, Omega_b0 = 37.5", 37.5), ("bar, Omega_b0 = 24", 24.0))


def integrate(pot, pot_axi, mass_history, r_half_kpc, n_samples=20, seed=42, t_gyr=10.0,
              dt_myr=0.25, n_out=401, t_today=8.0):
    """Backward leapfrog with Chandrasekhar friction for bound mass M(t_lookback); potential
    time coordinate t = t_today - t_lookback (bar defined on 0..8)."""
    ic = bm.present_day_samples(n_samples, seed)
    x = ic[:, :3].copy(); v = ic[:, 3:].copy()
    dt = -dt_myr * 1e-3; n = int(round(t_gyr / (dt_myr * 1e-3)))

    def accel(x, v, t):
        a = pot.force(x, t=t) / _KPC_PER_KMS_GYR
        m = mass_history(t_today - t)
        if m > 1e6:
            r = np.linalg.norm(x, axis=1); vm = np.linalg.norm(v, axis=1)
            rho = pot.density(x, t=t)
            R = np.hypot(x[:, 0], x[:, 1]); xr = np.column_stack([R, 0 * R, 0 * R])
            vc = np.sqrt(np.maximum(-R * pot_axi.force(xr)[:, 0], 1.0))
            X = vm / vc
            lnL = np.maximum(np.log(r / np.maximum(r_half_kpc, _G * m / vm**2)), 0.0)
            coeff = 4 * np.pi * _G**2 * m * rho * lnL * (erf(X) - 2 * X / np.sqrt(np.pi) * np.exp(-X**2))
            a = a - (coeff / vm**3)[:, None] * v / _KPC_PER_KMS_GYR
        return a

    every = max(1, n // (n_out - 1))
    t = t_today; ts, X = [t], [x.copy()]
    a = accel(x, v, t)
    for i in range(1, n + 1):
        vh = v + 0.5 * dt * a
        x = x + dt * vh * _KPC_PER_KMS_GYR
        t += dt
        a = accel(x, vh, t)
        v = vh + 0.5 * dt * a
        if i % every == 0 or i == n:
            ts.append(t); X.append(x.copy())
    ts = np.array(ts); X = np.stack(X)                     # (n_out, n_samples, 3), today first
    return t_today - ts, np.linalg.norm(X, axis=2), X       # look-back time, r, positions


def per_gyr(lb, r):
    """median over samples of the per-Gyr pericentre and apocentre; returns (t_mid, peri, apo)."""
    edges = np.arange(0, lb.max() + 1e-9, 1.0); out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (lb >= lo) & (lb < hi)
        out.append((0.5 * (lo + hi), np.median(r[m].min(0)), np.median(r[m].max(0))))
    return np.array(out).T


def main(n_samples=20):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from astropy.table import Table
    OUT.mkdir(parents=True, exist_ok=True); PLOTS.mkdir(exist_ok=True)
    pot_axi = bm.axisymmetric_potential()
    pots = {}
    for name, w in CASES:
        pots[name] = pot_axi if w is None else bm.slowing_bar_potential(bm.BarHistory(w))[0]
    rows, curves = [], {}
    for label, m_inf, tau, t_inf in CLASS2_SET:
        rh = 1.0 * (m_inf / 1e10) ** (1 / 3); hist = exponential_stripping(m_inf, tau, t_inf)
        for name, w in CASES:
            lb, r, X = integrate(pots[name], pot_axi, hist, rh, n_samples=n_samples, t_gyr=t_inf)
            tm, peri, apo = per_gyr(lb, r)
            curves[(label, name)] = (lb, r, tm, peri, apo)
            late = lb <= 1.0; inf = lb >= t_inf - 1.0
            row = dict(history=label, potential=name, peri_last_gyr=float(np.median(r[late].min(0))),
                       apo_last_gyr=float(np.median(r[late].max(0))), peri_5_6=float(peri[5]), apo_5_6=float(apo[5]),
                       apo_before_infall=float(np.median(r[inf].max(0))), apo_before_infall_16=float(np.quantile(r[inf].max(0), .16)),
                       apo_before_infall_84=float(np.quantile(r[inf].max(0), .84)), r_at_infall=float(np.median(r[-1])))
            rows.append(row)
            print("%-30s %-22s last Gyr peri/apo %.2f/%.1f | 5-6 Gyr %.2f/%.1f | apo before infall %5.1f [%5.1f, %5.1f] | r(t_inf) %5.1f" % (
                label, name, row["peri_last_gyr"], row["apo_last_gyr"], row["peri_5_6"], row["apo_5_6"],
                row["apo_before_infall"], row["apo_before_infall_16"], row["apo_before_infall_84"], row["r_at_infall"]), flush=True)
    tab = Table(rows=rows); tab.meta.update(n_samples=n_samples, seed=42, bar="Dillamore+2026 history, 8 Gyr old")
    tab.write(OUT / "class2_barred.ecsv", overwrite=True)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    ls = {"axisymmetric": "-", "bar, Omega_b0 = 37.5": "--", "bar, Omega_b0 = 24": ":"}
    for k, (label, *_r) in enumerate(CLASS2_SET):
        for name, w in CASES:
            lb, r, tm, peri, apo = curves[(label, name)]
            axes[0, k].plot(lb, r[:, 0], lw=0.5, ls=ls[name], label=name)          # nominal-ish sample 0
            axes[1, k].plot(tm, peri, "o", ls=ls[name], ms=3, color="C0"); axes[1, k].plot(tm, apo, "s", ls=ls[name], ms=3, color="C3")
        axes[0, k].set_yscale("log"); axes[0, k].set_ylim(0.3, 200); axes[0, k].set_title(label, fontsize=10); axes[0, k].grid(alpha=0.3)
        axes[1, k].set_yscale("log"); axes[1, k].set_xlabel("look-back time [Gyr]"); axes[1, k].grid(alpha=0.3)
        axes[1, k].axvline(8, color="grey", lw=0.6); axes[1, k].text(8.1, 0.5, "bar forms", fontsize=7, color="grey")
    axes[0, 0].set_ylabel("r [kpc] (one sample)"); axes[1, 0].set_ylabel("per-Gyr peri (blue) and apo (red), median of samples [kpc]")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Class 2 in the Hunter+2024 potential: axisymmetric vs growing, decelerating bar (solid / dashed 37.5 / dotted 24 km/s/kpc)", fontsize=11)
    fig.tight_layout(); fig.savefig(PLOTS / "class2_barred.png", dpi=150); plt.close(fig)
    return tab


if __name__ == "__main__":
    main()
