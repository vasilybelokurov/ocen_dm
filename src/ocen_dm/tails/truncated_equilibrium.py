"""Equilibrium of dark matter tidally truncated around a heavy nucleus (omega Cen).

Analytic pieces used in docs/dm_capture_constraints.tex:
* a tracer with rho ~ r^-gamma in a Kepler potential Phi = -GM/r has the isotropic Eddington
  DF f(E) ~ (-E)^(gamma - 3/2)  (gamma > 1/2);
* cutting the DF at E_t = -GM/r_J leaves at radius r < r_J the
  following fraction. This energy cut is stricter than an apocentre cut: it can
  also remove high-angular-momentum orbits confined inside r_J. The
  fraction  F_gamma(x) = 1 - I_x(gamma - 1/2, 3/2),  x = r / r_J,  I_x the regularised
  incomplete beta function  -- for gamma = 1: F = 1 - I_x(1/2, 3/2); gamma = 3/2: F = (1-x)^(3/2);
* impulsive tidal heating per pericentre passage (Spitzer 1958; Gnedin & Ostriker 1997) with the
  adiabatic correction of Gnedin, Hernquist & Ostriker (1999), A(x) = (1 + x^2)^-gamma_ad,
  x = omega * tau, for a particle at radius r around the nucleus.
"""
from __future__ import annotations
import numpy as np
from scipy.special import betainc

G_PC = 4.30091e-3          # pc (km/s)^2 / Msun
KMS_MYR_TO_PC = 1.02271    # pc travelled at 1 km/s in 1 Myr


def kepler_truncation_factor(x, gamma=1.0):
    """Fraction of the local density retained at r = x r_J after energy truncation at Phi(r_J),
    for rho ~ r^-gamma in a point-mass potential (isotropic)."""
    x = np.clip(np.asarray(x, float), 0, 1)
    return 1.0 - betainc(gamma - 0.5, 1.5, x)


def shock_heating(r_pc, m_host_enclosed, r_peri_pc, v_peri, m_nuc=3.55e6, gamma_ad=1.5):
    """Impulsive heating of a DM particle at radius r_pc around the nucleus by one pericentre passage.

    Tidal acceleration ~ T r with T = G M_host(<r_peri)/r_peri^3, acting for tau ~ r_peri / v_peri;
    Delta v = T r tau, compared with the local velocity scale sigma^2 ~ G M_nuc / r; the adiabatic
    correction (1 + (omega tau)^2)^-gamma_ad with omega = sqrt(G M_nuc / r^3) suppresses heating
    of orbits faster than the passage. Returns dict (Delta v, sigma, omega tau, dE/E per passage).
    """
    r = np.asarray(r_pc, float)
    T = G_PC * m_host_enclosed / r_peri_pc**3                       # (km/s)^2 / pc^2
    tau_myr = r_peri_pc / v_peri / KMS_MYR_TO_PC                    # Myr
    dv = T * r * tau_myr * KMS_MYR_TO_PC                            # km/s  (a [km^2/s^2/pc] * t)
    sigma = np.sqrt(G_PC * m_nuc / r)
    omega = np.sqrt(G_PC * m_nuc / r**3) * KMS_MYR_TO_PC            # 1/Myr
    x = omega * tau_myr
    A = (1 + x**2) ** (-gamma_ad)
    dE_E = (dv / sigma) ** 2 * A                                    # ~ 2/3 (dv^2/2) / (sigma^2/2), order unity factors dropped
    return dict(dv=dv, sigma=sigma, omega_tau=x, A=A, dE_E=dE_E, tau_myr=tau_myr, T=T)


if __name__ == "__main__":
    print("truncation factor F_gamma(x): x = r/r_J")
    for gamma in (1.0, 1.5):
        print("  gamma = %.1f: " % gamma + "  ".join("x=%.2f F=%.2f" % (x, kepler_truncation_factor(x, gamma)) for x in (0.043, 0.143, 0.286, 0.571)))
    print("  e^-x     : " + "  ".join("x=%.2f F=%.2f" % (x, np.exp(-x)) for x in (0.043, 0.143, 0.286, 0.571)))


def figure(path="plots/truncated_equilibrium.png"):
    """Left: retained fraction F_gamma(r/r_J) vs r/r_J against the e^-x stand-in. Right: heating per
    pericentre passage vs radius for the class-1 (1.57 kpc) and class-3 (0.42 kpc) pericentres."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    x = np.linspace(0, 1, 300)
    ax[0].plot(x, kepler_truncation_factor(x, 1.0), "C0", lw=2, label=r"energy truncation, $\rho\propto r^{-1}$: $1-I_x(1/2,3/2)$")
    ax[0].plot(x, kepler_truncation_factor(x, 1.5), "C3", lw=2, label=r"energy truncation, $\rho\propto r^{-3/2}$: $(1-x)^{3/2}$")
    ax[0].plot(x, np.exp(-x), "k--", lw=1.2, label=r"$e^{-x}$ (stand-in used in the profile figure)")
    for rj, c in ((35, "C6"), (70, "C0")):
        for r in (3, 10, 20):
            ax[0].axvline(r / rj, color=c, lw=0.6, ls=":")
    ax[0].text(0.62, 0.55, "dotted: 3, 10, 20 pc\nat $r_J$ = 35 (pink) / 70 pc (blue)", fontsize=7)
    ax[0].set_xlabel(r"$x = r / r_J$"); ax[0].set_ylabel("fraction of the initial density retained"); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    r = np.linspace(5, 100, 200)
    for (M, rp, vp, lab, c) in ((1.36e10, 1570.0, 387.0, "class 1-2 pericentre 1.57 kpc", "C0"), (2.22e9, 420.0, 507.0, "class 3 pericentre 0.42 kpc", "C6")):
        s = shock_heating(r, M, rp, vp)
        ax[1].plot(r, s["dE_E"], color=c, lw=2, label=lab)
        ax[1].plot(r, (s["dv"] / s["sigma"]) ** 2, color=c, lw=1, ls="--")
    for n, lab in ((1 / 114, "unity after 114 passages (class 1, 10 Gyr)"), (1 / 66, "unity after 66 passages (class 3, 8 Gyr)")):
        ax[1].axhline(n, color="grey", lw=0.8, ls=":"); ax[1].text(6, n * 1.2, lab, fontsize=7, color="grey")
    ax[1].set_yscale("log"); ax[1].set_ylim(1e-6, 1); ax[1].set_xlabel("radius around the nucleus [pc]")
    ax[1].set_ylabel(r"$\Delta E/|E|$ per pericentre passage"); ax[1].set_title("solid: with adiabatic correction; dashed: impulsive", fontsize=9)
    ax[1].legend(fontsize=7, loc="lower right"); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
