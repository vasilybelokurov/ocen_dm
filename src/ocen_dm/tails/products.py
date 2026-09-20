"""Build the plausible progenitor-orbit set: summary table, orbit tracks and figure.

Usage: PYTHONPATH=src python -m ocen_dm.tails.products
Outputs: results/tails/progenitor_orbits_summary.ecsv, results/tails/progenitor_orbit_tracks.npz,
         plots/progenitor_orbits.png
"""
from __future__ import annotations
import warnings
from pathlib import Path
import numpy as np


def build_products(out_dir: Path = Path("results/tails"), plot_dir: Path = Path("plots")):
    warnings.filterwarnings("ignore")
    from astropy.table import Table
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .progenitor_orbits import plausible_orbit_set, PRESENT

    out_dir.mkdir(parents=True, exist_ok=True); plot_dir.mkdir(parents=True, exist_ok=True)
    S = plausible_orbit_set()
    rows, tracks = [], {}
    nan = dict(r_infall_kpc=np.nan, m_infall=np.nan, tau_gyr=np.nan, peri_5_6=np.nan,
               apo_5_6=np.nan, peri_2_3=np.nan, apo_2_3=np.nan)
    for cls, orbs in S.items():
        for o in orbs:
            t = o.t_gyr
            extra = dict(nan)
            if cls == "class3":
                e = o.notes["early_0_1_gyr"]
                extra.update(peri_5_6=e["r_peri"], apo_5_6=e["r_apo"])   # pre-migration elements (7-8 Gyr ago)
            if cls == "class2":
                s = lambda lo, hi: o.r[(t <= -lo) & (t >= -hi)]
                extra.update(r_infall_kpc=o.notes["r_at_infall_kpc"], m_infall=o.notes["m_infall"],
                             tau_gyr=o.notes["tau_gyr"], peri_5_6=s(5, 6).min(), apo_5_6=s(5, 6).max(),
                             peri_2_3=s(2, 3).min(), apo_2_3=s(2, 3).max())
            rows.append(dict(cls=cls, label=o.label, potential=o.potential, r_peri_kpc=o.r_peri,
                             r_apo_kpc=o.r_apo, ecc=o.eccentricity, z_max_kpc=float(np.abs(o.z).max()),
                             **extra))
            key = "%s|%s" % (cls, o.label)
            for name, arr in (("t", t), ("r", o.r), ("R", o.R), ("z", o.z)):
                tracks[key + "|" + name] = arr
    tab = Table(rows=rows)
    tab.meta.update(present=PRESENT,
                    note="class2: exponential stripping since infall (peri/apo_5_6 = 5-6 Gyr ago); class3: "
                         "Dillamore+2026 back-integration, peri/apo_5_6 columns hold the pre-migration "
                         "elements 7-8 Gyr ago")
    tab.write(out_dir / "progenitor_orbits_summary.ecsv", overwrite=True)
    np.savez_compressed(out_dir / "progenitor_orbit_tracks.npz", **tracks)

    fig, axes = plt.subplots(3, 1, figsize=(9, 10.5))
    titles = ["Class 1: present orbit, static potentials (backwards 10 Gyr)",
              "Class 2: present orbit + dynamical friction backwards (satellite stripped since infall)",
              "Class 3: back-integrated through the growing, decelerating bar (Dillamore+2026 set-up, "
              r"Hunter+2024, $\Omega_{b,0}$ = 24)"]
    for ax, cls, title in zip(axes, ("class1", "class2", "class3"), titles):
        for o in S[cls]:
            ax.plot(-o.t_gyr, o.r, lw=0.8, label=o.label)
        ax.set_yscale("log"); ax.set_ylabel("r [kpc]"); ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.3)
        ax.axhline(7.0, color="k", ls=":", lw=0.7)
    axes[0].set_xlabel("look-back time [Gyr]"); axes[1].set_xlabel("look-back time [Gyr]")
    axes[2].set_xlabel("look-back time [Gyr]   (bar forms at 8, grows to 7, decelerates from 7; migration in the last ~2.5)")
    axes[1].axhline(100, color="grey", ls="--", lw=0.7)
    axes[1].text(5.0, 115, "~ virial radius of the young Milky Way", fontsize=8, color="grey")
    fig.tight_layout(); fig.savefig(plot_dir / "progenitor_orbits.png", dpi=150); plt.close(fig)
    return tab


if __name__ == "__main__":
    print(build_products())
