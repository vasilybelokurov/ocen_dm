"""Which bar pattern speeds put omega Cen's *present* orbit in resonance with the bar?

Orbital frequencies from AGAMA's action finder in the axisymmetrised Hunter et al. (2024)
potential; a resonance with a prograde bar requires
    l Omega_r + n Omega_z + m (Omega_phi - Omega_b) = 0
    =>  Omega_b = Omega_phi + (l Omega_r + n Omega_z) / m .
omega Cen is retrograde (Omega_phi < 0), so only combinations with enough radial/vertical
frequency give Omega_b > 0. Usage: PYTHONPATH=src python -m ocen_dm.tails.resonances
"""
from __future__ import annotations
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from astropy.table import Table
from . import bar_migration as bm

#: (name, l, n, m)
RESONANCES = (("corotation", 0, 0, 1), ("OLR-type (1,0,2)", 1, 0, 2), ("(3,0,4)", 3, 0, 4),
              ("vertical (0,1,1)", 0, 1, 1), ("(1,1,2)", 1, 1, 2), ("retrograde 1:1 (1,0,1)", 1, 0, 1),
              ("(5,0,4)", 5, 0, 4), ("(4,0,3)", 4, 0, 3), ("(3,0,2)", 3, 0, 2), ("(1,1,1)", 1, 1, 1))


def present_frequencies(n_samples: int = 1000, seed: int = 42):
    """(Omega_r, Omega_z, Omega_phi) [km/s/kpc] and actions for the present-day samples."""
    agama = bm._agama()
    af = agama.ActionFinder(bm.axisymmetric_potential())
    act, _, fr = af(bm.present_day_samples(n_samples, seed), angles=True)
    return fr, act


def resonance_table(n_samples: int = 1000, seed: int = 42) -> Table:
    fr, act = present_frequencies(n_samples, seed)
    Or, Oz, Op = fr[:, 0], fr[:, 1], fr[:, 2]
    rows = []
    for name, l, n, m in RESONANCES:
        ob = Op + (l * Or + n * Oz) / m
        p16, p50, p84 = np.percentile(ob, [16, 50, 84])
        rows.append(dict(resonance=name, l=l, n=n, m=m, omega_b=p50, omega_b_16=p16, omega_b_84=p84, physical=bool(p50 > 0)))
    t = Table(rows=rows)
    t.meta.update(Omega_r=float(np.median(Or)), Omega_z=float(np.median(Oz)), Omega_phi=float(np.median(Op)),
                  Jr=float(np.median(act[:, 0])), Jz=float(np.median(act[:, 1])), Lz=float(np.median(act[:, 2])),
                  potential="Hunter+2024 axisymmetrised", n_samples=n_samples)
    return t


if __name__ == "__main__":
    t = resonance_table()
    print("Omega_r %.1f Omega_z %.1f Omega_phi %.1f km/s/kpc" % (t.meta["Omega_r"], t.meta["Omega_z"], t.meta["Omega_phi"]))
    t.pprint(max_width=150)
    from pathlib import Path
    Path("results/tails").mkdir(parents=True, exist_ok=True); t.write("results/tails/resonances.ecsv", overwrite=True)
