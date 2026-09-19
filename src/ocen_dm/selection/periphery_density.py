"""Is there an overdensity of cluster-like stars at the radii where the spectroscopy reaches?

The kinematic plateau at 6-8 km/s beyond 60 pc (see :mod:`ocen_dm.kinematics.periphery`) is
only meaningful if the stars measured there are really omega Cen's. That is a **star-count**
question, not a kinematic one, and the project had never asked it beyond 42 arcmin.

The test uses the Pristine periphery catalogue (Kuzma 2025; Gaia proper motions, CaHK
metallicities, G0 < 16, uniform coverage to 5.1 deg -- the all-star surface density varies by
only 6 per cent between 0.67 and 5.1 deg, so counts divided by annulus area need no footprint
correction). Cluster-like means **both**:

* proper motion within ``window_masyr`` of the systemic motion, and
* ``[Fe/H] < feh_max``: the field median is -0.26, the cluster's inner stars -1.49.

The background is not assumed. It is measured by repeating the identical selection in eight
control windows of the same radius placed at the same distance from the field's own
proper-motion centroid, which samples the field at the same proper-motion amplitude.

Measured 2026-09-19 (``window_masyr=0.8``, ``feh_max=-1.2``): a decisive excess to about
100 pc, falling like r^-6, and nothing in a circular average beyond 130 pc.

Beyond 1.6 deg the 23 surviving candidates are **not** isotropic: an axial Rayleigh test
gives a preferred axis at position angle 148 deg with p = 0.016, and the along-axis excess is
0.34 +- 0.16 per square degree against 0.00 +- 0.16 across it. The spectroscopic members
beyond 1.2 deg sit on the same axis. That is what tidal tails look like, but the axis was
chosen after seeing the data and the significance is 2 sigma, so it is a lead and not a
detection.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from ..kinematics.outer_profile import OCEN_DEC, OCEN_RA
from ..paths import processed_dir

__all__ = ["MU_SYS", "TAIL_AXIS_PA_DEG", "load_periphery", "candidate_masks",
           "density_profile", "axial_rayleigh", "along_across"]

#: systemic proper motion used to centre the selection window (mas/yr)
MU_SYS = (-3.257, -6.730)
#: preferred axis of the outer candidates, position angle east of the x (east) axis
TAIL_AXIS_PA_DEG = 148.0


def load_periphery() -> tuple[Table, np.ndarray, np.ndarray]:
    """The catalogue, the radius in degrees and the position angle in degrees."""
    t = Table.read(processed_dir() / "tails" / "kuzma2025_periphery.ecsv")
    x = (np.asarray(t["ra"], float) - OCEN_RA) * np.cos(np.radians(OCEN_DEC))
    y = np.asarray(t["dec"], float) - OCEN_DEC
    return t, np.hypot(x, y), np.degrees(np.arctan2(y, x))


def candidate_masks(t: Table, window_masyr: float = 0.8, feh_max: float = -1.2,
                    n_control: int = 8) -> tuple[np.ndarray, list[np.ndarray]]:
    """Cluster-like selection, plus control selections that sample the field the same way."""
    pa, pd = np.asarray(t["pmra"], float), np.asarray(t["pmdec"], float)
    metal = np.asarray(t["feh"], float) < feh_max
    sig = (np.hypot(pa - MU_SYS[0], pd - MU_SYS[1]) < window_masyr) & metal
    centre = np.array([np.median(pa), np.median(pd)])
    v = np.array(MU_SYS) - centre
    radius, angle0 = np.hypot(*v), np.arctan2(v[1], v[0])
    step = 2 * np.pi / (n_control + 1)
    ctrl = []
    for i in range(1, n_control + 1):
        c = centre + radius * np.array([np.cos(angle0 + i * step), np.sin(angle0 + i * step)])
        ctrl.append((np.hypot(pa - c[0], pd - c[1]) < window_masyr) & metal)
    return sig, ctrl


def density_profile(edges_deg=(0.5, 0.7, 0.9, 1.2, 1.6, 2.1, 2.8, 3.6, 4.4, 5.1),
                    distance_kpc: float = 5.43, **kw) -> Table:
    """Background-subtracted surface density of cluster-like stars, per square degree."""
    t, r, _ = load_periphery()
    sig, ctrl = candidate_masks(t, **kw)
    rows = []
    for lo, hi in zip(edges_deg[:-1], edges_deg[1:]):
        m = (r >= lo) & (r < hi)
        area = np.pi * (hi ** 2 - lo ** 2)
        n = int((sig & m).sum())
        cs = np.array([float((c & m).sum()) for c in ctrl])
        bg, bg_sd = cs.mean() / area, cs.std(ddof=1) / area
        excess = n / area - bg
        err = float(np.hypot(np.sqrt(max(n, 1)) / area, bg_sd / np.sqrt(len(ctrl))))
        rows.append((lo, hi, float(0.5 * (lo + hi) * np.pi / 180 * distance_kpc * 1e3), area,
                     n, float(cs.mean()), n / area, bg, excess, err, excess / err))
    return Table(rows=rows, names=("r_lower_deg", "r_upper_deg", "r_pc", "area_deg2", "n_signal",
                                   "n_control_mean", "density", "background", "excess",
                                   "excess_err", "significance"))


def axial_rayleigh(angles_deg: np.ndarray) -> tuple[float, float, float]:
    """Rayleigh test on the doubled angle: returns (R, p, preferred axis in deg)."""
    a = np.radians(np.asarray(angles_deg, float)) * 2.0
    n = len(a)
    c, s = np.cos(a).mean(), np.sin(a).mean()
    rr = float(np.hypot(c, s))
    return rr, float(np.exp(-n * rr ** 2)), float(np.degrees(np.arctan2(s, c)) / 2.0 % 180.0)


def along_across(r_min: float = 1.6, r_max: float = 5.1, half_width_deg: float = 30.0,
                 axis_pa_deg: float = TAIL_AXIS_PA_DEG, **kw) -> Table:
    """Excess density inside a wedge about the axis, against the perpendicular wedge."""
    t, r, th = load_periphery()
    sig, ctrl = candidate_masks(t, **kw)
    dth = np.abs(((th - axis_pa_deg + 90.0) % 180.0) - 90.0)
    frac = half_width_deg / 90.0
    area = np.pi * (r_max ** 2 - r_min ** 2) * frac
    rows = []
    for label, wedge in (("along", dth < half_width_deg), ("across", dth >= 90.0 - half_width_deg)):
        m = (r >= r_min) & (r < r_max) & wedge
        n = int((sig & m).sum())
        cs = float(np.mean([(c & m).sum() for c in ctrl]))
        excess = (n - cs) / area
        err = float(np.sqrt(max(n, 1) + cs / len(ctrl)) / area)
        rows.append((label, n, cs, area, excess, err, excess / err))
    return Table(rows=rows, names=("wedge", "n_signal", "n_control_mean", "area_deg2",
                                   "excess", "excess_err", "significance"))
