"""Like-for-like literature presets: our pipeline under a published analysis's assumptions.

Each preset fixes what that paper fixed (datasets, distance, anisotropy family,
components switched off) so that the comparison tests the *pipeline*, not the
choice of assumptions. Published values carry their source; derived quantities
(M/L) state every input.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .fit import DarkMatterModel, NoDarkMatterModel, Prior

__all__ = ["Preset", "PRESETS", "build_preset", "v_band_luminosity"]

#: integrated V magnitude and reddening of omega Cen (Harris 1996, 2010 edition,
#: https://physwww.mcmaster.ca/~harris/mwgc.dat): V_t = 3.68, E(B-V) = 0.12.
HARRIS_V_T = 3.68
HARRIS_EBV = 0.12
M_V_SUN = 4.83


def v_band_luminosity(distance_kpc: float, v_t: float = HARRIS_V_T, ebv: float = HARRIS_EBV, r_v: float = 3.1) -> float:
    """Total V-band luminosity in Lsun at ``distance_kpc`` (A_V = R_V E(B-V))."""
    m_v = v_t - r_v * ebv - 5.0 * np.log10(distance_kpc * 1e3 / 10.0)
    return 10.0 ** (-0.4 * (m_v - M_V_SUN))


@dataclass(frozen=True)
class Preset:
    name: str
    reference: str
    url: str
    datasets: tuple[str, ...]
    family_kwargs: dict[str, Any]
    published: dict[str, tuple[float, float, str]]      # quantity -> (value, error, note)
    derived: tuple[str, ...] = ()                       # which of our derived quantities to report
    notes: str = ""

    def family(self):
        kw = dict(self.family_kwargs)
        gamma = kw.pop("gamma", None)
        return NoDarkMatterModel(**kw) if gamma is None else DarkMatterModel(gamma=gamma, **kw)


PRESETS: dict[str, Preset] = {
    "watkins2013": Preset(
        name="watkins2013",
        reference="Watkins, van de Ven, den Brok & van den Bosch 2013, MNRAS 436, 2598 (discrete axisymmetric Jeans on HST PMs)",
        url="https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2598W",
        datasets=("hst_pm_radial", "hst_pm_tangential"),
        family_kwargs=dict(tracer="trager", fix_distance=True, distance_kpc=4.59, instruments=(),
                           constant_beta=True, beta0_max=0.5, fixed={"M_rem": 1e4, "a_rem": 1.0, "M_bh": 1e2}),
        published={"ML_V": (2.71, 0.05, "M/L_V in Msun/Lsun at D = 4.59 kpc (their Table 3, isotropic-rotator model)"),
                   "beta_0": (0.10, 0.02, "global anisotropy; positive = radial")},
        derived=("ML_V",),
        notes="Our proper motions are the oMEGACat ones, not the Watkins+ 2013 sample; the distance, "
              "constant M/L, constant beta and light-weighted (Trager) tracer follow the paper. "
              "M/L_V uses Harris V_t = 3.68, E(B-V) = 0.12, R_V = 3.1 at 4.59 kpc.",
    ),
    "omegacat6": Preset(
        name="omegacat6",
        reference="Haeberle et al. 2025 (oMEGACat VI): kinematic distance from LOS vs PM dispersions",
        url="https://ui.adsabs.harvard.edu/abs/2025ApJ...983...95H",
        datasets=("hst_pm_radial", "hst_pm_tangential", "muse_los_dispersion"),
        family_kwargs=dict(tracer="composite", distance_prior=Prior("uniform", 4.5, 6.5), instruments=("MUSE",)),
        published={"distance": (5.494, 0.061, "kinematic distance, kpc")},
        notes="Flat distance prior so the number is ours, not the prior's; the MUSE/HST equipartition "
              "scale is left free as in the K1 fits (their analysis treats the samples differently).",
    ),
    "baumgardt2018": Preset(
        name="baumgardt2018",
        reference="Baumgardt & Hilker 2018, MNRAS 478, 1520 (N-body fits to dispersion profiles)",
        url="https://ui.adsabs.harvard.edu/abs/2018MNRAS.478.1520B",
        datasets=("hst_pm_radial", "hst_pm_tangential", "muse_los_dispersion", "gaia_dr2_pm"),
        family_kwargs=dict(tracer="composite", fix_distance=True, distance_kpc=5.24, instruments=("MUSE", "GaiaDR2")),
        published={"M_total": (3.55e6, 0.03e6, "total mass, Msun (docs/LITERATURE_BASELINES.md); the distance used "
                                                "for it is UNVERIFIED here -- 5.24 kpc assumed from the 2018 catalogue")},
        derived=("M_total",),
        notes="UNVERIFIED assumption: their 2018 catalogue distance taken as 5.24 kpc; check against the paper's "
              "Table 1 before quoting the comparison. Total mass = stars + remnants + BH.",
    ),
    "imbh_limit": Preset(
        name="imbh_limit",
        reference="van der Marel & Anderson 2010, ApJ 710, 1063 (HST PMs, anisotropic models: M_BH <~ 1.2e4 at 1 sigma, "
                  "1.8e4 at 3 sigma); Baumgardt et al. 2019, MNRAS 488, 5340 (no IMBH; 4.6 per cent of the mass in stellar BHs)",
        url="https://ui.adsabs.harvard.edu/abs/2010ApJ...710.1063V",
        datasets=("hst_pm_radial", "hst_pm_tangential", "muse_los_dispersion"),
        family_kwargs=dict(tracer="composite", instruments=("MUSE",)),
        published={"M_bh": (1.2e4, np.nan, "1-sigma upper limit on an IMBH mass, Msun (van der Marel & Anderson 2010)")},
        notes="Compared through the posterior 95 per cent upper limit on M_bh; remnants (M_rem, a_rem) stay free, "
              "as the dark-cluster alternative of both papers.",
    ),
}


def build_preset(name: str):
    """``(preset, family, datasets)`` for ``ocen fit --preset``."""
    if name not in PRESETS:
        raise KeyError(f"unknown preset {name!r}; known: {sorted(PRESETS)}")
    p = PRESETS[name]
    return p, p.family(), list(p.datasets)


def derived_quantities(preset: Preset, theta: dict[str, float], family) -> dict[str, float]:
    """Our values of the preset's published quantities for one parameter vector."""
    full = family.complete(theta)
    D = family.distance(full)
    out: dict[str, float] = {}
    if "ML_V" in preset.derived:
        out["ML_V"] = full["M_star"] / v_band_luminosity(D)
    if "M_total" in preset.derived:
        out["M_total"] = full["M_star"] + full["M_rem"] + full["M_bh"]
        if "M_dm_100" in full:
            out["M_total"] += full["M_dm_100"]
    for k in ("beta_0", "distance", "M_bh"):
        if k in preset.published and k in full:
            out[k] = float(full[k])
    return out
