"""Can we reproduce the published Gaia EDR3 dispersion profile from the published stars?

Vasiliev & Baumgardt (2021, MNRAS 505, 5978; https://arxiv.org/abs/2102.09568) release both
the member catalogue and the PM dispersion profile derived from it. Our own binned
measurement sits 4-6 per cent above their profile between 460 and 1000 arcsec, which is many
times their quoted uncertainty, so either they or we are wrong. This module answers it by
applying their documented prescriptions to their own stars in our annuli.

Two ingredients, both from their Section 3 and Table 1:

* **Error inflation.** Gaia's formal PM uncertainties are underestimated in crowded fields.
  They multiply them by ``eta = (1 + Sigma / Sigma_0) ** zeta`` with the local source density
  ``Sigma`` in stars per square arcmin, adopting for the PM the same coefficients they fit for
  the parallax: ``zeta = 0.04`` with ``Sigma_0 = 10`` for 5-parameter solutions and
  ``Sigma_0 = 5`` for 6-parameter ones (their Table 1, 'clean' rows). The released catalogue
  carries the **raw** Gaia errors -- verified star by star against ``gaia_edr3.gaia_source``
  on WSDB, median ratio 1.0000 for all 228055 members -- so the scaling is theirs to apply at
  fit time, and anyone deconvolving the catalogue errors directly under-removes the noise.
* **Low-error selection.** Only stars with ``err < kappa * sigma(R) / (eta - eta_0)``
  (``kappa = 0.2``, ``eta_0 = 0.9``) constrain the dispersion; the rest enter the astrometric
  fit convolved with a fixed profile. This drops roughly half the stars beyond 700 arcsec.

Measured outcome (2026-09-19): the inflation alone removes the discrepancy. Our annuli move
from +5.0, +4.6, +6.1, +4.8 per cent at 521, 635, 773 and 941 arcsec to +1.5, -0.6, 0.0 and
-1.0. The low-error cut changes little. **Our measurement was biased high because it
deconvolved raw Gaia errors.**

Remaining differences, which this module does *not* try to remove, are structural: their
dispersion is a cubic spline with 2-5 nodes over the whole radial range fitted jointly with a
Plummer cluster profile and a two-Gaussian field by MCMC, so its quoted 0.002 mas/yr
uncertainties are the posterior width of a stiff global model, not the precision of an
independent measurement in any annulus. Their rotation is one free amplitude on a fixed
functional form, their radial mean is fixed to perspective expansion, and their dispersion is
isotropic by construction (their Section 6 fits anisotropy separately).
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from ..paths import processed_dir
from .outer_profile import dispersion_ml, load_members

__all__ = ["ETA_COEFFS", "KAPPA", "ETA_FLOOR", "error_inflation", "low_error_mask",
           "published_profile", "replication_table"]

#: (Sigma_0 [stars/arcmin^2], zeta) for the 'clean' 5p and 6p subsets, their Table 1
ETA_COEFFS = {"5p": (10.0, 0.04), "6p": (5.0, 0.04)}
KAPPA = 0.2
ETA_FLOOR = 0.9


def error_inflation(source_density: np.ndarray, five_parameter: np.ndarray) -> np.ndarray:
    """``eta = (1 + Sigma / Sigma_0) ** zeta``, their Table 1 'clean' rows."""
    s5, z5 = ETA_COEFFS["5p"]; s6, z6 = ETA_COEFFS["6p"]
    return np.where(five_parameter, (1.0 + source_density / s5) ** z5,
                    (1.0 + source_density / s6) ** z6)


def low_error_mask(err: np.ndarray, sigma_at_r: np.ndarray, eta: np.ndarray) -> np.ndarray:
    """Their dispersion-constraining subset: ``err < kappa * sigma(R) / (eta - eta_0)``."""
    return err < KAPPA * sigma_at_r / np.maximum(eta - ETA_FLOOR, 1e-3)


def published_profile() -> tuple[np.ndarray, np.ndarray]:
    """Their tabulated radius (arcsec) and median PM dispersion (mas/yr)."""
    t = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    return np.asarray(t["r"], float), np.asarray(t["sigma_pm"], float)


def replication_table(edges: np.ndarray | None = None, prob_min: float = 0.9,
                      distance_kpc: float = 5.43) -> Table:
    """Our annulus dispersions with and without their prescriptions, against their profile.

    Columns: ``sigma_raw`` deconvolves the catalogue errors as released, ``sigma_eta``
    applies the inflation, ``sigma_eta_lowerr`` adds their low-error selection, and each
    ``*_ratio`` is that value divided by their profile interpolated to the same radius.
    """
    edges = np.geomspace(300.0, 2400.0, 11) if edges is None else np.asarray(edges, float)
    rp, sp = published_profile()
    cat = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    s = load_members(exact=True, distance_kpc=distance_kpc)
    qf = np.asarray(cat["quality_flag"], int)
    eta = error_inflation(np.asarray(cat["source_density"], float), (qf & 1) > 0)
    clean = (qf & 2) > 0
    err_mean = 0.5 * (s.err_r + s.err_t)
    lowerr = low_error_mask(err_mean, np.interp(s.r_arcsec, rp, sp), eta)
    ones = np.ones_like(eta)

    def sig(mask: np.ndarray, scale: np.ndarray) -> float:
        a = dispersion_ml(s.mu_r[mask], s.err_r[mask] * scale[mask])
        b = dispersion_ml(s.mu_t[mask], s.err_t[mask] * scale[mask])
        return float(np.sqrt(0.5 * (a[0] ** 2 + b[0] ** 2)))

    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = clean & (s.prob > prob_min) & (s.r_arcsec >= lo) & (s.r_arcsec < hi)
        if m.sum() < 20:
            continue
        rm = float(np.median(s.r_arcsec[m])); pv = float(np.interp(rm, rp, sp))
        m2 = m & lowerr
        raw, sc = sig(m, ones), sig(m, eta)
        lo_e = sig(m2, eta) if m2.sum() > 20 else np.nan
        rows.append((lo, rm, hi, int(m.sum()), int(m2.sum()), float(np.median(eta[m])), pv,
                     raw, raw / pv, sc, sc / pv, lo_e, lo_e / pv))
    return Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_stars", "n_low_error",
                                   "median_eta", "sigma_published", "sigma_raw", "raw_ratio",
                                   "sigma_eta", "eta_ratio", "sigma_eta_lowerr", "eta_lowerr_ratio"))
