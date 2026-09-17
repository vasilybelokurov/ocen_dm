"""Perspective and projection effects on proper motions across an extended field.

A cluster moving with a single 3-D velocity does not show a single proper motion across
a 1-degree field (Vasiliev & Belokurov 2020, MNRAS 497, 4162, section 5.1; van de Ven
et al. 2006 for omega Cen itself). For a star at sky position ``n_i`` and distance
``D (1 + zeta_i)``:

* the systemic velocity projects onto the star's **own** tangent basis, not the centre's
  -- the local east/north directions rotate across the field by ~ dRA sin(dec), and the
  line-of-sight component of the systemic velocity acquires a transverse part
  ``-(v_los / D) theta`` (perspective contraction/expansion);
* a star at depth ``zeta`` appears to move slower by ``mu_sys zeta``; the depth is not
  known star by star, so this term survives as an **apparent dispersion**
  ``|mu_sys| sigma_z(R) / D`` along the direction of the systemic proper motion.

The first item is removed exactly by :func:`systemic_pm_field`; the second is modelled in
the likelihood by :func:`depth_dispersion`, using the tracer's line-of-sight depth at each
projected radius. Both effects are of order 0.05-0.1 mas/yr at 2000 arcsec for omega Cen,
against a measured dispersion of 0.22 mas/yr.
"""

from __future__ import annotations

import numpy as np

from ..mass_models.base import MassComponent

__all__ = ["KMS_PER_MASYR_KPC", "unit_vectors", "systemic_velocity_vector", "systemic_pm_field",
           "depth_dispersion", "systemic_pa"]

KMS_PER_MASYR_KPC = 4.740470463533348
ARCSEC_PER_RAD = 206264.806


def unit_vectors(ra_deg, dec_deg):
    """Line-of-sight, east and north unit vectors at each position (shape (..., 3))."""
    a = np.radians(np.asarray(ra_deg, float)); d = np.radians(np.asarray(dec_deg, float))
    n = np.stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)], axis=-1)
    east = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
    north = np.stack([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)], axis=-1)
    return n, east, north


def systemic_velocity_vector(ra0, dec0, pmra_cosdec, pmdec, v_los, distance_kpc):
    """3-D velocity (km/s, ICRS cartesian) of a body at the centre with the given PM and v_los."""
    n, e, nth = unit_vectors(ra0, dec0)
    vt = KMS_PER_MASYR_KPC * distance_kpc
    return v_los * n + pmra_cosdec * vt * e + pmdec * vt * nth


def systemic_pm_field(ra_deg, dec_deg, v_sys, distance_kpc):
    """Proper motion (mas/yr; pmra*cosdec, pmdec) that a star at each position would show if it
    moved with exactly ``v_sys`` at distance ``distance_kpc``: the systemic velocity projected
    onto the star's own tangent basis. Exact in angle (no small-angle expansion); first order
    in the depth offset, which is treated statistically by :func:`depth_dispersion`."""
    n, e, nth = unit_vectors(ra_deg, dec_deg)
    vt = KMS_PER_MASYR_KPC * distance_kpc
    v = np.asarray(v_sys, float)
    return (e @ v) / vt, (nth @ v) / vt


def systemic_pa(pmra_cosdec, pmdec):
    """Position angle (rad, from north through east) of the systemic proper motion."""
    return float(np.arctan2(pmra_cosdec, pmdec))


def depth_dispersion(tracer: MassComponent, R_pc, mu_sys_masyr: float, distance_kpc: float,
                     z_max_pc: float = 500.0, n: int = 801) -> np.ndarray:
    """Apparent PM dispersion from the unknown line-of-sight depth, ``|mu_sys| sigma_z(R) / D``,
    where ``sigma_z^2(R) = int nu z^2 dz / int nu dz`` along the line of sight through the tracer
    at projected radius ``R``. Returned in mas/yr, along the direction of the systemic PM."""
    R = np.atleast_1d(np.asarray(R_pc, float))
    z = np.linspace(-z_max_pc, z_max_pc, n)
    r = np.sqrt(R[:, None] ** 2 + z[None, :] ** 2)
    nu = np.asarray(tracer.density(r.ravel()), float).reshape(r.shape)
    sz2 = np.trapezoid(nu * z[None, :] ** 2, z, axis=1) / np.trapezoid(nu, z, axis=1)
    return abs(mu_sys_masyr) * np.sqrt(sz2) / (distance_kpc * 1e3)
