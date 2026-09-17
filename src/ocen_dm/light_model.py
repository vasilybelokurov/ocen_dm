"""Luminous model of omega Cen: a multi-Gaussian expansion of the light profile.

Specification section 3.1 asks for the luminous component to come from a public
profile rather than a hard-coded Plummer sphere. v1 uses the Trager, King &
Djorgovski (1995) V-band surface-brightness profile (73 points, 10.5-2590 arcsec).

Method
------
A spherical MGE in projection is ``Sigma(R) = sum_j L_j / (2 pi s_j^2) exp(-R^2 / 2 s_j^2)``.
With the widths ``s_j`` fixed on a logarithmic grid the model is linear in the
amplitudes ``L_j``, so the fit is a non-negative least-squares problem and needs
no external MGE package. Deprojection is analytic and keeps the same widths:
``rho_j(r) = L_j / ((2 pi)^{3/2} s_j^3) exp(-r^2 / 2 s_j^2)``, which is exactly
:class:`ocen_dm.mass_models.MGE`.

Only the *shape* of the profile is used. The stellar mass-to-light ratio is a free
parameter of the mass model, so the photometric zero-point and the extinction
correction cancel; the amplitudes are normalised to unit total and scaled by the
total stellar mass at model-construction time.

Weighting: the fit minimises *fractional* residuals (equivalently, residuals in
magnitudes), scaled by the authors' per-point weights, so the 12-magnitude dynamic
range does not let the bright centre dominate the outskirts.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table
from scipy import optimize

from .mass_models import MGE
from .paths import processed_dir, raw_dir

__all__ = [
    "SurfaceBrightnessProfile",
    "MGEFit",
    "load_trager_profile",
    "fit_mge_projected",
    "projected_half_light_radius",
    "build_stellar_mge",
    "arcsec_to_pc",
    "star_count_profile", "composite_profile", "load_tracer_profile",
]

ARCSEC_PER_RAD = 206264.806


def arcsec_to_pc(arcsec: Any, distance_kpc: float) -> np.ndarray:
    """Convert an angle in arcsec to a length in pc at ``distance_kpc``."""
    return np.asarray(arcsec, dtype=float) * distance_kpc * 1e3 / ARCSEC_PER_RAD


@dataclass(frozen=True)
class SurfaceBrightnessProfile:
    """A measured surface-brightness profile.

    Attributes
    ----------
    r_arcsec : numpy.ndarray
        Radii in arcsec.
    mu : numpy.ndarray
        Surface brightness in mag/arcsec^2.
    weight : numpy.ndarray
        Relative weights in [0, 1].
    source : str
        Where the points came from.
    """

    r_arcsec: np.ndarray
    mu: np.ndarray
    weight: np.ndarray
    source: str

    @property
    def intensity(self) -> np.ndarray:
        """Relative surface intensity ``10**(-0.4 (mu - min mu))``; shape only."""
        return 10.0 ** (-0.4 * (self.mu - np.min(self.mu)))


def load_trager_profile(path: Path | None = None) -> SurfaceBrightnessProfile:
    """Load the processed Trager et al. (1995) omega Cen profile."""
    path = path or processed_dir() / "literature" / "trager1995_ocen_sbp.ecsv"
    table = Table.read(path)
    order = np.argsort(np.asarray(table["r"]))
    return SurfaceBrightnessProfile(
        r_arcsec=np.asarray(table["r"], dtype=float)[order],
        mu=np.asarray(table["mu_v"], dtype=float)[order],
        weight=np.asarray(table["weight"], dtype=float)[order],
        source="Trager, King & Djorgovski 1995, V band (VizieR J/AJ/109/218)",
    )


def star_count_profile(mag_cut: float = 19.0, r_min_arcsec: float = 1.0, r_max_arcsec: float = 120.0,
                       n_bins: int = 22, path: Path | None = None) -> SurfaceBrightnessProfile:
    """Number-density profile of oMEGACat stars, as a pseudo surface-brightness profile.

    The kinematic tracers are stars, not light: a surface-brightness profile is
    weighted by a handful of bright giants in the core, whereas the Jeans
    equation needs the *number* density of the tracer population. Inside ~20
    arcsec the HST counts lie 0.2-0.3 mag below the Trager MGE for every
    magnitude cut from F625W < 18 to < 20 (JOURNAL 2026-09-17), so this is not
    incompleteness.

    Uses the authors' ``selection_hq_f625w`` flag, the pixel-based centre
    (15000, 15000) at 0.04 arcsec/pixel (their notebook), Poisson errors. The
    output ``mu`` is ``-2.5 log10(counts / area)`` up to an arbitrary zero-point;
    ``weight`` is the Poisson weight relative to a bin of 100 stars, capped at 1.
    """
    path = path or raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits"
    t = Table.read(path)
    hq = np.asarray(t["selection_hq_f625w"]).astype(bool) & np.isfinite(np.asarray(t["f625w"], float))
    m = hq & (np.asarray(t["f625w"], float) < mag_cut)
    x = -0.04 * (np.asarray(t["x"], float)[m] - 15000.0)
    y = 0.04 * (np.asarray(t["y"], float)[m] - 15000.0)
    edges = np.geomspace(r_min_arcsec, r_max_arcsec, n_bins + 1)
    n, _ = np.histogram(np.hypot(x, y), edges)
    keep = n >= 5
    mid = np.sqrt(edges[1:] * edges[:-1])[keep]
    density = n[keep] / (np.pi * np.diff(edges**2))[keep]
    return SurfaceBrightnessProfile(
        r_arcsec=mid, mu=-2.5 * np.log10(density), weight=np.minimum(n[keep] / 100.0, 1.0),
        source=f"oMEGACat star counts, selection_hq_f625w, F625W < {mag_cut:g} (n={m.sum():,})",
    )


def composite_profile(outer: SurfaceBrightnessProfile, inner: SurfaceBrightnessProfile,
                      r_switch_arcsec: float = 25.0, anchor_arcsec: tuple[float, float] = (30.0, 100.0),
                      min_outer_weight: float = 0.5) -> SurfaceBrightnessProfile:
    """Splice a star-count profile inside ``r_switch`` onto a light profile outside.

    The inner profile's arbitrary zero-point is set by the weighted mean magnitude
    offset from the outer profile over ``anchor_arcsec``, where both are trusted
    (outside the giant-dominated core, inside the field of view). Only outer
    points with ``weight >= min_outer_weight`` enter the anchor: Trager's
    low-weight points sit 0.5-0.9 mag off the curve and, interpolated through,
    biased the zero-point by 0.2-0.3 mag (found 2026-09-17). Both profiles must
    have points in the anchor range.
    """
    lo, hi = anchor_arcsec
    a_in = (inner.r_arcsec >= lo) & (inner.r_arcsec <= hi)
    a_out = (outer.r_arcsec >= lo) & (outer.r_arcsec <= hi) & (outer.weight >= min_outer_weight)
    if a_in.sum() < 2 or a_out.sum() < 2:
        raise ValueError("both profiles need at least two points in the anchor range")
    # interpolate the outer profile (in mag) at the inner anchor radii
    mu_out_at_inner = np.interp(np.log(inner.r_arcsec[a_in]), np.log(outer.r_arcsec[a_out]), outer.mu[a_out])
    offset = np.average(mu_out_at_inner - inner.mu[a_in], weights=inner.weight[a_in])
    use_in = inner.r_arcsec < r_switch_arcsec
    use_out = outer.r_arcsec >= r_switch_arcsec
    r = np.concatenate([inner.r_arcsec[use_in], outer.r_arcsec[use_out]])
    mu = np.concatenate([inner.mu[use_in] + offset, outer.mu[use_out]])
    w = np.concatenate([inner.weight[use_in], outer.weight[use_out]])
    order = np.argsort(r)
    return SurfaceBrightnessProfile(
        r_arcsec=r[order], mu=mu[order], weight=w[order],
        source=f"{inner.source} inside {r_switch_arcsec:g} arcsec (zero-point offset {offset:+.3f} mag over "
               f"{lo:g}-{hi:g} arcsec) + {outer.source} outside",
    )


def load_tracer_profile(kind: str = "composite") -> SurfaceBrightnessProfile:
    """``'trager'``: the V-band light; ``'composite'``: HST star counts inside 25 arcsec + Trager outside."""
    trager = load_trager_profile()
    if kind == "trager":
        return trager
    if kind == "composite":
        return composite_profile(trager, star_count_profile())
    raise ValueError(f"unknown tracer profile kind {kind!r}")


@dataclass(frozen=True)
class MGEFit:
    """Result of a projected MGE fit.

    Attributes
    ----------
    sigmas_arcsec : numpy.ndarray
        Widths of the retained Gaussians, arcsec.
    fractions : numpy.ndarray
        Light fractions of the retained Gaussians; they sum to 1.
    rms_mag : float
        Weighted rms of the residuals in magnitudes.
    max_abs_resid_mag : float
        Largest absolute residual in magnitudes.
    n_grid : int
        Size of the width grid the fit started from.
    """

    sigmas_arcsec: np.ndarray
    fractions: np.ndarray
    rms_mag: float
    max_abs_resid_mag: float
    n_grid: int

    def surface_intensity(self, r_arcsec: Any) -> np.ndarray:
        """Reprojected relative surface intensity at ``r_arcsec``."""
        r = np.atleast_1d(np.asarray(r_arcsec, dtype=float))[:, None]
        s = self.sigmas_arcsec[None, :]
        return (self.fractions[None, :] / (2 * np.pi * s**2) * np.exp(-0.5 * (r / s) ** 2)).sum(axis=1)

    def surface_brightness(self, r_arcsec: Any, mu_ref: float) -> np.ndarray:
        """Reprojected surface brightness, given the zero-point that maps unit
        total light onto the data's magnitude scale."""
        return mu_ref - 2.5 * np.log10(self.surface_intensity(r_arcsec))

    @property
    def components(self) -> list[tuple[float, float]]:
        """``(sigma_arcsec, fraction)`` pairs, brightest-in-centre first."""
        return sorted(zip(self.sigmas_arcsec.tolist(), self.fractions.tolist()))


def _basis(r: np.ndarray, sigmas: np.ndarray) -> np.ndarray:
    """Projected unit-light Gaussians evaluated at ``r``: shape (n_r, n_sigma)."""
    rr, ss = r[:, None], sigmas[None, :]
    return np.exp(-0.5 * (rr / ss) ** 2) / (2 * np.pi * ss**2)


def fit_mge_projected(
    profile: SurfaceBrightnessProfile,
    n_grid: int = 64,
    sigma_range_arcsec: tuple[float, float] | None = None,
    min_fraction: float = 1e-4,
) -> MGEFit:
    """Fit a projected MGE to a surface-brightness profile by weighted NNLS.

    Parameters
    ----------
    profile : SurfaceBrightnessProfile
        The data.
    n_grid : int, optional
        Number of log-spaced trial widths. NNLS keeps only a handful, but the grid
        must be dense enough that some node sits close to each true width:
        measured on noiseless synthetic profiles, a single Gaussian is recovered
        to 0.36 / 0.10 / 0.027 / 0.007 mag rms with 16 / 32 / 64 / 128 nodes,
        while the real Trager data plateau at 0.184 mag from 32 nodes on.
    sigma_range_arcsec : tuple of float, optional
        Width grid limits; default half the innermost to twice the outermost radius.
    min_fraction : float, optional
        Components carrying less than this fraction of the light are dropped.

    Returns
    -------
    MGEFit
    """
    r, I, w = profile.r_arcsec, profile.intensity, profile.weight
    if sigma_range_arcsec is None:
        sigma_range_arcsec = (0.5 * r.min(), 2.0 * r.max())
    sigmas = np.geomspace(*sigma_range_arcsec, n_grid)

    # Fractional residuals: divide each row by the data value, so the fit is in
    # magnitudes to first order; the authors' weights then scale the rows.
    scale = (np.sqrt(np.clip(w, 1e-3, 1.0)) / I)[:, None]
    A = _basis(r, sigmas) * scale
    b = I * scale[:, 0]
    amplitudes, _ = optimize.nnls(A, b, maxiter=50 * n_grid)

    total = amplitudes.sum()
    if total <= 0:
        raise RuntimeError("MGE fit returned no light; check the profile")
    keep = amplitudes / total >= min_fraction
    sigmas, amplitudes = sigmas[keep], amplitudes[keep]
    fractions = amplitudes / amplitudes.sum()

    fit = MGEFit(sigmas, fractions, 0.0, 0.0, n_grid)
    model = fit.surface_intensity(r) * amplitudes.sum()
    resid_mag = -2.5 * np.log10(model / I)
    rms = float(np.sqrt(np.average(resid_mag**2, weights=w)))
    return MGEFit(sigmas, fractions, rms, float(np.max(np.abs(resid_mag))), n_grid)


def projected_half_light_radius(fit: MGEFit) -> float:
    """Radius enclosing half the projected light, in arcsec.

    For one Gaussian this is ``sigma * sqrt(2 ln 2) = 1.1774 sigma``; for a sum it
    is solved numerically from ``L(<R) = sum_j f_j [1 - exp(-R^2 / 2 s_j^2)]``.
    """
    def enclosed(R: float) -> float:
        return float(np.sum(fit.fractions * (1 - np.exp(-0.5 * (R / fit.sigmas_arcsec) ** 2))))

    lo, hi = fit.sigmas_arcsec.min() * 1e-3, fit.sigmas_arcsec.max() * 10
    return float(optimize.brentq(lambda R: enclosed(R) - 0.5, lo, hi))


def build_stellar_mge(fit: MGEFit, distance_kpc: float, total_mass: float = 1.0,
                      name: str = "stars") -> MGE:
    """Turn a projected fit into the spherical :class:`MGE` mass component.

    Parameters
    ----------
    fit : MGEFit
        Projected fit; widths are shared between projection and deprojection.
    distance_kpc : float
        Converts arcsec to pc.
    total_mass : float, optional
        Total stellar mass in Msun, i.e. the M/L nuisance times the total light.
    """
    return MGE(
        masses=(fit.fractions * total_mass).tolist(),
        sigmas=arcsec_to_pc(fit.sigmas_arcsec, distance_kpc).tolist(),
        name=name,
    )
