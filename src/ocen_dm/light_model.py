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
from .paths import processed_dir

__all__ = [
    "SurfaceBrightnessProfile",
    "MGEFit",
    "load_trager_profile",
    "fit_mge_projected",
    "projected_half_light_radius",
    "build_stellar_mge",
    "arcsec_to_pc",
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
