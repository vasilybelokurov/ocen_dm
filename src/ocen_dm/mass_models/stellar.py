"""Luminous-matter components: Plummer and a multi-Gaussian expansion.

The specification (section 3.1) asks for the light to be represented by an MGE or
another smooth basis that deprojects consistently, rather than a single hard-coded
Plummer sphere. Both are here: :class:`Plummer` for limiting cases and tests,
:class:`MGE` for the real light model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from scipy import special

from .base import G, MassComponent, as_array

__all__ = ["Plummer", "MGE"]


@dataclass(frozen=True)
class Plummer(MassComponent):
    """Plummer sphere.

    ``rho(r) = 3 M / (4 pi a^3) * (1 + r^2/a^2)^(-5/2)``

    Attributes
    ----------
    mass : float
        Total mass in Msun.
    scale : float
        Plummer scale radius ``a`` in pc.
    name : str
        Label used in reports.
    """

    mass: float
    scale: float
    name: str = "plummer"

    def __post_init__(self) -> None:
        if self.mass < 0:
            raise ValueError("Plummer mass must be non-negative")
        if self.scale <= 0:
            raise ValueError("Plummer scale radius must be positive")

    def density(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        return (
            3.0 * self.mass / (4.0 * np.pi * self.scale**3)
            * (1.0 + (radii / self.scale) ** 2) ** -2.5
        )

    def enclosed_mass(self, r: Any) -> np.ndarray:
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        out = np.where(
            np.isinf(radii),
            self.mass,
            self.mass * radii**3 / (radii**2 + self.scale**2) ** 1.5,
        )
        return np.asarray(out, dtype=float)

    def potential(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        return -G * self.mass / np.sqrt(radii**2 + self.scale**2)

    @property
    def total_mass(self) -> float:
        return float(self.mass)

    @property
    def half_mass_radius(self) -> float:
        """Radius containing half the mass, ``a / sqrt(2^(2/3) - 1)``."""
        return self.scale / np.sqrt(2.0 ** (2.0 / 3.0) - 1.0)


@dataclass(frozen=True)
class MGE(MassComponent):
    """Spherical multi-Gaussian expansion of a luminous profile.

    Each component is ``rho_j(r) = M_j / ((2 pi)^{3/2} s_j^3) exp(-r^2 / (2 s_j^2))``,
    the deprojection of a circular Gaussian in surface density, so a fit made in
    projection carries over without a separate deprojection step.

    Attributes
    ----------
    masses : sequence of float
        Component masses in Msun.
    sigmas : sequence of float
        Component dispersions (Gaussian widths) in pc.
    mass_to_light : float
        Multiplies every component mass. Kept separate so that an M/L nuisance
        parameter can be varied without rebuilding the expansion.
    """

    masses: Sequence[float]
    sigmas: Sequence[float]
    mass_to_light: float = 1.0
    name: str = "mge"

    def __post_init__(self) -> None:
        masses = np.asarray(self.masses, dtype=float)
        sigmas = np.asarray(self.sigmas, dtype=float)
        if masses.shape != sigmas.shape:
            raise ValueError("MGE masses and sigmas must have the same length")
        if masses.size == 0:
            raise ValueError("MGE needs at least one Gaussian")
        if np.any(masses < 0):
            raise ValueError("MGE component masses must be non-negative")
        if np.any(sigmas <= 0):
            raise ValueError("MGE component widths must be positive")
        if self.mass_to_light <= 0:
            raise ValueError("mass_to_light must be positive")

    @property
    def _masses(self) -> np.ndarray:
        return np.asarray(self.masses, dtype=float) * self.mass_to_light

    @property
    def _sigmas(self) -> np.ndarray:
        return np.asarray(self.sigmas, dtype=float)

    def density(self, r: Any) -> np.ndarray:
        radii = as_array(r)[:, None]
        masses, sigmas = self._masses[None, :], self._sigmas[None, :]
        terms = (
            masses / ((2.0 * np.pi) ** 1.5 * sigmas**3)
            * np.exp(-0.5 * (radii / sigmas) ** 2)
        )
        return terms.sum(axis=1)

    def enclosed_mass(self, r: Any) -> np.ndarray:
        radii = np.atleast_1d(np.asarray(r, dtype=float))[:, None]
        masses, sigmas = self._masses[None, :], self._sigmas[None, :]
        x = np.divide(radii, sigmas, out=np.full_like(radii * sigmas, np.inf),
                      where=np.isfinite(radii))
        # M_j(<r) = M_j [ erf(x/sqrt2) - sqrt(2/pi) x exp(-x^2/2) ]
        frac = special.erf(x / np.sqrt(2.0)) - np.sqrt(2.0 / np.pi) * x * np.exp(-0.5 * x**2)
        frac = np.where(np.isinf(x), 1.0, frac)
        return (masses * frac).sum(axis=1)

    def potential(self, r: Any) -> np.ndarray:
        radii = as_array(r)[:, None]
        masses, sigmas = self._masses[None, :], self._sigmas[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            phi = -G * masses * special.erf(radii / (np.sqrt(2.0) * sigmas)) / radii
            # r -> 0 limit of erf(r/(sqrt2 s))/r is sqrt(2/pi)/s
            phi = np.where(radii > 0, phi, -G * masses * np.sqrt(2.0 / np.pi) / sigmas)
        return phi.sum(axis=1)

    @property
    def total_mass(self) -> float:
        return float(self._masses.sum())
