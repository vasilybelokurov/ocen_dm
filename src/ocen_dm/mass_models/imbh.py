"""Intermediate-mass black hole: a point mass."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .base import G, MassComponent, as_array

__all__ = ["PointMass"]


@dataclass(frozen=True)
class PointMass(MassComponent):
    """Point mass at the centre.

    The density is a delta function, so :meth:`density` returns zero away from
    the origin and ``inf`` at it; the enclosed mass and potential are the
    meaningful quantities. ``verify()`` therefore skips the density-integral
    check for this component (see :class:`CompositeMassModel`).

    Attributes
    ----------
    mass : float
        Black-hole mass in Msun, ``>= 0``. Zero is allowed and is the no-IMBH case.
    """

    mass: float
    name: str = "imbh"

    def __post_init__(self) -> None:
        if self.mass < 0:
            raise ValueError("IMBH mass must be non-negative (0 means no IMBH)")

    def density(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        return np.where(radii > 0, 0.0, np.inf if self.mass > 0 else 0.0)

    def enclosed_mass(self, r: Any) -> np.ndarray:
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        return np.where(radii > 0, self.mass, 0.0)

    def potential(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        with np.errstate(divide="ignore"):
            return np.where(radii > 0, -G * self.mass / radii, -np.inf)

    @property
    def total_mass(self) -> float:
        return float(self.mass)
