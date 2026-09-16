"""Dark-matter components: a truncated generalized NFW and a Burkert core.

Specification section 3.1 asks for two families run separately -- a cusp and a
core -- with a smooth outer truncation, and warns against giving the first
inference more halo-shape freedom than the data support. So ``gamma`` is a
parameter with the two standard values as named constructors (``nfw``,
``cored``), and promoting it to a free parameter is a deliberate act.

Truncation
----------
``T(r) = (1 + (r/r_t)^2)^(-2)`` multiplies the density. It is smooth everywhere,
keeps the total mass finite, and falls as ``r^-4`` outside ``r_t`` in the manner
of a tidally truncated halo. The truncation makes the enclosed mass non-analytic,
so it is integrated numerically; the untruncated limit keeps its closed form and
is used to test the quadrature.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import integrate

from .base import G, MassComponent, as_array

__all__ = ["TruncatedGNFW", "Burkert"]


@dataclass(frozen=True)
class TruncatedGNFW(MassComponent):
    """Truncated generalized NFW halo.

    ``rho(r) = rho_s (r/r_s)^-gamma (1 + r/r_s)^(gamma-3) T(r; r_t)``

    Attributes
    ----------
    rho_s : float
        Scale density in Msun / pc^3.
    r_s : float
        Scale radius in pc.
    gamma : float
        Inner logarithmic slope: 1 is NFW, 0 is a core. Must lie in [0, 3).
    r_t : float
        Truncation radius in pc; ``inf`` disables the truncation.
    """

    rho_s: float
    r_s: float
    gamma: float = 1.0
    r_t: float = np.inf
    name: str = "dm_gnfw"

    def __post_init__(self) -> None:
        if self.rho_s < 0:
            raise ValueError("rho_s must be non-negative")
        if self.r_s <= 0:
            raise ValueError("r_s must be positive")
        if not 0.0 <= self.gamma < 3.0:
            raise ValueError("gamma must satisfy 0 <= gamma < 3")
        if self.r_t <= 0:
            raise ValueError("r_t must be positive (use inf for no truncation)")

    @classmethod
    def nfw(cls, rho_s: float, r_s: float, r_t: float = np.inf) -> "TruncatedGNFW":
        """Cuspy family, ``gamma = 1``."""
        return cls(rho_s=rho_s, r_s=r_s, gamma=1.0, r_t=r_t, name="dm_nfw")

    @classmethod
    def cored(cls, rho_s: float, r_s: float, r_t: float = np.inf) -> "TruncatedGNFW":
        """Cored family, ``gamma = 0``."""
        return cls(rho_s=rho_s, r_s=r_s, gamma=0.0, r_t=r_t, name="dm_cored_gnfw")

    @property
    def truncated(self) -> bool:
        """Whether an outer truncation is active."""
        return np.isfinite(self.r_t)

    def truncation(self, r: Any) -> np.ndarray:
        """Smooth truncation factor ``T(r)``, 1 when no truncation is set."""
        radii = as_array(r)
        if not self.truncated:
            return np.ones_like(radii)
        return (1.0 + (radii / self.r_t) ** 2) ** -2

    def density(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        x = radii / self.r_s
        with np.errstate(divide="ignore"):
            shape = x ** (-self.gamma) * (1.0 + x) ** (self.gamma - 3.0)
        return self.rho_s * shape * self.truncation(radii)

    def _analytic_nfw_mass(self, r: np.ndarray) -> np.ndarray:
        """Closed form for the untruncated ``gamma = 1`` case."""
        x = r / self.r_s
        return 4.0 * np.pi * self.rho_s * self.r_s**3 * (np.log1p(x) - x / (1.0 + x))

    def enclosed_mass(self, r: Any) -> np.ndarray:
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        if not self.truncated and np.isclose(self.gamma, 1.0):
            finite = np.isfinite(radii)
            out = np.full_like(radii, np.inf)
            out[finite] = self._analytic_nfw_mass(radii[finite])
            return out
        return MassComponent.enclosed_mass(self, radii)

    @property
    def total_mass(self) -> float:
        if not self.truncated:
            return float("inf")  # gNFW mass diverges logarithmically or worse
        return self.tables.total_mass


@dataclass(frozen=True)
class Burkert(MassComponent):
    """Burkert cored halo, optionally truncated.

    ``rho(r) = rho_0 r_0^3 / [ (r + r_0)(r^2 + r_0^2) ] * T(r; r_t)``

    Attributes
    ----------
    rho_0 : float
        Central density in Msun / pc^3.
    r_0 : float
        Core radius in pc.
    r_t : float
        Truncation radius in pc; ``inf`` disables it.
    """

    rho_0: float
    r_0: float
    r_t: float = np.inf
    name: str = "dm_burkert"

    def __post_init__(self) -> None:
        if self.rho_0 < 0:
            raise ValueError("rho_0 must be non-negative")
        if self.r_0 <= 0:
            raise ValueError("r_0 must be positive")
        if self.r_t <= 0:
            raise ValueError("r_t must be positive (use inf for no truncation)")

    @property
    def truncated(self) -> bool:
        """Whether an outer truncation is active."""
        return np.isfinite(self.r_t)

    def truncation(self, r: Any) -> np.ndarray:
        """Smooth truncation factor ``T(r)``."""
        radii = as_array(r)
        if not self.truncated:
            return np.ones_like(radii)
        return (1.0 + (radii / self.r_t) ** 2) ** -2

    def density(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        shape = (
            self.rho_0 * self.r_0**3
            / ((radii + self.r_0) * (radii**2 + self.r_0**2))
        )
        return shape * self.truncation(radii)

    def _analytic_mass(self, r: np.ndarray) -> np.ndarray:
        """Closed form for the untruncated Burkert profile."""
        x = r / self.r_0
        return (
            2.0 * np.pi * self.rho_0 * self.r_0**3
            * (np.log1p(x) + 0.5 * np.log1p(x**2) - np.arctan(x))
        )

    def enclosed_mass(self, r: Any) -> np.ndarray:
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        if not self.truncated:
            finite = np.isfinite(radii)
            out = np.full_like(radii, np.inf)
            out[finite] = self._analytic_mass(radii[finite])
            return out
        return MassComponent.enclosed_mass(self, radii)

    @property
    def total_mass(self) -> float:
        if not self.truncated:
            return float("inf")  # Burkert mass diverges logarithmically
        return self.tables.total_mass
