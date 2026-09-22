"""Single- and two-transition spherical anisotropy profiles.

``beta(r) = beta_0 + (beta_inf - beta_0) * r^2 / (r^2 + r_beta^2)``

This is JamPy's ``logistic`` anisotropy with ``alpha = 2`` and ``r_a = r_beta``,
so both Jeans engines can use the same anisotropy parametrisation. Their mass
representations and numerical projections still differ. The
Jeans equation needs ``exp(2 * int beta(t)/t dt)``, which this form gives in
closed form.
TurnoverAnisotropy adds a second ordered transition for the adopted rung 2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = ["Anisotropy", "TurnoverAnisotropy"]


@dataclass(frozen=True)
class Anisotropy:
    """Radial anisotropy ``beta = 1 - sigma_t^2 / sigma_r^2`` as a function of radius.

    Attributes
    ----------
    beta_0 : float
        Central anisotropy.
    beta_inf : float
        Anisotropy at large radius. Must be < 1 (beta = 1 is purely radial).
    r_beta : float
        Transition radius in pc.
    """

    beta_0: float = 0.0
    beta_inf: float = 0.0
    r_beta: float = 10.0

    def __post_init__(self) -> None:
        for name in ("beta_0", "beta_inf"):
            value = getattr(self, name)
            if not (-np.inf < value < 1.0):
                raise ValueError(f"{name} must be < 1 (got {value}); beta = 1 is purely radial")
        if self.r_beta <= 0:
            raise ValueError("r_beta must be positive")

    def beta(self, r: Any) -> np.ndarray:
        r = np.asarray(r, dtype=float)
        return self.beta_0 + (self.beta_inf - self.beta_0) * r**2 / (r**2 + self.r_beta**2)

    def log_integrating_factor(self, r: Any) -> np.ndarray:
        """``ln g(r)`` with ``g(r) = exp(2 int^r beta(t)/t dt)``, up to a constant.

        For this profile ``g(r) = r^(2 beta_0) (r^2 + r_beta^2)^(beta_inf - beta_0)``.
        Returned in the log so the Jeans integral can be done with differences of
        logs and never overflows for large radii or radial anisotropy.
        """
        r = np.asarray(r, dtype=float)
        with np.errstate(divide="ignore"):
            return 2.0 * self.beta_0 * np.log(r) + (self.beta_inf - self.beta_0) * np.log(r**2 + self.r_beta**2)

    @property
    def jampy_logistic(self) -> list[float]:
        """The equivalent JamPy ``beta`` array for ``logistic=True``: ``[r_a, beta_0, beta_inf, alpha=2]``.

        ``r_a`` must be handed to JamPy in the same units as its radii (arcsec), so
        callers convert ``r_beta`` from pc themselves.
        """
        return [self.r_beta, self.beta_0, self.beta_inf, 2.0]


@dataclass(frozen=True)
class TurnoverAnisotropy:
    """Two ordered smooth transitions, allowing a peak, trough or monotonic beta.

    With f_i=r²/(r²+r_i²), beta = beta_0 (1-f_1) + beta_mid (f_1-f_2)
    + beta_inf f_2. Since r_2>r_1, these are non-negative weights summing
    to one. The profile cannot overshoot its three control levels, so
    levels below one guarantee beta(r)<1. beta_mid is a control level,
    not necessarily the attained maximum or a flat intermediate plateau.

    This moment parametrisation does not certify a non-negative DF.
    Radii are in pc; delta_r_beta is the positive gap r_2-r_1.
    """

    beta_0: float = -0.75
    beta_mid: float = 0.25
    beta_inf: float = 0.0
    r_beta: float = 1.0
    delta_r_beta: float = 20.0

    def __post_init__(self) -> None:
        for name in ("beta_0", "beta_mid", "beta_inf"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) >= 1:
                raise ValueError(f"{name} must be finite and < 1")
        for name in ("r_beta", "delta_r_beta"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")

    @property
    def r_beta_outer(self) -> float:
        return self.r_beta + self.delta_r_beta

    def beta(self, r: Any) -> np.ndarray:
        r = np.asarray(r, dtype=float)
        f1 = r*r/(r*r+self.r_beta**2)
        f2 = r*r/(r*r+self.r_beta_outer**2)
        return self.beta_0*(1-f1)+self.beta_mid*(f1-f2)+self.beta_inf*f2

    def log_integrating_factor(self, r: Any) -> np.ndarray:
        r = np.asarray(r, dtype=float)
        with np.errstate(divide="ignore"):
            return (2*self.beta_0*np.log(r)
                    +(self.beta_mid-self.beta_0)*np.log(r*r+self.r_beta**2)
                    +(self.beta_inf-self.beta_mid)*np.log(r*r+self.r_beta_outer**2))
