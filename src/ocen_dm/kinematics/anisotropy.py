"""Single- and two-transition spherical anisotropy profiles.

``beta(r) = beta_0 + (beta_inf - beta_0) * r^2 / (r^2 + r_beta^2)``

This is JamPy's ``logistic`` anisotropy with ``alpha = 2`` and ``r_a = r_beta``,
so both Jeans engines can use the same anisotropy parametrisation. Their mass
representations and numerical projections still differ. The
Jeans equation needs ``exp(2 * int beta(t)/t dt)``, which this form gives in
closed form.
TurnoverAnisotropy adds a second ordered transition for the adopted rung 2.
SplineAnisotropy is beta linear in ln r through fixed nodes (the JAM spline, 2026-11-16).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = ["Anisotropy", "TurnoverAnisotropy", "SplineAnisotropy", "SPLINE_NODES_ARCSEC"]


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


#: spline nodes (arcsec, 3-D radius), as the JAM fits (bin/fit_jam_halo.py --bnodes 10 100 250 700 1200 2000)
SPLINE_NODES_ARCSEC = (10.0, 100.0, 250.0, 700.0, 1200.0, 2000.0)


@dataclass(frozen=True)
class SplineAnisotropy:
    """beta linear in ln r through ``(nodes_pc, values)``, constant inside the first and beyond the last node.

    Same form as the JAM spline (bin/jam_halo_models.SplineBeta), so the two engines share one anisotropy
    parametrisation. ``beta_0`` and ``beta_inf`` are the end values (the Jeans solver's tail uses beta_inf).
    With x = ln r, ln g = 2 int beta dx is piecewise quadratic in x and is evaluated exactly.
    """

    nodes_pc: tuple
    values: tuple

    def __post_init__(self) -> None:
        r, v = np.asarray(self.nodes_pc, float), np.asarray(self.values, float)
        if r.ndim != 1 or r.size < 2 or r.shape != v.shape or np.any(np.diff(r) <= 0) or np.any(r <= 0):
            raise ValueError("nodes_pc must be positive, increasing, and match values")
        if not np.all(np.isfinite(v)) or np.any(v >= 1):
            raise ValueError("spline beta values must be finite and < 1")

    @property
    def beta_0(self) -> float:
        return float(self.values[0])

    @property
    def beta_inf(self) -> float:
        return float(self.values[-1])

    def beta(self, r: Any) -> np.ndarray:
        x = np.log(np.maximum(np.asarray(r, float), 1e-300))
        return np.interp(x, np.log(np.asarray(self.nodes_pc, float)), np.asarray(self.values, float))

    def log_integrating_factor(self, r: Any) -> np.ndarray:
        """``ln g = 2 int_{x_0}^{x} beta dx`` (x = ln r), zero at the first node; exact for the piecewise-linear beta."""
        xk, vk = np.log(np.asarray(self.nodes_pc, float)), np.asarray(self.values, float)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (vk[1:] + vk[:-1]) * np.diff(xk))])   # int beta dx at the nodes
        with np.errstate(divide="ignore"):
            x = np.log(np.asarray(r, float))
        i = np.clip(np.searchsorted(xk, x, side="right") - 1, 0, len(xk) - 1)
        b_i = vk[i]
        b_x = self.beta(np.exp(x))
        inner = x < xk[0]
        # inside a segment (and beyond the last node, where beta is constant): trapezoid from node i to x
        val = cum[i] + 0.5 * (b_i + b_x) * (x - xk[i])
        val = np.where(inner, vk[0] * (x - xk[0]), val)
        return 2.0 * val
