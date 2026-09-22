"""Spherical Jeans solver for a tracer population in a composite potential.

Given a tracer density ``nu(r)`` (the deprojected light MGE), a total mass profile
``M(<r)`` (the composite model) and an anisotropy profile ``beta(r)``, the second
moment ``sigma_r^2`` follows from the spherical Jeans equation

``d(nu sigma_r^2)/dr + 2 beta nu sigma_r^2 / r = -nu G M(<r) / r^2``,

whose solution with the integrating factor ``g(r) = exp(2 int beta/t dt)`` is

``nu sigma_r^2 (r) = (1/g(r)) int_r^inf g(s) nu(s) G M(s) / s^2 ds``.

The projected second moments then follow (Binney & Mamon 1982; van der Marel 1994):

``Sigma sigma_los^2 (R) = 2 int_R^inf [1 - beta R^2/r^2]        nu sigma_r^2 r dr / sqrt(r^2 - R^2)``
``Sigma sigma_pmR^2 (R) = 2 int_R^inf [1 - beta + beta R^2/r^2] nu sigma_r^2 r dr / sqrt(r^2 - R^2)``
``Sigma sigma_pmT^2 (R) = 2 int_R^inf [1 - beta]                nu sigma_r^2 r dr / sqrt(r^2 - R^2)``

with ``Sigma(R) = 2 int_R^inf nu r dr / sqrt(r^2 - R^2)``. The substitution
``r = R cosh u`` removes the integrable singularity at ``r = R``: ``r dr / sqrt(r^2 - R^2)
= r du``, and every projection becomes a smooth integral over ``u`` that fixed
Gauss-Legendre nodes handle to high accuracy, vectorised over all requested ``R``.

Numerics: ``nu sigma_r^2`` is tabulated once per model on a 600-point log grid by
cumulative Simpson integration from the outside in (the same approach as the
mass tables) and splined in the log; a full set of projected profiles for ~70
radii costs a few milliseconds. Units are the project's: pc, Msun, km/s.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Any

import numpy as np
from scipy import integrate, interpolate

from ..mass_models.base import G, MassComponent
from .anisotropy import Anisotropy, TurnoverAnisotropy

__all__ = ["SphericalJeans", "KMS_PER_MASYR_KPC", "PROJECTIONS"]

#: v [km/s] = KMS_PER_MASYR_KPC * mu [mas/yr] * D [kpc]
KMS_PER_MASYR_KPC = 4.740470463533348
#: recognised projected components
PROJECTIONS = ("los", "pmr", "pmt")

_GRID_MIN, _GRID_MAX, _GRID_N = 1e-4, 1e5, 600
_U_NODES = 96   # Gauss-Legendre nodes in u for the projections


@dataclass(frozen=True)
class SphericalJeans:
    """Second moments of a spherical tracer in a spherical potential.

    Attributes
    ----------
    mass : MassComponent
        The total mass model (typically a :class:`CompositeMassModel`).
    tracer : MassComponent
        Whose ``density`` is the tracer number density ``nu(r)`` (normalisation
        irrelevant); usually the stellar MGE.
    anisotropy : Anisotropy
        The ``beta(r)`` profile.
    r_max : float
        Outer edge of the tracer in pc; ``nu`` is taken as zero beyond it. The
        Jeans integral formally runs to infinity, so this must be far enough out
        that ``nu`` is negligible there.
    """

    mass: MassComponent
    tracer: MassComponent
    anisotropy: Anisotropy | TurnoverAnisotropy = Anisotropy()
    r_max: float = 1e4

    # ------------------------------------------------------------ intrinsic ---
    @cached_property
    def _table(self) -> tuple[np.ndarray, interpolate.CubicSpline, interpolate.CubicSpline]:
        """Log grid, spline of ln(nu sigma_r^2), spline of ln(nu)."""
        u = np.linspace(np.log(_GRID_MIN), np.log(min(self.r_max, _GRID_MAX)), _GRID_N)
        r = np.exp(u)
        nu = np.asarray(self.tracer.density(r), dtype=float)
        # Truncate where the tracer has fallen 20 decades below its peak: a
        # Gaussian tracer underflows to zero far inside r_max, and a spline through
        # ln(0) develops kinks that the projection nodes then sample as garbage.
        alive = nu > 1e-20 * np.nanmax(nu)
        last = int(np.max(np.nonzero(alive)[0]))
        u, r, nu = u[: last + 1], r[: last + 1], nu[: last + 1]
        m = np.asarray(self.mass.enclosed_mass(r), dtype=float)
        ln_g = self.anisotropy.log_integrating_factor(r)

        # Integrand of the Jeans solution in ln s: g nu G M / s^2 * s. Work with the
        # log of g relative to its value at the current point to avoid overflow.
        ln_g_ref = ln_g[-1]
        integrand = np.exp(ln_g - ln_g_ref) * nu * G * m / r          # (g nu G M / s^2) * s
        # int_r^{r_max}, accumulated from the OUTSIDE inward. Subtracting a
        # cumulative integral from its total loses everything once the outer part
        # is below 1e-16 of the whole (the integrand spans 20 decades), which made
        # sigma_r collapse to zero at ~3 sigma_MGE, well inside r_cut (2026-09-18).
        outer = integrate.cumulative_simpson(integrand[::-1], x=-u[::-1], initial=0.0)[::-1]
        outer = np.maximum(outer, 0.0)
        # nu sigma_r^2 = exp(ln_g_ref - ln_g) * outer ; add the analytic tail beyond
        # r_max assuming nu ~ r^-s and M const there (tiny when r_max is generous).
        with np.errstate(divide="ignore", invalid="ignore"):
            slope = -(np.log(nu[-1]) - np.log(nu[-2])) / (u[-1] - u[-2])
        tail_power = slope + 1.0 - 2.0 * self.anisotropy.beta_inf     # exponent of s in g nu M/s^2 * s
        tail = (integrand[-1] / max(tail_power, 1e-3)) if slope > 1 else 0.0
        nu_sig2 = np.exp(ln_g_ref - ln_g) * (outer + tail)
        nu_sig2 = np.maximum(nu_sig2, 1e-300)
        ln_ns = interpolate.CubicSpline(u, np.log(nu_sig2), extrapolate=False)
        ln_nu = interpolate.CubicSpline(u, np.log(nu), extrapolate=False)
        object.__setattr__(self, "_r_cut", float(r[-1]))     # effective outer edge of the tracer
        return r, ln_ns, ln_nu

    @property
    def r_cut(self) -> float:
        """Outer edge actually used: where the tracer has fallen 20 decades below its peak."""
        self._table
        return self._r_cut

    def nu_sigma_r2(self, r: Any) -> np.ndarray:
        """``nu(r) sigma_r^2(r)`` in (tracer units) km^2/s^2."""
        _, ln_ns, _ = self._table
        rr = np.atleast_1d(np.asarray(r, dtype=float))
        with np.errstate(invalid="ignore"):
            out = np.exp(ln_ns(np.log(rr)))
        out[~np.isfinite(out) | (rr >= self.r_cut)] = 0.0
        return out

    def sigma_r(self, r: Any) -> np.ndarray:
        """Intrinsic radial velocity dispersion in km/s."""
        rr = np.atleast_1d(np.asarray(r, dtype=float))
        nu = np.asarray(self.tracer.density(rr), dtype=float)
        return np.sqrt(self.nu_sigma_r2(rr) / np.maximum(nu, 1e-300))

    def sigma_t(self, r: Any) -> np.ndarray:
        """Intrinsic tangential dispersion (one component) ``sqrt(1 - beta) sigma_r``."""
        rr = np.atleast_1d(np.asarray(r, dtype=float))
        return np.sqrt(np.clip(1.0 - self.anisotropy.beta(rr), 0.0, None)) * self.sigma_r(rr)

    # ------------------------------------------------------------ projected ---
    def _projection_nodes(self, R: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Radii ``r = R cosh u`` and weights for every R: shapes (nR, nU)."""
        x, w = np.polynomial.legendre.leggauss(_U_NODES)
        u_max = np.arccosh(np.maximum(self.r_cut / R, 1.0 + 1e-12))
        u = 0.5 * u_max[:, None] * (x[None, :] + 1.0)                 # map [-1,1] -> [0, u_max]
        wu = 0.5 * u_max[:, None] * w[None, :]
        r = R[:, None] * np.cosh(u)
        return r, wu, u

    def surface_density(self, R: Any) -> np.ndarray:
        """Projected tracer density ``Sigma(R)`` (tracer units / pc^2)."""
        RR = np.atleast_1d(np.asarray(R, dtype=float))
        r, wu, _ = self._projection_nodes(RR)
        _, _, ln_nu = self._table
        with np.errstate(invalid="ignore"):
            nu = np.exp(ln_nu(np.log(r)))
        nu[~np.isfinite(nu) | (r >= self.r_cut)] = 0.0
        return 2.0 * np.sum(nu * r * wu, axis=1)

    def projected_second_moment(self, R: Any, kind: str) -> np.ndarray:
        """``Sigma(R) * sigma_kind^2(R)`` for ``kind`` in ``{'los', 'pmr', 'pmt'}``."""
        if kind not in PROJECTIONS:
            raise ValueError(f"kind must be one of {PROJECTIONS}, got {kind!r}")
        RR = np.atleast_1d(np.asarray(R, dtype=float))
        r, wu, _ = self._projection_nodes(RR)
        ns = self.nu_sigma_r2(r.ravel()).reshape(r.shape)
        beta = self.anisotropy.beta(r)
        ratio2 = (RR[:, None] / r) ** 2
        if kind == "los":
            weight = 1.0 - beta * ratio2
        elif kind == "pmr":
            weight = 1.0 - beta + beta * ratio2
        else:
            weight = 1.0 - beta
        return 2.0 * np.sum(weight * ns * r * wu, axis=1)

    def _check_inside_tracer(self, RR: np.ndarray) -> None:
        """Refuse radii where the tracer has no stars: 0/0 there is not a dispersion."""
        if np.any(RR >= self.r_cut):
            raise ValueError(
                f"projected radius {RR.max():.3g} pc is at or beyond the tracer's edge "
                f"({self.r_cut:.3g} pc): no stars there, no dispersion to predict"
            )

    def sigma_projected(self, R: Any, kind: str) -> np.ndarray:
        """Projected velocity dispersion in km/s along ``kind``."""
        RR = np.atleast_1d(np.asarray(R, dtype=float))
        self._check_inside_tracer(RR)
        return np.sqrt(self.projected_second_moment(RR, kind) / self.surface_density(RR))

    def projected_moments(self, R: Any) -> dict[str, np.ndarray]:
        """``Sigma`` and the three second moments ``Sigma sigma_kind^2`` in one pass.

        Same integrals as :meth:`surface_density` and
        :meth:`projected_second_moment`, but the nodes, tracer density, ``nu
        sigma_r^2`` and ``beta`` are evaluated once and shared -- this is the
        routine the likelihood calls, so it is the hot path.
        """
        RR = np.atleast_1d(np.asarray(R, dtype=float))
        r, wu, _ = self._projection_nodes(RR)
        _, _, ln_nu = self._table
        with np.errstate(invalid="ignore"):
            nu = np.exp(ln_nu(np.log(r)))
        nu[~np.isfinite(nu) | (r >= self.r_cut)] = 0.0
        ns = self.nu_sigma_r2(r.ravel()).reshape(r.shape)
        beta = self.anisotropy.beta(r)
        ratio2 = (RR[:, None] / r) ** 2
        base = 2.0 * r * wu
        return {
            "Sigma": np.sum(nu * base, axis=1),
            "los": np.sum((1.0 - beta * ratio2) * ns * base, axis=1),
            "pmr": np.sum((1.0 - beta + beta * ratio2) * ns * base, axis=1),
            "pmt": np.sum((1.0 - beta) * ns * base, axis=1),
        }

    def dispersions_kms(self, R_pc: Any) -> dict[str, np.ndarray]:
        """All three projected dispersions in km/s at projected radii in pc.

        Raises
        ------
        ValueError
            If any radius lies at or beyond the tracer's edge ``r_cut``.
        """
        RR = np.atleast_1d(np.asarray(R_pc, dtype=float))
        self._check_inside_tracer(RR)
        m = self.projected_moments(RR)
        return {k: np.sqrt(m[k] / m["Sigma"]) for k in PROJECTIONS}

    def dispersions_observed(self, R_arcsec: Any, distance_kpc: float) -> dict[str, np.ndarray]:
        """Observables at projected radii in arcsec: ``sigma_los`` [km/s] and the
        proper-motion dispersions ``sigma_pmr``, ``sigma_pmt`` in mas/yr, for a
        cluster at ``distance_kpc``."""
        RR = np.atleast_1d(np.asarray(R_arcsec, dtype=float)) * distance_kpc * 1e3 / 206264.806
        kms = self.dispersions_kms(RR)
        to_masyr = 1.0 / (KMS_PER_MASYR_KPC * distance_kpc)
        return {"sigma_los": kms["los"], "sigma_pmr": kms["pmr"] * to_masyr,
                "sigma_pmt": kms["pmt"] * to_masyr}
