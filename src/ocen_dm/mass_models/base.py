"""Spherical mass components: the interface and its numerical safety net.

Unit system
-----------
Internal arrays are plain floats in **pc, solar masses, km/s**. Astropy units are
used at the boundaries (loaders, reports) and in the tests, not in the hot inner
loops that the inference will call millions of times. In this system

``G = 4.300917270036e-3  pc Msun^-1 (km/s)^2``

which is CODATA 2018 via ``astropy.constants``; ``test_mass_models.py`` asserts
the module constant against astropy so it cannot drift.

Contract
--------
Every component implements ``density``, ``enclosed_mass`` and ``potential``.
Whatever is given analytically, the base class supplies the rest numerically, and
``verify(...)`` checks the analytic forms against independent quadrature and
finite differences. A component whose analytic and numerical curves disagree is a
bug, not a tolerance to be widened.

Sign conventions: the potential is negative and tends to 0 at infinity; the
radial force ``force(r)`` is the inward (negative) radial acceleration
``-dPhi/dr``, in (km/s)^2 / pc.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Callable, Sequence

import numpy as np
from scipy import integrate, optimize

__all__ = ["G", "MassComponent", "VerificationReport", "as_array", "dphi_dr"]

#: Gravitational constant in pc Msun^-1 (km/s)^2 (CODATA 2018, via astropy).
G = 4.300917270036276e-3


def as_array(r: Any) -> np.ndarray:
    """Return ``r`` as a 1-D float array, rejecting negative radii.

    Parameters
    ----------
    r : array_like
        Radius or radii in pc.

    Returns
    -------
    numpy.ndarray
        At least 1-D, float64.

    Raises
    ------
    ValueError
        If any radius is negative or non-finite.
    """
    out = np.atleast_1d(np.asarray(r, dtype=float))
    if np.any(~np.isfinite(out)):
        raise ValueError("radii must be finite")
    if np.any(out < 0):
        raise ValueError("radii must be non-negative")
    return out


def dphi_dr(
    potential: Callable[[np.ndarray], np.ndarray],
    radii: np.ndarray,
    step: float = 1e-2,
) -> np.ndarray:
    """Differentiate a potential with a five-point stencil in ``ln r``.

    A naive central difference fails badly here, and the failure is arithmetic
    rather than physical. Deep inside a core the potential is a large near
    constant whose variation across a small step is only a few thousand times
    double-precision noise, so the difference is dominated by roundoff; and for
    components whose potential is itself computed by quadrature, the quadrature
    error is not smooth in ``r``, so differencing amplifies it.

    Both problems are solved by taking a *larger* step and recovering the
    accuracy with a higher-order stencil: the truncation error falls as
    ``step^4`` while the roundoff contribution falls as ``1/step``.

    Parameters
    ----------
    potential : callable
        Function of radius returning the potential in (km/s)^2.
    radii : numpy.ndarray
        Radii in pc.
    step : float, optional
        Fractional step in ``ln r``.

    Returns
    -------
    numpy.ndarray
        ``dPhi/dr`` in (km/s)^2 / pc.
    """
    shifts = (-2.0, -1.0, 1.0, 2.0)
    weights = (1.0 / 12.0, -2.0 / 3.0, 2.0 / 3.0, -1.0 / 12.0)
    derivative = np.zeros_like(radii, dtype=float)
    for shift, weight in zip(shifts, weights):
        derivative += weight * potential(radii * np.exp(shift * step))
    # d/dr = (1/r) d/dln r
    return derivative / (step * radii)


@dataclass(frozen=True)
class VerificationReport:
    """Outcome of checking a component's analytic forms against numerics."""

    name: str
    max_mass_error: float
    max_force_error: float
    finite_everywhere: bool
    monotonic_mass: bool
    passed: bool

    def __str__(self) -> str:  # pragma: no cover - display only
        return (
            f"{self.name}: mass err {self.max_mass_error:.2e}, "
            f"force err {self.max_force_error:.2e}, "
            f"finite={self.finite_everywhere}, monotonic={self.monotonic_mass}, "
            f"{'PASS' if self.passed else 'FAIL'}"
        )


class MassComponent(ABC):
    """A spherically symmetric mass distribution.

    Subclasses must implement :meth:`density`. They should override
    :meth:`enclosed_mass` and :meth:`potential` where closed forms exist; the
    numerical fallbacks here are correct but slow, and exist so that a new
    component is usable (and checkable) the moment its density is written.
    """

    #: short label used in reports and composite bookkeeping
    name: str = "component"

    @abstractmethod
    def density(self, r: Any) -> np.ndarray:
        """Mass density at radius ``r``.

        Parameters
        ----------
        r : array_like
            Radii in pc.

        Returns
        -------
        numpy.ndarray
            Density in Msun / pc^3.
        """

    @property
    def total_mass(self) -> float:
        """Total mass in Msun, ``inf`` when the profile does not converge."""
        return float(self.enclosed_mass(np.inf)[0])

    def enclosed_mass(self, r: Any) -> np.ndarray:
        """Mass enclosed within radius ``r``, by quadrature over the density.

        Parameters
        ----------
        r : array_like
            Radii in pc; ``inf`` is allowed and gives the total mass.

        Returns
        -------
        numpy.ndarray
            Enclosed mass in Msun.
        """
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        out = np.empty_like(radii)
        for i, radius in enumerate(radii):
            if radius == 0:
                out[i] = 0.0
                continue
            # Integrate in log radius: the integrand spans many decades and a
            # linear grid wastes its points at large radii.
            if np.isinf(radius):
                value, _ = integrate.quad(
                    lambda lu: 4 * np.pi * np.exp(3 * lu) * self.density(np.exp(lu))[0],
                    -30.0, 30.0, limit=400,
                )
            else:
                value, _ = integrate.quad(
                    lambda lu: 4 * np.pi * np.exp(3 * lu) * self.density(np.exp(lu))[0],
                    -30.0, np.log(radius), limit=400,
                )
            out[i] = value
        return out

    def potential(self, r: Any) -> np.ndarray:
        """Gravitational potential at ``r``, in (km/s)^2, zero at infinity.

        Uses the standard spherical result
        ``Phi(r) = -G[ M(<r)/r + 4*pi*int_r^inf rho(s) s ds ]``.
        """
        radii = as_array(r)
        inner = self.enclosed_mass(radii)
        out = np.empty_like(radii)
        for i, radius in enumerate(radii):
            outer, _ = integrate.quad(
                lambda lu: 4 * np.pi * np.exp(2 * lu) * self.density(np.exp(lu))[0],
                np.log(radius) if radius > 0 else -30.0, 30.0, limit=400,
            )
            term = inner[i] / radius if radius > 0 else 0.0
            out[i] = -G * (term + outer)
        return out

    def force(self, r: Any) -> np.ndarray:
        """Radial acceleration ``-dPhi/dr`` in (km/s)^2 / pc (negative inward).

        For a spherical system this is ``-G M(<r) / r^2``, which is what is used
        here: it is exact and avoids differentiating the potential.
        """
        radii = as_array(r)
        out = np.zeros_like(radii)
        nonzero = radii > 0
        out[nonzero] = -G * self.enclosed_mass(radii[nonzero]) / radii[nonzero] ** 2
        return out

    def circular_velocity(self, r: Any) -> np.ndarray:
        """Circular speed ``sqrt(G M(<r) / r)`` in km/s."""
        radii = as_array(r)
        out = np.zeros_like(radii)
        nonzero = radii > 0
        out[nonzero] = np.sqrt(G * self.enclosed_mass(radii[nonzero]) / radii[nonzero])
        return out

    # ------------------------------------------------------------------ checks
    def verify(
        self,
        r_min: float = 1e-3,
        r_max: float = 1e4,
        n: int = 40,
        mass_tol: float = 1e-6,
        force_tol: float = 1e-4,
        fd_step: float = 1e-2,
    ) -> VerificationReport:
        """Check the analytic forms against independent numerics.

        Three things are compared over a log grid: the analytic enclosed mass
        against quadrature of the density; the analytic force against a central
        finite difference of the potential; and the basic physical requirements
        that every value is finite and that enclosed mass does not decrease.

        Parameters
        ----------
        r_min, r_max : float, optional
            Radial range in pc.
        n : int, optional
            Number of log-spaced test radii.
        mass_tol, force_tol : float, optional
            Maximum acceptable relative error. The force tolerance is looser
            because it is limited by the finite-difference arithmetic, not by
            the model: see :func:`dphi_dr`.
        fd_step : float, optional
            Fractional step in ``ln r`` for the potential derivative.

        Returns
        -------
        VerificationReport
        """
        radii = np.geomspace(r_min, r_max, n)

        analytic_mass = self.enclosed_mass(radii)
        numeric_mass = MassComponent.enclosed_mass(self, radii)
        scale = np.maximum(np.abs(analytic_mass), 1e-30)
        mass_error = float(np.max(np.abs(analytic_mass - numeric_mass) / scale))

        derivative = dphi_dr(self.potential, radii, step=fd_step)
        analytic_force = self.force(radii)
        force_scale = np.maximum(np.abs(analytic_force), 1e-30)
        force_error = float(np.max(np.abs(analytic_force + derivative) / force_scale))

        finite = bool(
            np.all(np.isfinite(analytic_mass))
            and np.all(np.isfinite(self.density(radii)))
            and np.all(np.isfinite(self.potential(radii)))
        )
        monotonic = bool(np.all(np.diff(analytic_mass) >= -1e-9 * np.max(analytic_mass)))

        return VerificationReport(
            name=self.name,
            max_mass_error=mass_error,
            max_force_error=force_error,
            finite_everywhere=finite,
            monotonic_mass=monotonic,
            passed=bool(
                mass_error < mass_tol and force_error < force_tol and finite and monotonic
            ),
        )

    def radius_enclosing(self, mass: float, bracket: tuple[float, float] = (1e-6, 1e6)) -> float:
        """Return the radius containing ``mass`` Msun.

        Raises
        ------
        ValueError
            If the requested mass is not reached inside ``bracket``.
        """
        lo, hi = bracket
        if self.enclosed_mass(hi)[0] < mass:
            raise ValueError(
                f"{self.name}: encloses only {self.enclosed_mass(hi)[0]:.3e} Msun "
                f"within {hi} pc, less than the requested {mass:.3e}"
            )
        return float(
            optimize.brentq(lambda x: self.enclosed_mass(x)[0] - mass, lo, hi, xtol=1e-10)
        )
