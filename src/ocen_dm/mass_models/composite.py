"""The composite cluster model: ``Phi_tot = Phi_star + Phi_rem + Phi_IMBH + Phi_DM``.

This is the object the kinematic likelihood and the stream models both consume,
and the place where the specification's headline quantities are defined:
``M_DM(<r)``, ``f_DM(<r)``, ``v_esc(r)`` and the Jacobi radius (sections 1 and 6).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy import optimize

from .base import G, MassComponent, VerificationReport, as_array, dphi_dr
from .imbh import PointMass

__all__ = ["CompositeMassModel"]


@dataclass(frozen=True)
class CompositeMassModel(MassComponent):
    """Sum of spherical components, with the derived quantities of the project.

    Attributes
    ----------
    components : sequence of MassComponent
        The parts of the model. Names should be unique; they key the per-component
        outputs. A component may be omitted entirely (no IMBH, no DM) rather than
        being included with zero mass, and both spellings behave identically.
    """

    components: Sequence[MassComponent]
    name: str = "composite"

    def __post_init__(self) -> None:
        if not self.components:
            raise ValueError("a composite model needs at least one component")
        names = [c.name for c in self.components]
        if len(set(names)) != len(names):
            raise ValueError(f"component names must be unique, got {names}")

    # ------------------------------------------------------------------ access
    def __getitem__(self, name: str) -> MassComponent:
        for component in self.components:
            if component.name == name:
                return component
        raise KeyError(f"no component named {name!r}; have {[c.name for c in self.components]}")

    def get(self, name: str) -> MassComponent | None:
        """Return the named component, or ``None`` when it is absent."""
        try:
            return self[name]
        except KeyError:
            return None

    @property
    def component_names(self) -> tuple[str, ...]:
        """Names of the components, in order."""
        return tuple(c.name for c in self.components)

    # ------------------------------------------------------------- the profile
    def density(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        # A point mass has zero density away from the origin; including its
        # formal delta function here would make the sum meaningless.
        return sum(
            (c.density(radii) for c in self.components if not isinstance(c, PointMass)),
            start=np.zeros_like(radii),
        )

    def enclosed_mass(self, r: Any) -> np.ndarray:
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        return sum(
            (c.enclosed_mass(radii) for c in self.components),
            start=np.zeros_like(radii),
        )

    def potential(self, r: Any) -> np.ndarray:
        radii = as_array(r)
        return sum(
            (c.potential(radii) for c in self.components),
            start=np.zeros_like(radii),
        )

    @property
    def total_mass(self) -> float:
        return float(sum(c.total_mass for c in self.components))

    # -------------------------------------------------------- derived products
    def mass_breakdown(self, r: Any) -> dict[str, np.ndarray]:
        """Enclosed mass of every component plus the total, at ``r``.

        Returns
        -------
        dict
            ``{component name: M(<r)}`` with an extra ``'total'`` entry. These are
            the per-sample quantities the specification asks to store for every
            posterior draw (section 6, Experiment K2).
        """
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        out = {c.name: c.enclosed_mass(radii) for c in self.components}
        out["total"] = sum(out.values(), start=np.zeros_like(radii))
        return out

    def dark_matter_mass(self, r: Any, prefix: str = "dm") -> np.ndarray:
        """Enclosed mass of the dark-matter components only.

        Parameters
        ----------
        r : array_like
            Radii in pc.
        prefix : str, optional
            Components whose name starts with this are treated as dark matter.

        Returns
        -------
        numpy.ndarray
            ``M_DM(<r)`` in Msun; zeros when the model has no DM component.
        """
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        parts = [c.enclosed_mass(radii) for c in self.components if c.name.startswith(prefix)]
        return sum(parts, start=np.zeros_like(radii))

    def dark_matter_fraction(self, r: Any, prefix: str = "dm") -> np.ndarray:
        """``f_DM(<r) = M_DM(<r) / M_total(<r)``, zero where the total mass is zero."""
        radii = np.atleast_1d(np.asarray(r, dtype=float))
        total = self.enclosed_mass(radii)
        dark = self.dark_matter_mass(radii, prefix=prefix)
        return np.divide(dark, total, out=np.zeros_like(total), where=total > 0)

    def escape_speed(self, r: Any, r_esc: float | None = None) -> np.ndarray:
        """Escape speed in km/s.

        Parameters
        ----------
        r : array_like
            Radii in pc.
        r_esc : float, optional
            Radius at which a star is considered to have escaped. With ``None``
            the reference is infinity, giving ``sqrt(-2 Phi(r))``. In a tidal
            field the physically relevant choice is the Jacobi radius, and
            passing it gives ``sqrt(2[Phi(r_esc) - Phi(r)])`` -- a smaller,
            more honest escape speed. Which one is used matters for the tails,
            so it is explicit rather than assumed.

        Returns
        -------
        numpy.ndarray
            Escape speed in km/s, zero where the star is already unbound.
        """
        radii = as_array(r)
        reference = 0.0 if r_esc is None else float(self.potential(r_esc)[0])
        delta = 2.0 * (reference - self.potential(radii))
        return np.sqrt(np.clip(delta, 0.0, None))

    def jacobi_radius(
        self,
        galactocentric_radius: float,
        enclosed_galactic_mass: float,
        bracket: tuple[float, float] = (1e-3, 1e4),
    ) -> float:
        """Solve ``r_J = R (M_c(<r_J) / (3 M_G(<R)))^(1/3)`` for ``r_J``.

        The cluster mass inside the Jacobi radius depends on the radius itself,
        so this is implicit and solved numerically rather than evaluated with the
        total mass.

        Parameters
        ----------
        galactocentric_radius : float
            ``R`` in pc.
        enclosed_galactic_mass : float
            ``M_G(<R)`` in Msun.
        bracket : tuple of float, optional
            Search bracket in pc.

        Returns
        -------
        float
            Jacobi radius in pc.

        Raises
        ------
        ValueError
            If no solution exists inside ``bracket``.
        """
        if galactocentric_radius <= 0 or enclosed_galactic_mass <= 0:
            raise ValueError("galactocentric radius and enclosed galactic mass must be positive")

        def residual(radius: float) -> float:
            cluster = self.enclosed_mass(radius)[0]
            return radius - galactocentric_radius * (
                cluster / (3.0 * enclosed_galactic_mass)
            ) ** (1.0 / 3.0)

        lo, hi = bracket
        if residual(lo) * residual(hi) > 0:
            raise ValueError(
                f"no Jacobi radius in [{lo}, {hi}] pc for R={galactocentric_radius} pc "
                f"and M_G={enclosed_galactic_mass:.3e} Msun"
            )
        return float(optimize.brentq(residual, lo, hi, xtol=1e-10))

    def radial_profile(self, radii: Any, prefix: str = "dm") -> dict[str, np.ndarray]:
        """All the per-sample derived quantities on one radial grid.

        Returns
        -------
        dict
            ``r``, the per-component and total enclosed masses, ``M_dm``,
            ``f_dm``, ``rho_total``, ``v_circ`` and ``v_esc``. This is the record
            the specification asks to store for every posterior sample
            (section 6, Experiment K2).
        """
        grid = np.atleast_1d(np.asarray(radii, dtype=float))
        out: dict[str, np.ndarray] = {"r": grid}
        for key, value in self.mass_breakdown(grid).items():
            out[f"M_{key}"] = value
        out["M_dm"] = self.dark_matter_mass(grid, prefix=prefix)
        out["f_dm"] = self.dark_matter_fraction(grid, prefix=prefix)
        out["rho_total"] = self.density(grid)
        out["v_circ"] = self.circular_velocity(grid)
        out["v_esc"] = self.escape_speed(grid)
        return out

    # ------------------------------------------------------------------ checks
    def verify_all(self, **kwargs: Any) -> list[VerificationReport]:
        """Verify every component and the composite itself.

        A :class:`PointMass` is skipped for the density-integral check, because
        its density is a delta function; its enclosed mass and potential are
        exercised through the composite.
        """
        reports = [
            c.verify(**kwargs) for c in self.components if not isinstance(c, PointMass)
        ]
        reports.append(self._verify_composite(**kwargs))
        return reports

    def _verify_composite(self, **kwargs: Any) -> VerificationReport:
        """Check the composite's own force and potential against its mass profile.

        The density-integral test does not apply when a point mass is present,
        since that mass is deliberately absent from ``density``; the force test
        below covers it, because the force uses the enclosed mass.
        """
        r_min = kwargs.get("r_min", 1e-3)
        r_max = kwargs.get("r_max", 1e4)
        n = kwargs.get("n", 40)
        radii = np.geomspace(r_min, r_max, n)

        has_point_mass = any(isinstance(c, PointMass) for c in self.components)
        if has_point_mass:
            mass_error = 0.0
        else:
            analytic = self.enclosed_mass(radii)
            numeric = MassComponent.enclosed_mass(self, radii)
            mass_error = float(
                np.max(np.abs(analytic - numeric) / np.maximum(np.abs(analytic), 1e-30))
            )

        derivative = dphi_dr(self.potential, radii, step=kwargs.get("fd_step", 1e-2))
        force = self.force(radii)
        force_error = float(
            np.max(np.abs(force + derivative) / np.maximum(np.abs(force), 1e-30))
        )

        masses = self.enclosed_mass(radii)
        finite = bool(np.all(np.isfinite(masses)) and np.all(np.isfinite(self.potential(radii))))
        monotonic = bool(np.all(np.diff(masses) >= -1e-9 * np.max(masses)))
        return VerificationReport(
            name=self.name,
            max_mass_error=mass_error,
            max_force_error=force_error,
            finite_everywhere=finite,
            monotonic_mass=monotonic,
            passed=bool(
                mass_error < kwargs.get("mass_tol", 1e-6)
                and force_error < kwargs.get("force_tol", 1e-4)
                and finite and monotonic
            ),
        )
