"""Compact dark-remnant population.

Specification section 3.1: the remnants are a separate compact component with a
free total mass and scale radius, not absorbed into the luminous profile. v1 is
phenomenological -- a Plummer sphere -- and may later be replaced by a multi-mass
collisional prediction. It is a distinct class so that swap is a one-line change
and so the remnant mass appears under its own name in every posterior.
"""

from __future__ import annotations

from dataclasses import dataclass

from .stellar import Plummer

__all__ = ["RemnantPlummer"]


@dataclass(frozen=True)
class RemnantPlummer(Plummer):
    """Plummer sphere standing in for the dark-remnant population.

    Attributes
    ----------
    mass : float
        Total remnant mass ``M_rem`` in Msun.
    scale : float
        Remnant scale radius ``a_rem`` in pc, typically well inside the
        luminous scale since remnants are mass-segregated.
    """

    name: str = "remnants"
