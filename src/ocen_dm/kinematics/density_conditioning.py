"""Condition an archived halo model on its local density without changing priors.

For log-uniform M100, rho = M100 / A(rs) implies
p(rho, rs) = p(rs) / [rho ln(Mmax/Mmin)] inside the mass support.
Thus fixing rho preserves p(rs) if every rs remains inside that support.
This adapter checks that condition and refuses cases requiring a truncated
scale-radius prior, rather than silently implementing a different conditional.
"""
from __future__ import annotations

import numpy as np

from .fit import DarkMatterModel
from ..mass_models import TruncatedGNFW


class DensityConditionedModel:
    def __init__(self, parent: DarkMatterModel, density: float, radius_pc: float = 20.):
        if type(parent) is not DarkMatterModel:
            raise TypeError("density conditioning requires a DarkMatterModel parent")
        if not np.isfinite(density) or density <= 0 or not 0 < radius_pc < 100:
            raise ValueError("need positive finite density and a radius inside 100 pc")
        priors = {p.name:p.prior for p in parent.parameters}
        if "M_dm_100" not in priors or "r_s" not in priors:
            raise ValueError("both halo mass and scale radius must be sampled in the parent")
        if priors["M_dm_100"].kind != "loguniform" or priors["r_s"].kind != "loguniform":
            raise ValueError("conditional transform requires log-uniform halo mass and radius")
        self.parent = parent
        self.density = float(density)
        self.radius_pc = float(radius_pc)
        self.mass_prior = priors["M_dm_100"]
        self.parameters = tuple(p for p in parent.parameters if p.name != "M_dm_100")
        self.names = tuple(p.name for p in self.parameters)
        self.label = f"{parent.label}_rho{radius_pc:g}_{density:g}"

        # Normalised halo shapes are monotonic with rs separately inside and
        # outside radius_pc. Mixing endpoint integrals gives conservative bounds
        # on M100 for EVERY rs, including possible interior extrema.
        masses = []
        for rs in (priors["r_s"].lo, priors["r_s"].hi):
            h = self.halo({"r_s":rs})
            inner, total = h.enclosed_mass([radius_pc, 100.])
            masses.append((inner, total-inner))
        lower = min(m[0] for m in masses)+min(m[1] for m in masses)
        upper = max(m[0] for m in masses)+max(m[1] for m in masses)
        self.mass_support_envelope = [float(lower), float(upper)]
        if lower <= self.mass_prior.lo or upper >= self.mass_prior.hi:
            raise ValueError("fixed density may truncate the parent mass support; "
                             "this requires a separately normalised conditional radius prior")

    def __getattr__(self, name):
        return getattr(self.parent, name)

    def transform(self, u):
        u = np.asarray(u, float)
        return np.array([p.prior.transform(u[..., i]) for i,p in enumerate(self.parameters)]).T

    def to_dict(self, x):
        return dict(zip(self.names, np.asarray(x, float)))

    def halo(self, theta):
        unit = TruncatedGNFW(rho_s=1., r_s=theta["r_s"], gamma=self.parent.gamma,
                             r_t=self.parent.r_t)
        normalisation = self.density/float(unit.density(self.radius_pc).item())
        return TruncatedGNFW(rho_s=normalisation, r_s=theta["r_s"], gamma=self.parent.gamma,
                             r_t=self.parent.r_t, name="dm_cored" if self.parent.gamma==0 else "dm_nfw")

    def complete(self, theta):
        full = self.parent.complete(theta)
        mass = float(self.halo(full).enclosed_mass(100.).item())
        if not self.mass_prior.lo <= mass <= self.mass_prior.hi:
            raise ValueError("derived halo mass lies outside parent prior support")
        full["M_dm_100"] = mass
        return full

    def build(self, theta):
        return self.parent.build(self.complete(theta))

    @property
    def prior_density_at_constraint(self):
        """Parent marginal prior PDF per unit density, valid after support check."""
        return 1/(self.density*np.log(self.mass_prior.hi/self.mass_prior.lo))
