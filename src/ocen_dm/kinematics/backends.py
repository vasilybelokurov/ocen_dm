"""Alternative dynamical engines behind the same binned-profile likelihood.

Every back-end exposes ``projected_moments(R_pc) -> {'Sigma', 'los', 'pmr', 'pmt'}``
(tracer surface density and ``Sigma * sigma^2`` for the three projections, in
tracer units / pc^2 and km^2/s^2) and ``_check_inside_tracer(R)``, exactly like
:class:`~ocen_dm.kinematics.jeans.SphericalJeans`, so
:class:`~ocen_dm.kinematics.likelihood.ProfileLikelihood` accepts any of them.

* :class:`JamBackend` -- JamPy's spherical anisotropic Jeans solution
  (Cappellari 2008 eq. 50; PMs from Cappellari 2020 App. B3). Same physics as our
  solver, independent implementation. Mass components that are not Gaussian
  (remnant Plummer, dark haloes) are projected numerically and fitted with a
  projected MGE by NNLS before being handed over. JamPy is imported, never
  vendored (non-commercial licence).
* :class:`AgamaDFBackend` -- a positive distribution function of the
  Cuddeford-Osipkov-Merritt family built by AGAMA for the tracer density in the
  total potential (``beta(r) = (beta0 r_a^2 + r^2) / (r_a^2 + r^2)``), with the
  projected moments from ``GalaxyModel.moments``. Different anisotropy family
  from the Jeans models (beta -> 1 at large r unless r_a = inf), and it exists
  only where a positive DF does.
"""

from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from functools import cached_property
from typing import Any

import numpy as np
from scipy import integrate, optimize

from ..mass_models import MGE, CompositeMassModel, PointMass
from ..mass_models.base import MassComponent
from .anisotropy import Anisotropy

__all__ = ["JamBackend", "AgamaDFBackend", "project_density", "mge_from_component"]

ARCSEC_PER_RAD = 206264.806


# ------------------------------------------------------------- projection ---
def project_density(component: MassComponent, R: np.ndarray, r_max: float = 1e4) -> np.ndarray:
    """Surface density ``Sigma(R) = 2 int_R^inf rho(r) r dr / sqrt(r^2 - R^2)`` (Msun/pc^2)."""
    R = np.atleast_1d(np.asarray(R, float))
    out = np.empty_like(R)
    for i, Ri in enumerate(R):
        # substitute r = R cosh u to remove the endpoint singularity
        # r = R cosh u, dr / sqrt(r^2 - R^2) = du  ->  Sigma = 2 int rho(R cosh u) R cosh u du
        f = lambda u: float(np.asarray(component.density(Ri * np.cosh(u))).reshape(-1)[0]) * Ri * np.cosh(u)
        u_max = np.arccosh(max(r_max / Ri, 1.0 + 1e-9))
        out[i] = 2.0 * integrate.quad(f, 0.0, u_max, limit=200)[0]
    return out


def mge_from_component(component: MassComponent, sigmas_pc: np.ndarray, r_fit: np.ndarray | None = None) -> MGE:
    """Projected NNLS MGE of a spherical mass component, so JamPy can use it.

    Fits ``Sigma(R)`` on ``r_fit`` (default 100 log points 0.05-2000 pc) with
    Gaussians of the given widths; the returned MGE's total mass is renormalised
    to the component's mass within ``r_fit.max()``.
    """
    r_fit = np.geomspace(0.05, 2000.0, 100) if r_fit is None else np.asarray(r_fit, float)
    Sigma = project_density(component, r_fit)
    basis = np.exp(-0.5 * (r_fit[:, None] / sigmas_pc[None, :]) ** 2) / (2 * np.pi * sigmas_pc[None, :] ** 2)
    # weight by 1/Sigma so the fit is in relative terms over the dynamic range
    w = 1.0 / np.maximum(Sigma, 1e-300)
    masses, _ = optimize.nnls(basis * w[:, None], Sigma * w, maxiter=200 * len(sigmas_pc))
    keep = masses > 0
    mge = MGE(masses=masses[keep].tolist(), sigmas=sigmas_pc[keep].tolist(), name=f"mge_of_{component.name}")
    target = float(np.asarray(component.enclosed_mass(r_fit.max())).reshape(-1)[0])
    got = float(np.asarray(mge.enclosed_mass(r_fit.max())).reshape(-1)[0])
    return MGE(masses=(np.array(mge.masses) * target / got).tolist(), sigmas=mge.sigmas, name=mge.name)


def _split_mass_model(mass: MassComponent) -> tuple[list[MGE], float, list[MassComponent]]:
    """(Gaussian components, point mass, everything else)."""
    comps = list(mass.components) if isinstance(mass, CompositeMassModel) else [mass]
    gaussians, mbh, other = [], 0.0, []
    for c in comps:
        if isinstance(c, MGE):
            gaussians.append(c)
        elif isinstance(c, PointMass):
            mbh += float(c.mass)
        else:
            other.append(c)
    return gaussians, mbh, other


# ------------------------------------------------------------------ JamPy ---
@dataclass(frozen=True)
class JamBackend:
    """JamPy spherical Jeans back-end (see module docstring)."""

    mass: MassComponent
    tracer: MGE
    anisotropy: Anisotropy = Anisotropy()
    distance_kpc: float = 5.43
    epsrel: float = 1e-4
    mge_sigmas_pc: tuple[float, ...] = tuple(np.geomspace(0.05, 500.0, 40))

    @cached_property
    def _pc_per_arcsec(self) -> float:
        return self.distance_kpc * 1e3 / ARCSEC_PER_RAD

    @cached_property
    def _mass_mge(self) -> tuple[np.ndarray, np.ndarray, float]:
        """(surf_pot [Msun/pc^2], sigma_pot [arcsec], mbh) for the total potential."""
        gaussians, mbh, other = _split_mass_model(self.mass)
        for c in other:
            gaussians.append(mge_from_component(c, np.asarray(self.mge_sigmas_pc)))
        m = np.concatenate([np.asarray(g.masses) * g.mass_to_light for g in gaussians])
        s = np.concatenate([np.asarray(g.sigmas) for g in gaussians])
        return m / (2 * np.pi * s**2), s / self._pc_per_arcsec, mbh

    @property
    def r_cut(self) -> float:
        return 50.0 * float(np.max(self.tracer.sigmas))

    def _check_inside_tracer(self, RR: np.ndarray) -> None:
        if np.any(RR >= self.r_cut):
            raise ValueError("projected radius beyond the tracer")

    def surface_density(self, R: Any) -> np.ndarray:
        R = np.atleast_1d(np.asarray(R, float))[:, None]
        m = np.asarray(self.tracer.masses)[None, :] * self.tracer.mass_to_light
        s = np.asarray(self.tracer.sigmas)[None, :]
        return np.sum(m / (2 * np.pi * s**2) * np.exp(-0.5 * (R / s) ** 2), axis=1)

    def projected_moments(self, R: Any) -> dict[str, np.ndarray]:
        from jampy.sph.jam_sph_proj import jam_sph_proj

        RR = np.atleast_1d(np.asarray(R, float))
        rad_arcsec = RR / self._pc_per_arcsec
        surf_pot, sig_pot, mbh = self._mass_mge
        m = np.asarray(self.tracer.masses) * self.tracer.mass_to_light
        s = np.asarray(self.tracer.sigmas)
        surf_lum, sig_lum = m / (2 * np.pi * s**2), s / self._pc_per_arcsec
        beta = [self.anisotropy.r_beta / self._pc_per_arcsec, self.anisotropy.beta_0, self.anisotropy.beta_inf, 2.0]
        Sigma = self.surface_density(RR)
        out = {"Sigma": Sigma}
        # JamPy sorts/needs monotonic radii; pass unique sorted radii and map back
        uniq, inv = np.unique(rad_arcsec, return_inverse=True)
        for tensor in ("los", "pmr", "pmt"):
            with contextlib.redirect_stdout(io.StringIO()):
                res = jam_sph_proj(surf_lum, sig_lum, surf_pot, sig_pot, mbh, self.distance_kpc / 1e3, uniq,
                                   beta=beta, logistic=True, tensor=tensor, epsrel=self.epsrel, plot=False, quiet=True)
            sigma = np.asarray(res.model, float)[inv]
            out[tensor] = Sigma * sigma**2
        return out

    def dispersions_kms(self, R_pc: Any) -> dict[str, np.ndarray]:
        m = self.projected_moments(R_pc)
        return {k: np.sqrt(m[k] / m["Sigma"]) for k in ("los", "pmr", "pmt")}


# ------------------------------------------------------------------ AGAMA ---
@dataclass(frozen=True)
class AgamaDFBackend:
    """AGAMA QuasiSpherical (Cuddeford-Osipkov-Merritt) DF back-end.

    Parameters
    ----------
    beta0 : float
        Central anisotropy (<= 1/2 for a cored tracer to admit a positive DF).
    r_a : float
        Anisotropy radius in pc; ``inf`` keeps ``beta = beta0`` everywhere.
    """

    mass: MassComponent
    tracer: MassComponent
    beta0: float = 0.0
    r_a: float = np.inf
    distance_kpc: float = 5.43
    r_min: float = 1e-3
    r_max: float = 1e4

    @cached_property
    def _model(self):
        import agama
        agama.setUnits(mass=1, length=1e-3, velocity=1)          # Msun, pc, km/s
        gaussians, mbh, other = _split_mass_model(self.mass)
        smooth = CompositeMassModel(gaussians + other) if (gaussians or other) else None
        parts = []
        if smooth is not None:
            dens = lambda xyz: np.asarray(smooth.density(np.sqrt(np.sum(np.atleast_2d(xyz) ** 2, axis=1))), float)
            parts.append(agama.Potential(type="Multipole", density=dens, symmetry="s", lmax=0,
                                         rmin=self.r_min, rmax=self.r_max, gridSizeR=60))
        if mbh > 0:
            parts.append(agama.Potential(type="Plummer", mass=mbh, scaleRadius=1e-4))
        pot = agama.Potential(*parts) if len(parts) > 1 else parts[0]
        tracer = self.tracer
        tdens = agama.Density(lambda xyz: np.asarray(tracer.density(np.sqrt(np.sum(np.atleast_2d(xyz) ** 2, axis=1))), float),
                              symmetry="s") if not hasattr(tracer, "_agama_density") else tracer._agama_density
        kwargs = dict(type="QuasiSpherical", potential=pot, density=tdens, beta0=self.beta0)
        if np.isfinite(self.r_a):
            kwargs["r_a"] = self.r_a
        df = agama.DistributionFunction(**kwargs)
        return agama.GalaxyModel(pot, df)

    @property
    def r_cut(self) -> float:
        return self.r_max

    def _check_inside_tracer(self, RR: np.ndarray) -> None:
        if np.any(RR >= self.r_cut):
            raise ValueError("projected radius beyond the tracer")

    def beta(self, r: Any) -> np.ndarray:
        r = np.asarray(r, float)
        return (self.beta0 * self.r_a**2 + r**2) / (self.r_a**2 + r**2) if np.isfinite(self.r_a) else np.full_like(r, self.beta0)

    def density_check(self, r: Any) -> np.ndarray:
        """Realised / input tracer density: 1 where a positive DF reproduces the tracer.

        Departures mark where the requested anisotropy has no positive
        distribution function for this tracer in this potential (An & Evans
        2006); the Jeans equation still has a solution there, the DF does not.
        """
        rr = np.atleast_1d(np.asarray(r, float))
        rho, _ = self._model.moments(np.column_stack([rr, 0 * rr, 0 * rr]), dens=True, vel=False, vel2=True)
        return np.asarray(rho, float) / np.asarray(self.tracer.density(rr), float)

    def projected_moments(self, R: Any) -> dict[str, np.ndarray]:
        RR = np.atleast_1d(np.asarray(R, float))
        gm = self._model
        pts = np.column_stack([RR, np.zeros_like(RR)])            # observed X = R, Y = 0 -> XX radial, YY tangential, ZZ los
        Sigma, v2 = gm.moments(pts, dens=True, vel=False, vel2=True)
        Sigma = np.asarray(Sigma, float); v2 = np.asarray(v2, float)
        # moments() returns the surface density and the density-NORMALISED second
        # moments <v_i v_j> integrated along Z; our convention is Sigma * <v^2>
        return {"Sigma": Sigma, "pmr": Sigma * v2[:, 0], "pmt": Sigma * v2[:, 1], "los": Sigma * v2[:, 2]}

    def dispersions_kms(self, R_pc: Any) -> dict[str, np.ndarray]:
        m = self.projected_moments(R_pc)
        return {k: np.sqrt(m[k] / m["Sigma"]) for k in ("los", "pmr", "pmt")}
