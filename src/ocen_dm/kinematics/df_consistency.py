"""Necessary DF conditions; these tests never certify existence of a positive DF.

The central black-hole condition applies to regular spherical equilibrium DFs.
The finite-radius conditions assume a *separable augmented density*
nu_tilde(Psi, r^2) = P(Psi) R(r^2). Jeans fits do not impose this extra assumption.
See An & Evans (2006), and An, Van Hese & Baes (2012), sections 4.1--4.2.
All quantities concern the stellar tracer, not an unspecified DM/remnant DF.
"""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp

from .anisotropy import Anisotropy, TurnoverAnisotropy
from ..mass_models.stellar import MGE


def mge_log_slopes(tracer: MGE, r):
    """Return gamma=-d ln(nu)/d ln(r) and d gamma/d ln(r), without differencing."""
    r = np.asarray(r, float)
    if r.ndim != 1 or np.any(~np.isfinite(r)) or np.any(r <= 0):
        raise ValueError("r must be a finite positive vector")
    mass, sigma = tracer._masses, tracer._sigmas
    positive = mass > 0
    if not np.any(positive):
        raise ValueError("tracer must have nonzero mass")
    z = (r[:, None] / sigma[positive])**2
    logw = np.log(mass[positive]/sigma[positive]**3) - z/2
    w = np.exp(logw-logsumexp(logw, axis=1)[:, None])
    gamma = (w*z).sum(axis=1)
    variance = (w*(z-gamma[:, None])**2).sum(axis=1)
    return gamma, 2*gamma-variance


def beta_log_derivative(anisotropy, r):
    """Analytic d beta/d ln(r) for both fitted anisotropy families."""
    r = np.asarray(r, float)
    if isinstance(anisotropy, TurnoverAnisotropy):
        transitions = [(anisotropy.r_beta, anisotropy.beta_mid-anisotropy.beta_0),
                       (anisotropy.r_beta_outer, anisotropy.beta_inf-anisotropy.beta_mid)]
    elif isinstance(anisotropy, Anisotropy):
        transitions = [(anisotropy.r_beta, anisotropy.beta_inf-anisotropy.beta_0)]
    else:
        raise TypeError("unsupported anisotropy family")
    out = np.zeros_like(r)
    for scale, amplitude in transitions:
        fraction = r*r/(r*r+scale*scale)
        out += 2*amplitude*fraction*(1-fraction)
    return out


def necessary_profiles(model, r):
    """Return dimensionless quantities with signs of P', P'', and R_2.

    P=nu*g, g=exp(2 integral beta/r dr), F=GM/r^2, q=d ln M/d ln r.
    P'=P/(r F) B1; P''=P/(r^2 F^2) B2. Positive prefactors are removed.
    P' >= 0 is necessary only for beta0 <= 1/2; P'' for beta0 <= -1/2.
    R_2/R=(1-beta)(2-beta)-(d beta/d ln r)/2 is always necessary within
    the separable class. Passing these finite-order tests is not sufficient.
    """
    r = np.asarray(r, float)
    gamma, gamma_dot = mge_log_slopes(model.tracer, r)
    beta = model.anisotropy.beta(r)
    beta_dot = beta_log_derivative(model.anisotropy, r)
    mass = model.mass.enclosed_mass(r)
    rho = model.mass.density(r)  # excludes the point-mass delta at the origin
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0) or np.any(~np.isfinite(rho)):
        raise ValueError("invalid mass profile")
    q = 4*np.pi*r**3*rho/mass
    s = -gamma+2*beta
    return {"r_pc": r, "gamma": gamma, "beta": beta,
            "B1": -s, "B2": s*s+s-gamma_dot+2*beta_dot-s*q,
            "R2_over_R": (1-beta)*(2-beta)-beta_dot/2,
            "rho_total": rho, "rho_stars": model.tracer.density(r)}


def sampled_negative_intervals(r, values, tolerance=1e-8):
    """Each contiguous run of negative grid points, not an interpolated root."""
    r, values = np.asarray(r), np.asarray(values)
    if np.any(~np.isfinite(values)):
        raise ValueError("nonfinite diagnostic")
    indices = np.flatnonzero(values < -tolerance)
    if not len(indices):
        return []
    groups = np.split(indices, np.flatnonzero(np.diff(indices) != 1)+1)
    return [[float(r[g[0]]), float(r[g[-1]])] for g in groups]


def summarize_profiles(profiles, beta0, tolerance=1e-8):
    applicability = {"B1": beta0 <= .5, "B2": beta0 <= -.5,
                     "R2_over_R": True}
    out = {}
    for key, applicable in applicability.items():
        values = profiles[key]
        idx = int(np.argmin(values))
        intervals = sampled_negative_intervals(profiles["r_pc"], values, tolerance)
        out[key] = {"necessary_for_separable_df": bool(applicable),
                    "minimum": float(values[idx]),
                    "radius_at_minimum_pc": float(profiles["r_pc"][idx]),
                    "negative_intervals_pc": intervals,
                    "fails_necessary_condition": bool(applicable and intervals)}
    out["separable_df_ruled_out"] = any(v["fails_necessary_condition"] for v in out.values())
    out["df_nonnegative_certified"] = False
    return out
