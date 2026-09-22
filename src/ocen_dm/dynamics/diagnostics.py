"""Aperture, shell, binding and instantaneous tidal diagnostics.

The iterative binding diagnostic uses a spherical, self-excluded particle
potential about the prescribed nucleus. It is not a conserved Jacobi energy
or an exact bound/unbound classification in an eccentric, barred host.
All particles, including diagnosed escapers, remain in the simulation.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from . import G


def effective_number(mass):
    mass = np.asarray(mass)
    return float(mass.sum()**2/np.sum(mass**2)) if len(mass) and np.any(mass) else 0.


def weighted_median(values, mass):
    if not len(mass):
        return np.nan
    order = np.argsort(values)
    weights = np.asarray(mass)[order]
    cumulative = np.cumsum(weights)-.5*weights
    return float(np.interp(.5*weights.sum(), cumulative, np.asarray(values)[order]))


def nucleus_potential(r, nucleus):
    if nucleus.live:
        return np.zeros_like(np.asarray(r), dtype=float)
    return -G*nucleus.mass_msun/np.sqrt(np.asarray(r)**2+(nucleus.scale_pc/1000)**2)


def nucleus_mass(r, nucleus):
    r = np.asarray(r)
    if nucleus.live:
        return np.zeros_like(r, dtype=float)
    return nucleus.mass_msun*r**3/(r*r+(nucleus.scale_pc/1000)**2)**1.5


def spherical_particle_potential(radius, mass, sources):
    """Monopole potential at each particle, excluding its own contribution."""
    r = np.maximum(np.asarray(radius), 1e-15)
    order = np.argsort(r)
    rr = r[order]
    mm = np.where(sources, mass, 0.)[order]
    inner = np.cumsum(mm)
    outer = np.cumsum((mm/rr)[::-1])[::-1] - mm/rr
    result = np.empty_like(r)
    result[order] = -G*((inner-mm)/rr+outer)
    return result


def bound_particles(xv, mass, nucleus, iterations=100):
    radius = np.linalg.norm(xv[:, :3], axis=1)
    kinetic = .5*np.sum(xv[:, 3:]**2, axis=1)
    external = nucleus_potential(radius, nucleus)
    keep = np.ones(len(radius), dtype=bool)
    for _ in range(iterations):
        energy = kinetic + external + spherical_particle_potential(radius, mass, keep)
        next_keep = keep & (energy < 0)
        if np.array_equal(keep, next_keep):
            return keep, energy
        keep = next_keep
    raise RuntimeError("Iterative spherical unbinding did not converge")


def tidal_tensor(host, position, time=0., step_kpc=None):
    x = np.asarray(position, float)
    step = step_kpc or max(np.linalg.norm(x)*1e-5, 1e-6)
    columns = []
    for axis in np.eye(3)*step:
        columns.append((host.force(x+axis, t=time)-host.force(x-axis, t=time))/(2*step))
    tensor = np.column_stack(columns)
    return (tensor+tensor.T)/2


def tidal_scale(host, centre, time, mass_at_radius):
    """Largest-eigenvalue scale with instantaneous omega=x cross v / |x|^2.

    This convention includes the centrifugal matrix, but neither the Euler
    term nor Coriolis velocity terms in the spatial eigenvalue diagnostic.
    It recovers the circular spherical denominator 3-dlnM/dlnR.
    """
    if host is None:
        return np.inf, 0.
    x, v = centre[:3], centre[3:]
    rr = np.dot(x, x)
    if rr == 0:
        return np.nan, np.nan
    omega = np.cross(x, v)/rr
    tensor = tidal_tensor(host, x, time) + np.dot(omega, omega)*np.eye(3)-np.outer(omega, omega)
    lam = float(np.linalg.eigvalsh(tensor)[-1])
    if lam <= 0:
        return np.inf, lam
    # Select the outermost crossing; a cored satellite may have no crossing.
    radii = np.geomspace(1e-8, max(np.sqrt(rr), 1.), 400)
    values = np.array([G*mass_at_radius(r)/r**3-lam for r in radii])
    crossings = np.where((values[:-1] >= 0) & (values[1:] < 0))[0]
    if not len(crossings):
        return (0. if values[0] <= 0 else np.nan), lam
    i = crossings[-1]
    root = brentq(lambda r: G*mass_at_radius(r)/r**3-lam, radii[i], radii[i+1], xtol=1e-12)
    return float(root), lam


def velocity_moments(xv, mass):
    if len(mass) < 3 or effective_number(mass) < 3:
        return np.nan, np.nan
    x, v = xv[:, :3], xv[:, 3:]
    radius = np.linalg.norm(x, axis=1)
    er = x/np.maximum(radius[:, None], 1e-30)
    phi = np.arctan2(x[:, 1], x[:, 0])
    ep = np.column_stack((-np.sin(phi), np.cos(phi), np.zeros(len(phi))))
    et = np.cross(ep, er)
    spherical = np.column_stack((np.sum(v*er, axis=1), np.sum(v*et, axis=1), np.sum(v*ep, axis=1)))
    mean = np.average(spherical, axis=0, weights=mass)
    variance = np.average((spherical-mean)**2, axis=0, weights=mass)
    beta = 1-(variance[1]+variance[2])/(2*variance[0]) if variance[0] > 0 else np.nan
    return float(np.sqrt(variance[0])), float(beta)


def stellar_centre(xv, particles, config, reference):
    """Shrinking-sphere centre of the live nucleus; outer escapers do not set it.

    This locates a stellar concentration; survival is assessed separately from
    bound mass and size. It is not evidence of survival by itself.
    """
    if not config.nucleus.live:
        return np.asarray(reference)
    number = len(config.components)
    member = particles["component"] == number
    stars, mass = np.asarray(xv)[member], particles["mass"][member]
    if len(stars) < 16:
        raise ValueError("Too few nucleus particles to locate the stellar centre")
    centre = np.array([weighted_median(stars[:, k], mass) for k in range(3)])
    radius = np.linalg.norm(stars[:, :3]-centre, axis=1)
    keep = radius <= np.quantile(radius, .8)
    minimum = min(len(stars), max(16, int(.02*len(stars))))
    while keep.sum() > minimum:
        centre = np.average(stars[keep, :3], weights=mass[keep], axis=0)
        radius = np.linalg.norm(stars[:, :3]-centre, axis=1)
        new_keep = radius <= np.quantile(radius[keep], .8)
        if new_keep.sum() < minimum or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    centre = np.average(stars[keep, :3], weights=mass[keep], axis=0)
    # More stars than in the positional core reduce bulk-velocity noise.
    radius = np.linalg.norm(stars[:, :3]-centre, axis=1)
    inner = radius <= np.quantile(radius, .5)
    velocity = np.average(stars[inner, 3:], weights=mass[inner], axis=0)
    return np.r_[centre, velocity]


def mass_radius(radius, mass, fraction=.5):
    if not len(mass) or np.sum(mass) <= 0:
        return np.nan
    order = np.argsort(radius)
    cumulative = np.cumsum(mass[order])
    return float(np.interp(fraction*cumulative[-1], cumulative, radius[order]))


def snapshot_diagnostics(xv, particles, config, centre, time, host=None):
    prescribed = np.asarray(centre)
    centre = stellar_centre(xv, particles, config, prescribed)
    local = xv-centre
    r = np.linalg.norm(local[:, :3], axis=1)
    mass = particles["mass"]
    bound, _ = bound_particles(local, mass, config.nucleus)
    order = np.argsort(r[bound])
    bound_r, cumulative = r[bound][order], np.cumsum(mass[bound][order])

    def enclosed(radius):
        i = np.searchsorted(bound_r, radius, side="right")
        return float(nucleus_mass(radius, config.nucleus)+(cumulative[i-1] if i else 0))

    rj, lam = tidal_scale(host, centre, time, enclosed)
    result = dict(tidal_scale_pc=rj*1000, disruptive_eigenvalue=lam,
                  bound_particle_mass_msun=float(mass[bound].sum()),
                  nucleus_mass_msun=float(mass[particles["component"] == len(config.components)].sum()) if config.nucleus.live else config.nucleus.mass_msun,
                  centre_offset_from_track_pc=float(np.linalg.norm(centre[:3]-prescribed[:3])*1000),
                  centre_velocity_offset_from_track_kms=float(np.linalg.norm(centre[3:]-prescribed[3:])),
                  centre_radius_kpc=float(np.linalg.norm(centre[:3])))
    result.update({f"centre_{name}": float(value) for name, value in zip(
        ["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"], centre)})
    for number, component in enumerate(config.particle_components):
        member = particles["component"] == number
        prefix = component.name
        result[f"{prefix}_total_mass_msun"] = float(mass[member].sum())
        result[f"{prefix}_bound_mass_msun"] = float(mass[member & bound].sum())
        result[f"{prefix}_bound_neff"] = effective_number(mass[member & bound])
        chosen = member & bound
        result[f"{prefix}_bound_fraction"] = float(mass[chosen].sum()/mass[member].sum())
        result[f"{prefix}_half_mass_radius_pc"] = 1000*mass_radius(r[chosen], mass[chosen])
        for axis, label in enumerate(["x", "y", "z"]):
            projected = np.sqrt(np.maximum(r*r-local[:, axis]**2, 0))
            result[f"{prefix}_projected_half_mass_{label}_pc"] = 1000*mass_radius(projected[chosen], mass[chosen])
            if np.any(chosen):
                v = local[chosen, axis+3]
                mean = np.average(v, weights=mass[chosen])
                result[f"{prefix}_sigma_{label}_kms"] = float(np.sqrt(np.average((v-mean)**2, weights=mass[chosen])))
            else:
                result[f"{prefix}_sigma_{label}_kms"] = np.nan
        if component.refinement_radius_pc is not None:
            inner = member & (r < component.refinement_radius_pc/1000)
            heavy = member & (mass > 10*np.min(mass[member]))
            result[f"{prefix}_heavy_mass_fraction_inside_refinement_radius"] = float(
                mass[inner & heavy].sum()/mass[inner].sum()) if np.any(inner) else np.nan
        for radius in config.radii_pc:
            aperture = member & (r <= radius/1000)
            shell = member & (r >= radius/1000*np.exp(-.1)) & (r < radius/1000*np.exp(.1))
            volume_pc3 = 4*np.pi/3*radius**3*(np.exp(.3)-np.exp(-.3))
            tag = f"{prefix}_{radius:g}pc"
            result[f"{tag}_mass_msun"] = float(mass[aperture].sum())
            result[f"{tag}_bound_mass_msun"] = float(mass[aperture & bound].sum())
            result[f"{tag}_rho_msun_pc3"] = float(mass[shell].sum()/volume_pc3)
            result[f"{tag}_bound_rho_msun_pc3"] = float(mass[shell & bound].sum()/volume_pc3)
            result[f"{tag}_aperture_neff"] = effective_number(mass[aperture])
            result[f"{tag}_shell_neff"] = effective_number(mass[shell])
            sigma, beta = velocity_moments(local[shell], mass[shell])
            result[f"{tag}_sigma_r_kms"], result[f"{tag}_beta"] = sigma, beta
    return result
