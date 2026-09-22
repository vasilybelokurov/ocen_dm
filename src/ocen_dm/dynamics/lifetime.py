"""Physical heating estimates and resolved stellar-observable comparisons."""
from __future__ import annotations

import numpy as np

from . import G, TIME_UNIT_MYR
from .diagnostics import bound_particles, effective_number, stellar_centre


def heating_time_gyr(vrms_kms, density_msun_pc3, effective_mass_msun, coulomb_log):
    """Bertone & Fairbairn (2008), eq. 11; local equal-speed approximation.

    Effective mass <m^2>/<m> represents a physical stellar mass spectrum.
    This estimate is not a collisional evolution operator or a depletion law.
    """
    values = np.broadcast_arrays(vrms_kms, density_msun_pc3, effective_mass_msun, coulomb_log)
    if any(np.any(~np.isfinite(v) | (v <= 0)) for v in values):
        raise ValueError('Heating inputs must be positive and finite')
    velocity, density, mass, logarithm = values
    return .814*velocity**3/(G**2*mass*density*1e9*logarithm)*TIME_UNIT_MYR/1000


def projected_stellar_profiles(xv, particles, config, reference, data, distance_kpc=5.34,
                               minimum_neff=30):
    """Bin the bound nucleus like the measured kinematics, along three axes.

    Constant stellar M/L and complete annular selection are assumed. Cartesian
    views bracket orientation; these are diagnostic residuals, not a calibrated
    joint likelihood or a survival acceptance cut. Low-count bins stay missing.
    """
    if not config.nucleus.live:
        raise ValueError('Stellar survival requires a live nucleus')
    local = xv-stellar_centre(xv, particles, config, reference)
    bound, _ = bound_particles(local, particles['mass'], config.nucleus)
    chosen = bound & (particles['component'] == len(config.components))
    position, velocity, weight = local[chosen, :3], local[chosen, 3:], particles['mass'][chosen]
    rows = []
    for axis, label in enumerate('xyz'):
        plane = [k for k in range(3) if k != axis]
        xy, vxy = position[:, plane], velocity[:, plane]
        radius = np.linalg.norm(xy, axis=1)
        radial = np.sum(xy*vxy, axis=1)/np.maximum(radius, 1e-30)
        tangential = (xy[:, 0]*vxy[:, 1]-xy[:, 1]*vxy[:, 0])/np.maximum(radius, 1e-30)
        arcsec = radius/distance_kpc*206264.806247
        for profile in data.profiles:
            values = {'los': velocity[:, axis], 'pmr': radial/(4.74047*distance_kpc),
                      'pmt': tangential/(4.74047*distance_kpc)}[profile.kind]
            for index, (lo, hi) in enumerate(zip(profile.r_lower, profile.r_upper)):
                select = (arcsec >= lo) & (arcsec < hi)
                neff = effective_number(weight[select])
                predicted = np.nan
                if neff >= minimum_neff:
                    mean = np.average(values[select], weights=weight[select])
                    predicted = np.sqrt(np.average((values[select]-mean)**2, weights=weight[select]))
                observed = profile.value[index]
                error = profile.err_hi[index] if predicted > observed else profile.err_lo[index]
                rows.append(dict(view=label, dataset=profile.name, kind=profile.kind,
                                 radius_arcsec=float(profile.r[index]), observed=float(observed),
                                 err_lo=float(profile.err_lo[index]), err_hi=float(profile.err_hi[index]),
                                 prediction=float(predicted), neff=neff,
                                 residual_sigma=float((predicted-observed)/error)))
    return rows
