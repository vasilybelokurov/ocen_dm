"""Joint spherical equilibria, with an independent signed DF check.

AGAMA's QuasiSpherical constructor clips negative DF values. We therefore
check the *signed* constant-beta inversion before constructing/sampling it,
and also require recovery of the requested density and anisotropy.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import gamma as gamma_function, roots_jacobi

from . import G


def agama_module():
    import agama
    agama.setUnits(length=1, velocity=1, mass=1)
    return agama


def xyz(r):
    r = np.atleast_1d(r)
    return np.column_stack((r, np.zeros_like(r), np.zeros_like(r)))


def density_parameters(component, norm=1.):
    if component.profile == "plummer":
        return dict(type="Plummer", mass=norm, scaleRadius=component.scale_pc/1000)
    return dict(type="Spheroid", densityNorm=norm, scaleRadius=component.scale_pc / 1000,
                gamma=component.gamma, beta=component.outer_slope, alpha=1,
                outerCutoffRadius=component.cutoff_pc / 1000, cutoffStrength=2)


def density_log_derivatives(component, r):
    """d ln rho/d ln r and its logarithmic derivative, before anisotropy."""
    x = np.asarray(r) / (component.scale_pc / 1000)
    if component.profile == "plummer":
        return -5*x*x/(1+x*x), -10*x*x/(1+x*x)**2
    y = np.asarray(r) / (component.cutoff_pc / 1000)
    slope = -component.gamma + (component.gamma - component.outer_slope) * x / (1 + x)
    curvature = (component.gamma - component.outer_slope) * x / (1 + x)**2
    return slope - 2*y*y, curvature - 4*y*y


def pericentre_selection(potential, radius_pc, floor):
    """Importance sampling on orbital pericentre, not instantaneous radius.

    p=1 inside the protection radius, declining as r_peri^-2 outside it to
    a positive floor. Inverse-p particle weights recover the original DF.
    These are fixed masses; later tidal evolution can bring heavy particles
    inward, which must be monitored and tested with finer refinements.
    """
    def probability(xv):
        # AGAMA also evaluates its selection function at trial points outside
        # the bound DF's support, where Rperiapo is intentionally undefined.
        energy = potential.potential(xv[:, :3])+.5*np.sum(xv[:, 3:]**2, axis=1)
        valid = np.isfinite(energy) & (energy < 0)
        result = np.zeros(len(xv))
        if not np.any(valid):
            return result
        peri = potential.Rperiapo(xv[valid])[:, 0]
        if not np.all(np.isfinite(peri)):
            raise ValueError("Invalid pericentre in phase-space refinement")
        result[valid] = np.maximum(floor, np.minimum(1., (radius_pc/1000/np.maximum(peri, 1e-30))**2))
        return result
    return probability


def signed_df_grid(potential, density, slopes, beta, rmin, rmax, n_energy=160,
                   n_grid=8192, quadrature_order=192):
    """Unclipped Cuddeford inversion for constant -1/2 < beta < 1/2.

    Returns the energy factor in f(E,L)=f_E(E)*L**(-2*beta). The outer
    boundary term vanishes for the exponentially tapered densities used here.
    Gauss-Jacobi quadrature integrates the endpoint weight exactly. This is a
    finite-grid diagnostic, not a proof of positivity at every energy.
    """
    if not -.5 < beta < .5:
        raise ValueError("Signed inversion requires -1/2 < beta < 1/2")
    r = np.geomspace(rmin / 10, rmax * 10, n_grid)
    pos = xyz(r)
    psi = -potential.potential(pos)
    force = -potential.force(pos)[:, 0]
    if np.any(force <= 0) or np.any(np.diff(psi) >= 0):
        raise ValueError("The spherical relative potential must decrease strictly with radius")
    rho = density.density(pos)
    s, ds = slopes(r)
    s = s + 2*beta
    augmented = rho * r**(2*beta)
    # Phi'' = 4 pi G rho_total - 2 Phi'/r in spherical symmetry.
    second_phi = 4*np.pi*G*potential.density(pos) - 2*force/r
    second = augmented / r**2 / force**2 * (s*s - s + ds - s*r*second_phi/force)
    spline = PchipInterpolator(np.log(psi[::-1]), second[::-1], extrapolate=False)
    energy_r = np.geomspace(rmin, rmax, n_energy)
    binding = -potential.potential(xyz(energy_r))
    nodes, weights = roots_jacobi(quadrature_order, 0, beta-.5)
    fraction = (nodes + 1)/2
    weights = weights / 2**(beta+.5)
    integrand = spline(np.log(binding[:, None] * (1-fraction)))
    # Only the negligible far tail can lie below the interpolation domain.
    integrand = np.where(binding[:, None]*(1-fraction) < psi[-1], 0., integrand)
    if not np.all(np.isfinite(integrand)):
        raise ValueError("Signed DF integration left its tabulated potential domain")
    coefficient = 2**beta / ((2*np.pi)**1.5 * gamma_function(1-beta) * gamma_function(.5+beta))
    factor = coefficient * binding**(beta+.5)
    signed = factor * (integrand @ weights)
    absolute = factor * (np.abs(integrand) @ weights)
    cancellation = np.divide(signed, absolute, out=np.zeros_like(signed), where=absolute > 0)
    return dict(radius_kpc=energy_r, binding_energy_kms2=binding, f_energy=signed,
                signed_over_absolute=cancellation)


def build_equilibrium(config, directory):
    """Build all DFs in the same nucleus+stars+DM potential, then sample.

    The rigid nucleus is included exactly once. The finite-density profiles
    already specify a self-consistent spherical mass model; inverting them in
    its *combined* potential needs no iterative adjustment of their masses.
    """
    agama = agama_module()
    nucleus = agama.Potential(type="Plummer", mass=config.nucleus.mass_msun,
                              scaleRadius=config.nucleus.scale_pc/1000)
    # AGAMA export() does not serialize parameters of analytic Plummer models.
    (directory / "nucleus.ini").write_text(
        f"[Potential]\ntype=Plummer\nmass={config.nucleus.mass_msun:.17g}\n"
        f"scaleRadius={config.nucleus.scale_pc/1000:.17g}\n")
    densities, potentials, normalizations = [], [], []
    for component in config.particle_components:
        factory = agama.Potential if component.profile == "plummer" else agama.Density
        unit = factory(**density_parameters(component))
        reference = (unit.totalMass() if component.mass_msun is not None else
                     unit.enclosedMass(component.radius_pc / 1000))
        target = component.mass_msun if component.mass_msun is not None else component.mass_within_msun
        norm = target / reference
        density = factory(**density_parameters(component, norm))
        potential = (density if component.profile == "plummer" else
                     agama.Potential(type="Multipole", density=density, symmetry="s",
                                    lmax=0, gridSizeR=240,
                                    rmin=min(component.scale_pc, config.nucleus.scale_pc)*1e-7,
                                    rmax=component.cutoff_pc*.01))
        if component.profile == "plummer":
            (directory/f"{component.name}_potential.ini").write_text(
                f"[Potential]\ntype=Plummer\nmass={norm:.17g}\nscaleRadius={component.scale_pc/1000:.17g}\n")
        else:
            potential.export(str(directory / f"{component.name}_potential.ini"))
        densities.append(density)
        potentials.append(potential)
        normalizations.append(norm)
    original = agama.Potential(*potentials) if config.nucleus.live else agama.Potential(nucleus, *potentials)
    sections = [] if config.nucleus.live else ["[Potential nucleus]\nfile=nucleus.ini\n"]
    sections.extend(f"[Potential {c.name}]\nfile={c.name}_potential.ini\n" for c in config.particle_components)
    (directory / "satellite_initial.ini").write_text("\n".join(sections))
    total = agama.Potential(str(directory / "satellite_initial.ini"))
    probes = xyz(np.geomspace(config.nucleus.scale_pc*1e-5,
                              max(c.cutoff_pc for c in config.components)*.003, 100))
    if not np.allclose(total.force(probes), original.force(probes), rtol=1e-6, atol=1e-8):
        raise ValueError("Saved satellite potential does not reproduce its original force")
    positions, masses, softenings, labels = [], [], [], []
    report = {"method": "constant-beta DF in combined spherical potential",
              "nucleus": {"representation": "live Plummer" if config.nucleus.live else "rigid Plummer", "mass_msun": config.nucleus.mass_msun,
                          "scale_pc": config.nucleus.scale_pc}, "components": {}}
    for number, (component, density, norm) in enumerate(zip(config.particle_components, densities, normalizations)):
        print(f"Validating {component.name} equilibrium (beta={component.beta:g})", flush=True)
        rmin = min(component.scale_pc, config.nucleus.scale_pc)*1e-6
        rmax = component.cutoff_pc*.005
        signed = signed_df_grid(total, density, lambda r: density_log_derivatives(component, r),
                                component.beta, rmin, rmax)
        np.savez_compressed(directory / f"{component.name}_signed_df.npz", **signed)
        if np.any(signed["signed_over_absolute"] < -.005):
            worst = np.argmin(signed["signed_over_absolute"])
            raise ValueError(f"{component.name}: negative distribution function near energy "
                             f"{ -signed['binding_energy_kms2'][worst]:.6g} (km/s)^2; "
                             "change the density/beta model instead of clipping the DF")
        df = agama.DistributionFunction(type="QuasiSpherical", potential=total,
                                        density=density, beta0=component.beta)
        galaxy = agama.GalaxyModel(total, df)
        check_r = np.geomspace(min(component.scale_pc, config.nucleus.scale_pc)*1e-4,
                               component.cutoff_pc*.002, 40)
        desired = density.density(xyz(check_r))
        recovered, moments = galaxy.moments(xyz(check_r), dens=True, vel=False, vel2=True)
        beta_recovered = 1 - (moments[:, 1] + moments[:, 2])/(2*moments[:, 0])
        error = np.abs(recovered/desired - 1)
        np.savez_compressed(directory / f"{component.name}_density_check.npz",
                            radius_kpc=check_r, requested_density=desired,
                            recovered_density=recovered, recovered_beta=beta_recovered)
        if not np.all(np.isfinite(error)) or np.max(error) > .05:
            raise ValueError(f"{component.name}: DF density recovery failed (max error {np.max(error):.3g})")
        if not np.all(np.isfinite(beta_recovered)) or np.max(np.abs(beta_recovered-component.beta)) > .03:
            raise ValueError(f"{component.name}: DF anisotropy recovery failed")
        # The public AGAMA RNG is independent of numpy. Single-threaded sampling
        # and a separately derived component seed make runs reproducible.
        agama.setRandomSeed(config.seed + number)
        selection = (pericentre_selection(total, component.refinement_radius_pc, component.refinement_floor)
                     if component.refinement_radius_pc is not None else None)
        sampling_model = agama.GalaxyModel(total, df, sf=selection) if selection else galaxy
        with agama.setNumThreads(1):
            xv, weight = sampling_model.sample(component.n_particles)
        if selection:
            weight = weight/selection(xv)
        if not np.all(np.isfinite(weight)) or np.any(weight <= 0):
            raise ValueError(f"{component.name}: invalid particle masses")
        expected_mass = float(density.totalMass())
        sampled_mass = float(np.sum(weight))
        neff = weight.sum()**2/np.sum(weight**2)
        if abs(sampled_mass/expected_mass-1) > max(.03, 5/np.sqrt(neff)):
            raise ValueError(f"{component.name}: sampled DF normalization differs from requested mass")
        weight = weight * (expected_mass/sampled_mass)
        if not np.all(np.isfinite(xv)) or np.any(total.potential(xv[:, :3]) + .5*np.sum(xv[:, 3:]**2, axis=1) >= 0):
            raise ValueError(f"{component.name}: invalid or unbound initial particles")
        # Do not move a sampled cusp to its noisy finite-N centre of mass:
        # the density and DF were constructed about the analytic nucleus.
        positions.append(xv)
        masses.append(weight)
        softenings.append(component.softening_pc/1000*(weight/weight.min())**(1/3))
        labels.append(np.full(len(weight), number, dtype=np.int16))
        report["components"][component.name] = {
            "profile": component.profile,
            "density_norm_msun_kpc3": float(norm), "total_mass_msun": expected_mass,
            "sampled_df_mass_msun": sampled_mass,
            "particle_mass_range_msun": [float(weight.min()), float(weight.max())],
            "n_particles": len(weight), "effective_particle_number": float(neff),
            "sampling": "pericentre importance with inverse-probability masses" if selection else "equal mass",
            "refinement_radius_pc": component.refinement_radius_pc,
            "refinement_floor": component.refinement_floor if selection else None,
            "max_density_fractional_error": float(np.max(error)),
            "max_beta_error": float(np.max(np.abs(beta_recovered-component.beta))),
            "minimum_signed_df_ratio": float(np.min(signed["signed_over_absolute"])),
            "signed_df_radius_range_pc": [rmin*1000, rmax*1000],
            "mass_inside_pc": {str(r): float(density.enclosedMass(r/1000)) for r in config.radii_pc},
            "sampled_aperture_neff": {
                str(r): float(weight[np.linalg.norm(xv[:, :3], axis=1) <= r/1000].sum()**2 /
                              max(np.sum(weight[np.linalg.norm(xv[:, :3], axis=1) <= r/1000]**2), 1e-300))
                for r in config.radii_pc},
            "dm_to_nucleus_mass_ratio" if component.name == "dm" else "stars_to_nucleus_mass_ratio":
                {str(r): float(density.enclosedMass(r/1000)/nucleus.enclosedMass(r/1000))
                 for r in config.radii_pc},
        }
    xv = np.vstack(positions)
    particles = dict(xv=xv, mass=np.concatenate(masses), softening=np.concatenate(softenings),
                     component=np.concatenate(labels), particle_id=np.arange(len(xv), dtype=np.int64))
    if config.nucleus.live:
        bulk = np.average(particles["xv"][:, 3:], axis=0, weights=particles["mass"])
        particles["xv"][:, 3:] -= bulk
        report["removed_sample_bulk_velocity_kms"] = bulk.tolist()
    np.savez_compressed(directory / "initial_local.npz", **particles)
    return particles, report
