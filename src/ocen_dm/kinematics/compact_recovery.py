"""Noiseless compact-DF recovery utilities; all stars share one light/mass DF.

Observed radii, selection nodes and error scales are retained, but mock values
come from an equilibrium and streaming is zero. Errors weight a diagnostic
objective; this module does not compute evidences or posterior uncertainties.
"""
from dataclasses import replace

import numpy as np
from scipy.special import roots_legendre

from .df_fit import DFJointProblem, PhotometricData
from .likelihood import ARCSEC_PER_RAD, KinematicData, ProfileLikelihood


def mock_problem(template, photometry, truth):
    """Preserve bin selection and asymmetry, replace values with noiseless means."""
    data = KinematicData(tuple(replace(p, streaming2=None) for p in template.profiles))
    predictions = ProfileLikelihood(data).predict(truth, truth.config.distance_kpc)
    data = KinematicData(tuple(replace(p, value=predictions[p.name].copy(),
                                      note="Noiseless compact DF mock; no streaming")
                               for p in data.profiles))
    radius = photometry.r_arcsec*truth.config.distance_kpc*1000/ARCSEC_PER_RAD
    sigma = truth.projected_moments(radius)["Sigma"]
    photo = PhotometricData(photometry.r_arcsec.copy(), -2.5*np.log10(sigma),
                           photometry.sigma_mag.copy(), "Noiseless compact DF mock",
                           photometry.error_model)
    return DFJointProblem(data, photo)


def residual_vector(problem, evaluation):
    """Signed residuals whose squared norm equals the existing joint objective."""
    rows = []
    for p in problem.data.profiles:
        prediction = evaluation["predictions"][p.name]
        error = np.where(prediction > p.value, p.err_hi, p.err_lo)
        rows.append((prediction-p.value)/error)
    rows.append((evaluation["photometry"]["prediction"]-problem.photometry.mu)
                / problem.photometry.sigma_mag)
    residual = np.concatenate(rows)
    if not np.all(np.isfinite(residual)):
        raise ValueError("nonfinite mock residual")
    return residual


def refined_config(config):
    """Double potential, velocity, projection AND action/contour resolutions."""
    n = config.numerics
    kw = dict(numerics=replace(n, potential_nodes=2*n.potential_nodes,
                              velocity_nodes=2*n.velocity_nodes,
                              moment_nodes=2*n.moment_nodes,
                              projection_nodes=2*n.projection_nodes,
                              iteration_tolerance=n.iteration_tolerance/2))
    if hasattr(config, "contours"):
        c = config.contours
        kw["contours"] = replace(c, action_nodes=2*c.action_nodes,
                                 circularity_nodes=2*c.circularity_nodes-1,
                                 normalization_order=2*c.normalization_order,
                                 ode_rtol=c.ode_rtol/2)
    return replace(config, **kw)


def mass_profiles(model, radii):
    """Enclosed stellar, remnant and halo masses on the same physical radii."""
    radii = np.asarray(radii, float)
    if np.any(radii <= 0) or np.any(~np.isfinite(radii)):
        raise ValueError("mass radii must be finite and positive")
    xyz = np.column_stack((radii, radii*0, radii*0))
    from .positive_df import agama_pc
    stars = -model.stellar_potential.force(xyz)[:, 0]*radii**2/agama_pc().G
    matter = model.config.matter
    remnants = matter.M_rem*radii**3/(radii**2+matter.a_rem**2)**1.5
    z, w = roots_legendre(128)
    r = radii[:, None]*(z+1)/2
    halo = 2*np.pi*radii*np.sum(w*r*r*matter.halo_density(r), axis=1)
    return dict(r_pc=radii.tolist(), stars=stars.tolist(), remnants=remnants.tolist(),
                halo=halo.tolist(), total=(stars+remnants+halo).tolist(),
                M_star=float(model.diagnostics["stellar_mass"]), rho20=matter.rho20)


def recovery_metrics(truth, recovered):
    """Physical accuracy distinct from goodness of fit; zero halo uses abs error."""
    if truth["r_pc"] != recovered["r_pc"]:
        raise ValueError("recovery profiles must use identical radii")
    star = abs(recovered["M_star"]/truth["M_star"]-1)
    total = float(np.max(abs(np.array(recovered["total"])/truth["total"]-1)))
    rho = abs(recovered["rho20"]-truth["rho20"])
    return dict(stellar_mass_fractional_error=star, total_mass_max_fractional_error=total,
                rho20_absolute_error=rho,
                rho20_relative_error=rho/truth["rho20"] if truth["rho20"] else None,
                physical_accuracy_passed=bool(star < .05 and total < .1 and
                    (rho/truth["rho20"] < .1 if truth["rho20"] else rho < .1)))


def recovery_gates(optimizer_success, numerical_passed, residual, metrics):
    residual = np.asarray(residual)
    rms = float(np.sqrt(np.mean(residual**2)))
    maximum = float(np.max(abs(residual)))
    observable = bool(rms < .1 and maximum < .3)
    return dict(optimizer_terminated=bool(optimizer_success),
                numerical_passed=bool(numerical_passed), rms_residual_sigma=rms,
                max_residual_sigma=maximum, observable_accuracy_passed=observable,
                physical_accuracy_passed=metrics["physical_accuracy_passed"],
                passed=bool(optimizer_success and numerical_passed and observable
                            and metrics["physical_accuracy_passed"]))
