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


class Converged(Exception):
    """Raised by the fitter when the absolute objective-change rule is met."""


class AbsoluteStop:
    """Stop when accepted iterates improve the objective by < delta over `window` steps.

    The trust-region optimizer evaluates the Jacobian only at accepted iterates,
    so `update` is called there with the objective sum(residual**2). An absolute
    rule is used because relative tolerances are unreachable near a noiseless
    minimum, where evaluation noise is comparable to the objective itself.
    """

    def __init__(self, delta=0.01, window=2):
        if delta <= 0 or window < 1:
            raise ValueError("delta must be positive and window at least one")
        self.delta, self.window, self.history = float(delta), int(window), []

    def update(self, objective):
        self.history.append(float(objective))
        if len(self.history) > self.window and \
                self.history[-1-self.window]-self.history[-1] < self.delta:
            raise Converged(f"objective improved by less than {self.delta:g} over "
                            f"{self.window} accepted steps")


def data_fit_gates(optimizer_success, numerical_passed, evaluation, thresholds):
    """Goodness of fit to observed data; says nothing about the mass decomposition.

    Photometry is judged by chi2/N with its adopted errors. The observed splice
    (oMEGACat counts plus the Trager et al. 1995 compilation) has errors of
    0.1-0.58 mag and 0.13 mag scatter between duplicate full-weight points, so an
    unweighted magnitude RMS is not a meaningful criterion. The error-weighted
    RMS in mag is reported for reference only. Legacy threshold sets that name
    ``photometric_rms_mag`` apply it to the weighted RMS.
    """
    terms = evaluation["terms"]
    n_kin = sum(t["n"] for t in terms.values())
    chi2_kin = sum(t["chi2"] for t in terms.values())
    photo = evaluation["photometry"]
    residual = np.asarray(photo["prediction"])-np.asarray(evaluation["photometry_mu"])
    sigma = np.asarray(evaluation.get("photometry_sigma", np.ones_like(residual)))
    weighted_rms = float(np.sqrt(np.average(residual**2, weights=sigma**-2)))
    n_phot = residual.size
    per_term = {k: t["chi2"]/t["n"] for k, t in terms.items()}
    phot_ok = (photo["chi2"]/n_phot < thresholds["chi2_phot_per_point"]
               if "chi2_phot_per_point" in thresholds
               else weighted_rms < thresholds["photometric_rms_mag"])
    fit_ok = bool(chi2_kin/n_kin < thresholds["chi2_kin_per_point"] and
                  max(per_term.values()) < thresholds["chi2_per_point_any_dataset"] and phot_ok)
    return dict(optimizer_terminated=bool(optimizer_success), numerical_passed=bool(numerical_passed),
                chi2_kinematic=chi2_kin, n_kinematic=n_kin, chi2_kin_per_point=chi2_kin/n_kin,
                chi2_per_point_by_dataset=per_term, chi2_photometric=photo["chi2"],
                n_photometric=n_phot, chi2_phot_per_point=photo["chi2"]/n_phot,
                photometric_weighted_rms_mag=weighted_rms, data_fit_passed=fit_ok,
                passed=bool(optimizer_success and numerical_passed and fit_ok))


def data_radius_pc(problem, distance_kpc):
    """Outermost projected radius (pc) used by any kinematic bin edge or photometric point."""
    arcsec = [problem.photometry.r_arcsec.max()]
    for p in problem.data.profiles:
        arcsec.append((p.r_upper if p.r_upper is not None else p.r).max())
        if p.r_nodes is not None:
            arcsec.append(np.max(p.r_nodes))
    return float(max(arcsec)*distance_kpc*1000/ARCSEC_PER_RAD)
