"""Our own proper-motion dispersion profile from the Vasiliev & Baumgardt (2021) members.

The published EDR3 profile is a smooth fitted model whose percentile errors are ~1 per cent
and which no bound-tracer Jeans model reproduces (JOURNAL 2026-09-18). This module measures
the dispersion directly from the member catalogue so that the outer kinematics can be
checked rather than taken on trust:

* per-star radial and tangential proper motions about the cluster centre, with the error
  **covariance projected onto the same directions** (``pmra_error``, ``pmdec_error`` and
  ``pmra_pmdec_corr`` -- ignoring the correlation biases the projected errors by up to the
  correlation itself);
* a one-dimensional maximum-likelihood dispersion per bin with the per-star errors
  deconvolved, fitting the mean simultaneously (the mean radial PM is not zero at large
  radius: the perspective contraction of a cluster with v_los = 232.6 km/s reaches
  0.1 mas/yr at 2400 arcsec);
* the diagnostics that decide whether an outer dispersion is real: its dependence on the
  membership cut and on magnitude, and the expected field contamination ``sum(1 - P)``.

Nothing here enters a fit; it is a measurement and its diagnostics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from astropy.table import Table
from scipy import optimize

from ..paths import processed_dir
from .perspective import depth_dispersion, systemic_pm_field, systemic_velocity_vector

__all__ = ["MemberSample", "load_members", "systemic_pm", "dispersion_ml", "binned_dispersion",
           "mixture_dispersion", "mixture_dispersion_free", "field_pdf_factory", "contamination_audit",
           "OCEN_RA", "OCEN_DEC", "OCEN_VSYS_KMS", "KMS_PER_MASYR_KPC"]

OCEN_RA = 201.696833          # Baumgardt catalogue centre, deg
OCEN_DEC = -47.476583
OCEN_VSYS_KMS = 232.6         # systemic line-of-sight velocity (oMEGACat VI)
KMS_PER_MASYR_KPC = 4.740470463533348


@dataclass(frozen=True)
class MemberSample:
    """Members with radial/tangential PMs and projected per-star errors (mas/yr).

    The systemic proper motion has been subtracted before the projection: projecting the
    absolute PM onto the radial/tangential directions turns the systemic motion into an
    azimuthal pattern of amplitude ``|mu_sys|``, which reads as a spurious dispersion of
    ``|mu_sys| / sqrt(2)`` ~ 5.3 mas/yr for omega Cen (found 2026-09-18).
    """

    r_arcsec: np.ndarray
    mu_r: np.ndarray
    mu_t: np.ndarray
    err_r: np.ndarray
    err_t: np.ndarray
    prob: np.ndarray
    g_mag: np.ndarray
    quality_flag: np.ndarray
    phi: np.ndarray                       # position angle of each star, rad, north through east
    mu_sys: tuple[float, float] = (0.0, 0.0)
    exact: bool = True                    # systemic field projected star by star, or a constant

    def __len__(self) -> int:
        return len(self.r_arcsec)

    def select(self, mask: np.ndarray) -> "MemberSample":
        return MemberSample(*(np.asarray(getattr(self, f))[mask] for f in
                              ("r_arcsec", "mu_r", "mu_t", "err_r", "err_t", "prob", "g_mag", "quality_flag", "phi")),
                            mu_sys=self.mu_sys, exact=self.exact)


def systemic_pm(table: Table, r_arcsec: np.ndarray, prob_min: float = 0.9,
                r_max: float = 600.0) -> tuple[float, float]:
    """Error-weighted mean PM of the secure members inside ``r_max`` (mas/yr)."""
    m = (np.asarray(table["membership_prob"], float) >= prob_min) & (r_arcsec < r_max)
    out = []
    for col, err in (("pmra", "pmra_error"), ("pmdec", "pmdec_error")):
        w = 1.0 / np.asarray(table[err], float)[m] ** 2
        out.append(float(np.sum(np.asarray(table[col], float)[m] * w) / np.sum(w)))
    return out[0], out[1]


def load_members(path: Any = None, ra0: float = OCEN_RA, dec0: float = OCEN_DEC,
                 mu_sys: tuple[float, float] | None = None, v_los: float = OCEN_VSYS_KMS,
                 distance_kpc: float = 5.43, exact: bool = True) -> MemberSample:
    """Read the processed member table and project PMs and their covariance to (R, T).

    The systemic motion is removed **star by star**: the cluster's 3-D velocity (from
    ``mu_sys`` at the centre, ``v_los`` and ``distance_kpc``) is projected onto each star's
    own tangent basis (:func:`~ocen_dm.kinematics.perspective.systemic_pm_field`), which
    carries the perspective contraction and the rotation of the equatorial basis across the
    field exactly. ``exact=False`` subtracts the constant ``mu_sys`` instead (the naive
    treatment, kept for comparison). ``mu_sys`` defaults to the error-weighted mean of the
    secure members inside 600 arcsec, iterated once so that the field's own perspective
    term does not bias it.
    """
    table = Table.read(path or processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    ra = np.asarray(table["ra"], float); dec = np.asarray(table["dec"], float)
    x = (ra - ra0) * np.cos(np.radians(dec0)) * 3600.0          # arcsec, east positive
    y = (dec - dec0) * 3600.0
    r = np.hypot(x, y)
    safe = np.maximum(r, 1e-9)
    cos_p, sin_p = x / safe, y / safe                            # radial unit vector in (pmra*, pmdec)
    if mu_sys is None:
        mu_sys = systemic_pm(table, r)
        if exact:                                   # remove the field once, re-estimate the centre value
            v0 = systemic_velocity_vector(ra0, dec0, mu_sys[0], mu_sys[1], v_los, distance_kpc)
            fa, fd = systemic_pm_field(ra, dec, v0, distance_kpc)
            fa0, fd0 = systemic_pm_field(np.array([ra0]), np.array([dec0]), v0, distance_kpc)
            corr = Table({"pmra": np.asarray(table["pmra"], float) - (fa - fa0[0]),
                          "pmdec": np.asarray(table["pmdec"], float) - (fd - fd0[0]),
                          "pmra_error": table["pmra_error"], "pmdec_error": table["pmdec_error"],
                          "membership_prob": table["membership_prob"]})
            mu_sys = systemic_pm(corr, r)
    if exact:
        v_sys = systemic_velocity_vector(ra0, dec0, mu_sys[0], mu_sys[1], v_los, distance_kpc)
        exp_a, exp_d = systemic_pm_field(ra, dec, v_sys, distance_kpc)
    else:
        exp_a, exp_d = mu_sys[0], mu_sys[1]
    pmra = np.asarray(table["pmra"], float) - exp_a
    pmdec = np.asarray(table["pmdec"], float) - exp_d
    mu_r = pmra * cos_p + pmdec * sin_p
    mu_t = -pmra * sin_p + pmdec * cos_p
    ea = np.asarray(table["pmra_error"], float); ed = np.asarray(table["pmdec_error"], float)
    rho = np.asarray(table["pmra_pmdec_corr"], float)
    cov_ad = rho * ea * ed
    err_r = np.sqrt(np.maximum((ea * cos_p) ** 2 + (ed * sin_p) ** 2 + 2 * cov_ad * cos_p * sin_p, 0.0))
    err_t = np.sqrt(np.maximum((ea * sin_p) ** 2 + (ed * cos_p) ** 2 - 2 * cov_ad * cos_p * sin_p, 0.0))
    return MemberSample(r, mu_r, mu_t, err_r, err_t,
                        np.asarray(table["membership_prob"], float), np.asarray(table["g_mag"], float),
                        np.asarray(table["quality_flag"], int), np.arctan2(x, y),
                        mu_sys=(float(mu_sys[0]), float(mu_sys[1])), exact=exact)


def dispersion_ml(values: np.ndarray, errors: np.ndarray, err_scale: float = 1.0,
                  sigma_sys: float = 0.0, extra_var: np.ndarray | float = 0.0) -> tuple[float, float, float]:
    """Maximum-likelihood mean and intrinsic dispersion with per-star errors deconvolved.

    Model: ``value_i ~ N(mean, sigma^2 + (err_scale * err_i)^2 + sigma_sys^2 + extra_var_i)``;
    ``extra_var`` carries known per-star apparent-dispersion terms such as the line-of-sight
    depth effect of the systemic proper motion.

    Parameters
    ----------
    err_scale : float
        Multiplies every quoted error, to test their calibration.
    sigma_sys : float
        A floor added in quadrature (e.g. Gaia's spatially correlated systematics).

    Returns
    -------
    (sigma, sigma_error, mean)
        ``sigma`` is the intrinsic dispersion in the units of ``values``; ``sigma_error``
        comes from the curvature of the profile likelihood.
    """
    v = np.asarray(values, float); e2 = (err_scale * np.asarray(errors, float)) ** 2 + sigma_sys**2 + np.asarray(extra_var, float)
    n = len(v)
    if n < 5:
        return np.nan, np.nan, np.nan

    def nll(theta: np.ndarray) -> float:
        mean, ln_s2 = theta[0], theta[1]
        var = np.exp(ln_s2) + e2
        return 0.5 * float(np.sum((v - mean) ** 2 / var + np.log(var)))

    var0 = max(np.var(v) - np.mean(e2), 1e-6 * np.mean(e2))
    res = optimize.minimize(nll, [np.mean(v), np.log(var0)], method="Nelder-Mead",
                            options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 2000})
    mean, s2 = res.x[0], float(np.exp(res.x[1]))
    sigma = np.sqrt(max(s2, 0.0))
    # error on sigma from the second derivative of the profile likelihood in s2
    h = max(0.02 * s2, 1e-8)
    f0, fp, fm = nll(res.x), nll([mean, np.log(s2 + h)]), nll([mean, np.log(max(s2 - h, 1e-12))])
    d2 = (fp - 2 * f0 + fm) / (np.log((s2 + h) / max(s2 - h, 1e-12)) / 2) ** 2
    sig_ln_s2 = 1.0 / np.sqrt(max(d2, 1e-12))
    return sigma, 0.5 * sigma * sig_ln_s2, mean


def binned_dispersion(sample: MemberSample, edges: np.ndarray, prob_min: float = 0.9,
                      g_range: tuple[float, float] = (-np.inf, np.inf), err_scale: float = 1.0,
                      sigma_sys: float = 0.0, quality_mask: int | None = None,
                      depth_tracer: Any = None, distance_kpc: float = 5.43) -> Table:
    """Radial and tangential PM dispersion in each annulus, with diagnostics.

    ``depth_tracer`` (a mass component whose density is the tracer's) switches on the
    line-of-sight depth term: the systemic PM times the depth spread appears as a
    dispersion ``|mu_sys| sigma_z(R) / D`` along the systemic-PM direction, i.e. a per-star
    variance ``sigma_depth^2 cos^2(phi - phi_sys)`` in the radial and ``sin^2`` in the
    tangential component. It is 0.01-0.03 mas/yr for omega Cen and is removed here so the
    quoted dispersion is internal.

    Columns: ``r_lower``, ``r_median``, ``r_upper``, ``n_stars``, ``sigma_pmr``,
    ``sigma_pmr_err``, ``mean_pmr``, ``sigma_pmt``, ``sigma_pmt_err``, ``mean_pmt``,
    ``sigma_pm`` (the 1-D combined value), ``median_err``, ``expected_contaminants``
    (``sum(1 - P)``) and ``median_g``.
    """
    keep = (sample.prob >= prob_min) & (sample.g_mag >= g_range[0]) & (sample.g_mag <= g_range[1])
    if quality_mask is not None:
        keep &= (sample.quality_flag & quality_mask) > 0
    s = sample.select(keep)
    phi_sys = np.arctan2(s.mu_sys[0], s.mu_sys[1])
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (s.r_arcsec >= lo) & (s.r_arcsec < hi)
        if m.sum() < 5:
            continue
        xr = xt = 0.0
        if depth_tracer is not None:
            sd = depth_dispersion(depth_tracer, s.r_arcsec[m] * distance_kpc * 1e3 / 206264.806,
                                  float(np.hypot(*s.mu_sys)), distance_kpc)
            dphi = s.phi[m] - phi_sys
            xr, xt = sd**2 * np.cos(dphi) ** 2, sd**2 * np.sin(dphi) ** 2
        sr, esr, mr = dispersion_ml(s.mu_r[m], s.err_r[m], err_scale, sigma_sys, xr)
        st, est, mt = dispersion_ml(s.mu_t[m], s.err_t[m], err_scale, sigma_sys, xt)
        comb = np.sqrt(0.5 * (sr**2 + st**2))
        rows.append((lo, float(np.median(s.r_arcsec[m])), hi, int(m.sum()), sr, esr, mr, st, est, mt,
                     comb, float(np.median(0.5 * (s.err_r[m] + s.err_t[m]))),
                     float(np.sum(1.0 - s.prob[m])), float(np.median(s.g_mag[m]))))
    return Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_stars", "sigma_pmr", "sigma_pmr_err",
                                   "mean_pmr", "sigma_pmt", "sigma_pmt_err", "mean_pmt", "sigma_pm",
                                   "median_err", "expected_contaminants", "median_g"))


# ------------------------------------------------------------ contamination ---
def _two_gaussian_em(v: np.ndarray, n_iter: int = 60) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit a 2-component 1-D Gaussian mixture (weights, means, sigmas) by EM."""
    w = np.array([0.6, 0.4]); mu = np.array([np.median(v), np.median(v)]) + np.array([-1.0, 1.0]) * np.std(v) * 0.5
    sd = np.array([np.std(v), 2 * np.std(v)])
    for _ in range(n_iter):
        pdf = np.stack([w[k] * np.exp(-0.5 * ((v - mu[k]) / sd[k]) ** 2) / (sd[k] * np.sqrt(2 * np.pi)) for k in range(2)])
        resp = pdf / np.maximum(pdf.sum(axis=0), 1e-300)
        nk = resp.sum(axis=1)
        w = nk / len(v); mu = (resp * v).sum(axis=1) / nk
        sd = np.sqrt(np.maximum((resp * (v - mu[:, None]) ** 2).sum(axis=1) / nk, 1e-8))
    return w, mu, sd


def field_pdf_factory(field_values: np.ndarray):
    """Empirical PM distribution of the field stars in an annulus, as a callable pdf."""
    w, mu, sd = _two_gaussian_em(np.asarray(field_values, float))

    def pdf(v: np.ndarray) -> np.ndarray:
        v = np.asarray(v, float)
        return sum(w[k] * np.exp(-0.5 * ((v - mu[k]) / sd[k]) ** 2) / (sd[k] * np.sqrt(2 * np.pi)) for k in range(2))

    pdf.params = (w, mu, sd)
    return pdf


def mixture_dispersion(values: np.ndarray, errors: np.ndarray, field_pdf, extra_var=0.0
                       ) -> tuple[float, float, float, float, float]:
    """Cluster dispersion with contamination modelled instead of cut away.

    ``value_i ~ (1 - f) N(mean, sigma^2 + err_i^2 + extra_var_i) + f field_pdf(value_i)``, fitted to
    *all* stars (no membership cut), so neither field stars nor the truncated wings of the
    member distribution bias ``sigma``. Returns ``(sigma, sigma_err, mean, f, f_err)``.
    """
    v = np.asarray(values, float); e2 = np.asarray(errors, float) ** 2 + np.asarray(extra_var, float)
    fp = np.maximum(field_pdf(v), 1e-300)

    def nll(theta):
        mean, ln_s2, logit_f = theta
        f = 1.0 / (1.0 + np.exp(-logit_f))
        var = np.exp(ln_s2) + e2
        pc = np.exp(-0.5 * (v - mean) ** 2 / var) / np.sqrt(2 * np.pi * var)
        return -float(np.sum(np.log((1 - f) * pc + f * fp)))

    var0 = max(np.var(v[np.abs(v - np.median(v)) < 1.0]) - np.mean(e2), 1e-4)
    res = optimize.minimize(nll, [np.median(v), np.log(var0), np.log(0.02 / 0.98)], method="Nelder-Mead",
                            options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 4000})
    mean, s2, lf = res.x[0], float(np.exp(res.x[1])), res.x[2]
    f = 1.0 / (1.0 + np.exp(-lf))
    # errors from the diagonal curvature of the profile likelihood
    def curv(i, h):
        x0 = res.x.copy(); xp = x0.copy(); xm = x0.copy(); xp[i] += h; xm[i] -= h
        return (nll(xp) - 2 * nll(x0) + nll(xm)) / h**2
    sig_ln_s2 = 1.0 / np.sqrt(max(curv(1, 0.05), 1e-12)); sig_lf = 1.0 / np.sqrt(max(curv(2, 0.1), 1e-12))
    sigma = np.sqrt(s2)
    return sigma, 0.5 * sigma * sig_ln_s2, mean, f, f * (1 - f) * sig_lf


def contamination_audit(sample: MemberSample, edges: np.ndarray, quality_mask: int | None = 2,
                        field_prob_max: float = 0.05, depth_tracer: Any = None, distance_kpc: float = 5.43) -> Table:
    """Per annulus: the P-cut dispersions (P>0.9, P>0.99), the mixture-model dispersion fitted to
    all quality stars with the field distribution taken from the P<``field_prob_max`` stars, the
    fitted contaminant fraction, and the published-probability estimate ``sum(1-P)/N``."""
    q = np.ones(len(sample), bool) if quality_mask is None else (sample.quality_flag & quality_mask) > 0
    s = sample.select(q)
    phi_sys = np.arctan2(s.mu_sys[0], s.mu_sys[1])
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (s.r_arcsec >= lo) & (s.r_arcsec < hi)
        mem9 = m & (s.prob >= 0.9); mem99 = m & (s.prob >= 0.99); fld = m & (s.prob < field_prob_max)
        if mem9.sum() < 20 or fld.sum() < 50:
            continue
        xr = xt = 0.0
        if depth_tracer is not None:
            sd = depth_dispersion(depth_tracer, s.r_arcsec[m] * distance_kpc * 1e3 / 206264.806, float(np.hypot(*s.mu_sys)), distance_kpc)
            dphi = s.phi[m] - phi_sys; xr, xt = sd**2 * np.cos(dphi) ** 2, sd**2 * np.sin(dphi) ** 2
        out = {}
        for comp, val, err, xv in (("r", s.mu_r, s.err_r, xr), ("t", s.mu_t, s.err_t, xt)):
            fpdf = field_pdf_factory(val[fld])
            out[comp] = mixture_dispersion(val[m], err[m], fpdf, xv)
            out[comp + "_field"] = fpdf.params
        s9r, _, _ = dispersion_ml(s.mu_r[mem9], s.err_r[mem9]); s9t, _, _ = dispersion_ml(s.mu_t[mem9], s.err_t[mem9])
        s99r, _, _ = dispersion_ml(s.mu_r[mem99], s.err_r[mem99]); s99t, _, _ = dispersion_ml(s.mu_t[mem99], s.err_t[mem99])
        (sr, esr, _, fr, efr), (st, est, _, ft, eft) = out["r"], out["t"]
        rows.append((lo, float(np.median(s.r_arcsec[m])), hi, int(m.sum()), int(mem9.sum()), int(fld.sum()),
                     np.sqrt(0.5 * (s9r**2 + s9t**2)), np.sqrt(0.5 * (s99r**2 + s99t**2)),
                     np.sqrt(0.5 * (sr**2 + st**2)), 0.5 * np.hypot(esr, est),
                     0.5 * (fr + ft), 0.5 * np.hypot(efr, eft),
                     float(np.sum(1 - s.prob[mem9]) / mem9.sum()),
                     float(np.sqrt(np.sum(out["r_field"][0] * (out["r_field"][2] ** 2 + out["r_field"][1] ** 2))))))
    return Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_all", "n_p90", "n_field",
                                   "sigma_p90", "sigma_p99", "sigma_mixture", "sigma_mixture_err",
                                   "f_contam", "f_contam_err", "f_from_probabilities", "field_rms"))


def mixture_dispersion_free(values: np.ndarray, errors: np.ndarray, extra_var=0.0, n_iter: int = 200,
                            field_sigma_min: float = 1.5) -> dict[str, float]:
    """Cluster dispersion from a 3-component mixture with the field fitted freely.

    ``value_i ~ (1-f) N(mean, sigma^2 + err_i^2 + extra_var_i) + f [w N(m1, s1^2) + (1-w) N(m2, s2^2)]``
    by EM over *all* stars in the annulus: no membership probability is used anywhere, so the
    field density at the cluster's proper motion is an interpolation of the smooth broad
    components rather than a template with a hole where the members were removed. The field
    widths are floored at ``field_sigma_min`` (mas/yr): the Galactic field has sigma ~ 5-7
    mas/yr here, and without the floor EM places a narrow "field" component on the wings of
    the member distribution and biases ``sigma`` low (seen 2026-09-18). Returns
    ``sigma``, its error (profile-likelihood curvature), ``mean``, ``f`` (field fraction of the
    annulus) and the field components.
    """
    v = np.asarray(values, float); e2 = np.asarray(errors, float) ** 2 + np.asarray(extra_var, float)
    med = np.median(v); core = np.abs(v - med) < 1.0
    mean, s2 = float(np.median(v[core])), max(float(np.var(v[core]) - np.mean(e2[core])), 1e-4)
    f = 0.3; w = 0.6; m = np.array([med, med]); s = np.array([3.0, 8.0])
    for _ in range(n_iter):
        var = s2 + e2
        pc = (1 - f) * np.exp(-0.5 * (v - mean) ** 2 / var) / np.sqrt(2 * np.pi * var)
        p1 = f * w * np.exp(-0.5 * ((v - m[0]) / s[0]) ** 2) / (s[0] * np.sqrt(2 * np.pi))
        p2 = f * (1 - w) * np.exp(-0.5 * ((v - m[1]) / s[1]) ** 2) / (s[1] * np.sqrt(2 * np.pi))
        tot = np.maximum(pc + p1 + p2, 1e-300)
        rc, r1, r2 = pc / tot, p1 / tot, p2 / tot
        nc, n1, n2 = rc.sum(), r1.sum(), r2.sum()
        f = (n1 + n2) / len(v); w = n1 / max(n1 + n2, 1e-12)
        m = np.array([(r1 * v).sum() / n1, (r2 * v).sum() / n2])
        s = np.sqrt(np.maximum([(r1 * (v - m[0]) ** 2).sum() / n1, (r2 * (v - m[1]) ** 2).sum() / n2],
                               field_sigma_min**2))
        # cluster component: heteroscedastic, so its M-step is a weighted ML in (mean, s2)
        wts = rc / var
        mean = float((wts * v).sum() / wts.sum())
        # one Newton-like step on s2 from the weighted score, kept positive
        score = 0.5 * np.sum(rc * ((v - mean) ** 2 / var**2 - 1.0 / var))
        fisher = 0.5 * np.sum(rc / var**2)
        s2 = max(s2 + score / fisher, 1e-6)
    # curvature error on sigma at the converged responsibilities
    def nll_s2(ls2):
        var = np.exp(ls2) + e2
        return -float(np.sum(rc * (-0.5 * (v - mean) ** 2 / var - 0.5 * np.log(var))))
    h = 0.05; l0 = np.log(s2)
    d2 = (nll_s2(l0 + h) - 2 * nll_s2(l0) + nll_s2(l0 - h)) / h**2
    sigma = np.sqrt(s2)
    return {"sigma": sigma, "sigma_err": 0.5 * sigma / np.sqrt(max(d2, 1e-12)), "mean": mean, "f": f,
            "field": (w, m.copy(), s.copy()), "n_cluster": float(nc)}
