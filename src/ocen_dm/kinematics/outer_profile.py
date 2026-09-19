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
           "mixture_dispersion", "mixture_dispersion_free", "dispersion_with_field_template", "dispersion_2d",
           "field_pdf_factory", "contamination_audit",
           "OCEN_RA", "OCEN_DEC", "OCEN_VSYS_KMS", "KMS_PER_MASYR_KPC"]

#: Reproduce the two estimator defects fixed on 2026-09-19, for the before/after audit only.
#: ``"signs"`` restores the wrong sign on the covariance cross-terms in the mean update;
#: ``"interval"`` restores the fixed 2 per cent step that floored every uncertainty. Never
#: set these outside :mod:`ocen_dm.kinematics.estimator_audit`.
LEGACY: set[str] = set()

from ..cluster import OCEN_DEC, OCEN_RA  # single definition, derived from the data
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
    mu_a: np.ndarray = None               # equatorial residuals after the systemic field
    mu_d: np.ndarray = None
    err_a: np.ndarray = None              # per-star error covariance in the equatorial frame
    err_d: np.ndarray = None
    err_corr: np.ndarray = None
    sys_a: np.ndarray = None              # the systemic field subtracted at each star's position
    sys_d: np.ndarray = None
    mu_sys: tuple[float, float] = (0.0, 0.0)
    exact: bool = True                    # systemic field projected star by star, or a constant

    def __len__(self) -> int:
        return len(self.r_arcsec)

    @property
    def absolute_pm(self) -> tuple[np.ndarray, np.ndarray]:
        """Observed (absolute) proper motions: residual plus the systemic field."""
        return self.mu_a + self.sys_a, self.mu_d + self.sys_d

    def scale_errors(self, factor: np.ndarray | float) -> "MemberSample":
        """Copy with every per-star uncertainty multiplied by ``factor``.

        Used to switch between the raw Gaia uncertainties as released and the
        density-dependent inflation of Vasiliev & Baumgardt (2021); see
        :mod:`ocen_dm.kinematics.vb2021_replication`.
        """
        f = np.asarray(factor, float)
        fields = ("r_arcsec", "mu_r", "mu_t", "err_r", "err_t", "prob", "g_mag", "quality_flag", "phi",
                  "mu_a", "mu_d", "err_a", "err_d", "err_corr", "sys_a", "sys_d")
        scaled = {"err_r", "err_t", "err_a", "err_d"}
        vals = [None if getattr(self, n) is None else
                (np.asarray(getattr(self, n)) * f if n in scaled else np.asarray(getattr(self, n)))
                for n in fields]
        return MemberSample(*vals, mu_sys=self.mu_sys, exact=self.exact)

    def select(self, mask: np.ndarray) -> "MemberSample":
        fields = ("r_arcsec", "mu_r", "mu_t", "err_r", "err_t", "prob", "g_mag", "quality_flag", "phi",
                  "mu_a", "mu_d", "err_a", "err_d", "err_corr", "sys_a", "sys_d")
        vals = [None if getattr(self, f) is None else np.asarray(getattr(self, f))[mask] for f in fields]
        return MemberSample(*vals, mu_sys=self.mu_sys, exact=self.exact)


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
                        mu_a=pmra, mu_d=pmdec, err_a=ea, err_d=ed, err_corr=rho,
                        sys_a=np.broadcast_to(exp_a, pmra.shape).copy(), sys_d=np.broadcast_to(exp_d, pmdec.shape).copy(),
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


def dispersion_with_field_template(values: np.ndarray, errors: np.ndarray, field_pdf,
                                   extra_var=0.0, sigma_grid: np.ndarray | None = None) -> dict[str, float]:
    """Cluster dispersion with the field taken from an **independent** sample.

    ``value_i ~ (1-f) N(mean, sigma^2 + err_i^2 + extra_var_i) + f field_pdf(value_i)``, where
    ``field_pdf`` is measured outside the cluster (see
    :mod:`ocen_dm.selection.field_template`) and only its normalisation ``f`` is free. This
    is the treatment to trust where the field dominates: a field fitted to the same stars
    can absorb the cluster peak, and does (JOURNAL 2026-09-18).

    Solved as a profile likelihood on a grid in ``sigma`` with ``(mean, f)`` maximised by EM
    at each node -- robust where the cluster is a per-cent-level component and a simplex
    search on three parameters wanders (found 2026-09-18 in the outermost annulus).
    """
    v = np.asarray(values, float); e2 = np.asarray(errors, float) ** 2 + np.asarray(extra_var, float)
    fp = np.maximum(np.asarray(field_pdf(v), float), 1e-300)
    coarse = np.geomspace(0.03, 1.2, 120) if sigma_grid is None else np.asarray(sigma_grid, float)

    def profile(sigma: float) -> tuple[float, float, float]:
        """Maximised log-likelihood at fixed sigma, with (mean, f) by EM."""
        var = sigma**2 + e2
        mean, f = float(np.median(v)), 0.5
        ll = -np.inf
        for _ in range(200):
            pc = (1 - f) * np.exp(-0.5 * (v - mean) ** 2 / var) / np.sqrt(2 * np.pi * var)
            pf = f * fp
            tot = np.maximum(pc + pf, 1e-300)
            r = pc / tot
            new_ll = float(np.sum(np.log(tot)))
            f = 1.0 - r.sum() / len(v)
            w = r / var
            mean = float((w * v).sum() / max(w.sum(), 1e-300))
            if abs(new_ll - ll) < 1e-8:
                ll = new_ll
                break
            ll = new_ll
        return ll, mean, f

    ll_coarse = np.array([profile(sg)[0] for sg in coarse])
    j0 = int(np.argmax(ll_coarse))
    if sigma_grid is None:
        # refine around the maximum so that the 1-sigma interval is not quantised by the
        # coarse grid step (3 per cent in sigma, which is larger than the statistical error
        # for 10^4 stars -- found 2026-09-18)
        lo = coarse[max(j0 - 2, 0)]; hi = coarse[min(j0 + 2, len(coarse) - 1)]
        grid = np.geomspace(lo, hi, 161)
    else:
        grid = coarse
    # two stages: locate the maximum on a coarse grid, then refine, so that the 1-sigma
    # interval is not quantised by the 3 per cent grid step (larger than the statistical
    # error for 10^4 stars -- found 2026-09-18)
    ll_coarse = np.array([profile(sg)[0] for sg in coarse])
    j0 = int(np.argmax(ll_coarse))
    grid = (np.geomspace(coarse[max(j0 - 3, 0)], coarse[min(j0 + 3, len(coarse) - 1)], 241)
            if sigma_grid is None else coarse)
    out = np.array([profile(sg) for sg in grid])
    ll = out[:, 0]
    i = int(np.argmax(ll))
    sigma, mean, f = float(grid[i]), float(out[i, 1]), float(out[i, 2])
    # 1-sigma interval from the profile likelihood (Delta ln L = 1/2), interpolated
    target = ll[i] - 0.5
    def cross(side):
        idx = range(i, len(grid) - 1) if side > 0 else range(i, 0, -1)
        for j in idx:
            k = j + side
            if ll[k] < target:
                t = (target - ll[j]) / (ll[k] - ll[j]) if ll[k] != ll[j] else 0.0
                return float(grid[j] + t * (grid[k] - grid[j]))
        return float(grid[-1] if side > 0 else grid[0])
    lo_s, hi_s = cross(-1), cross(+1)
    if lo_s <= grid[0] * 1.001 or hi_s >= grid[-1] * 0.999:      # interval wider than the refined window
        out2 = np.array([profile(sg) for sg in coarse]); ll = out2[:, 0]; grid = coarse
        i = int(np.argmax(ll)); target = ll[i] - 0.5; lo_s, hi_s = cross(-1), cross(+1)
    return {"sigma": sigma, "sigma_err": 0.5 * (hi_s - lo_s), "sigma_lo": lo_s, "sigma_hi": hi_s,
            "mean": mean, "f": f, "n_cluster": float((1 - f) * len(v)), "lnl": float(ll[i]),
            "at_grid_edge": bool(i == 0 or i == len(grid) - 1)}


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


# ---------------------------------------------------------------- 2-D fit ---
def dispersion_2d(sample: MemberSample, mask: np.ndarray, field_density, depth_var: np.ndarray | float = 0.0,
                  sigma_grid: np.ndarray | None = None, n_iter: int = 300, sigma_max: float = 1.2,
                  field_at: tuple[np.ndarray, np.ndarray] | None = None) -> dict[str, float]:
    """Cluster dispersion from the **two-dimensional** proper-motion distribution.

    The field is not projected: each star is scored against a two-dimensional empirical
    density ``field_density(mu_alpha*, mu_delta)`` measured outside the cluster, in the same
    equatorial frame in which the systemic motion has been removed. Projecting onto the
    radial direction mixes the field's two unequal widths (3.1 and 2.0 mas/yr for omega Cen)
    around the annulus and manufactures non-Gaussian structure that no one-dimensional model
    fits; in two dimensions that problem does not arise.

    The cluster is a two-dimensional Gaussian whose covariance is
    ``R(phi_i) diag(sigma_R^2 + depth_R^2, sigma_T^2 + depth_T^2) R(phi_i)^T`` plus the
    star's own error covariance -- so the anisotropy and the rotation are fitted in the
    frame where they are defined, star by star.

    Returns ``sigma_r``, ``sigma_t``, their 1-sigma profile-likelihood intervals, the mean
    radial and tangential motions, the field fraction ``f`` and the cluster count.
    """
    s = sample.select(mask)
    a, d = np.asarray(s.mu_a, float), np.asarray(s.mu_d, float)
    ea, ed, rho = np.asarray(s.err_a, float), np.asarray(s.err_d, float), np.asarray(s.err_corr, float)
    phi = np.asarray(s.phi, float)
    cos_p, sin_p = np.sin(phi), np.cos(phi)          # radial unit vector in (alpha*, delta)
    # The cluster is modelled in residual proper motion (systemic field removed), but the
    # FIELD has no systemic motion of its own: its distribution is position-independent in
    # *absolute* proper motion. ``field_at`` therefore gives the absolute proper motions at
    # which to score the field template, which is itself built in absolute proper motion.
    fa, fd = (a, d) if field_at is None else (np.asarray(field_at[0], float)[mask],
                                              np.asarray(field_at[1], float)[mask])
    fp = np.maximum(np.asarray(field_density(fa, fd), float), 1e-300)
    dv = np.broadcast_to(np.asarray(depth_var, float), a.shape)

    # the line-of-sight depth term acts along the systemic proper motion, not along the
    # radial direction: it enters as a rank-1 covariance along that fixed direction
    mu_hat = np.asarray(sample.mu_sys, float)
    mu_hat = mu_hat / max(float(np.hypot(*mu_hat)), 1e-12)

    def loglike(sr2: float, st2: float) -> tuple[float, float, float, float]:
        """Maximised log-likelihood at fixed (sigma_R^2, sigma_T^2); means and f by EM."""
        # cluster covariance per star, in the equatorial frame
        caa = sr2 * cos_p**2 + st2 * sin_p**2 + ea**2 + dv * mu_hat[0] ** 2
        cdd = sr2 * sin_p**2 + st2 * cos_p**2 + ed**2 + dv * mu_hat[1] ** 2
        cad = (sr2 - st2) * cos_p * sin_p + rho * ea * ed + dv * mu_hat[0] * mu_hat[1]
        det = caa * cdd - cad**2
        mr, mt, f = 0.0, 0.0, 0.5
        ll = -np.inf
        for _ in range(n_iter):
            ma = mr * cos_p + mt * (-sin_p); md = mr * sin_p + mt * cos_p
            da, dd = a - ma, d - md
            q = (cdd * da**2 - 2 * cad * da * dd + caa * dd**2) / det
            pc = (1 - f) * np.exp(-0.5 * q) / (2 * np.pi * np.sqrt(det))
            tot = np.maximum(pc + f * fp, 1e-300)
            r = pc / tot
            new_ll = float(np.sum(np.log(tot)))
            f = 1.0 - r.sum() / len(a)
            # Weighted least squares for the mean in the rotated frame. With
            # J = [[cos, -sin], [sin, cos]] and Sigma^-1 = (1/det) [[cdd, -cad], [-cad, caa]],
            # the normal equations are (sum r J^T Sigma^-1 J) u = sum r J^T Sigma^-1 x. The
            # off-diagonal cad terms in A12 and b2 carry a MINUS sign; they were positive
            # until 2026-09-19, which biased the fitted tangential mean whenever the cluster
            # was anisotropic (synthetic test: -0.125 recovered for a true -0.200).
            w = r / det
            A11 = np.sum(w * (cdd * cos_p**2 - 2 * cad * cos_p * sin_p + caa * sin_p**2))
            A22 = np.sum(w * (cdd * sin_p**2 + 2 * cad * cos_p * sin_p + caa * cos_p**2))
            sgn = +1.0 if "signs" in LEGACY else -1.0
            A12 = np.sum(w * (-cdd * cos_p * sin_p + sgn * cad * (cos_p**2 - sin_p**2) + caa * sin_p * cos_p))
            b1 = np.sum(w * (cdd * a * cos_p - cad * (a * sin_p + d * cos_p) + caa * d * sin_p))
            b2 = np.sum(w * (-cdd * a * sin_p + sgn * cad * (a * cos_p - d * sin_p) + caa * d * cos_p))
            det_A = A11 * A22 - A12**2
            if abs(det_A) > 1e-30:
                mr = (b1 * A22 - b2 * A12) / det_A
                mt = (b2 * A11 - b1 * A12) / det_A
            if abs(new_ll - ll) < 1e-7:
                ll = new_ll
                break
            ll = new_ll
        return ll, mr, mt, f

    # maximise over (sigma_R, sigma_T) with a simplex on the log variances -- a grid over both
    # would cost thousands of EM solves on 3 x 10^4 stars -- then scan each axis for the interval
    from scipy.optimize import minimize

    if sigma_grid is not None:
        grid = np.asarray(sigma_grid, float)
        tab = np.array([[loglike(sr**2, st**2)[0] for st in grid] for sr in grid])
        i, j = np.unravel_index(np.argmax(tab), tab.shape)
        sr_hat, st_hat = float(grid[i]), float(grid[j])
    else:
        # multi-start, bounded to cluster-like dispersions: started from the sample variance the
        # simplex runs away to the field solution once the cluster is a few per cent of the
        # stars (the outermost annulus, found 2026-09-18)
        best = (-np.inf, sigma_max / 3, sigma_max / 3)
        for s0 in (0.15, 0.25, 0.40, 0.60):
            res = minimize(lambda p: -loglike(min(np.exp(p[0]), sigma_max**2), min(np.exp(p[1]), sigma_max**2))[0],
                           [np.log(s0**2)] * 2, method="Nelder-Mead",
                           options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 300})
            sr, st = float(np.sqrt(min(np.exp(res.x[0]), sigma_max**2))), float(np.sqrt(min(np.exp(res.x[1]), sigma_max**2)))
            if -res.fun > best[0] and max(sr, st) < 0.98 * sigma_max:
                best = (-res.fun, sr, st)
        sr_hat, st_hat = best[1], best[2]
    ll_max, mr, mt, f = loglike(sr_hat**2, st_hat**2)

    def interval(which: str) -> tuple[float, float]:
        """1-sigma profile-likelihood interval on one dispersion, the other re-maximised.

        The crossing of ``ln L = ln L_max - 1/2`` is bracketed with a geometrically growing
        step and then **interpolated**. Stepping by a fixed 2 per cent of sigma and returning
        the first point past the crossing, as this did until 2026-09-19, imposed a floor of
        2 per cent per component (1.4 per cent combined) on every reported uncertainty
        regardless of sample size: every HST bin with more than 10^5 stars sat exactly on it.
        """
        def profile(x: float) -> float:
            if which == "r":
                return -minimize(lambda p: -loglike(x**2, np.exp(p[0]))[0], [np.log(st_hat**2)],
                                 method="Nelder-Mead",
                                 options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 60}).fun
            return -minimize(lambda p: -loglike(np.exp(p[0]), x**2)[0], [np.log(sr_hat**2)],
                             method="Nelder-Mead",
                             options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 60}).fun

        x0 = sr_hat if which == "r" else st_hat
        if "interval" in LEGACY:                 # the floored version, audit only
            out = []
            for side in (-1, +1):
                x = x0
                for _ in range(60):
                    x = x + 0.02 * x0 * side
                    if x <= 0.01 or ll_max - profile(x) >= 0.5:
                        break
                out.append(x)
            return out[0], out[1]
        out = []
        for side in (-1, +1):
            step = 0.002 * x0
            x_in, d_in = x0, 0.0                      # last point inside the interval
            x_out = None
            for _ in range(80):
                x = x_in + side * step
                if x <= 0.01:
                    x_out, d_out = 0.01, ll_max - profile(0.01)
                    break
                d = ll_max - profile(x)
                if d >= 0.5:
                    x_out, d_out = x, d
                    break
                x_in, d_in, step = x, d, step * 1.6
            if x_out is None:                          # never crossed: report the last point
                out.append(x_in)
                continue
            # linear interpolation in the log-likelihood drop between the bracketing points
            if d_out > d_in:
                out.append(x_in + (x_out - x_in) * (0.5 - d_in) / (d_out - d_in))
            else:
                out.append(x_out)
        return out[0], out[1]

    lo_r, hi_r = interval("r"); lo_t, hi_t = interval("t")
    sigma = np.sqrt(0.5 * (sr_hat**2 + st_hat**2))
    return {"sigma_r": sr_hat, "sigma_r_err": 0.5 * (hi_r - lo_r), "sigma_t": st_hat,
            "sigma_t_err": 0.5 * (hi_t - lo_t), "sigma": sigma,
            "sigma_err": float(0.5 * np.hypot(0.5 * (hi_r - lo_r), 0.5 * (hi_t - lo_t))),
            "mean_r": mr, "mean_t": mt, "f": f, "n_cluster": float((1 - f) * len(a)), "lnl": ll_max,
            "n_stars": int(len(a)), "at_bound": bool(max(sr_hat, st_hat) > 0.9 * sigma_max)}
