"""Poisson likelihood for binned star counts as the tracer-density constraint.

A count profile holds, per annulus, the number of selected stars ``N_i`` and the
effective area ``A_i`` (footprint coverage times annulus area, arcsec^2). The
model predicts the tracer surface density ``Sigma(R)`` up to an amplitude, so

    mu_i = A_i * (a * Sigma_bar_i + b),      N_i ~ Poisson(mu_i),

with ``Sigma_bar_i`` the area average of the model over the bin's effective
area, ``sum_j w_ij Sigma(R_ij)`` with nodes ``R_ij`` and covered-area weights
``w_ij`` (``node_weights``; rows sum to 1). Products without weights (schema 1,
equal-count quantile nodes of the stars' radii with an equal-weight mean) are
read for replay only: that mean is star-weighted, i.e. it estimates
int Sigma^2 dA / int Sigma dA instead of the area mean, and is biased high by
~3% per bin for Sigma ~ R^-2 at the HST bin width. ``a >= 0`` the number of stars per
unit model surface density (profiled analytically at every evaluation, or held
at a given value when the mass per counted star is an external assumption) and
``b >= 0`` an optional uniform field density (profiled when ``fit_field``).
The residual vector uses signed deviance residuals, whose squared sum is the
Poisson deviance ``D = 2 sum_i [N_i ln(N_i/mu_i) - (N_i - mu_i)]``
(= -2 ln L + const), so a least-squares fit maximizes the Poisson likelihood.
No error bars are adopted anywhere.

Completeness. A product may carry a per-bin detection completeness ``c_i``
(measured externally, e.g. Gaia against HST): mu_i = A_i * (a * c_i * Sigma_bar_i + b),
so ``a`` counts stars of a complete sample. The field term is not scaled.

Tied amplitudes. Two products of the same population (``CountTie``) can share one
amplitude: a_den = a_num / r, with r the measured ratio of complete-sample counts and
a Gaussian prior r ~ N(ratio, sigma). ``tied_compare`` maximizes the joint Poisson
likelihood plus that prior over (a_num, r, field of each product); the prior adds one
residual (r - ratio) / sigma to the objective.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["CountProfile", "CountTie", "poisson_profile_fit", "deviance_residuals", "tied_compare"]


def poisson_profile_fit(counts, area, shape, fit_field=True, iterations=60, amplitude=None):
    """Maximize the Poisson likelihood over a >= 0 (and b >= 0) for mu = A (a shape + b).

    Damped Newton iterations on the concave log-likelihood; the field term is
    fixed at zero unless ``fit_field``. With ``amplitude`` given, ``a`` is held
    at that value and only ``b`` is profiled (1D Newton, bounded below by 0).
    Returns (a, b, mu).
    """
    n = np.asarray(counts, float)
    A = np.asarray(area, float)
    s = np.asarray(shape, float)
    if np.any(s < 0) or np.any(A <= 0) or np.any(n < 0):
        raise ValueError("counts and areas must be non-negative, shape non-negative")
    if amplitude is not None:
        if not np.isfinite(amplitude) or amplitude <= 0:
            raise ValueError("fixed count amplitude must be positive")
        a = float(amplitude)
        b = max((n.sum()-a*(A*s).sum())/A.sum(), 0.) if fit_field else 0.
        if fit_field:
            for _ in range(iterations):
                mu = np.maximum(A*(a*s+b), 1e-300)
                gb = np.sum((n/mu-1.)*A)
                Hbb = np.sum(n/mu**2*A*A)
                b_new = max(b+gb/max(Hbb, 1e-300), 0.)
                if abs(b_new-b) <= 1e-12*max(abs(b), 1e-300):
                    b = b_new
                    break
                b = b_new
        return a, float(b), A*(a*s+b)
    a = max(n.sum()/max((A*s).sum(), 1e-300), 1e-300)
    b = 0.
    for _ in range(iterations):
        mu = A*(a*s+b)
        if np.any(mu <= 0):
            mu = np.maximum(mu, 1e-300)
        r = n/mu-1.                       # d lnL / d mu_i
        w = n/mu**2                       # -d2 lnL / d mu_i^2
        ga = np.sum(r*A*s)
        Haa = np.sum(w*(A*s)**2)
        if fit_field:
            gb = np.sum(r*A)
            Hab = np.sum(w*A*A*s)
            Hbb = np.sum(w*A*A)
            H = np.array([[Haa, Hab], [Hab, Hbb]])
            g = np.array([ga, gb])
            try:
                step = np.linalg.solve(H+1e-12*np.eye(2), g)
            except np.linalg.LinAlgError:
                step = g/np.maximum(np.diag(H), 1e-300)
            a_new, b_new = a+step[0], b+step[1]
            if b_new < 0:                 # active bound: refit a alone at b = 0
                b_new = 0.
                a_new = a+ga/max(Haa, 1e-300) if step[0] != 0 else a
        else:
            a_new, b_new = a+ga/max(Haa, 1e-300), 0.
        a_new = max(a_new, 1e-300*a if a > 0 else 1e-300)
        if abs(a_new-a) <= 1e-12*abs(a) and abs(b_new-b) <= 1e-12*max(abs(b), 1e-300):
            a, b = a_new, b_new
            break
        a, b = a_new, b_new
    mu = A*(a*s+b)
    return float(a), float(b), mu


def deviance_residuals(counts, mu):
    """Signed Poisson deviance residuals; their squared sum is the deviance."""
    n = np.asarray(counts, float)
    mu = np.asarray(mu, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(n > 0, n*np.log(n/mu), 0.)-(n-mu)
    term = np.maximum(term, 0.)           # numerical noise near n = mu
    return np.sign(n-mu)*np.sqrt(2.*term)


@dataclass(frozen=True)
class CountProfile:
    """Binned star counts of one tracer selection with per-bin effective areas."""
    name: str
    r_lower: np.ndarray            # arcsec
    r_upper: np.ndarray
    r_nodes: np.ndarray            # (n_bins, k) quadrature radii of each bin, arcsec
    counts: np.ndarray             # integers
    area_arcsec2: np.ndarray       # coverage-corrected effective area per bin
    source: str
    selection: str
    fit_field: bool = True
    note: str = ""
    meta: dict = field(default_factory=dict)
    node_weights: np.ndarray | None = None   # (n_bins, k) covered-area fractions; None = legacy equal weights
    completeness: np.ndarray | None = None   # per-bin detection completeness in (0, 1]; None = 1

    def __post_init__(self):
        arrays = (self.r_lower, self.r_upper, self.counts, self.area_arcsec2)
        n = len(self.counts)
        if n < 2 or any(len(np.asarray(a)) != n for a in arrays) or np.asarray(self.r_nodes).shape[0] != n:
            raise ValueError(f"{self.name}: inconsistent count-profile arrays")
        if np.any(np.asarray(self.r_upper) <= np.asarray(self.r_lower)) or np.any(np.asarray(self.r_lower) < 0):
            raise ValueError(f"{self.name}: bin edges must be positive and increasing")
        if np.any(np.asarray(self.counts) < 0) or np.any(np.asarray(self.counts) != np.round(self.counts)):
            raise ValueError(f"{self.name}: counts must be non-negative integers")
        if np.any(np.asarray(self.area_arcsec2) <= 0) or np.any(~np.isfinite(self.r_nodes)):
            raise ValueError(f"{self.name}: areas must be positive and nodes finite")
        if self.node_weights is not None:
            w = np.asarray(self.node_weights, float)
            if w.shape != np.asarray(self.r_nodes).shape or np.any(~np.isfinite(w)) or np.any(w < 0) \
                    or not np.allclose(w.sum(axis=1), 1., rtol=0, atol=1e-10):
                raise ValueError(f"{self.name}: node weights must match the nodes, be >= 0 and sum to 1 per bin")
        if self.completeness is not None:
            c = np.asarray(self.completeness, float)
            if c.shape != (n,) or np.any(~np.isfinite(c)) or np.any(c <= 0) or np.any(c > 1):
                raise ValueError(f"{self.name}: completeness must be one value in (0, 1] per bin")

    @property
    def n(self):
        return len(self.counts)

    @property
    def r_median(self):
        return np.median(self.r_nodes, axis=1)

    def bin_average(self, sigma_nodes):
        """Model surface density averaged over each bin (area-weighted if weights are stored)."""
        sigma_nodes = np.asarray(sigma_nodes, float).reshape(self.r_nodes.shape)
        if self.node_weights is None:
            return sigma_nodes.mean(axis=1)
        return np.sum(np.asarray(self.node_weights, float)*sigma_nodes, axis=1)

    def model_shape(self, sigma_nodes):
        """Bin-averaged model surface density times the detection completeness."""
        sigma_nodes = np.asarray(sigma_nodes, float).reshape(self.r_nodes.shape)
        if np.any(~np.isfinite(sigma_nodes)) or np.any(sigma_nodes <= 0):
            raise ValueError(f"{self.name}: invalid projected surface density")
        shape = self.bin_average(sigma_nodes)
        return shape if self.completeness is None else shape*np.asarray(self.completeness, float)

    def terms(self, mu, a, b, amplitude_fixed):
        """Deviance residuals and log-likelihood for the predicted counts ``mu``."""
        res = deviance_residuals(self.counts, mu)
        n = np.asarray(self.counts, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            loglike = float(np.sum(np.where(n > 0, n*np.log(mu), 0.)-mu-_lgamma1(n)))
        return dict(deviance=float(res @ res), residual=res, mu=mu, amplitude=a, field=b,
                    amplitude_fixed=amplitude_fixed,
                    field_density_per_arcmin2=b*3600., loglike=loglike, n=self.n)

    def compare(self, sigma_nodes, amplitude=None):
        """Poisson fit given the model surface density at ``r_nodes`` (same shape).

        ``amplitude`` (stars per unit surface density per arcsec^2) holds ``a``
        fixed instead of profiling it; the field is still profiled if enabled.
        """
        shape = self.model_shape(sigma_nodes)
        a, b, mu = poisson_profile_fit(self.counts, self.area_arcsec2, shape, self.fit_field, amplitude=amplitude)
        return self.terms(mu, a, b, amplitude is not None)

    def to_dict(self):
        row = dict(schema_version=1, name=self.name, r_lower=self.r_lower.tolist(),
                   r_upper=self.r_upper.tolist(), r_nodes=np.asarray(self.r_nodes).tolist(),
                   counts=[int(c) for c in self.counts], area_arcsec2=self.area_arcsec2.tolist(),
                   source=self.source, selection=self.selection, fit_field=self.fit_field,
                   note=self.note, meta=dict(self.meta))
        if self.node_weights is not None:
            row.update(schema_version=2, node_weights=np.asarray(self.node_weights).tolist())
        if self.completeness is not None:
            row["completeness"] = np.asarray(self.completeness, float).tolist()
        return row

    @classmethod
    def from_dict(cls, row):
        row = dict(row)
        version = row.pop("schema_version")
        if version not in (1, 2):
            raise ValueError("unsupported count-profile snapshot")
        keys = ("r_lower", "r_upper", "r_nodes", "area_arcsec2")+(("node_weights",) if version == 2 else ())
        arrays = {k: np.array(row.pop(k), float) for k in keys}
        counts = np.array(row.pop("counts"), int)
        if row.get("completeness") is not None:
            row["completeness"] = np.array(row["completeness"], float)
        return cls(counts=counts, **arrays, **row)


@dataclass(frozen=True)
class CountTie:
    """a_denominator = a_numerator / r for complete samples, with r ~ N(ratio, sigma)."""
    numerator: str
    denominator: str
    ratio: float
    sigma: float
    note: str = ""

    def __post_init__(self):
        if self.numerator == self.denominator or not (self.ratio > 0 and self.sigma > 0):
            raise ValueError("count tie needs two different products, ratio > 0 and sigma > 0")

    def to_dict(self):
        return dict(numerator=self.numerator, denominator=self.denominator, ratio=self.ratio,
                    sigma=self.sigma, note=self.note)

    @classmethod
    def from_dict(cls, row):
        return cls(**row)


def _joint_ab(n, A, S, field_rows, iterations=80):
    """Maximize the Poisson likelihood of mu = A * (a S + b * field_rows) over a > 0, b >= 0."""
    a = max(n.sum()/max((A*S).sum(), 1e-300), 1e-300)
    b = 0.
    F = field_rows.astype(float)
    use_b = bool(F.any())
    for _ in range(iterations):
        mu = np.maximum(A*(a*S+b*F), 1e-300)
        r, w = n/mu-1., n/mu**2
        ga, Haa = np.sum(r*A*S), np.sum(w*(A*S)**2)
        if use_b:
            gb, Hab, Hbb = np.sum(r*A*F), np.sum(w*A*A*S*F), np.sum(w*(A*F)**2)
            step = np.linalg.solve(np.array([[Haa, Hab], [Hab, Hbb]])+1e-12*np.eye(2), np.array([ga, gb]))
            a_new, b_new = a+step[0], b+step[1]
            if b_new < 0:
                b_new, a_new = 0., a+ga/max(Haa, 1e-300)
        else:
            a_new, b_new = a+ga/max(Haa, 1e-300), 0.
        a_new = max(a_new, 1e-3*a)
        done = abs(a_new-a) <= 1e-12*a and abs(b_new-b) <= 1e-12*max(b, 1e-300)
        a, b = a_new, b_new
        if done:
            break
    return a, b


def tied_compare(numerator: CountProfile, denominator: CountProfile, sigma_num, sigma_den, tie: CountTie):
    """Joint Poisson fit of two products sharing one amplitude (see module docstring).

    The numerator product must have no field term. Returns (terms_num, terms_den, tie_term);
    tie_term holds the fitted ratio r, its prior residual and chi2 = residual^2.
    """
    from scipy.optimize import minimize_scalar
    if numerator.fit_field:
        raise ValueError("the numerator product of a count tie must not fit a field term")
    s_num, s_den = numerator.model_shape(sigma_num), denominator.model_shape(sigma_den)
    n = np.concatenate([numerator.counts, denominator.counts]).astype(float)
    A = np.concatenate([numerator.area_arcsec2, denominator.area_arcsec2]).astype(float)
    F = np.concatenate([np.zeros(numerator.n, bool), np.full(denominator.n, denominator.fit_field)])

    def solve(r):
        S = np.concatenate([s_num, s_den/r])
        a, b = _joint_ab(n, A, S, F)
        mu = A*(a*S+b*F)
        dev = float(np.sum(deviance_residuals(n, mu)**2))
        return dev+((r-tie.ratio)/tie.sigma)**2, a, b, mu

    lo, hi = max(tie.ratio-8*tie.sigma, 1e-6*tie.ratio), tie.ratio+8*tie.sigma
    best = minimize_scalar(lambda r: solve(r)[0], bounds=(lo, hi), method="bounded",
                           options=dict(xatol=1e-7*tie.ratio))
    r = float(best.x)
    _, a, b, mu = solve(r)
    t_num = numerator.terms(mu[:numerator.n], a, 0., False)
    t_den = denominator.terms(mu[numerator.n:], a/r, b, False)
    res = (r-tie.ratio)/tie.sigma
    tie_term = dict(ratio=r, prior_ratio=tie.ratio, prior_sigma=tie.sigma, residual=np.array([res]),
                    chi2=float(res**2), numerator=tie.numerator, denominator=tie.denominator,
                    at_bound=bool(min(r-lo, hi-r) < 1e-3*tie.sigma))
    t_num["tied"] = t_den["tied"] = True
    return t_num, t_den, tie_term


def _lgamma1(n):
    from scipy.special import gammaln
    return gammaln(np.asarray(n, float)+1.)
