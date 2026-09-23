"""Poisson likelihood for binned star counts as the tracer-density constraint.

A count profile holds, per annulus, the number of selected stars ``N_i`` and the
effective area ``A_i`` (footprint coverage times annulus area, arcsec^2). The
model predicts the tracer surface density ``Sigma(R)`` up to an amplitude, so

    mu_i = A_i * (a * Sigma_bar_i + b),      N_i ~ Poisson(mu_i),

with ``Sigma_bar_i`` the model averaged over the stars' radii (equal-count
quantile nodes, as for the kinematic bins), ``a >= 0`` the number of stars per
unit model surface density (profiled analytically at every evaluation) and
``b >= 0`` an optional uniform field density (profiled when ``fit_field``).
The residual vector uses signed deviance residuals, whose squared sum is the
Poisson deviance ``D = 2 sum_i [N_i ln(N_i/mu_i) - (N_i - mu_i)]``
(= -2 ln L + const), so a least-squares fit maximizes the Poisson likelihood.
No error bars are adopted anywhere.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["CountProfile", "poisson_profile_fit", "deviance_residuals"]


def poisson_profile_fit(counts, area, shape, fit_field=True, iterations=60):
    """Maximize the Poisson likelihood over a >= 0 (and b >= 0) for mu = A (a shape + b).

    Damped Newton iterations on the concave log-likelihood; the field term is
    fixed at zero unless ``fit_field``. Returns (a, b, mu).
    """
    n = np.asarray(counts, float)
    A = np.asarray(area, float)
    s = np.asarray(shape, float)
    if np.any(s < 0) or np.any(A <= 0) or np.any(n < 0):
        raise ValueError("counts and areas must be non-negative, shape non-negative")
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
    r_nodes: np.ndarray            # (n_bins, k) equal-count quantile radii of the counted stars, arcsec
    counts: np.ndarray             # integers
    area_arcsec2: np.ndarray       # coverage-corrected effective area per bin
    source: str
    selection: str
    fit_field: bool = True
    note: str = ""
    meta: dict = field(default_factory=dict)

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

    @property
    def n(self):
        return len(self.counts)

    @property
    def r_median(self):
        return np.median(self.r_nodes, axis=1)

    def compare(self, sigma_nodes):
        """Poisson fit given the model surface density at ``r_nodes`` (same shape)."""
        sigma_nodes = np.asarray(sigma_nodes, float).reshape(self.r_nodes.shape)
        if np.any(~np.isfinite(sigma_nodes)) or np.any(sigma_nodes <= 0):
            raise ValueError(f"{self.name}: invalid projected surface density")
        shape = sigma_nodes.mean(axis=1)
        a, b, mu = poisson_profile_fit(self.counts, self.area_arcsec2, shape, self.fit_field)
        res = deviance_residuals(self.counts, mu)
        n = np.asarray(self.counts, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            loglike = float(np.sum(np.where(n > 0, n*np.log(mu), 0.)-mu-_lgamma1(n)))
        return dict(deviance=float(res @ res), residual=res, mu=mu, amplitude=a, field=b,
                    field_density_per_arcmin2=b*3600., loglike=loglike, n=self.n)

    def to_dict(self):
        return dict(schema_version=1, name=self.name, r_lower=self.r_lower.tolist(),
                    r_upper=self.r_upper.tolist(), r_nodes=np.asarray(self.r_nodes).tolist(),
                    counts=[int(c) for c in self.counts], area_arcsec2=self.area_arcsec2.tolist(),
                    source=self.source, selection=self.selection, fit_field=self.fit_field,
                    note=self.note, meta=dict(self.meta))

    @classmethod
    def from_dict(cls, row):
        row = dict(row)
        if row.pop("schema_version") != 1:
            raise ValueError("unsupported count-profile snapshot")
        arrays = {k: np.array(row.pop(k), float) for k in ("r_lower", "r_upper", "r_nodes", "area_arcsec2")}
        counts = np.array(row.pop("counts"), int)
        return cls(counts=counts, **arrays, **row)


def _lgamma1(n):
    from scipy.special import gammaln
    return gammaln(np.asarray(n, float)+1.)
