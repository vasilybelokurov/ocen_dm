"""Before and after the two estimator defects fixed on 2026-09-19.

Both were found by an external code review and verified independently before anything was
changed (JOURNAL 2026-09-19):

* **signs** -- the covariance cross-terms in the mean-update normal equations carried a plus
  where the derivation gives a minus. They vanish for an isotropic cluster with uncorrelated
  errors, which every test mock was, so nothing caught it. It biased the fitted streaming and
  the anisotropy whenever the cluster was anisotropic.
* **interval** -- the profile-likelihood interval stepped by a fixed 2 per cent of sigma and
  returned the first point past the crossing, flooring every reported uncertainty at 2 per
  cent per component regardless of sample size.

This module reruns the measurement with either or both restored, purely so the effect can be
plotted. Nothing here is used by the science path.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from . import outer_profile as _op

__all__ = ["with_legacy", "gaia_audit", "hst_audit"]


class with_legacy:
    """Context manager restoring the named defects inside its block."""

    def __init__(self, *names: str):
        self.names = set(names)

    def __enter__(self):
        self._saved = set(_op.LEGACY)
        _op.LEGACY |= self.names
        return self

    def __exit__(self, *exc):
        _op.LEGACY.clear()
        _op.LEGACY |= self._saved
        return False


def gaia_audit(edges_arcsec=None) -> Table:
    """The Gaia EDR3 profile measured four ways: fixed, and with each defect restored."""
    from .outer_gaia import DEFAULT_EDGES, build_edr3_profile
    edges = DEFAULT_EDGES if edges_arcsec is None else np.asarray(edges_arcsec, float)
    out = {}
    for label, names in (("fixed", ()), ("legacy_signs", ("signs",)),
                         ("legacy_interval", ("interval",)), ("legacy_both", ("signs", "interval"))):
        with with_legacy(*names):
            out[label] = Table.read(build_edr3_profile(edges=edges))
    base = out["fixed"]
    t = Table({"r_median": base["r_median"], "n_stars": base["n_stars"]})
    for label, tab in out.items():
        t[f"sigma_{label}"] = tab["sigma_pm"]
        t[f"err_{label}"] = tab["sigma_pm_err"]
        t[f"ratio_{label}"] = np.asarray(tab["sigma_pmt"]) / np.asarray(tab["sigma_pmr"])
        t[f"meant_{label}"] = np.abs(np.asarray(tab["mean_pmt"]))
    # leave the product on disk in its corrected form
    build_edr3_profile(edges=edges)
    return t


def hst_audit(edges_arcsec=(150.0, 200.0, 250.0, 300.0, 340.0)) -> Table:
    """The HST profile with and without the floored uncertainty."""
    from .hst_profile import hst_profile
    fixed = hst_profile(edges_arcsec=edges_arcsec)
    with with_legacy("interval", "signs"):
        legacy = hst_profile(edges_arcsec=edges_arcsec)
    return Table({"r_median": fixed["r_median"], "n_stars": fixed["n_stars"],
                  "sigma_fixed": fixed["sigma_pm"], "err_fixed": fixed["sigma_pm_err"],
                  "sigma_legacy": legacy["sigma_pm"], "err_legacy": legacy["sigma_pm_err"]})
