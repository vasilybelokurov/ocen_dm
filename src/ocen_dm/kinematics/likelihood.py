"""Likelihood of the binned kinematic profiles given a spherical Jeans model.

Three rules, each enforced in code rather than left to the user:

1. **Bins, not bin centres.** A published dispersion is a statistic of the stars
   between ``r_lower`` and ``r_upper``; the model is compared through the
   tracer-weighted mean of ``sigma^2`` over the bin,
   ``<sigma^2> = int Sigma(R) sigma^2(R) R dR / int Sigma(R) R dR``.
   Profiles published without bin edges (a mean radius only) are evaluated at
   that radius, and say so in their record.
2. **Asymmetric errors** are used as published: a split normal with the upper
   error when the model lies above the datum and the lower one below.
3. **Independence.** The oMEGACat combined PM profile is built from the same
   stars as its radial and tangential profiles; including it with either is
   refused. Instrument-level nuisance scales exist for the Gaia profiles because
   they disagree with HST by 20-25 per cent where they overlap (JOURNAL, 2026-09-16).

Units: data radii in arcsec, LOS dispersions in km/s, PM dispersions in mas/yr;
the model is converted at the model's distance, which is a parameter.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from astropy.table import Table

from ..paths import processed_dir
from .jeans import KMS_PER_MASYR_KPC, SphericalJeans

__all__ = ["BinnedProfile", "KinematicData", "ProfileLikelihood", "load_profile", "DATASETS"]

ARCSEC_PER_RAD = 206264.806
_BIN_NODES = 6   # Gauss-Legendre nodes per bin for the tracer-weighted average

#: kind -> how the model dispersion is formed from the three projections
#:   'los'  : line of sight, km/s
#:   'pmr'  : PM radial component, mas/yr
#:   'pmt'  : PM tangential component, mas/yr
#:   'pmc'  : combined 1-D PM dispersion sqrt((pmr^2 + pmt^2)/2), mas/yr
KINDS = ("los", "pmr", "pmt", "pmc")


@dataclass(frozen=True)
class BinnedProfile:
    """One published dispersion profile.

    Attributes
    ----------
    name : str
        Dataset key.
    kind : str
        One of :data:`KINDS`.
    r : numpy.ndarray
        Representative radius per bin, arcsec.
    r_lower, r_upper : numpy.ndarray or None
        Bin edges in arcsec; ``None`` when the source publishes none.
    value, err_lo, err_hi : numpy.ndarray
        Dispersion and its lower/upper 1-sigma errors, km/s (los) or mas/yr (pm).
    instrument : str
        Label the nuisance scale is keyed on.
    shares_stars_with : tuple of str
        Other dataset keys built from the same stars.
    streaming2 : numpy.ndarray or None
        Mean-square streaming velocity ``<vbar^2>`` in each bin, in the square of
        the value's unit. Published dispersions are measured about a fitted mean
        (a rotation curve for MUSE), whereas the spherical Jeans model predicts
        the full second moment ``sigma^2 + vbar^2``; the model dispersion compared
        with the datum is therefore ``sqrt(sigma_model^2 - streaming2)``.
    """

    name: str
    kind: str
    r: np.ndarray
    r_lower: np.ndarray | None
    r_upper: np.ndarray | None
    value: np.ndarray
    err_lo: np.ndarray
    err_hi: np.ndarray
    instrument: str
    shares_stars_with: tuple[str, ...] = ()
    note: str = ""
    streaming2: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.streaming2 is not None and (len(self.streaming2) != len(self.r) or np.any(self.streaming2 < 0)):
            raise ValueError(f"{self.name}: streaming2 must be non-negative with one entry per bin")
        if self.kind not in KINDS:
            raise ValueError(f"kind must be one of {KINDS}")
        n = len(self.r)
        for a in (self.value, self.err_lo, self.err_hi):
            if len(a) != n:
                raise ValueError(f"{self.name}: column lengths differ")
        if np.any(self.err_lo <= 0) or np.any(self.err_hi <= 0):
            raise ValueError(f"{self.name}: errors must be positive")
        if (self.r_lower is None) != (self.r_upper is None):
            raise ValueError(f"{self.name}: give both bin edges or neither")
        if self.r_lower is not None and np.any(self.r_upper <= self.r_lower):
            raise ValueError(f"{self.name}: r_upper must exceed r_lower")

    @property
    def has_edges(self) -> bool:
        return self.r_lower is not None

    @property
    def n(self) -> int:
        return len(self.r)


# ----------------------------------------------------------------- loading ---
def pm_rotation_curve(r_arcsec: np.ndarray) -> np.ndarray:
    """Mean tangential proper motion of the rotating cluster at projected radius ``r``.

    Interpolated (linearly, in mas/yr) from the Vasiliev & Baumgardt (2021) EDR3
    profile of "the rotational PM component" (their readme), i.e. the azimuthally
    averaged mean tangential PM per annulus: 0 at the centre, peaking at
    0.25 mas/yr near 430 arcsec, 0.02 mas/yr by 2000 arcsec. Used as ``<mu_T>``;
    the streaming term ``<mu_T^2>`` >= ``<mu_T>^2`` is therefore a lower bound
    when the rotation axis is not in the plane of the sky.
    """
    t = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    return np.interp(np.asarray(r_arcsec, float), np.asarray(t["r"], float), np.asarray(t["vrot_pm"], float))


def _asym(t: Table, base: str) -> tuple[np.ndarray, np.ndarray]:
    return np.asarray(t[f"{base}_err_lo"], float), np.asarray(t[f"{base}_err_hi"], float)


def _omegacat(name: str, value: str, kind: str, shares: tuple[str, ...]) -> BinnedProfile:
    t = Table.read(processed_dir() / "kinematics" / f"omegacat_vi_{name}.ecsv")
    lo, hi = _asym(t, value)
    streaming2 = None
    note = "oMEGACat VI adaptive log bins"
    if kind == "los":
        # sigma_los was fitted jointly with vbar(phi) = v_rot sin(phi - theta_0) in each
        # annulus (columns vlos, theta_0 of los_profile.fits, our los_rotation product,
        # same bins); averaged over the annulus <vbar^2> = v_rot^2 / 2. This is the
        # term that turns the model's second moment into the published dispersion.
        rot = Table.read(processed_dir() / "kinematics" / "omegacat_vi_los_rotation.ecsv")
        if len(rot) != len(t) or not np.allclose(rot["r_median"], t["r_median"]):
            raise ValueError("los_rotation and los_dispersion products are not on the same bins")
        streaming2 = 0.5 * np.asarray(rot["v_rot"], float) ** 2
        note += "; dispersion about the fitted rotation curve, <vbar^2> = v_rot^2/2 added to the model"
    elif kind == "pmt":
        # The catalogue PMs are 'locally corrected' (local mean removed; measured mean
        # mu_T = 0 +- 0.005 mas/yr per annulus, 2026-09-17), so sigma_pmt is a dispersion
        # about the rotating mean. The rotation itself comes from the Vasiliev &
        # Baumgardt (2021) EDR3 rotation curve, evaluated at the bin median radius.
        streaming2 = pm_rotation_curve(np.asarray(t["r_median"], float)) ** 2
        note += "; <mu_T>^2 from the Vasiliev & Baumgardt 2021 PM rotation curve added to the model"
    elif kind == "pmr":
        note += "; radial PM carries no rotation term (rotation is tangential)"
    return BinnedProfile(
        name=f"hst_{name}" if kind != "los" else "muse_los_dispersion", kind=kind,
        r=np.asarray(t["r_median"], float), r_lower=np.asarray(t["r_lower"], float),
        r_upper=np.asarray(t["r_upper"], float), value=np.asarray(t[value], float),
        err_lo=lo, err_hi=hi, instrument="HST" if kind != "los" else "MUSE",
        shares_stars_with=shares, note=note, streaming2=streaming2,
    )


def _baumgardt2019() -> BinnedProfile:
    t = Table.read(processed_dir() / "kinematics" / "baumgardt2019_ocen_pm_dispersion.ecsv")
    return BinnedProfile(
        name="gaia_dr2_pm", kind="pmc", r=np.asarray(t["r"], float), r_lower=None, r_upper=None,
        value=np.asarray(t["sigma_pm"], float), err_lo=np.asarray(t["sigma_pm_err_lo"], float),
        err_hi=np.asarray(t["sigma_pm_err_hi"], float), instrument="GaiaDR2",
        shares_stars_with=("gaia_edr3_pm",),
        note="mean radius only, no bin edges; ASSUMED to be the 1-D PM dispersion (fitted scale "
             "~0.9-1.0, not 1.4, supports this); whether it was fitted about a rotating mean is "
             "unknown from the local files -- NO rotation term applied; if the published values "
             "include rotation this is correct, otherwise the model is biased high by up to "
             "vrot^2/2 (~10 per cent of sigma^2 at 700 arcsec)",
    )


def _vasiliev2021(r_min_arcsec: float = 0.0, n_max: int | None = 8) -> BinnedProfile:
    """Vasiliev & Baumgardt (2021) EDR3 PM dispersion profile.

    The authors publish the profile on a fine radial grid; neighbouring points are
    derived from overlapping stellar samples and are not independent. Using all
    of them would give this one dataset several times the weight of the HST
    profiles, so by default the grid is thinned to ``n_max`` points spread evenly
    in log radius. Set ``n_max=None`` to keep every point (only for plotting).
    """
    t = Table.read(processed_dir() / "kinematics" / "vasiliev2021_ocen_pm_profiles.ecsv")
    r = np.asarray(t["r"], float); s50 = np.asarray(t["sigma_pm"], float)
    lo = s50 - np.asarray(t["sigma_pm_p16"], float); hi = np.asarray(t["sigma_pm_p84"], float) - s50
    keep = np.flatnonzero((r > max(r_min_arcsec, 0.0)) & (s50 > 0) & (lo > 0) & (hi > 0))
    if n_max is not None and len(keep) > n_max:
        targets = np.geomspace(r[keep][0], r[keep][-1], n_max)
        keep = np.unique(keep[np.abs(np.log(r[keep])[:, None] - np.log(targets)[None, :]).argmin(axis=0)])
    thinned = n_max is not None
    vrot = np.asarray(t["vrot_pm"], float)[keep]
    return BinnedProfile(
        name="gaia_edr3_pm", kind="pmc", r=r[keep], r_lower=None, r_upper=None,
        value=s50[keep], err_lo=lo[keep], err_hi=hi[keep], instrument="GaiaEDR3",
        shares_stars_with=("gaia_dr2_pm",),
        note=("thinned to %d log-spaced points because " % len(keep) if thinned else "ALL points kept although ")
             + "adjacent grid points come from overlapping samples and are not independent; "
             "percentile errors; dispersion fitted jointly with rotation (readme: 'PM dispersion' "
             "and 'PM rotation' profiles), so <mu_T>^2/2 is added to the 1-D model; ASSUMED 1-D "
             "PM dispersion -- the fitted scale s_GaiaEDR3 ~ 1.05 (not 1.4) supports this",
        streaming2=0.5 * vrot ** 2,
    )


def _edr3_ours() -> BinnedProfile:
    """Our own Gaia EDR3 dispersion profile, measured beyond 460 arcsec.

    Replaces ``gaia_edr3_pm`` (a 2-5 node spline whose inner points are an inward
    continuation over a region with no usable Gaia star). Built from quality-flagged stars
    whose errors are small next to the signal, so the answer does not depend on the error
    model; see :mod:`ocen_dm.kinematics.outer_gaia`. Bins are independent, the quoted error
    already carries the error-model systematic, and rotation is returned as ``streaming2``.
    """
    from .outer_gaia import load_edr3_profile
    t = load_edr3_profile()
    err = np.asarray(t["sigma_pm_err"], float)
    return BinnedProfile(
        name="gaia_edr3_ours", kind="pmc", r=np.asarray(t["r_median"], float),
        r_lower=np.asarray(t["r_lower"], float), r_upper=np.asarray(t["r_upper"], float),
        value=np.asarray(t["sigma_pm"], float), err_lo=err, err_hi=err, instrument="GaiaEDR3",
        shares_stars_with=("gaia_dr2_pm", "gaia_edr3_pm"),
        note="our measurement: quality-flagged EDR3 stars with err < 0.4 sigma(R), 2-D empirical "
             "field template, exact perspective and depth removal, independent annuli; the value "
             "is the midpoint of the raw-error and density-inflated-error fits and half their "
             "separation is included in the error as a systematic",
        streaming2=np.asarray(t["streaming2"], float),
    )


def _hst_ours(component: str) -> BinnedProfile:
    """Our own HST measurement, flagged stars only, reaching 360 arcsec.

    Replaces the published oMEGACat profile, which stops at 300 arcsec and so never
    overlapped Gaia. Same stars the survey vouches for, our estimator, 60 arcsec further
    out. The tangential dataset carries the external rotation term for the same reason the
    published one does: HST proper motions are locally corrected, so the measured dispersion
    is about a rotation-free mean while the Jeans model predicts the full second moment.
    """
    from .hst_profile import load_hst_product
    t = load_hst_product()
    col = "sigma_pmr" if component == "pmr" else "sigma_pmt"
    r = np.asarray(t["r_median"], float)
    err = np.asarray(t[col + "_err"], float)
    other = "hst_pm_tangential_ours" if component == "pmr" else "hst_pm_radial_ours"
    note = ("our measurement from selection_hq_astrometry stars, cluster+field mixture, "
            "reaching 360 arcsec where the published profile stops at 300")
    return BinnedProfile(
        name="hst_pm_radial_ours" if component == "pmr" else "hst_pm_tangential_ours",
        kind=component, r=r, r_lower=np.asarray(t["r_lower"], float),
        r_upper=np.asarray(t["r_upper"], float), value=np.asarray(t[col], float),
        err_lo=err, err_hi=err, instrument="HST",
        shares_stars_with=(other, "hst_pm_radial", "hst_pm_tangential", "hst_pm_combined"),
        note=note + ("; radial PM carries no rotation term" if component == "pmr" else
                     "; <mu_T>^2 from the Vasiliev & Baumgardt 2021 curve added to the model"),
        streaming2=None if component == "pmr" else pm_rotation_curve(r) ** 2)


def _edr3_ours_component(component: str) -> BinnedProfile:
    """Our Gaia EDR3 measurement split into its radial and tangential parts.

    The two are measured jointly but are very nearly uncorrelated: over 60 independent
    synthetic realisations with uniform position angles the correlation between the fitted
    sigma_R and sigma_T is -0.08, so treating them as two datasets costs a 0.6 per cent
    error in the joint chi2 and buys the anisotropy at large radius (2026-09-19).
    """
    from .outer_gaia import load_edr3_profile
    t = load_edr3_profile()
    col = "sigma_pmr" if component == "pmr" else "sigma_pmt"
    err = np.asarray(t[col + "_err"], float)
    mean = np.asarray(t["mean_pmr" if component == "pmr" else "mean_pmt"], float)
    other = "gaia_edr3_ours_tangential" if component == "pmr" else "gaia_edr3_ours_radial"
    return BinnedProfile(
        name="gaia_edr3_ours_radial" if component == "pmr" else "gaia_edr3_ours_tangential",
        kind=component, r=np.asarray(t["r_median"], float),
        r_lower=np.asarray(t["r_lower"], float), r_upper=np.asarray(t["r_upper"], float),
        value=np.asarray(t[col], float), err_lo=err, err_hi=err, instrument="GaiaEDR3",
        shares_stars_with=(other, "gaia_edr3_ours", "gaia_edr3_pm", "gaia_dr2_pm"),
        note="our measurement, one component; the pair is nearly uncorrelated (rho = -0.08)",
        streaming2=mean ** 2)


DATASETS: dict[str, Any] = {
    "hst_pm_radial": lambda: _omegacat("pm_radial", "sigma_pmr", "pmr", ("hst_pm_tangential", "hst_pm_combined")),
    "hst_pm_tangential": lambda: _omegacat("pm_tangential", "sigma_pmt", "pmt", ("hst_pm_radial", "hst_pm_combined")),
    "hst_pm_combined": lambda: _omegacat("pm_combined", "sigma_pmc", "pmc", ("hst_pm_radial", "hst_pm_tangential")),
    "muse_los_dispersion": lambda: _omegacat("los_dispersion", "sigma_los", "los", ()),
    "gaia_dr2_pm": _baumgardt2019,
    "gaia_edr3_pm": _vasiliev2021,
    "gaia_edr3_ours": _edr3_ours,
    "gaia_edr3_ours_radial": lambda: _edr3_ours_component("pmr"),
    "gaia_edr3_ours_tangential": lambda: _edr3_ours_component("pmt"),
    "hst_pm_radial_ours": lambda: _hst_ours("pmr"),
    "hst_pm_tangential_ours": lambda: _hst_ours("pmt"),
}

#: combinations that would count the same stars twice
_FORBIDDEN_TOGETHER = (("hst_pm_combined", "hst_pm_radial"), ("hst_pm_combined", "hst_pm_tangential"),
                       ("gaia_edr3_ours", "gaia_edr3_pm"),
                       ("gaia_edr3_ours", "gaia_edr3_ours_radial"),
                       ("gaia_edr3_ours", "gaia_edr3_ours_tangential"),
                       ("gaia_edr3_ours_radial", "gaia_edr3_pm"),
                       ("gaia_edr3_ours_tangential", "gaia_edr3_pm"),
                       ("hst_pm_radial", "hst_pm_radial_ours"),
                       ("hst_pm_tangential", "hst_pm_tangential_ours"),
                       ("hst_pm_combined", "hst_pm_radial_ours"),
                       ("hst_pm_combined", "hst_pm_tangential_ours"))


def load_profile(name: str, **kwargs: Any) -> BinnedProfile:
    """Load one dataset by key."""
    if name not in DATASETS:
        raise KeyError(f"unknown dataset {name!r}; known: {sorted(DATASETS)}")
    return DATASETS[name](**kwargs) if kwargs else DATASETS[name]()


@dataclass(frozen=True)
class KinematicData:
    """The set of profiles entering one likelihood."""

    profiles: tuple[BinnedProfile, ...]

    def __post_init__(self) -> None:
        names = [p.name for p in self.profiles]
        if len(set(names)) != len(names):
            raise ValueError("duplicate datasets in the likelihood")
        for a, b in _FORBIDDEN_TOGETHER:
            if a in names and b in names:
                raise ValueError(
                    f"{a!r} and {b!r} are built from the same stars; use radial+tangential "
                    "OR combined, never both"
                )

    @classmethod
    def load(cls, names: Iterable[str], **per_dataset: dict[str, Any]) -> "KinematicData":
        return cls(tuple(load_profile(n, **per_dataset.get(n, {})) for n in names))

    @property
    def instruments(self) -> tuple[str, ...]:
        return tuple(sorted({p.instrument for p in self.profiles}))

    @property
    def n_points(self) -> int:
        return sum(p.n for p in self.profiles)


# -------------------------------------------------------------- likelihood ---
class ProfileLikelihood:
    """Log-likelihood of binned dispersion profiles under a Jeans model.

    Parameters
    ----------
    data : KinematicData
    """

    def __init__(self, data: KinematicData) -> None:
        self.data = data
        self._nodes_gl, self._weights_gl = np.polynomial.legendre.leggauss(_BIN_NODES)

    # ---- model prediction --------------------------------------------------
    def _nodes(self, profile: BinnedProfile, pc_per_arcsec: float) -> tuple[np.ndarray, np.ndarray | None]:
        """Projected radii (pc) at which the model is needed for ``profile`` and the
        quadrature weights (``None`` for point-evaluated profiles)."""
        if not profile.has_edges:
            return profile.r * pc_per_arcsec, None
        lo = profile.r_lower * pc_per_arcsec
        hi = profile.r_upper * pc_per_arcsec
        lo = np.maximum(lo, 1e-3 * hi)                                  # r_lower = 0 in the innermost bin
        half = 0.5 * (hi - lo)
        R = lo[:, None] + half[:, None] * (self._nodes_gl[None, :] + 1.0)   # (n_bins, n_nodes)
        return R.ravel(), (half[:, None] * self._weights_gl[None, :]).ravel()

    @staticmethod
    def _sigma2(m: dict[str, np.ndarray], kind: str) -> np.ndarray:
        if kind == "pmc":
            return 0.5 * (m["pmr"] + m["pmt"]) / m["Sigma"]
        return m[kind] / m["Sigma"]

    def _reduce(self, profile: BinnedProfile, m: dict[str, np.ndarray], R: np.ndarray,
                w: np.ndarray | None, distance_kpc: float) -> np.ndarray:
        """Bin-average ``sigma^2`` with tracer weight ``Sigma R`` and convert units."""
        if w is None:
            value2 = self._sigma2(m, profile.kind)
        else:
            shape = (profile.n, _BIN_NODES)
            weight = (m["Sigma"] * R * w).reshape(shape)
            s2 = self._sigma2(m, profile.kind).reshape(shape)
            value2 = np.sum(weight * s2, axis=1) / np.sum(weight, axis=1)
        if profile.kind != "los":
            value2 = value2 / (KMS_PER_MASYR_KPC * distance_kpc) ** 2  # -> (mas/yr)^2
        if profile.streaming2 is not None:
            value2 = np.maximum(value2 - profile.streaming2, 1e-12)     # second moment -> dispersion about the mean
        return np.sqrt(value2)

    def predict_profile(self, jeans: SphericalJeans, distance_kpc: float,
                        profile: BinnedProfile) -> np.ndarray:
        """Model dispersion for each bin of ``profile``, in the profile's units."""
        R, w = self._nodes(profile, distance_kpc * 1e3 / ARCSEC_PER_RAD)
        jeans._check_inside_tracer(R)
        return self._reduce(profile, jeans.projected_moments(R), R, w, distance_kpc)

    def predict(self, jeans: SphericalJeans, distance_kpc: float,
                scales: dict[str, float] | None = None) -> dict[str, np.ndarray]:
        """Model dispersions for every profile, with instrument scales applied.

        All profiles' radii go through one call to the solver's projection.
        """
        scales = scales or {}
        pc_per_arcsec = distance_kpc * 1e3 / ARCSEC_PER_RAD
        nodes = [self._nodes(p, pc_per_arcsec) for p in self.data.profiles]
        R_all = np.concatenate([R for R, _ in nodes])
        jeans._check_inside_tracer(R_all)
        m_all = jeans.projected_moments(R_all)
        out, start = {}, 0
        for p, (R, w) in zip(self.data.profiles, nodes):
            sl = slice(start, start + len(R)); start += len(R)
            m = {k: v[sl] for k, v in m_all.items()}
            out[p.name] = self._reduce(p, m, R, w, distance_kpc) * scales.get(p.instrument, 1.0)
        return out

    # ---- likelihood ----------------------------------------------------------
    @staticmethod
    def _split_normal_lnlike(model: np.ndarray, p: BinnedProfile) -> np.ndarray:
        err = np.where(model > p.value, p.err_hi, p.err_lo)
        return -0.5 * ((model - p.value) / err) ** 2 - np.log(err) - 0.5 * np.log(2 * np.pi)

    def lnlike_terms(self, jeans: SphericalJeans, distance_kpc: float,
                     scales: dict[str, float] | None = None) -> dict[str, float]:
        """Per-dataset log-likelihoods."""
        pred = self.predict(jeans, distance_kpc, scales)
        return {p.name: float(np.sum(self._split_normal_lnlike(pred[p.name], p))) for p in self.data.profiles}

    def lnlike(self, jeans: SphericalJeans, distance_kpc: float,
               scales: dict[str, float] | None = None) -> float:
        """Total log-likelihood; ``-inf`` if the model cannot be evaluated."""
        try:
            return float(sum(self.lnlike_terms(jeans, distance_kpc, scales).values()))
        except (ValueError, FloatingPointError):
            return -np.inf

    def chi2_terms(self, jeans: SphericalJeans, distance_kpc: float,
                   scales: dict[str, float] | None = None) -> dict[str, tuple[float, int]]:
        """``(chi^2, n_points)`` per dataset, for goodness-of-fit reporting."""
        pred = self.predict(jeans, distance_kpc, scales)
        out = {}
        for p in self.data.profiles:
            err = np.where(pred[p.name] > p.value, p.err_hi, p.err_lo)
            out[p.name] = (float(np.sum(((pred[p.name] - p.value) / err) ** 2)), p.n)
        return out
