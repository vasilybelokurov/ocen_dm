"""Measure the HST proper-motion dispersion ourselves, in the annuli Gaia uses.

The published oMEGACat profile stops at **300 arcsec**, and the catalogue's high-quality
astrometry flag stops at about **340 arcsec**. Neither is where the data stop: the catalogue
carries proper motions for 1475096 stars out to **466 arcsec**, including 23730 between 340
and 380 and 7725 between 380 and 420. That was the whole reason the HST and Gaia profiles
appeared not to overlap -- a choice of published product, not a property of the data.

Taking the raw stars beyond the flag does not work on its own: their error-deconvolved
dispersion is 25-36 per cent above the flagged stars at the same radius and it *rises*
outwards (0.80, 0.88, 1.78 mas/yr at 340-380, 380-420, 420-466), which no cluster does. The
cause is field contamination, and it is easy to see once stated: HST proper motions are
**relative to the cluster**, so cluster members sit at the origin while field stars sit near
minus the systemic motion, 7.5 mas/yr away. Only 1-2 per cent of stars lie beyond 3 mas/yr,
but in a raw variance those few dominate.

So this module applies the same treatment as everywhere else in the project: a cluster
two-dimensional Gaussian with its axes along each star's radial direction, convolved with
that star's error covariance, plus the **same empirical Gaia DR3 field template** used for
the Gaia and Pristine profiles, scored in absolute proper motion (relative plus systemic).
No proper-motion window is imposed, so the answer is not a readout of a cut.

**What the flag is.** ``selection_hq_astrometry`` is oMEGACat's own published high-quality
selection (Haeberle et al. 2025, ApJ, arXiv:2503.04903, their Section 2.1). It keeps 48 per
cent of the stars with proper motions and requires, all together: a temporal baseline longer
than 10 years; more than 75 per cent of the individual measurements surviving the clipping
stage; reduced chi-square below 5 in both proper-motion components; a proper-motion error
inside the lower 95 per cent of the error distribution in its own 0.5-magnitude bin; and
reliable photometry in **both** F625W and F814W (unsaturated, point-spread-function fit
quality above the 85th percentile of its magnitude bin, neighbour flux inside the fit
aperture below half the star's own). Stars fainter than F625W = 24 are cut outright, because
there the error limit reaches 0.3 mas/yr, comparable to half the outer velocity dispersion.

Their stated principle is the one this project arrived at independently for Gaia: **a star
whose measurement error approaches the dispersion being measured makes the answer a function
of the error model rather than of the cluster.** The two-filter photometric requirement is
not cosmetic either -- F625W and F814W span 2002 to 2022, so demanding both is how they
verify the astrometry is supported across the whole temporal baseline.

**The unflagged stars are not as good as the flagged ones, and it shows.** Running the same
mixture on each subset where both exist gives an unflagged-to-flagged dispersion ratio of
1.092, 1.088, 1.079 and 1.052 at 150-200, 200-250, 250-300 and 300-340 arcsec, a weighted
mean of **1.077 +- 0.011**. Their errors are underestimated by about that much: the same
pathology as unflagged Gaia, at a fifth of the amplitude. Outside 340 arcsec every star is
unflagged, so the raw measurement there is biased high by the full 8 per cent.

Those outer stars fail for a specific reason: **beyond 380 arcsec not one has F625W or F814W
photometry**, so the flag's two-filter requirement rejects them automatically. That could
have meant their astrometry was fine and the correction inappropriate, so it was tested.
Inner unflagged stars that lack photometry give dispersions of 1.174, 1.066 and 1.044 times
the flagged value at 150-250, 250-300 and 300-340 arcsec, against 1.090, 1.080 and 1.066 for
unflagged stars that have it. Missing photometry is therefore not a free pass: those stars
are inflated as much as the rest, which is what the paper's own reasoning predicts. The
adopted 1.077 sits inside that range, though the value nearest the boundary is 1.044, so the
outer correction could be up to 3 per cent too large.
``hst_profile`` therefore divides each bin by ``1 + f_unflagged * (ratio - 1)`` and carries
the ratio's uncertainty as a systematic. Validation of the whole chain: with the flag
required, our measurement sits 1.1-2.3 per cent above the published oMEGACat profile in its
own bins; without the flag and without this correction it sits 4.7-6.1 per cent above it.

Two approximations, both small over this footprint and both flagged: the systemic motion is
removed as a constant rather than projected star by star (the perspective term is below
0.001 mas/yr inside 8 arcmin), and the field template's magnitude selection is Gaia's, not
HST's, which matters at the per-cent level when the field fraction is itself 1-2 per cent.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from pathlib import Path

from ..paths import raw_dir
from .outer_profile import KMS_PER_MASYR_KPC, MemberSample, OCEN_DEC, OCEN_RA, dispersion_2d

__all__ = ["MU_SYS", "HST_R_MAX_ARCSEC", "HST_FLAG_MAX_ARCSEC", "UNFLAGGED_BIAS",
           "DEFAULT_EDGES", "PRODUCT", "load_hst_sample", "unflagged_bias", "hst_profile",
           "build_hst_profile", "load_hst_product"]

#: where the quality flag runs out: 66 stars in 340-360 arcsec, none beyond
HST_FLAG_MAX_ARCSEC = 360.0
PRODUCT = "ocen_pm_dispersion_hst_ours"
#: log-spaced inside, then the three annuli that carry the overlap with Gaia
DEFAULT_EDGES = np.concatenate([np.geomspace(2.0, 250.0, 21), [300.0, 340.0, 360.0]])

#: measured ratio of the unflagged to the flagged dispersion, and its uncertainty
UNFLAGGED_BIAS = (1.077, 0.011)

#: systemic proper motion, used only to place the field template in the relative frame
MU_SYS = (-3.257, -6.730)
#: the catalogue's own outer edge
HST_R_MAX_ARCSEC = 466.0


def load_hst_sample(require_flag: bool = False) -> MemberSample:
    """oMEGACat proper motions as a :class:`MemberSample` in the cluster's relative frame."""
    t = Table.read(raw_dir() / "omegacat_vi_kinematics" / "catalog_and_selections.fits")
    ra = np.asarray(t["RA"], float); dec = np.asarray(t["DEC"], float)
    a = np.asarray(t["pmra_corrected"], float); d = np.asarray(t["pmdec_corrected"], float)
    ea = np.asarray(t["pmra_corrected_err"], float); ed = np.asarray(t["pmdec_corrected_err"], float)
    keep = (np.isfinite(ra) & np.isfinite(dec) & np.isfinite(a) & np.isfinite(d)
            & np.isfinite(ea) & np.isfinite(ed) & (ea > 0) & (ed > 0))
    if require_flag:
        keep &= np.asarray(t["selection_hq_astrometry"], float) > 0
    ra, dec, a, d, ea, ed = ra[keep], dec[keep], a[keep], d[keep], ea[keep], ed[keep]
    gcol = next((c for c in ("f625w", "F625W", "f814w", "mag_F625W") if c in t.colnames), None)
    g = np.asarray(t[gcol], float)[keep] if gcol else np.full(int(keep.sum()), np.nan)
    flag = (np.asarray(t["selection_hq_astrometry"], float)[keep] > 0).astype(int) * 2

    x = (ra - OCEN_RA) * np.cos(np.radians(OCEN_DEC)) * 3600.0
    y = (dec - OCEN_DEC) * 3600.0
    r = np.hypot(x, y); safe = np.maximum(r, 1e-9)
    cos_p, sin_p = x / safe, y / safe
    # the catalogue's proper motions are already relative to the cluster; re-centre on the
    # flagged members so the cluster sits exactly at the origin
    ref = flag > 0
    a = a - np.median(a[ref]); d = d - np.median(d[ref])
    return MemberSample(
        r, a * cos_p + d * sin_p, -a * sin_p + d * cos_p,
        np.sqrt((ea * cos_p) ** 2 + (ed * sin_p) ** 2),
        np.sqrt((ea * sin_p) ** 2 + (ed * cos_p) ** 2),
        np.ones(len(r)), g, flag, np.arctan2(x, y),
        mu_a=a, mu_d=d, err_a=ea, err_d=ed, err_corr=np.zeros(len(r)),
        sys_a=np.full(len(r), MU_SYS[0]), sys_d=np.full(len(r), MU_SYS[1]),
        mu_sys=MU_SYS, exact=False)


def unflagged_bias(edges_arcsec=(150.0, 200.0, 250.0, 300.0, 340.0)) -> tuple[float, float]:
    """Remeasure the unflagged-to-flagged dispersion ratio where both subsets exist."""
    from ..selection.field_template import field_density_2d
    s = load_hst_sample()
    dens = field_density_2d()
    flag = s.quality_flag > 0
    edges = np.asarray(edges_arcsec, float)

    def mix(mask):
        o = dispersion_2d(s, mask, dens, depth_var=0.0, field_at=s.absolute_pm, sigma_max=1.5)
        return (float(np.sqrt(0.5 * (o["sigma_r"] ** 2 + o["sigma_t"] ** 2))),
                float(0.5 * np.hypot(o["sigma_r_err"], o["sigma_t_err"])))

    ratios, errs = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        base = (s.r_arcsec >= lo) & (s.r_arcsec < hi)
        if (base & flag).sum() < 1000 or (base & ~flag).sum() < 1000:
            continue
        a, ea = mix(base & flag); b, eb = mix(base & ~flag)
        ratios.append(b / a); errs.append((b / a) * float(np.hypot(ea / a, eb / b)))
    r, e = np.asarray(ratios), np.asarray(errs)
    w = float(np.sum(r / e ** 2) / np.sum(1 / e ** 2))
    return w, float(1 / np.sqrt(np.sum(1 / e ** 2)))


def hst_profile(edges_arcsec=(300.0, 340.0), require_flag: bool = True,
                distance_kpc: float = 5.43, min_stars: int = 100,
                correct_unflagged: bool = False) -> Table:
    """Cluster dispersion per annulus from the mixture fit, with no proper-motion window.

    **The default is flagged stars only.** If a star's astrometry is not trusted, the answer
    is to leave it out, not to model it: inside 340 arcsec there are tens of thousands of
    flagged stars and nothing is gained by adding the rest.

    ``correct_unflagged`` exists for one case only. Beyond 340 arcsec the catalogue contains
    **zero** flagged stars, so a flagged dispersion cannot be computed there at all; the
    option divides each bin by ``1 + f_unflagged * (k - 1)`` with ``k`` from
    :data:`UNFLAGGED_BIAS` and carries the ratio's uncertainty. Those points are a
    cross-check, never a measurement, and never enter a fit.
    """
    from ..selection.field_template import field_density_2d
    s = load_hst_sample(require_flag=require_flag)
    dens = field_density_2d()
    k = KMS_PER_MASYR_KPC * distance_kpc
    edges = np.asarray(edges_arcsec, float)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (s.r_arcsec >= lo) & (s.r_arcsec < hi)
        if m.sum() < min_stars:
            continue
        o = dispersion_2d(s, m, dens, depth_var=0.0, field_at=s.absolute_pm, sigma_max=1.5)
        sig = float(np.sqrt(0.5 * (o["sigma_r"] ** 2 + o["sigma_t"] ** 2)))
        err = float(0.5 * np.hypot(o["sigma_r_err"], o["sigma_t_err"]))
        n_flag = int((m & (s.quality_flag > 0)).sum())
        f_un = 1.0 - n_flag / m.sum()
        corr = 1.0 + f_un * (UNFLAGGED_BIAS[0] - 1.0) if correct_unflagged else 1.0
        sys_err = sig / corr * f_un * UNFLAGGED_BIAS[1] if correct_unflagged else 0.0
        sig_c = sig / corr
        err_c = float(np.hypot(err / corr, sys_err))
        rows.append((lo, float(np.median(s.r_arcsec[m])), hi, int(m.sum()), n_flag, f_un,
                     o["f"], o["sigma_r"] / corr, o["sigma_r_err"] / corr,
                     o["sigma_t"] / corr, o["sigma_t_err"] / corr,
                     sig_c, err_c, sig, sys_err, sig_c * k, err_c * k,
                     float(np.median(s.r_arcsec[m]) * distance_kpc * 1e3 / 206264.806)))
    return Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_stars", "n_flagged",
                                   "f_unflagged", "f_field", "sigma_pmr", "sigma_pmr_err",
                                   "sigma_pmt", "sigma_pmt_err", "sigma_pm", "sigma_pm_err",
                                   "sigma_raw", "sigma_sys", "sigma_kms", "sigma_kms_err",
                                   "r_pc"),
                 meta={"require_flag": require_flag, "correct_unflagged": correct_unflagged,
                       "unflagged_bias": UNFLAGGED_BIAS,
                       "note": "HST PMs are relative; only dispersions are meaningful. Field is "
                               "the Gaia DR3 empirical template scored in absolute PM."})


def build_hst_profile(edges_arcsec=None, distance_kpc: float = 5.43,
                      min_stars: int = 60) -> "Path":
    """Measure and write ``ocen_pm_dispersion_hst_ours.ecsv``, flagged stars only.

    This is the product that replaces the published oMEGACat profile in the likelihood: it
    uses the same stars the survey vouches for, measured with the project's own estimator,
    and it reaches **360 arcsec** rather than 300, which is what gives a genuine overlap
    with Gaia (JOURNAL 2026-09-19).
    """
    from datetime import datetime, timezone
    from ..paths import processed_dir
    t = hst_profile(edges_arcsec=DEFAULT_EDGES if edges_arcsec is None else edges_arcsec,
                    require_flag=True, correct_unflagged=False,
                    distance_kpc=distance_kpc, min_stars=min_stars)
    t.meta.update({
        "product": PRODUCT, "flag": "selection_hq_astrometry",
        "r_max_arcsec": HST_FLAG_MAX_ARCSEC,
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": "flagged stars only, no unflagged-star correction; HST proper motions are "
                "locally corrected, so the dispersion is about a rotation-free local mean and "
                "the tangential dataset carries the external rotation term",
    })
    path = processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def load_hst_product(path=None) -> Table:
    from ..paths import processed_dir
    path = path or processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    if not path.exists():
        return Table.read(build_hst_profile())
    return Table.read(path)
