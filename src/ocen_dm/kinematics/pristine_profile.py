"""Extend the dispersion profile outwards with Pristine stars and our own mixture model.

The proper-motion dispersion measured from the Pristine periphery catalogue by selecting on
proper motion is circular: the catalogue's membership uses the proper motions, and a fixed
window is itself part of the answer (6.5 km/s at a 0.8 mas/yr window against 12.8 at 1.2 --
see :mod:`ocen_dm.kinematics.periphery`). This module removes that circularity by applying
the same **cluster-plus-field mixture** we use for Gaia: no proper-motion cut at all.

* the sample is every Pristine star with ``[Fe/H] < feh_max``, a chemical cut that is
  independent of the kinematics (field median -0.26, cluster -1.49);
* the cluster is a two-dimensional Gaussian whose axes follow each star's radial direction,
  convolved with that star's full error covariance, the correlation coefficient coming from
  the matched Gaia record;
* the field is an **empirical two-dimensional proper-motion density** measured from the same
  catalogue, under the same metallicity cut, beyond ``field_r_min_deg`` where the star-count
  excess is consistent with zero (:mod:`ocen_dm.selection.periphery_density`);
* the cluster fraction, the two dispersions and the two mean motions are fitted per annulus.

The line-of-sight depth term is omitted: it is at most 0.003 mas/yr, an order below the
statistical errors here, and the tracer model it needs is only calibrated inside 66 pc.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from ..paths import processed_dir
from .outer_profile import (KMS_PER_MASYR_KPC, MemberSample, OCEN_DEC, OCEN_RA, OCEN_VSYS_KMS,
                            dispersion_2d)
from .perspective import systemic_pm_field, systemic_velocity_vector

__all__ = ["MU_SYS", "load_pristine_sample", "pristine_field_density", "pristine_profile"]

#: systemic proper motion of the cluster (mas/yr), Vasiliev & Baumgardt 2021
MU_SYS = (-3.257, -6.730)


def load_pristine_sample(feh_max: float = -1.2, distance_kpc: float = 5.43) -> MemberSample:
    """Pristine metal-poor stars as a :class:`MemberSample`, systemic field removed exactly."""
    t = Table.read(processed_dir() / "tails" / "kuzma2025_periphery.ecsv")
    cov = Table.read(processed_dir() / "tails" / "kuzma2025_periphery_gaia_covariance.ecsv")
    order = {int(s): i for i, s in enumerate(np.asarray(cov["source_id"], np.int64))}
    idx = np.array([order[int(s)] for s in np.asarray(t["source_id"], np.int64)])
    rho = np.asarray(cov["pmra_pmdec_corr"], float)[idx]

    keep = np.asarray(t["feh"], float) < feh_max
    t = t[keep]; rho = rho[keep]
    ra = np.asarray(t["ra"], float); dec = np.asarray(t["dec"], float)
    x = (ra - OCEN_RA) * np.cos(np.radians(OCEN_DEC)) * 3600.0
    y = (dec - OCEN_DEC) * 3600.0
    r = np.hypot(x, y); safe = np.maximum(r, 1e-9)
    cos_p, sin_p = x / safe, y / safe

    v_sys = systemic_velocity_vector(OCEN_RA, OCEN_DEC, MU_SYS[0], MU_SYS[1],
                                     OCEN_VSYS_KMS, distance_kpc)
    exp_a, exp_d = systemic_pm_field(ra, dec, v_sys, distance_kpc)
    pmra = np.asarray(t["pmra"], float) - exp_a
    pmdec = np.asarray(t["pmdec"], float) - exp_d
    ea = np.asarray(t["pmra_error"], float); ed = np.asarray(t["pmdec_error"], float)
    cov_ad = rho * ea * ed
    return MemberSample(
        r, pmra * cos_p + pmdec * sin_p, -pmra * sin_p + pmdec * cos_p,
        np.sqrt(np.maximum((ea * cos_p) ** 2 + (ed * sin_p) ** 2 + 2 * cov_ad * cos_p * sin_p, 0.0)),
        np.sqrt(np.maximum((ea * sin_p) ** 2 + (ed * cos_p) ** 2 - 2 * cov_ad * cos_p * sin_p, 0.0)),
        np.asarray(t["membership_prob"], float), np.asarray(t["g0_mag"], float),
        np.ones(len(t), int) * 2, np.arctan2(x, y),
        mu_a=pmra, mu_d=pmdec, err_a=ea, err_d=ed, err_corr=rho,
        sys_a=exp_a, sys_d=exp_d, mu_sys=MU_SYS, exact=True)


def pristine_field_density(sample: MemberSample, field_r_min_deg: float = 3.6,
                           bw: float = 0.25, step: float = 0.1, vmax: float = 30.0):
    """Empirical 2-D field density from the sample's own outer stars, in absolute PM."""
    from ..selection.field_template import field_density_2d
    a, d = sample.absolute_pm
    t = Table({"r_arcsec": sample.r_arcsec, "mu_a": sample.mu_a, "mu_d": sample.mu_d,
               "sys_a": sample.sys_a, "sys_d": sample.sys_d})
    return field_density_2d(t, bw=bw, step=step, vmax=vmax,
                            r_min_arcsec=field_r_min_deg * 3600.0, absolute=True)


def pristine_profile(edges_deg=(0.30, 0.50, 0.70, 0.90, 1.20, 1.60, 2.10, 2.80),
                     feh_max: float = -1.2, field_r_min_deg: float = 3.6,
                     distance_kpc: float = 5.43, min_stars: int = 30,
                     min_cluster: float = 10.0) -> Table:
    """Cluster dispersion per annulus from the mixture fit, with no proper-motion selection."""
    s = load_pristine_sample(feh_max=feh_max, distance_kpc=distance_kpc)
    dens = pristine_field_density(s, field_r_min_deg=field_r_min_deg)
    k = KMS_PER_MASYR_KPC * distance_kpc
    rows = []
    for lo, hi in zip(edges_deg[:-1], edges_deg[1:]):
        m = (s.r_arcsec >= lo * 3600.0) & (s.r_arcsec < hi * 3600.0)
        if m.sum() < min_stars:
            continue
        o = dispersion_2d(s, m, dens, depth_var=0.0, field_at=s.absolute_pm)
        sig = float(np.sqrt(0.5 * (o["sigma_r"] ** 2 + o["sigma_t"] ** 2)))
        err = float(0.5 * np.hypot(o["sigma_r_err"], o["sigma_t_err"]))
        rpc = float(np.median(s.r_arcsec[m]) * distance_kpc * 1e3 / 206264.806)
        n_cl = (1.0 - o["f"]) * m.sum()
        # below ~10 effective cluster stars the mixture is degenerate: it can put the
        # cluster fraction to zero and the dispersion follows it (the 136 pc bin collapses
        # to exactly 0). Such bins are reported but flagged, never plotted as measurements.
        reliable = bool(n_cl >= min_cluster and sig > 0.0 and err > 0.0)
        rows.append((lo, hi, rpc, int(m.sum()), o["f"], n_cl,
                     o["sigma_r"], o["sigma_t"], sig, err, sig * k, err * k,
                     o["mean_r"], o["mean_t"], reliable))
    return Table(rows=rows, names=("r_lower_deg", "r_upper_deg", "r_pc", "n_stars", "f_field",
                                   "n_cluster", "sigma_pmr", "sigma_pmt", "sigma_pm",
                                   "sigma_pm_err", "sigma_kms", "sigma_kms_err",
                                   "mean_pmr", "mean_pmt", "reliable"),
                 meta={"feh_max": feh_max, "field_r_min_deg": field_r_min_deg,
                       "note": "no proper-motion selection; field is empirical beyond "
                               "field_r_min_deg under the same metallicity cut"})
