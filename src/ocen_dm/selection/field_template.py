"""An independent Gaia DR3 field sample, as the contamination template for the outskirts.

At 1800-2400 arcsec the Vasiliev & Baumgardt member catalogue is 85-90 per cent
Galactic field, so the dispersion of the cluster there depends entirely on how the
field is modelled. Fitting a two-Gaussian field to the same stars does not work:
the flexible component absorbs part of the cluster peak and biases sigma low
(JOURNAL 2026-09-18). This module fetches the field **from an annulus outside the
cluster** -- 0.75 to 1.6 degrees from the centre, where the cluster contributes
nothing -- and the measurement then uses its observed PM distribution as a fixed
shape with only the normalisation free.

The template is transformed into the same frame as the cluster stars (systemic
velocity projected star by star, then radial/tangential about the centre), so it
can be used component by component.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table

from ..paths import processed_dir

__all__ = ["QUERY", "PRODUCT", "fetch_field", "build_product", "load_field_template", "field_kde",
           "field_density_2d"]

#: centre and annulus of the field sample (deg)
from ..cluster import OCEN_DEC, OCEN_RA  # single definition, derived from the data
R_INNER_DEG, R_OUTER_DEG = 0.75, 1.6

QUERY = """
SELECT ra, dec, pmra, pmdec, pmra_error, pmdec_error, pmra_pmdec_corr,
       phot_g_mean_mag, ruwe, astrometric_excess_noise,
       q3c_dist(ra, dec, {ra0}, {dec0}) * 3600 AS r_arcsec
FROM gaia_dr3.gaia_source
WHERE q3c_radial_query(ra, dec, {ra0}, {dec0}, {r_out})
  AND q3c_dist(ra, dec, {ra0}, {dec0}) > {r_in}
  AND pmra IS NOT NULL AND pmra <> 'NaN'::float8
  AND phot_g_mean_mag < {g_max} AND parallax_over_error < 5
"""

PRODUCT = "ocen_field_template_dr3"


def fetch_field(g_max: float = 20.5) -> dict[str, np.ndarray]:
    """One WSDB cone query for the field annulus (measured: 226k rows in 3.3 s)."""
    import sqlutilpy as sqlutil

    q = QUERY.format(ra0=OCEN_RA, dec0=OCEN_DEC, r_in=R_INNER_DEG, r_out=R_OUTER_DEG, g_max=g_max)
    return sqlutil.get(q, asDict=True)


def build_product(g_max: float = 20.5, ruwe_max: float = 1.4, excess_noise_max: float = 1.0,
                  pm_error_max: float = 1.5, distance_kpc: float = 5.43, v_los: float = 232.6,
                  mu_sys: tuple[float, float] | None = None) -> Path:
    """Fetch the field, put it in the cluster's radial/tangential frame, and write the product.

    The quality cuts mirror what a member catalogue keeps (``ruwe``, astrometric excess
    noise, a PM-error ceiling); the systemic-motion subtraction is the same exact
    projection used for the members, so the template lives in the same coordinates.
    """
    from ..kinematics.outer_profile import load_members
    from ..kinematics.perspective import systemic_pm_field, systemic_velocity_vector

    d = fetch_field(g_max=g_max)
    if mu_sys is None:
        mu_sys = load_members(distance_kpc=distance_kpc).mu_sys
    v_sys = systemic_velocity_vector(OCEN_RA, OCEN_DEC, mu_sys[0], mu_sys[1], v_los, distance_kpc)
    exp_a, exp_d = systemic_pm_field(d["ra"], d["dec"], v_sys, distance_kpc)
    x = (d["ra"] - OCEN_RA) * np.cos(np.radians(OCEN_DEC)) * 3600.0
    y = (d["dec"] - OCEN_DEC) * 3600.0
    r = np.hypot(x, y); cos_p, sin_p = x / r, y / r
    pmra = d["pmra"] - exp_a; pmdec = d["pmdec"] - exp_d
    keep = ((d["ruwe"] < ruwe_max) & (d["astrometric_excess_noise"] < excess_noise_max)
            & (0.5 * (d["pmra_error"] + d["pmdec_error"]) < pm_error_max))
    t = Table({
        "r_arcsec": r[keep],
        "sys_a": exp_a[keep],
        "sys_d": exp_d[keep],
        "mu_a": pmra[keep],
        "mu_d": pmdec[keep],
        "mu_r": (pmra * cos_p + pmdec * sin_p)[keep],
        "mu_t": (-pmra * sin_p + pmdec * cos_p)[keep],
        "pm_error": (0.5 * (d["pmra_error"] + d["pmdec_error"]))[keep],
        "g_mag": d["phot_g_mean_mag"][keep],
    })
    area = np.pi * ((R_OUTER_DEG * 60) ** 2 - (R_INNER_DEG * 60) ** 2)
    t.meta.update({
        "product": PRODUCT, "source": "gaia_dr3.gaia_source via WSDB (q3c; Koposov & Bartunov 2006)",
        "query": QUERY.format(ra0=OCEN_RA, dec0=OCEN_DEC, r_in=R_INNER_DEG, r_out=R_OUTER_DEG, g_max=g_max),
        "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_fetched": int(len(r)), "n_kept": int(keep.sum()),
        "cuts": {"g_max": g_max, "ruwe_max": ruwe_max, "excess_noise_max": excess_noise_max,
                 "pm_error_max": pm_error_max, "parallax_over_error_max": 5},
        "annulus_deg": [R_INNER_DEG, R_OUTER_DEG], "area_arcmin2": float(area),
        "surface_density_per_arcmin2": float(keep.sum() / area),
        "frame": "systemic 3-D velocity projected star by star and subtracted, then radial/tangential "
                 "about the cluster centre; mu_sys = (%.4f, %.4f) mas/yr, v_los = %.1f km/s, D = %.2f kpc"
                 % (mu_sys[0], mu_sys[1], v_los, distance_kpc),
    })
    path = processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def load_field_template(path: Path | None = None) -> Table:
    """The field template; raises with the build command if it has not been fetched."""
    path = path or processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} missing -- run `ocen fetch-field-template` (needs WSDB access)")
    return Table.read(path)


def field_density_2d(table: Table | None = None, bw: float = 0.15, step: float = 0.05,
                     vmax: float = 30.0, r_min_arcsec: float = 3600.0, absolute: bool = True,
                     err_max: float | None = None, g_range: tuple[float, float] | None = None,
                     min_template: int = 3000):
    """Two-dimensional empirical density of the field proper motions, per (mas/yr)^2.

    Built in the equatorial frame after the systemic field has been removed -- the frame in
    which the field distribution does not depend on where the star sits in the annulus --
    by binning and Gaussian smoothing, and returned as a bilinear interpolator. The field
    is measured beyond ``r_min_arcsec`` so that the cluster's own outskirts and tails do not
    enter the template (inside 1 degree the core fraction is 20 per cent higher, which is
    the cluster, JOURNAL 2026-09-18).

    The two-dimensional form matters: the field's widths are unequal (3.1 mas/yr in
    alpha*, 2.0 in delta) and strongly non-Gaussian (kurtosis 3-17 even in 0.5 degree
    cells), so projecting it onto each star's radial direction mixes the two widths around
    the annulus and produces shoulders that no one-dimensional model reproduces.
    """
    from scipy.ndimage import gaussian_filter

    t = table if table is not None else load_field_template()
    keep = np.asarray(t["r_arcsec"], float) >= r_min_arcsec
    # Selection matching. The template's own magnitude and error distribution differ from
    # the sample being fitted -- its median PM error is 0.28 mas/yr against a science
    # ceiling of 0.09 in the outermost Gaia bin -- and that changes both the error
    # convolution and, more importantly, which Galactic population the template describes.
    # Restricting it to the same errors and magnitudes is the fix (Codex review, 2026-09-19).
    matched = keep.copy()
    if err_max is not None and "pm_error" in t.colnames:
        matched &= np.asarray(t["pm_error"], float) <= err_max
    if g_range is not None and "g_mag" in t.colnames:
        g = np.asarray(t["g_mag"], float)
        matched &= (g >= g_range[0]) & (g <= g_range[1])
    n_matched = int(matched.sum())
    if n_matched >= min_template:
        keep = matched
    selection_matched = bool(n_matched >= min_template)
    a = np.asarray(t["mu_a"], float)[keep]; d = np.asarray(t["mu_d"], float)[keep]
    if absolute:
        # The field has no systemic motion of its own: its distribution is position
        # independent in ABSOLUTE proper motion, so the template is built there and scored
        # at each target star's absolute proper motion. (Measured difference from doing it
        # in the residual frame: <= 0.1 per cent in sigma, because the perspective term is
        # radial and cancels around an annulus -- JOURNAL 2026-09-18.)
        a = a + np.asarray(t["sys_a"], float)[keep]
        d = d + np.asarray(t["sys_d"], float)[keep]
    edges = np.arange(-vmax, vmax + step, step)
    centres = 0.5 * (edges[1:] + edges[:-1])
    h, _, _ = np.histogram2d(a, d, bins=[edges, edges])
    dens = gaussian_filter(h, bw / step, mode="constant")
    dens /= max(dens.sum() * step * step, 1e-300)
    floor = 1e-12

    def density(x, y):
        from scipy.ndimage import map_coordinates
        xi = (np.asarray(x, float) - centres[0]) / step
        yi = (np.asarray(y, float) - centres[0]) / step
        v = map_coordinates(dens, np.vstack([xi.ravel(), yi.ravel()]), order=1, mode="constant", cval=floor)
        return np.maximum(v.reshape(np.shape(x)), floor)

    density.n_template = int(keep.sum())
    density.bandwidth = bw
    density.selection_matched = selection_matched
    density.n_matched = n_matched
    density.err_max = err_max
    density.g_range = g_range
    return density


def field_kde(component: str = "r", table: Table | None = None, bw_method: float = 0.05,
              g_range: tuple[float, float] = (-np.inf, np.inf), grid_step: float = 0.01,
              grid_max: float = 60.0):
    """Kernel density of the field PM distribution in the cluster frame (per mas/yr).

    Evaluated once on a regular grid by binning and Gaussian smoothing -- exactly a
    Gaussian KDE up to the bin width -- and returned as an interpolating callable, so
    scoring 30,000 stars costs microseconds instead of the minutes a direct
    ``gaussian_kde`` would take against a 192,000-star sample.

    The bandwidth is in units of the sample's standard deviation; 0.05 of ~6 mas/yr is
    0.3 mas/yr, comparable to the per-star errors, and the result is insensitive to it
    (sigma moves by 0.003 between 0.03 and 0.10).
    """
    from scipy.ndimage import gaussian_filter1d

    t = table if table is not None else load_field_template()
    v = np.asarray(t[f"mu_{component}"], float)
    g = np.asarray(t["g_mag"], float)
    v = v[(g >= g_range[0]) & (g <= g_range[1])]
    bw = bw_method * float(np.std(v))                      # scipy's scott/silverman convention
    edges = np.arange(-grid_max, grid_max + grid_step, grid_step)
    centres = 0.5 * (edges[1:] + edges[:-1])
    counts, _ = np.histogram(v, edges)
    dens = gaussian_filter1d(counts.astype(float), bw / grid_step, mode="constant")
    dens /= max(dens.sum() * grid_step, 1e-300)
    floor = 1e-12

    def pdf(x: np.ndarray) -> np.ndarray:
        return np.maximum(np.interp(np.asarray(x, float), centres, dens, left=floor, right=floor), floor)

    pdf.bandwidth = bw
    pdf.n_template = len(v)
    return pdf
