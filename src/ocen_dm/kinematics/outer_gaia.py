"""Measure Gaia EDR3 proper-motion dispersions with reduced error-model sensitivity.

Quality-flagged stars are selected with mean component error below 0.4 times a
pilot dispersion. This reduces, but does not eliminate, calibration sensitivity.
The inner edge is 300 arcsec; the profile uses a two-dimensional empirical field
model, full error covariance, and perspective/depth corrections.

The default values and statistical errors are from the raw-error fit. The
eta-inflated fit is stored separately. Half their separation is a diagnostic
column, not an uncertainty added to the default errors and not a demonstrated
bracket on the truth. The component likelihood selects the rotation prescription
separately. Its default is the published curve.

``build_edr3_profile`` writes the product. ``load_edr3_profile`` builds it if
missing. Exact defaults and current counts are in ``docs/data_analysis.tex``.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table

from ..paths import processed_dir
from .outer_profile import dispersion_2d, load_members
from .vb2021_replication import error_inflation, published_profile

__all__ = ["PRODUCT", "R_MIN_ARCSEC", "DEFAULT_EDGES", "build_edr3_profile", "load_edr3_profile"]


#: how many equal-count quantiles of the selected stars' radii to store per bin, so the
#: likelihood can average the model over the radii actually observed rather than over a
#: complete annulus. With coverage that changes inside a bin the two differ: in HST's
#: 300-340 arcsec bin the selection-weighted mean radius is 310.8 arcsec against 319.5 for
#: full-annulus tracer weighting, worth 1.3 statistical errors (Codex review, 2026-09-20).
R_NODES = 8


def _radial_nodes(r: "np.ndarray") -> "np.ndarray":
    """Equal-count quantile midpoints of the selected radii: an unweighted average over
    these reproduces an average over the stars themselves."""
    q = (np.arange(R_NODES) + 0.5) / R_NODES
    return np.quantile(np.asarray(r, float), q)

PRODUCT = "ocen_pm_dispersion_edr3_ours"
#: inside this radius the quality flag passes 5 stars in total: nothing to measure
R_MIN_ARCSEC = 300.0
#: two narrow bins carry the sparse but very precise 300-460 arcsec stars, then log spacing
DEFAULT_EDGES = np.concatenate([[R_MIN_ARCSEC, 380.0], np.geomspace(460.0, 2400.0, 8)])
QUALITY_BIT = 2


def _composite_tracer(distance_kpc: float):
    from ..light_model import build_stellar_mge, fit_mge_projected, load_tracer_profile
    fit = fit_mge_projected(load_tracer_profile("composite"), sigma_range_arcsec=(7.0, 3000.0))
    return build_stellar_mge(fit, distance_kpc, 1.0)


def build_edr3_profile(edges: np.ndarray | None = None, err_max_frac: float = 0.4,
                       distance_kpc: float = 5.43) -> Path:
    """Measure the profile and write ``ocen_pm_dispersion_edr3_ours.ecsv``."""
    from ..selection.field_template import field_density_2d
    from .perspective import depth_dispersion

    edges = DEFAULT_EDGES if edges is None else np.asarray(edges, float)
    cat = Table.read(processed_dir() / "tails" / "vasiliev2021_ocen_members.ecsv")
    qf = np.asarray(cat["quality_flag"], int)
    eta = error_inflation(np.asarray(cat["source_density"], float), (qf & 1) > 0)
    rp, sp = published_profile()

    base = load_members(exact=True, distance_kpc=distance_kpc)
    err = 0.5 * (base.err_r + base.err_t)
    keep = ((qf & QUALITY_BIT) > 0) & (err < err_max_frac * np.interp(base.r_arcsec, rp, sp))
    samples = {"raw": base, "eta": base.scale_errors(eta)}
    tracer = _composite_tracer(distance_kpc)

    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = keep & (base.r_arcsec >= lo) & (base.r_arcsec < hi)
        if m.sum() < 25:
            continue
        # a field template matched to THIS bin's selection: same error ceiling, same
        # magnitude range. Without it the template describes a fainter, worse-measured
        # population than the stars being fitted.
        gm = base.g_mag[m]
        dens2d = field_density_2d(err_max=float(np.max(err[m])),
                                  g_range=(float(np.min(gm)) - 0.5, float(np.max(gm)) + 0.5))
        r_pc = base.r_arcsec[m] * distance_kpc * 1e3 / 206264.806
        sd = depth_dispersion(tracer, r_pc, float(np.hypot(*base.mu_sys)), distance_kpc)
        out = {k: dispersion_2d(s, m, dens2d, depth_var=sd ** 2, field_at=s.absolute_pm)
               for k, s in samples.items()}
        sig = {k: float(np.sqrt(0.5 * (o["sigma_r"] ** 2 + o["sigma_t"] ** 2))) for k, o in out.items()}
        value = sig["raw"]
        sys_err = 0.5 * abs(sig["raw"] - sig["eta"])      # recorded, no longer added
        stat = float(0.5 * np.hypot(out["raw"]["sigma_r_err"], out["raw"]["sigma_t_err"]))
        mean_r = 0.5 * (out["raw"]["mean_r"] + out["eta"]["mean_r"])
        mean_t = 0.5 * (out["raw"]["mean_t"] + out["eta"]["mean_t"])
        # Store each error-model fit separately; the sensitivity is not added in quadrature.
        # Both error models are stored per component. The DEFAULT is the raw one: the
        # catalogue's own uncertainties, unmodified. The density-dependent inflation is kept
        # as an alternative to be refitted and compared, rather than folded into a midpoint
        # with half the gap as a systematic, which is what this did until 2026-09-20. The
        # comparison of two fits answers the question better than a fudged error bar
        # (user's call, 2026-09-20).
        comp = {}
        for c in ("r", "t"):
            comp[c] = (out["raw"][f"sigma_{c}"], out["raw"][f"sigma_{c}_err"],
                       out["eta"][f"sigma_{c}"], out["eta"][f"sigma_{c}_err"])
        nodes = _radial_nodes(base.r_arcsec[m])
        rows.append((lo, float(np.median(base.r_arcsec[m])), hi, int(m.sum()), value, stat, sys_err,
                     stat, sig["raw"], sig["eta"],
                     comp["r"][0], comp["r"][1], comp["t"][0], comp["t"][1],
                     comp["r"][2], comp["r"][3], comp["t"][2], comp["t"][3],
                     mean_r, mean_t, 0.5 * (mean_r ** 2 + mean_t ** 2),
                     0.5 * (out["raw"]["f"] + out["eta"]["f"]), float(np.median(eta[m])),
                     float(np.median(base.g_mag[m])), *nodes))
    t = Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_stars", "sigma_pm", "sigma_stat",
                                "sigma_sys", "sigma_pm_err", "sigma_raw", "sigma_eta", "sigma_pmr",
                                "sigma_pmr_err", "sigma_pmt", "sigma_pmt_err",
                                "sigma_pmr_eta", "sigma_pmr_eta_err", "sigma_pmt_eta",
                                "sigma_pmt_eta_err", "mean_pmr",
                                "mean_pmt", "streaming2", "f_field", "median_eta", "median_g",
                                *[f"r_node{i}" for i in range(R_NODES)]))
    t.meta.update({
        "product": PRODUCT, "err_max_frac": err_max_frac, "r_min_arcsec": float(edges[0]),
        "distance_kpc": distance_kpc, "quality_bit": QUALITY_BIT,
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": "sigma_pm and sigma_pm* are the RAW-error fit (the catalogue's own uncertainties). "
                "The eta-inflated alternative is stored in sigma_*_eta for a comparison run; "
                "sigma_sys records half their separation but is NOT added to the error. "
                "streaming2 = (mean_pmr^2 + mean_pmt^2)/2, to be subtracted from the model second moment.",
        "source": "Vasiliev & Baumgardt 2021 EDR3 member catalogue, arXiv:2102.09568",
    })
    path = processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    t.write(path, format="ascii.ecsv", overwrite=True)
    return path


def load_edr3_profile(path: Path | None = None) -> Table:
    path = path or processed_dir() / "kinematics" / f"{PRODUCT}.ecsv"
    if not path.exists():
        return Table.read(build_edr3_profile())
    return Table.read(path)
