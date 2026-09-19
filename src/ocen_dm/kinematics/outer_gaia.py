"""Our own Gaia EDR3 proper-motion dispersion profile, built to be error-model independent.

Why this exists, in one paragraph. The published EDR3 profile (Vasiliev & Baumgardt 2021,
https://arxiv.org/abs/2102.09568) is a cubic spline with 2-5 nodes fitted over 0-2400 arcsec,
so five of the eight points our likelihood took from it sit inside 380 arcsec where their own
quality flag passes no star inside 200: an inward continuation, not a measurement. Our first
replacement fixed that but deconvolved the **raw** Gaia uncertainties the catalogue carries,
which are underestimated in crowded fields, and so ran 4-6 per cent high between 460 and 1000
arcsec. Their density-dependent inflation removes that offset but over-corrects the faint end
(:func:`~ocen_dm.kinematics.vb2021_replication.magnitude_consistency`). Neither error model is
right, so this profile avoids depending on one:

* **only stars whose errors are small next to the signal** (``err < err_max_frac * sigma(R)``,
  default 0.4). A rescaling by ``eta`` can then move the deconvolved dispersion by at most
  ``err_max_frac**2 (eta**2 - 1) / 2``, about 2 per cent, by construction;
* the value quoted is the **midpoint** of the raw-error and inflated-error fits, and **half
  their separation is carried as a systematic** added in quadrature to the statistical error.
  The magnitude test shows the truth lies between the two;
* **the inner edge is 300 arcsec**, where the quality flag first leaves anything at all: it
  passes 5 stars inside 300, then 53 in 300-380 and 441 in 380-460. Those are few but they
  are not nothing, and they are the *cleanest* stars in the whole Gaia sample -- errors of
  0.025-0.044 mas/yr against a dispersion near 0.47, a ratio of 5-9 per cent, so the error
  model is irrelevant there by a wide margin. They also fill the gap between HST's last
  point at 311 arcsec and the bulk of the Gaia sample;
* contamination is the two-dimensional empirical field template, and the fitted mean radial
  and tangential motions absorb rotation, which is returned as ``streaming2`` so the Jeans
  model can be compared with the full second moment.
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
        value = 0.5 * (sig["raw"] + sig["eta"])
        sys_err = 0.5 * abs(sig["raw"] - sig["eta"])
        stat = float(0.5 * np.hypot(out["raw"]["sigma_r_err"], out["raw"]["sigma_t_err"]))
        mean_r = 0.5 * (out["raw"]["mean_r"] + out["eta"]["mean_r"])
        mean_t = 0.5 * (out["raw"]["mean_t"] + out["eta"]["mean_t"])
        # per-component values and errors: the midpoint of the two error models, with half
        # their separation added in quadrature exactly as for the combined dispersion
        comp = {}
        for c in ("r", "t"):
            v = 0.5 * (out["raw"][f"sigma_{c}"] + out["eta"][f"sigma_{c}"])
            sy = 0.5 * abs(out["raw"][f"sigma_{c}"] - out["eta"][f"sigma_{c}"])
            comp[c] = (v, float(np.hypot(out["raw"][f"sigma_{c}_err"], sy)))
        rows.append((lo, float(np.median(base.r_arcsec[m])), hi, int(m.sum()), value, stat, sys_err,
                     float(np.hypot(stat, sys_err)), sig["raw"], sig["eta"],
                     comp["r"][0], comp["r"][1], comp["t"][0], comp["t"][1],
                     mean_r, mean_t, 0.5 * (mean_r ** 2 + mean_t ** 2),
                     0.5 * (out["raw"]["f"] + out["eta"]["f"]), float(np.median(eta[m])),
                     float(np.median(base.g_mag[m]))))
    t = Table(rows=rows, names=("r_lower", "r_median", "r_upper", "n_stars", "sigma_pm", "sigma_stat",
                                "sigma_sys", "sigma_pm_err", "sigma_raw", "sigma_eta", "sigma_pmr",
                                "sigma_pmr_err", "sigma_pmt", "sigma_pmt_err", "mean_pmr",
                                "mean_pmt", "streaming2", "f_field", "median_eta", "median_g"))
    t.meta.update({
        "product": PRODUCT, "err_max_frac": err_max_frac, "r_min_arcsec": float(edges[0]),
        "distance_kpc": distance_kpc, "quality_bit": QUALITY_BIT,
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "note": "sigma_pm is the midpoint of the raw-error and eta-inflated fits; sigma_sys is half "
                "their separation; sigma_pm_err adds it in quadrature to the statistical error. "
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
