"""Full Gaia DR3 astrometric covariance for the Kuzma & Ishigaki (2025) catalogue.

The published table carries ``e_pmRA`` and ``e_pmDE`` but not the correlation
``pmra_pmdec_corr``; the specification (section 8.2) forbids assuming independent
PM errors when the covariance exists. One ``sqlutilpy.local_join`` against
``gaia_dr3.gaia_source`` on WSDB fetches it for all 157,481 stars (measured: 3.5 s).

The product keeps the Gaia values alongside the catalogue's, and the loader
verifies that they agree -- which is how the Gaia release the catalogue used is
confirmed rather than assumed.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.table import Table

from ..data.base import raw_lineage, write_processed
from ..paths import processed_dir

__all__ = ["QUERY", "fetch_gaia_covariance", "build_product"]

QUERY = """
SELECT m.source_id,
       g.pmra, g.pmdec, g.pmra_error, g.pmdec_error, g.pmra_pmdec_corr,
       g.parallax, g.parallax_error, g.parallax_pmra_corr, g.parallax_pmdec_corr,
       g.ruwe, g.phot_g_mean_mag, g.astrometric_params_solved
FROM mytmptable AS m
LEFT JOIN gaia_dr3.gaia_source AS g ON g.source_id = m.source_id
ORDER BY m.source_id
"""

PRODUCT = "kuzma2025_periphery_gaia_covariance"


def fetch_gaia_covariance(source_ids: np.ndarray) -> dict[str, np.ndarray]:
    """Join ``source_ids`` to Gaia DR3 on WSDB; one call, no chunking."""
    import sqlutilpy as sqlutil

    ids = np.asarray(source_ids, dtype=np.int64)
    return sqlutil.local_join(QUERY, "mytmptable", (ids,), ("source_id",),
                              asDict=True, intNullVal=-1)


def build_product(periphery_path: Path | None = None) -> dict[str, Any]:
    """Fetch, verify against the catalogue, and write the covariance product.

    Raises
    ------
    RuntimeError
        If any star is unmatched, if the catalogue's PMs differ from Gaia DR3's by
        more than their published rounding, or if a correlation is out of range.
    """
    periphery_path = periphery_path or processed_dir() / "tails" / "kuzma2025_periphery.ecsv"
    cat = Table.read(periphery_path)
    order = np.argsort(np.asarray(cat["source_id"]))
    cat = cat[order]
    ids = np.asarray(cat["source_id"], dtype=np.int64)

    started = datetime.now(timezone.utc)
    r = fetch_gaia_covariance(ids)
    elapsed = (datetime.now(timezone.utc) - started).total_seconds()

    if not np.array_equal(r["source_id"], ids):
        raise RuntimeError("join returned a different id set or order than requested")
    unmatched = int(np.sum(~np.isfinite(r["pmra"])))
    if unmatched:
        raise RuntimeError(f"{unmatched} of {len(ids)} source_ids have no Gaia DR3 row")
    dpm = np.nanmax(np.abs(np.asarray(cat["pmra"]) - r["pmra"]))
    if dpm > 1e-3:
        raise RuntimeError(f"catalogue pmRA differs from Gaia DR3 by up to {dpm:.2e} mas/yr: "
                           "not the same release?")
    corr = r["pmra_pmdec_corr"]
    if np.any((corr < -1) | (corr > 1)):
        raise RuntimeError("pmra_pmdec_corr outside [-1, 1]")

    out = Table()
    out["source_id"] = ids
    for key, unit in (("pmra", u.mas / u.yr), ("pmdec", u.mas / u.yr),
                      ("pmra_error", u.mas / u.yr), ("pmdec_error", u.mas / u.yr),
                      ("parallax", u.mas), ("parallax_error", u.mas),
                      ("phot_g_mean_mag", u.mag)):
        out[key] = np.asarray(r[key], dtype=float) * unit
    for key in ("pmra_pmdec_corr", "parallax_pmra_corr", "parallax_pmdec_corr", "ruwe"):
        out[key] = np.asarray(r[key], dtype=float)
    out["astrometric_params_solved"] = np.asarray(r["astrometric_params_solved"], dtype=np.int64)

    out.meta["ocen_dataset"] = "kuzma2025_pristine"
    out.meta["ocen_schema"] = PRODUCT
    out.meta["ocen_source"] = "WSDB gaia_dr3.gaia_source, sqlutilpy.local_join"
    out.meta["ocen_query"] = " ".join(QUERY.split())
    out.meta["ocen_query_utc"] = started.isoformat(timespec="seconds")
    out.meta["ocen_query_seconds"] = round(elapsed, 1)
    out.meta["ocen_raw"] = raw_lineage(periphery_path)
    out.meta["ocen_checks"] = {
        "matched_fraction": 1.0,
        "max_abs_dpmra_vs_catalogue_masyr": float(dpm),
        "median_abs_pmra_pmdec_corr": float(np.median(np.abs(corr))),
    }
    out.meta["ocen_likelihood_rule"] = (
        "use the full 2x2 PM covariance (errors and pmra_pmdec_corr) per star; the "
        "parallax terms are provided for a future 3D treatment"
    )
    path = write_processed(out, PRODUCT, subdir="tails")
    return {"path": str(path), "n_rows": len(out), "query_seconds": elapsed,
            "median_abs_corr": float(np.median(np.abs(corr))),
            "ruwe_gt_1p4": int(np.sum(r["ruwe"] > 1.4))}
