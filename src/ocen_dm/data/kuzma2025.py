"""Loader for Kuzma & Ishigaki (2025), the Pristine+Gaia Omega Cen periphery.

Data: https://doi.org/10.5281/zenodo.14791430 (``wCen_table.fits``).

The real column names of that file have not been seen (Zenodo returned HTTP 504
throughout 2026-09-16), so the candidate lists below are Gaia archive
conventions only -- https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/
-- and an unresolved role fails loudly rather than guessing.

Likelihood rule (spec section 2.1A): the published membership probability may
be used for diagnostics and sample exploration, but not as a likelihood weight
while the same proper-motion/photometric information also enters the
scientific likelihood.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._product import build_product
from .base import find_raw, read_table, standardize
from .schema import validate_table
from .schema import ColumnRole, TableSchema

__all__ = ["DATASET", "SCHEMA", "raw_path", "load", "preprocess"]

DATASET = "kuzma2025_pristine"

LIKELIHOOD_RULE = (
    "membership_prob is diagnostic only; do not use it as a likelihood weight "
    "alongside the same PM/photometric information"
)

SCHEMA = TableSchema(
    name="kuzma2025_periphery",
    unique_id="source_id",
    notes=(
        "Wide-field Omega Cen periphery catalogue. Membership probability is "
        "diagnostic only; see spec section 2.1A."
    ),
    roles=(
        # 'ID' is deliberately NOT a candidate: a target number must not be able
        # to masquerade as a Gaia identifier.
        ColumnRole(
            "source_id", None, ("source_id", "SOURCE_ID", "sourceid", "gaia_source_id"),
            kind="integer", is_identifier=True, description="Gaia DR3 source_id",
        ),
        ColumnRole("ra", "deg", ("ra", "RA", "RAdeg", "RA_ICRS", "ra_deg"),
                   valid_range=(0.0, 360.0)),
        ColumnRole("dec", "deg", ("dec", "DEC", "DEdeg", "DE_ICRS", "dec_deg", "de"),
                   valid_range=(-90.0, 90.0)),
        ColumnRole("pmra", "mas / yr", ("pmra", "PMRA", "pmRA"),
                   required=False, allow_missing=True),
        ColumnRole("pmdec", "mas / yr", ("pmdec", "PMDEC", "pmDE"),
                   required=False, allow_missing=True),
        ColumnRole("pmra_error", "mas / yr", ("pmra_error", "e_pmRA", "pmra_err"),
                   required=False, allow_missing=True, valid_range=(0.0, None)),
        ColumnRole("pmdec_error", "mas / yr", ("pmdec_error", "e_pmDE", "pmdec_err"),
                   required=False, allow_missing=True, valid_range=(0.0, None)),
        ColumnRole("pmra_pmdec_corr", None, ("pmra_pmdec_corr", "pmra_pmdec_correlation"),
                   required=False, allow_missing=True, valid_range=(-1.0, 1.0),
                   dimensionless=True,
                   description="astrometric correlation; needed for a correct PM likelihood"),
        ColumnRole("parallax", "mas", ("parallax", "plx", "Plx"),
                   required=False, allow_missing=True),
        ColumnRole("parallax_error", "mas", ("parallax_error", "e_Plx"),
                   required=False, allow_missing=True, valid_range=(0.0, None)),
        # The published catalogue carries DEREDDENED magnitudes (G_0, BP_0, RP_0,
        # CaHK_0), so the roles are named for what they are. Calling G_0
        # "phot_g_mean_mag" would mislabel a corrected quantity as a raw one.
        ColumnRole("g0_mag", "mag", ("G_0",), required=False, allow_missing=True,
                   description="dereddened Gaia G"),
        ColumnRole("bp0_mag", "mag", ("BP_0",), required=False, allow_missing=True,
                   description="dereddened Gaia BP"),
        ColumnRole("rp0_mag", "mag", ("RP_0",), required=False, allow_missing=True,
                   description="dereddened Gaia RP"),
        ColumnRole("cahk0_mag", "mag", ("CaHK_0",), required=False, allow_missing=True,
                   description="dereddened Pristine CaHK"),
        ColumnRole("membership_prob", None,
                   ("prob", "p_member", "membership", "memb_prob", "pmemb"),
                   required=False, allow_missing=True, valid_range=(0.0, 1.0),
                   dimensionless=True,
                   description="published membership probability -- DIAGNOSTIC ONLY"),
        ColumnRole("feh", None, ("feh", "FeH", "[Fe/H]", "feh_pristine"),
                   required=False, allow_missing=True, dimensionless=True,
                   description="photometric [Fe/H] from synthetic CaHK, dex"),
        ColumnRole("feh_lo", None, (), required=False, allow_missing=True,
                   dimensionless=True, description="16th percentile of [Fe/H]"),
        ColumnRole("feh_hi", None, (), required=False, allow_missing=True,
                   dimensionless=True, description="84th percentile of [Fe/H]"),
    ),
)


def raw_path() -> Path:
    """Return the local path of ``wCen_table.fits``."""
    return find_raw(DATASET, "wCen_table.fits")


def load(path: Path | str | None = None, hdu: int | str | None = 1,
         *, validate: bool = True) -> Any:
    """Load, standardize and validate the periphery catalogue.

    Set ``validate=False`` only to inspect a table that is known not to
    conform; nothing in the pipeline may consume an unvalidated product.
    """
    table = read_table(path or raw_path(), hdu=hdu)
    out = standardize(table, SCHEMA, keep_extra_columns=True)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_likelihood_rule"] = LIKELIHOOD_RULE
    if validate:
        validate_table(out, SCHEMA)
    return out


def preprocess() -> dict[str, Any]:
    """Standardize, validate, and write the processed periphery product."""
    return build_product(
        dataset=DATASET,
        schema=SCHEMA,
        path=raw_path(),
        hdu=1,
        subdir="tails",
        likelihood_rule=LIKELIHOOD_RULE,
    )
