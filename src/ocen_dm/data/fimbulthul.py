"""Loader for the Fimbulthul stream candidate members.

Paper: Ibata et al. (2019), https://doi.org/10.1038/s41550-019-0751-x.

Retrieval: CDS/VizieR catalogue ``J/other/NatAs/3.667/tables1``, "Candidate
members of the Fimbulthul stellar stream" (309 rows). This is a CDS-curated
representation of Supplementary Table 1, not the publisher's own file.

Column names below were read from the live VizieR table description on
2026-09-16, not guessed. Note that the catalogue has **no Gaia source_id
column**: 306 of 309 rows carry a ``Gaia DR3 <id>`` string in the CDS-added
``SimbadName``, and three carry CRTS or UCAC4 designations instead. The
identifier role is therefore a string, not an integer.

Likelihood rule (spec section 2.1C): STREAMFINDER selection is complex, so
Fimbulthul enters v1 as a stream-track constraint only, never as an absolute
surface density.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._product import build_product
from .base import find_raw, read_table, standardize
from .schema import validate_table
from .schema import ColumnRole, TableSchema

__all__ = ["DATASET", "SCHEMA", "raw_path", "load", "preprocess"]

DATASET = "fimbulthul_ibata2019"

LIKELIHOOD_RULE = (
    "stream-track constraint only; no absolute density likelihood without a "
    "STREAMFINDER selection model"
)

SCHEMA = TableSchema(
    name="fimbulthul_members",
    unique_id="name",
    notes=(
        "STREAMFINDER candidates from VizieR J/other/NatAs/3.667/tables1. "
        "Track constraint only, no density likelihood."
    ),
    roles=(
        ColumnRole(
            "name", None, ("SimbadName", "Name"), kind="string", is_identifier=True,
            description="CDS cross-identification; 'Gaia DR3 <source_id>' for most rows",
        ),
        ColumnRole("ra", "deg", ("RAJ2000", "ra", "RA_ICRS"), valid_range=(0.0, 360.0),
                   description="right ascension, J2000"),
        ColumnRole("dec", "deg", ("DEJ2000", "dec", "DE_ICRS"), valid_range=(-90.0, 90.0),
                   description="declination, J2000"),
        ColumnRole("pmra", "mas / yr", ("pmRA",), required=False, allow_missing=True,
                   description="Gaia proper motion in RA, pmRA*cos(dec)"),
        ColumnRole("pmra_error", "mas / yr", ("e_pmRA",), required=False,
                   allow_missing=True, valid_range=(0.0, None)),
        ColumnRole("pmdec", "mas / yr", ("pmDE",), required=False, allow_missing=True),
        ColumnRole("pmdec_error", "mas / yr", ("e_pmDE",), required=False,
                   allow_missing=True, valid_range=(0.0, None)),
        ColumnRole("parallax", "mas", ("plx",), required=False, allow_missing=True),
        ColumnRole("parallax_error", "mas", ("e_plx",), required=False,
                   allow_missing=True, valid_range=(0.0, None)),
        ColumnRole("g0_mag", "mag", ("G0mag",), required=False, allow_missing=True,
                   description="dereddened Gaia G"),
        ColumnRole("bp_rp0", "mag", ("__BP-RP_0", "(BP-RP)0"), required=False,
                   allow_missing=True, description="dereddened Gaia BP-RP"),
    ),
)


def raw_path() -> Path:
    """Return the local path of the Fimbulthul candidate table."""
    return find_raw(DATASET, "fimbulthul_supp_table1.*")


def load(path: Path | str | None = None, *, validate: bool = True) -> Any:
    """Load, standardize and validate the Fimbulthul candidate list.

    ``validate=False`` is for inspection only.
    """
    table = read_table(path or raw_path())
    out = standardize(table, SCHEMA, keep_extra_columns=True)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_likelihood_rule"] = LIKELIHOOD_RULE
    if validate:
        validate_table(out, SCHEMA)
    return out


def preprocess() -> dict[str, Any]:
    """Standardize, validate, and write the processed Fimbulthul product."""
    return build_product(
        dataset=DATASET,
        schema=SCHEMA,
        path=raw_path(),
        subdir="tails",
        likelihood_rule=LIKELIHOOD_RULE,
        extra_meta={
            "ocen_vizier_catalogue": "J/other/NatAs/3.667/tables1",
            "ocen_representation": "cds_curated_reproduction",
        },
    )
