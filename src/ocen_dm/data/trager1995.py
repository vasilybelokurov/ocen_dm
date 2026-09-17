"""Trager, King & Djorgovski (1995) V-band surface-brightness profile of omega Cen.

VizieR J/AJ/109/218/tables, retrieved in full and filtered to ``Name == 'ngc5139'``.
Units are carried by the VOTable itself (read on 2026-09-17): ``logr`` in
log10(arcsec), ``muV`` / ``muVf`` / ``Resid`` in mag/arcsec^2.

The product keeps the measured surface brightness, the authors' Chebyshev fit, the
residual, the per-point weight and the data-set label. It adds ``r`` in arcsec,
derived from ``logr`` -- a unit conversion, not new information.

Only the *shape* of this profile enters the mass model: the stellar mass-to-light
ratio is a free parameter (specification section 3.1), so the photometric
zero-point and the extinction correction cancel in the MGE normalisation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.table import Table

from ._product import build_product
from .base import find_raw
from .schema import ColumnRole, TableSchema

__all__ = ["DATASET", "CLUSTER_KEY", "SCHEMA", "raw_path", "load", "preprocess"]

DATASET = "trager1995_sbp"
CLUSTER_KEY = "ngc5139"

LIKELIHOOD_RULE = (
    "surface-brightness SHAPE only: the M/L nuisance parameter absorbs the zero-point; "
    "heterogeneous data sets -- use the authors' weights, not equal weighting"
)

SCHEMA = TableSchema(
    name="trager1995_ocen_sbp",
    notes="V-band surface-brightness profile of omega Cen, 73 points.",
    roles=(
        ColumnRole("cluster", None, ("Name",), kind="string", required=False),
        ColumnRole("r", "arcsec", ("r_arcsec",), valid_range=(0.0, None),
                   description="radius, 10**logr; the catalogue publishes log10(r/arcsec)"),
        ColumnRole("mu_v", "mag / arcsec2", ("muV",), description="measured V surface brightness"),
        ColumnRole("mu_v_fit", "mag / arcsec2", ("muVf",), required=False, allow_missing=True,
                   description="authors' Chebyshev fit"),
        ColumnRole("residual", "mag / arcsec2", ("Resid",), required=False, allow_missing=True),
        ColumnRole("weight", None, ("Weight",), dimensionless=True, valid_range=(0.0, 1.0)),
        ColumnRole("data_set", None, ("DataSet",), kind="string", required=False),
    ),
)


def raw_path() -> Path:
    return find_raw(DATASET, "trager1995_tables.vot")


def _is_ocen(table: Table) -> np.ndarray:
    return np.char.strip(np.asarray(table["Name"]).astype(str)) == CLUSTER_KEY


def _linear_radius(table: Table) -> Table:
    """Add ``r_arcsec = 10**logr``: a re-expression of the published column, no new information."""
    table = table.copy()
    table["r_arcsec"] = (10.0 ** np.asarray(table["logr"], dtype=float)) * u.arcsec
    return table


TRANSFORM_NOTE = "r_arcsec = 10**logr (catalogue gives log10 of the radius in arcsec)"


def load(validate: bool = True) -> Any:
    from .base import read_table, standardize
    from .schema import validate_table

    table = read_table(raw_path())
    out = standardize(_linear_radius(table[_is_ocen(table)]), SCHEMA, keep_extra_columns=False)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_likelihood_rule"] = LIKELIHOOD_RULE
    if validate:
        validate_table(out, SCHEMA)
    return out


def preprocess() -> dict[str, Any]:
    return build_product(
        dataset=DATASET, schema=SCHEMA, path=raw_path(), subdir="literature",
        likelihood_rule=LIKELIHOOD_RULE, row_filter=_is_ocen, keep_extra_columns=False,
        transform=_linear_radius, transform_note=TRANSFORM_NOTE,
        extra_meta={"ocen_vizier_catalogue": "J/AJ/109/218/tables", "ocen_cluster_key": CLUSTER_KEY},
    )
