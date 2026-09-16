"""Systemic parameters of omega Cen from the Baumgardt & Hilker database compilation.

Working file: ``~/data/catalogues/gc_catalog_updated.fits``, the user's compilation
(mtime 2025-05-16) of https://people.smp.uq.edu.au/HolgerBaumgardt/globular/. The
compilation records neither the database version nor the column units, so units
are declared in ``configs/column_maps.yaml`` from the database's conventions and
the whole product is flagged UNVERIFIED against upstream until re-checked.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table

from ..provenance import load_registry
from ._product import build_product
from .schema import ColumnRole, TableSchema

__all__ = ["DATASET", "SCHEMA", "catalogue_path", "load", "preprocess"]

DATASET = "baumgardt_gc_catalogue"

SCHEMA = TableSchema(
    name="baumgardt_ocen_parameters",
    notes="One row: omega Cen's systemic and structural parameters.",
    roles=(
        ColumnRole("name", None, ("NAME",), kind="string"),
        ColumnRole("ra", "deg", ("RA",), valid_range=(0.0, 360.0)),
        ColumnRole("dec", "deg", ("DEC",), valid_range=(-90.0, 90.0)),
        ColumnRole("pmra", "mas / yr", ("PMRA",)),
        ColumnRole("pmdec", "mas / yr", ("PMDEC",)),
        ColumnRole("pmra_error", "mas / yr", ("PMRA_ERR",), valid_range=(0.0, None)),
        ColumnRole("pmdec_error", "mas / yr", ("PMDEC_ERR",), valid_range=(0.0, None)),
        ColumnRole("distance", "kpc", ("DIST",), valid_range=(0.0, None)),
        ColumnRole("distance_error", "kpc", ("DIST_ERR",), valid_range=(0.0, None)),
        ColumnRole("rv", "km / s", ("RV",)),
        ColumnRole("rv_error", "km / s", ("RV_ERR",), valid_range=(0.0, None)),
        ColumnRole("mass", "Msun", ("MASS",), valid_range=(0.0, None)),
        ColumnRole("r_h", "pc", ("RH",), valid_range=(0.0, None),
                   description="RH column; projected half-light or 3D half-mass is NOT "
                               "recorded in the compilation -- verify against the database"),
        ColumnRole("r_c", "pc", ("RC",), valid_range=(0.0, None), required=False,
                   description="RC column; core radius, unit assumed pc"),
        ColumnRole("r_peri", "kpc", ("RPERI",), required=False, allow_missing=True),
        ColumnRole("r_apo", "kpc", ("RAPO",), required=False, allow_missing=True),
        ColumnRole("feh", None, ("FEH",), dimensionless=True, required=False, allow_missing=True),
    ),
)


def catalogue_path() -> Path:
    return load_registry()[DATASET].files[0].path


def _is_ocen(table: Table) -> np.ndarray:
    names = np.char.upper(np.char.strip(np.asarray(table["NAME"]).astype(str)))
    return np.char.find(names, "5139") >= 0


def load(validate: bool = True) -> Any:
    from .base import read_table, standardize
    from .schema import validate_table

    table = read_table(catalogue_path())
    table = table[_is_ocen(table)]
    out = standardize(table, SCHEMA, keep_extra_columns=True)
    out.meta["ocen_dataset"] = DATASET
    if validate:
        validate_table(out, SCHEMA)
    return out


def preprocess() -> dict[str, Any]:
    return build_product(
        dataset=DATASET,
        schema=SCHEMA,
        path=catalogue_path(),
        subdir="literature",
        likelihood_rule=(
            "compiled literature values, UNVERIFIED against the live database; use as "
            "priors and cross-checks, not as data in the likelihood"
        ),
        row_filter=_is_ocen,
        keep_extra_columns=True,
    )
