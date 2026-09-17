"""Baumgardt et al. (2019) Gaia DR2 proper-motion dispersion profile of omega Cen.

VizieR J/MNRAS/482/5138/table4, retrieved in full and filtered to ``Name == 'NGC 5139'``.
Nine bins from 179 to 1747 arcsec (4.7-46 pc), i.e. the radii beyond the HST field
where the specification's dark-matter question lives. Errors are asymmetric
(``E_sigma`` upper, ``e_sigma`` lower); units are carried by the VOTable.

Independence: these stars overlap the Vasiliev & Baumgardt (2021) EDR3 sample and,
at the inner bins, the Kuzma (2025) members. A likelihood must not treat this
profile as independent of a per-star likelihood built on either.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from astropy.table import Table

from ._product import build_product
from .base import find_raw
from .schema import ColumnRole, TableSchema

__all__ = ["DATASET", "CLUSTER_KEY", "SCHEMA", "raw_path", "load", "preprocess"]

DATASET = "baumgardt2019_pm_profiles"
CLUSTER_KEY = "NGC 5139"

LIKELIHOOD_RULE = (
    "binned DR2 PM dispersion; stars overlap the EDR3 member sample and the Kuzma 2025 "
    "members -- not independent of likelihoods built on those"
)

SCHEMA = TableSchema(
    name="baumgardt2019_ocen_pm_dispersion",
    notes="Gaia DR2 PM dispersion of omega Cen versus radius, 9 bins, asymmetric errors.",
    roles=(
        ColumnRole("cluster", None, ("Name",), kind="string", required=False),
        ColumnRole("r", "arcsec", ("r",), valid_range=(0.0, None),
                   description="average distance of the bin's stars from the centre"),
        ColumnRole("n_stars", None, ("NPM",), kind="integer", valid_range=(1, None)),
        ColumnRole("sigma_pm", "mas / yr", ("sigma",), valid_range=(0.0, None)),
        ColumnRole("sigma_pm_err_hi", "mas / yr", ("E_sigma",), valid_range=(0.0, None)),
        ColumnRole("sigma_pm_err_lo", "mas / yr", ("e_sigma",), valid_range=(0.0, None)),
    ),
)


def raw_path() -> Path:
    return find_raw(DATASET, "baumgardt2019_table4.vot")


def _is_ocen(table: Table) -> np.ndarray:
    return np.char.strip(np.asarray(table["Name"]).astype(str)) == CLUSTER_KEY


def load(validate: bool = True) -> Any:
    from .base import read_table, standardize
    from .schema import validate_table

    table = read_table(raw_path())
    out = standardize(table[_is_ocen(table)], SCHEMA, keep_extra_columns=False)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_likelihood_rule"] = LIKELIHOOD_RULE
    if validate:
        validate_table(out, SCHEMA)
    return out


def preprocess() -> dict[str, Any]:
    return build_product(
        dataset=DATASET, schema=SCHEMA, path=raw_path(), subdir="kinematics",
        likelihood_rule=LIKELIHOOD_RULE, row_filter=_is_ocen, keep_extra_columns=False,
        extra_meta={"ocen_vizier_catalogue": "J/MNRAS/482/5138/table4", "ocen_cluster_key": CLUSTER_KEY},
    )
