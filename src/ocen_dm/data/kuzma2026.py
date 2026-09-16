"""Loader for Kuzma et al. (2026) periphery/tail spectroscopy.

Paper: https://arxiv.org/abs/2605.23474, https://doi.org/10.1093/mnras/stag1147.
The star list is journal supplementary material and must be retrieved by hand
(``docs/MANUAL_DOWNLOADS.md``); ESO programme ``108.22MM.001`` (PI Kuzma) holds
the underlying FLAMES data, which v1 does not re-reduce. The table is not in
VizieR as of 2026-09-16 (checked: MNRAS volume 550 holds two unrelated
catalogues).

Likelihood rule (spec section 2.1B): the sample is target-selected, so
velocities and metallicities are used *conditional on the observed target
positions*. Spectroscopic counts are not a tail surface density unless the
target-selection function is reconstructed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._product import build_product
from .base import find_raw, read_table, standardize
from .schema import validate_table
from .schema import ColumnRole, TableSchema

__all__ = ["DATASET", "SCHEMA", "raw_path", "load", "preprocess"]

DATASET = "kuzma2026_spectroscopy"

LIKELIHOOD_RULE = (
    "targeted sample: condition on target positions; do not treat detection "
    "counts as an unbiased tail surface density"
)

SCHEMA = TableSchema(
    name="kuzma2026_spectroscopy",
    unique_id="source_id",
    notes="Targeted spectroscopy; use conditional on target positions only.",
    roles=(
        # Generic 'ID'/'id' is not a candidate: see kuzma2025.SCHEMA.
        ColumnRole("star_id", None, ("Star_ID", "star_id"), kind="string",
                   required=False, description="observing target name, e.g. OT_T_325"),
        ColumnRole("source_id", None,
                   ("source_id", "SOURCE_ID", "gaia_id", "GaiaDR3", "DR3_Source_ID"),
                   kind="integer", is_identifier=True),
        ColumnRole("ra", "deg", ("ra", "RA", "RAdeg", "RA_ICRS"), valid_range=(0.0, 360.0)),
        ColumnRole("dec", "deg", ("dec", "DEC", "DEdeg", "DE_ICRS", "de"),
                   valid_range=(-90.0, 90.0)),
        ColumnRole("vlos", "km / s", ("vlos", "rv", "RV", "vrad", "HRV", "v_helio", "vhelio"),
                   description="line-of-sight velocity; frame recorded in meta"),
        ColumnRole("vlos_error", "km / s",
                   ("vlos_error", "e_RV", "erv", "rv_error", "e_HRV", "vrad_err"),
                   valid_range=(0.0, None)),
        ColumnRole("feh", None, ("feh", "FeH", "[Fe/H]", "met"),
                   required=False, allow_missing=True, dimensionless=True),
        ColumnRole("feh_error", None, ("feh_error", "e_FeH", "e_feh"),
                   required=False, allow_missing=True, valid_range=(0.0, None),
                   dimensionless=True),
        ColumnRole("pmra", "mas / yr", ("pmra", "pmRA", "PMRA"),
                   required=False, allow_missing=True),
        ColumnRole("pmdec", "mas / yr", ("pmdec", "pmDE", "PMDEC"),
                   required=False, allow_missing=True),
        ColumnRole("g_mag", "mag", ("G", "Gmag", "phot_g_mean_mag"),
                   required=False, allow_missing=True,
                   description="Gaia DR3 G (not dereddened in this table)"),
        ColumnRole("membership_prob", None, ("P_mem", "p_mem"), required=False,
                   allow_missing=True, valid_range=(0.0, 1.0), dimensionless=True,
                   description="membership probability -- see the likelihood rule"),
        ColumnRole("member_flag", None,
                   ("member", "member_flag", "is_member", "memb", "Member"),
                   kind="any", required=False, allow_missing=True),
    ),
)


def read_settings() -> dict[str, Any]:
    """Return the reader settings for this product from ``configs/data.yaml``.

    The delivered file is whitespace-separated ASCII whose column-name line is
    the LAST comment line, which the default reader does not find.
    """
    import yaml

    from ..paths import configs_dir

    path = configs_dir() / "data.yaml"
    if not path.is_file():
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        doc = yaml.safe_load(fh) or {}
    entry = (doc.get(DATASET) or {})
    settings: dict[str, Any] = {}
    if entry.get("format"):
        settings["format"] = str(entry["format"])
    settings.update(dict(entry.get("reader_kwargs") or {}))
    return settings


def raw_path() -> Path:
    """Return the local path of the supplementary member table."""
    return find_raw(DATASET, "kuzma2026_members.*")


def load(path: Path | str | None = None, *, validate: bool = True) -> Any:
    """Load, standardize and validate the spectroscopic member table.

    ``validate=False`` is for inspection only.
    """
    table = read_table(path or raw_path(), **read_settings())
    out = standardize(table, SCHEMA, keep_extra_columns=True)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_likelihood_rule"] = LIKELIHOOD_RULE
    if validate:
        validate_table(out, SCHEMA)
    return out


def preprocess() -> dict[str, Any]:
    """Standardize, validate, and write the processed spectroscopic product."""
    return build_product(
        dataset=DATASET,
        schema=SCHEMA,
        path=raw_path(),
        subdir="tails",
        likelihood_rule=LIKELIHOOD_RULE,
        read_kwargs=read_settings(),
    )
