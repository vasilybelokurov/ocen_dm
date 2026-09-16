"""One shared, fail-closed path from a raw file to a processed product."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .base import (
    check_against_manifest,
    raw_lineage,
    read_table,
    schema_report,
    standardize,
    write_processed,
)
from .schema import TableSchema

__all__ = ["build_product"]


def build_product(
    *,
    dataset: str,
    schema: TableSchema,
    path: Path,
    hdu: int | str | None = None,
    subdir: str,
    likelihood_rule: str,
    extra_meta: Mapping[str, Any] | None = None,
    keep_extra_columns: bool = True,
    read_kwargs: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Standardize, validate and write one product, in that order.

    Validation happens **before** the product is written: an invalid table
    raises and nothing is published, so a previously valid product is never
    replaced by a bad one.

    Parameters
    ----------
    dataset : str
        Registry key, recorded in the product metadata.
    schema : TableSchema
        Schema to map onto and enforce.
    path : pathlib.Path
        Raw input file.
    hdu : int or str, optional
        FITS extension to read.
    subdir : str
        Subdirectory of ``data/processed`` to write into.
    likelihood_rule : str
        The specification's usage restriction for these data.
    extra_meta : mapping, optional
        Additional metadata to attach.
    keep_extra_columns : bool, optional
        Carry unclaimed source columns through as ``src_*``.
    read_kwargs : mapping, optional
        Reader settings (e.g. ``format``, ``header_start``) recorded in the
        product lineage, so the exact read is reproducible.

    Returns
    -------
    dict
        Schema report plus ``path``, ``dataset`` and ``lineage``.

    Raises
    ------
    SchemaError
        If the raw file changed since it was verified, if a role cannot be
        resolved, or if validation fails.
    """
    lineage = raw_lineage(path, hdu)
    lineage["manifest"] = check_against_manifest(dataset, path, lineage["sha256"])

    read_kwargs = dict(read_kwargs or {})
    lineage["read_kwargs"] = read_kwargs
    table = read_table(path, hdu=hdu, **read_kwargs)
    out = standardize(table, schema, keep_extra_columns=keep_extra_columns)
    out.meta["ocen_dataset"] = dataset
    out.meta["ocen_likelihood_rule"] = likelihood_rule
    out.meta["ocen_raw"] = lineage
    for key, value in dict(extra_meta or {}).items():
        out.meta[key] = value

    report = schema_report(out, schema)          # raises on any violation
    report["path"] = str(write_processed(out, schema.name, subdir=subdir))
    report["dataset"] = dataset
    report["lineage"] = lineage
    return report
