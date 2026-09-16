"""Shared loader machinery: reading raw tables and writing processed products.

Every processed product carries the lineage needed to regenerate it: the raw
file path, its SHA-256, the HDU read, the resolved column map with source
units, and the package version. Products are written atomically, so an
interrupted run cannot leave a truncated table where a valid one was.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from .. import __version__
from ..paths import configs_dir, processed_dir, provenance_dir, raw_dir
from .schema import (
    ColumnBinding,
    SchemaError,
    TableSchema,
    describe_table,
    load_column_map,
    resolve_columns,
    validate_table,
)

__all__ = [
    "column_map_path",
    "read_table",
    "find_raw",
    "standardize",
    "write_processed",
    "schema_report",
    "raw_lineage",
    "check_against_manifest",
    "atomic_write_text",
]

_FORMAT_BY_SUFFIX = {
    ".fits": "fits",
    ".fit": "fits",
    ".fz": "fits",
    ".csv": "ascii.csv",
    ".tsv": "ascii.tab",
    ".dat": "ascii",
    ".txt": "ascii",
    ".ecsv": "ascii.ecsv",
    ".vot": "votable",
    ".xml": "votable",
    ".parquet": "parquet",
}


def column_map_path() -> Path:
    """Return ``<root>/configs/column_maps.yaml``, resolved on every call.

    Resolved dynamically so that redirecting ``OCEN_DM_ROOT`` cannot leave the
    pipeline reading data from one project and column maps from another.
    """
    return configs_dir() / "column_maps.yaml"


def atomic_write_text(path: Path, text: str) -> Path:
    """Write ``text`` to ``path`` atomically via a temporary file in the same directory."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
    return path


def read_table(path: Path | str, hdu: int | str | None = None,
               format: str | None = None, **kwargs: Any):
    """Read a catalogue file into an :class:`astropy.table.Table`.

    Parameters
    ----------
    path : path-like
        File to read.
    hdu : int or str, optional
        FITS extension; defaults to the first table extension.
    format : str, optional
        Explicit astropy format, overriding the suffix guess. Needed for ASCII
        products whose header is not where the default reader expects it.
    **kwargs
        Passed to :meth:`astropy.table.Table.read`.

    Returns
    -------
    astropy.table.Table

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    from astropy.table import Table

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing; run `ocen fetch-data` first")

    fmt = format or _FORMAT_BY_SUFFIX.get(path.suffix.lower())
    if fmt == "fits" and hdu is not None:
        return Table.read(path, format="fits", hdu=hdu, **kwargs)
    if fmt is not None:
        return Table.read(path, format=fmt, **kwargs)
    return Table.read(path, **kwargs)


def find_raw(dataset: str, pattern: str) -> Path:
    """Return the single raw file of ``dataset`` matching ``pattern``.

    Raises
    ------
    FileNotFoundError
        If nothing matches.
    SchemaError
        If more than one file matches, which would make the product ambiguous.
    """
    directory = raw_dir() / dataset
    matches = sorted(directory.glob(pattern)) if directory.is_dir() else []
    if not matches:
        raise FileNotFoundError(
            f"no file matching {pattern!r} in {directory}; "
            "run `ocen fetch-data` (or see docs/MANUAL_DOWNLOADS.md)"
        )
    if len(matches) > 1:
        raise SchemaError(
            f"{directory}: {pattern!r} matches {[m.name for m in matches]}; "
            "leave exactly one raw file per product"
        )
    return matches[0]


def raw_lineage(path: Path | str, hdu: int | str | None = None) -> dict[str, Any]:
    """Return the identity record of a raw input file.

    Parameters
    ----------
    path : path-like
        Raw file that was read.
    hdu : int or str, optional
        Extension read from it.

    Returns
    -------
    dict
        ``path``, ``bytes``, ``sha256``, ``hdu``, ``read_utc``.
    """
    from ..provenance import file_digest

    path = Path(path)
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": file_digest(path, "sha256"),
        "hdu": hdu,
        "read_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def check_against_manifest(dataset: str, path: Path | str, sha256: str) -> str:
    """Compare a raw file's hash with the one recorded at download time.

    Parameters
    ----------
    dataset : str
        Dataset key.
    path : path-like
        Raw file being read.
    sha256 : str
        Hash just computed for it.

    Returns
    -------
    str
        ``'match'``, or ``'unregistered'`` when the manifest has no entry for
        this file (e.g. a manual download that predates a fetch).

    Raises
    ------
    SchemaError
        If the manifest holds a different hash for the same file: the raw data
        changed after it was verified, so the product would not be reproducible.
    """
    manifest = provenance_dir() / "manifest.json"
    if not manifest.is_file():
        return "unregistered"
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return "unregistered"

    name = Path(path).name
    for entry in payload.get("files", []):
        if entry.get("dataset") != dataset:
            continue
        if Path(entry.get("path", "")).name != name and entry.get("name") != name:
            continue

        recorded = entry.get("sha256")
        anchor = entry.get("last_verified") or {}
        verified_hash = anchor.get("sha256") if anchor else (
            recorded if entry.get("ok") else None
        )

        if recorded and recorded != sha256 and verified_hash != sha256:
            raise SchemaError(
                f"{dataset}/{name}: the raw file has changed since it was "
                f"verified (manifest sha256 {recorded[:12]}..., file now "
                f"{sha256[:12]}...). Re-run `ocen fetch-data` rather than "
                "processing unverified bytes."
            )

        # A recorded hash is not a recorded *pass*: a file whose checksum did
        # not match the publisher must not become trusted input just because it
        # is unchanged since that failure.
        if entry.get("checksum_ok") is False and verified_hash != sha256:
            raise SchemaError(
                f"{dataset}/{name}: this file FAILED checksum verification "
                f"(published {str(entry.get('published_checksum'))[:12]}..., "
                f"observed {str(entry.get('md5'))[:12]}...). Re-fetch it; do not "
                "process a file the publisher's checksum rejects."
            )
        if entry.get("status") == "failed" and not entry.get("ok") and not anchor:
            raise SchemaError(
                f"{dataset}/{name}: the last retrieval attempt failed and no "
                "verified copy is recorded. Re-run `ocen fetch-data`."
            )
        if verified_hash == sha256:
            return "verified"
        if recorded == sha256:
            return "match"
    return "unregistered"


def standardize(
    table: Any,
    schema: TableSchema,
    *,
    overrides: Mapping[str, Any] | None = None,
    keep_extra_columns: bool = False,
):
    """Rename a raw table's columns to schema roles and convert units.

    Unit handling is deliberate: a source column that carries a unit is
    **converted** to the schema's target unit; a source column with no unit
    requires the input unit to be declared in ``configs/column_maps.yaml``.
    Assigning a unit to unlabelled numbers would silently reinterpret them.

    Parameters
    ----------
    table : astropy.table.Table
        Raw table.
    schema : TableSchema
        Target schema.
    overrides : mapping, optional
        Explicit role map; defaults to the entry for ``schema.name`` in
        ``configs/column_maps.yaml``.
    keep_extra_columns : bool, optional
        Also carry through the source columns that no role claimed, prefixed
        with ``src_``, so that information is not discarded before later
        likelihoods need it.

    Returns
    -------
    astropy.table.Table
        Table of role columns in schema units, with the mapping, the source
        units and any omitted columns recorded in ``meta``.

    Raises
    ------
    SchemaError
        If a required role is unresolved, or an input unit is neither carried by
        the column nor declared in the column map.
    """
    import numpy as np
    from astropy import units as u
    from astropy.table import Column, MaskedColumn, Table

    def _convert(column, target):
        """Convert a column to ``target``, preserving any mask.

        ``Column.to()`` drops the mask, which would silently turn a missing
        measurement into a real-looking number.
        """
        mask = getattr(column, "mask", None)
        converted = column.to(target)
        values = np.asarray(converted)
        if mask is not None:
            return MaskedColumn(
                values, mask=np.asarray(mask, dtype=bool), unit=target,
                name=column.info.name, description=column.info.description,
            )
        return Column(values, unit=target, name=column.info.name,
                      description=column.info.description)

    if overrides is None:
        overrides = load_column_map(column_map_path(), schema.name)

    bindings: dict[str, ColumnBinding] = resolve_columns(schema, table.colnames, overrides)

    out = Table()
    source_units: dict[str, str | None] = {}
    for role_name, binding in bindings.items():
        role = binding.role
        column = table[binding.source].copy()

        if role.unit is not None:
            target = u.Unit(role.unit)
            if column.unit is not None:
                source_units[role_name] = str(column.unit)
                if column.unit != target:
                    if not column.unit.is_equivalent(target):
                        raise SchemaError(
                            f"{schema.name}.{role_name}: source column "
                            f"{binding.source!r} has unit {column.unit}, which is "
                            f"not convertible to {role.unit}"
                        )
                    column = _convert(column, target)
            elif binding.declared_unit:
                source_units[role_name] = binding.declared_unit
                declared = u.Unit(binding.declared_unit)
                if not declared.is_equivalent(target):
                    raise SchemaError(
                        f"{schema.name}.{role_name}: declared input unit "
                        f"{binding.declared_unit!r} is not convertible to {role.unit}"
                    )
                column.unit = declared
                if declared != target:
                    column = _convert(column, target)
            else:
                raise SchemaError(
                    f"{schema.name}.{role_name}: source column {binding.source!r} "
                    f"carries no unit and the schema needs {role.unit!r}. Declare "
                    f"the input unit in configs/column_maps.yaml, e.g.\n"
                    f"  {schema.name}:\n"
                    f"    {role_name}: {{column: {binding.source}, unit: <the unit "
                    "printed in the paper or ReadMe>}}\n"
                    "Assigning the target unit to unlabelled numbers would "
                    "silently reinterpret them."
                )
        elif role.dimensionless:
            # A dimensionless role (probability, correlation, ratio) must be a
            # pure number: 0.5 percent is 0.005, not 0.5.
            source_units[role_name] = str(column.unit) if column.unit is not None else None
            if column.unit is not None and column.unit != u.dimensionless_unscaled:
                if not column.unit.is_equivalent(u.dimensionless_unscaled):
                    raise SchemaError(
                        f"{schema.name}.{role_name}: source column "
                        f"{binding.source!r} has unit {column.unit}, but this role "
                        "is a pure number"
                    )
                column = _convert(column, u.dimensionless_unscaled)
            elif binding.declared_unit:
                declared = u.Unit(binding.declared_unit)
                if not declared.is_equivalent(u.dimensionless_unscaled):
                    raise SchemaError(
                        f"{schema.name}.{role_name}: declared unit "
                        f"{binding.declared_unit!r} is not dimensionless"
                    )
                column.unit = declared
                column = _convert(column, u.dimensionless_unscaled)
        else:
            source_units[role_name] = str(column.unit) if column.unit is not None else None

        out[role_name] = column

    claimed = {b.source for b in bindings.values()}
    omitted = [c for c in table.colnames if c not in claimed]
    if keep_extra_columns:
        for name in omitted:
            out[f"src_{name}"] = table[name].copy()

    out.meta.update(dict(table.meta))
    out.meta["ocen_schema"] = schema.name
    out.meta["ocen_package_version"] = __version__
    out.meta["ocen_created_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    out.meta["ocen_source_columns"] = {k: b.source for k, b in bindings.items()}
    out.meta["ocen_source_units"] = source_units
    out.meta["ocen_omitted_columns"] = omitted if not keep_extra_columns else []
    return out


def write_processed(table: Any, name: str, subdir: str | None = None) -> Path:
    """Write a processed product as ECSV, atomically.

    Parameters
    ----------
    table : astropy.table.Table
        Standardized, already validated table.
    name : str
        File stem, e.g. ``'kuzma2025_periphery'``.
    subdir : str, optional
        Subdirectory under ``data/processed``.

    Returns
    -------
    pathlib.Path
        Path written. An existing product is replaced only once the new file is
        complete, so an interrupted write cannot destroy a valid product.
    """
    directory = processed_dir() / subdir if subdir else processed_dir()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}.ecsv"
    tmp = directory / f".{name}.ecsv.tmp{os.getpid()}"
    table.write(tmp, format="ascii.ecsv", overwrite=True)
    os.replace(tmp, path)
    return path


def schema_report(table: Any, schema: TableSchema | None = None) -> dict[str, Any]:
    """Describe a table and validate it against ``schema``.

    Unlike the first implementation, a validation failure is **not** downgraded
    to a report field: it propagates, so an invalid table can never be written
    as a processed product.

    Parameters
    ----------
    table : astropy.table.Table
        Table to describe.
    schema : TableSchema, optional
        If given, the table is validated against it.

    Returns
    -------
    dict
        ``columns`` plus, when a schema is supplied, ``validation``.

    Raises
    ------
    SchemaError
        If validation fails.
    """
    report: dict[str, Any] = {
        "n_rows": len(table),
        "columns": describe_table(table),
    }
    if schema is not None:
        report["validation"] = validate_table(table, schema)
    return report
