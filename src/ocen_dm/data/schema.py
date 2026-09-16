"""Column-role resolution and table schema validation.

Design rules, tightened after the 2026-09-16 review:

* **No invented metadata.** Roles resolve only against documented candidate
  names (e.g. the Gaia archive data model,
  https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/).
  An unresolved or ambiguous role raises and prints the real columns.
* **Units are converted, never asserted.** A source column carrying a unit is
  converted to the schema's target unit. A source column carrying no unit
  requires an explicitly declared input unit in ``configs/column_maps.yaml``;
  otherwise loading fails. Assigning a unit to unlabelled numbers silently
  reinterprets the values and is not allowed.
* **Identifiers are exact.** A role marked ``is_identifier`` must be integer or
  string. A float column is rejected: IEEE doubles cannot hold Gaia
  ``source_id`` values above 2**53 exactly.
* **Missing data are explicit.** Masked entries and NaNs in a role that does not
  set ``allow_missing`` invalidate the table.
* **Configuration mistakes fail.** Unknown override keys, and two roles mapped
  onto the same source column, are errors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml

__all__ = [
    "ColumnRole",
    "ColumnBinding",
    "TableSchema",
    "SchemaError",
    "resolve_columns",
    "validate_table",
    "describe_table",
    "load_column_map",
]


class SchemaError(ValueError):
    """Raised when a table cannot be mapped onto, or does not satisfy, a schema."""


#: numpy dtype kinds accepted for each declared role kind
_KIND_DTYPES = {
    "float": "fiu",      # integers are acceptable where a float is expected
    "integer": "iu",
    "string": "USO",
    "any": None,
}


@dataclass(frozen=True)
class ColumnRole:
    """A quantity the pipeline needs from a catalogue.

    Attributes
    ----------
    name : str
        Canonical role name used downstream (e.g. ``'ra'``).
    unit : str or None
        Target unit, astropy-parsable. Source values are **converted** to it.
        ``None`` for dimensionless, identifier or string columns.
    candidates : tuple of str
        Documented source-column names, matched case-insensitively.
    required : bool
        Whether resolution failure is an error.
    kind : {'float', 'integer', 'string', 'any'}
        Required dtype class of the source column.
    is_identifier : bool
        Marks a row/object identifier: floats are rejected outright.
    valid_range : tuple of (float or None, float or None) or None
        Inclusive physical bounds enforced after unit conversion.
    allow_missing : bool
        Whether masked entries or NaNs are acceptable in this role.
    dimensionless : bool
        Marks a pure number (probability, correlation, ratio, [Fe/H]). A source
        column carrying a scaled unit such as ``percent`` is converted, not
        accepted at face value.
    description : str
        What the quantity means. Recorded in the product metadata.
    """

    name: str
    unit: str | None = None
    candidates: tuple[str, ...] = ()
    required: bool = True
    kind: str = "float"
    is_identifier: bool = False
    valid_range: tuple[float | None, float | None] | None = None
    allow_missing: bool = False
    dimensionless: bool = False
    description: str = ""


@dataclass(frozen=True)
class ColumnBinding:
    """Resolution of one role onto a source column.

    Attributes
    ----------
    role : ColumnRole
        The role being satisfied.
    source : str
        Source column name, exactly as it appears in the file.
    declared_unit : str or None
        Unit declared in the column map for a source column that carries none.
    """

    role: ColumnRole
    source: str
    declared_unit: str | None = None


@dataclass(frozen=True)
class TableSchema:
    """The set of roles a standardized product must provide.

    Attributes
    ----------
    ordered_triples : tuple of tuple of str
        Role-name triples ``(lower, middle, upper)`` that must satisfy
        ``lower <= middle <= upper`` row by row, e.g. radial bin edges.
    """

    name: str
    roles: tuple[ColumnRole, ...]
    unique_id: str | None = None
    notes: str = ""
    ordered_triples: tuple[tuple[str, str, str], ...] = ()

    @property
    def required_roles(self) -> tuple[str, ...]:
        """Names of the roles that must resolve."""
        return tuple(r.name for r in self.roles if r.required)

    @property
    def role_names(self) -> tuple[str, ...]:
        """Names of every role in the schema."""
        return tuple(r.name for r in self.roles)

    def role(self, name: str) -> ColumnRole:
        """Return the role called ``name``."""
        for candidate in self.roles:
            if candidate.name == name:
                return candidate
        raise KeyError(f"schema {self.name!r} has no role {name!r}")


def load_column_map(path: Path | str, table: str) -> dict[str, dict[str, Any]]:
    """Load an explicit ``role -> {column, unit}`` map from a YAML config.

    Two entry forms are accepted::

        ra: RA_ICRS                        # source column only
        r:  {column: R, unit: arcmin}      # source column plus its input unit

    Parameters
    ----------
    path : path-like
        YAML file structured as ``{schema_name: {role: entry}}``.
    table : str
        Schema name to read.

    Returns
    -------
    dict
        ``role -> {'column': str, 'unit': str or None}``. Empty when the file
        or the key is absent.

    Raises
    ------
    SchemaError
        If an entry is neither a string nor a mapping with a ``column`` key.
    """
    path = Path(path)
    if not path.is_file():
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        doc = yaml.safe_load(fh) or {}
    entry = (doc.get(table) or {}) if isinstance(doc, dict) else {}

    out: dict[str, dict[str, Any]] = {}
    for role, value in entry.items():
        if value is None:
            continue
        if isinstance(value, str):
            out[str(role)] = {"column": value, "unit": None}
        elif isinstance(value, Mapping) and "column" in value:
            out[str(role)] = {
                "column": str(value["column"]),
                "unit": str(value["unit"]) if value.get("unit") else None,
            }
        else:
            raise SchemaError(
                f"{path}: entry for {table}.{role} must be a column name or a "
                f"mapping with a 'column' key, got {value!r}"
            )
    return out


def resolve_columns(
    schema: TableSchema,
    columns: Sequence[str],
    overrides: Mapping[str, Any] | None = None,
) -> dict[str, ColumnBinding]:
    """Map schema roles onto the actual column names of a table.

    Parameters
    ----------
    schema : TableSchema
        Roles to resolve.
    columns : sequence of str
        Column names present in the source table.
    overrides : mapping, optional
        Explicit ``role -> column`` or ``role -> {'column':…, 'unit':…}`` pairs,
        taking precedence over candidate matching.

    Returns
    -------
    dict
        ``role name -> ColumnBinding`` for every role that resolved.

    Raises
    ------
    SchemaError
        If a required role cannot be resolved; if a role matches more than one
        distinct column; if an override names an unknown role or a missing
        column; or if two roles resolve onto the same source column.
    """
    normalized: dict[str, dict[str, Any]] = {}
    for role_name, value in dict(overrides or {}).items():
        if isinstance(value, str):
            normalized[role_name] = {"column": value, "unit": None}
        elif isinstance(value, Mapping) and "column" in value:
            normalized[role_name] = {
                "column": str(value["column"]),
                "unit": str(value["unit"]) if value.get("unit") else None,
            }
        else:
            raise SchemaError(
                f"{schema.name}: override for role {role_name!r} must be a "
                f"column name or a mapping with a 'column' key, got {value!r}"
            )

    unknown = sorted(set(normalized) - set(schema.role_names))
    if unknown:
        raise SchemaError(
            f"{schema.name}: column map names unknown role(s) {unknown}; "
            f"schema roles are {sorted(schema.role_names)}. A misspelled key "
            "would otherwise be ignored silently."
        )

    lookup = {str(c).lower(): str(c) for c in columns}
    if len(lookup) != len(list(columns)):
        raise SchemaError(f"{schema.name}: source table has case-colliding column names")

    bindings: dict[str, ColumnBinding] = {}
    unresolved: list[str] = []

    for role in schema.roles:
        if role.name in normalized:
            wanted = normalized[role.name]["column"]
            if wanted.lower() not in lookup:
                raise SchemaError(
                    f"{schema.name}: column map sends role {role.name!r} to "
                    f"{wanted!r}, which is not in the table"
                )
            bindings[role.name] = ColumnBinding(
                role, lookup[wanted.lower()], normalized[role.name]["unit"]
            )
            continue

        # Candidate lists carry case variants of one name, so deduplicate by the
        # resolved source column: only genuinely distinct columns are ambiguous.
        hits: list[str] = []
        for candidate in role.candidates:
            match = lookup.get(candidate.lower())
            if match is not None and match not in hits:
                hits.append(match)

        if len(hits) == 1:
            bindings[role.name] = ColumnBinding(role, hits[0], None)
        elif len(hits) > 1:
            raise SchemaError(
                f"{schema.name}: role {role.name!r} is ambiguous -- candidates "
                f"{hits} all present; resolve it in configs/column_maps.yaml"
            )
        elif role.required:
            unresolved.append(role.name)

    if unresolved:
        raise SchemaError(
            f"{schema.name}: could not resolve required role(s) {unresolved}. "
            f"Table columns are {sorted(map(str, columns))}. Add explicit entries "
            "to configs/column_maps.yaml -- do not guess."
        )

    reused: dict[str, list[str]] = {}
    for role_name, binding in bindings.items():
        reused.setdefault(binding.source, []).append(role_name)
    collisions = {src: names for src, names in reused.items() if len(names) > 1}
    if collisions:
        raise SchemaError(
            f"{schema.name}: one source column cannot serve several roles: "
            f"{collisions}. Two roles reading the same numbers is a mapping "
            "mistake, not a valid configuration."
        )

    return bindings


def _mask_of(column: Any) -> Any:
    """Return a boolean mask array for ``column``, or ``None`` when unmasked."""
    import numpy as np

    mask = getattr(column, "mask", None)
    if mask is None:
        return None
    mask = np.asarray(mask, dtype=bool)
    return mask if mask.any() else None


def validate_table(table: Any, schema: TableSchema) -> dict[str, Any]:
    """Check a standardized table against its schema.

    Parameters
    ----------
    table : astropy.table.Table
        Table whose columns are already named by role and in schema units.
    schema : TableSchema
        Schema to enforce.

    Returns
    -------
    dict
        Report with ``n_rows``, ``units``, ``dtypes``, ``unique_id_ok``,
        ``n_nonfinite``, ``n_masked``.

    Raises
    ------
    SchemaError
        On a missing required role, an empty table, a wrong unit or dtype, a
        float identifier, a value outside ``valid_range``, missing data in a
        role that does not allow it, or duplicate identifiers.
    """
    import numpy as np
    from astropy import units as u

    missing = [r for r in schema.required_roles if r not in table.colnames]
    if missing:
        raise SchemaError(f"{schema.name}: standardized table is missing {missing}")

    if len(table) == 0:
        raise SchemaError(
            f"{schema.name}: table has no rows; an empty product is never a "
            "valid scientific input"
        )

    units: dict[str, str | None] = {}
    dtypes: dict[str, str] = {}
    nonfinite: dict[str, int] = {}
    masked_counts: dict[str, int] = {}

    for role in schema.roles:
        if role.name not in table.colnames:
            continue
        column = table[role.name]
        units[role.name] = str(column.unit) if column.unit is not None else None
        dtypes[role.name] = str(column.dtype)

        if role.unit is not None:
            if column.unit is None:
                raise SchemaError(
                    f"{schema.name}.{role.name}: expected unit {role.unit!r}, column has none"
                )
            if column.unit != u.Unit(role.unit):
                raise SchemaError(
                    f"{schema.name}.{role.name}: unit is {column.unit}, expected "
                    f"{role.unit}. standardize() converts to the schema unit, so a "
                    "mismatch here means the product was not built by this pipeline."
                )

        allowed = _KIND_DTYPES.get(role.kind)
        if allowed is not None and column.dtype.kind not in allowed:
            raise SchemaError(
                f"{schema.name}.{role.name}: dtype kind {column.dtype.kind!r} is "
                f"not valid for a {role.kind!r} role"
            )
        if role.is_identifier:
            if column.dtype.kind == "f":
                raise SchemaError(
                    f"{schema.name}.{role.name}: identifiers must not be floating "
                    "point -- doubles cannot represent Gaia source_id above 2**53 "
                    "exactly. Read the column as int64 or as a string."
                )
            if column.dtype.kind == "O":
                # dtype 'O' hides the real content, including floats.
                offenders = [
                    v for v in np.asarray(column).ravel()[:1000]
                    if not isinstance(v, (str, bytes, int, np.integer))
                ]
                if offenders:
                    raise SchemaError(
                        f"{schema.name}.{role.name}: object-dtype identifier "
                        f"contains non-identifier values, e.g. {offenders[0]!r}"
                    )
            if column.dtype.kind in "USO":
                blank = [
                    v for v in np.asarray(column).ravel()
                    if isinstance(v, (str, bytes)) and not str(v).strip()
                ]
                if blank:
                    raise SchemaError(
                        f"{schema.name}.{role.name}: {len(blank)} blank identifier(s)"
                    )

        mask = _mask_of(column)
        masked_counts[role.name] = int(mask.sum()) if mask is not None else 0

        values = np.asarray(getattr(column, "data", column))
        if column.dtype.kind == "f":
            with np.errstate(invalid="ignore"):
                finite = np.isfinite(values)
            nonfinite[role.name] = int(np.count_nonzero(~finite))
        else:
            finite = np.ones(len(values), dtype=bool)

        # A row is missing if it is masked OR non-finite; the two can overlap,
        # so count the union rather than the sum.
        invalid = ~finite
        if mask is not None:
            invalid = invalid | mask
        bad_missing = int(np.count_nonzero(invalid))
        if not role.allow_missing and bad_missing:
            raise SchemaError(
                f"{schema.name}.{role.name}: {bad_missing} missing value(s) "
                f"(masked or non-finite) in a role that requires complete data"
            )

        if role.valid_range is not None and column.dtype.kind in "fiu":
            low, high = role.valid_range
            usable = ~invalid
            sample = values[usable]
            if low is not None and sample.size and np.nanmin(sample) < low:
                raise SchemaError(
                    f"{schema.name}.{role.name}: minimum {np.nanmin(sample)} is "
                    f"below the physical bound {low} {role.unit or ''}".rstrip()
                )
            if high is not None and sample.size and np.nanmax(sample) > high:
                raise SchemaError(
                    f"{schema.name}.{role.name}: maximum {np.nanmax(sample)} is "
                    f"above the physical bound {high} {role.unit or ''}".rstrip()
                )

    for lower, middle, upper in schema.ordered_triples:
        if not {lower, middle, upper} <= set(table.colnames):
            continue
        low = np.asarray(table[lower])
        mid = np.asarray(table[middle])
        high = np.asarray(table[upper])
        bad = int(np.count_nonzero(~((low <= mid) & (mid <= high))))
        if bad:
            raise SchemaError(
                f"{schema.name}: {bad} row(s) violate {lower} <= {middle} <= "
                f"{upper}; the bin geometry is inconsistent"
            )
        order = np.argsort(low)
        if np.any(high[order][:-1] - low[order][1:] > 1e-9 * np.abs(high[order][:-1])):
            raise SchemaError(
                f"{schema.name}: radial bins defined by {lower}/{upper} overlap"
            )

    unique_ok: bool | None = None
    if schema.unique_id:
        if schema.unique_id not in table.colnames:
            raise SchemaError(
                f"{schema.name}: identifier column {schema.unique_id!r} is absent"
            )
        column = table[schema.unique_id]
        if _mask_of(column) is not None:
            raise SchemaError(f"{schema.name}: identifier {schema.unique_id!r} has masked rows")
        values = np.asarray(column)
        unique_ok = len(np.unique(values)) == len(values)
        if not unique_ok:
            raise SchemaError(
                f"{schema.name}: id column {schema.unique_id!r} has duplicates "
                f"({len(values) - len(np.unique(values))} repeats)"
            )

    return {
        "schema": schema.name,
        "n_rows": len(table),
        "roles": list(schema.required_roles),
        "units": units,
        "dtypes": dtypes,
        "unique_id_ok": unique_ok,
        "n_nonfinite": nonfinite,
        "n_masked": masked_counts,
    }


def describe_table(table: Any, max_columns: int | None = None) -> list[dict[str, Any]]:
    """Summarize a table's columns for a schema report.

    Parameters
    ----------
    table : astropy.table.Table
        Any table.
    max_columns : int, optional
        Truncate the report after this many columns.

    Returns
    -------
    list of dict
        One entry per column: ``name``, ``dtype``, ``unit``, ``description``.
    """
    rows = []
    for name in table.colnames[:max_columns]:
        column = table[name]
        rows.append(
            {
                "name": name,
                "dtype": str(column.dtype),
                "unit": str(column.unit) if getattr(column, "unit", None) else None,
                "description": str(getattr(column, "description", "") or ""),
            }
        )
    return rows
