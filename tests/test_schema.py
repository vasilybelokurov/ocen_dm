"""Tests for column-role resolution and schema validation."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table

from ocen_dm.data.schema import (
    ColumnRole,
    SchemaError,
    TableSchema,
    describe_table,
    load_column_map,
    resolve_columns,
    validate_table,
)

SCHEMA = TableSchema(
    name="demo",
    unique_id="source_id",
    roles=(
        ColumnRole("source_id", None, ("source_id", "ID"), kind="integer"),
        ColumnRole("ra", "deg", ("ra", "RAdeg")),
        ColumnRole("dec", "deg", ("dec", "DEdeg")),
        ColumnRole("vlos", "km / s", ("rv", "HRV"), required=False),
    ),
)


def test_resolution_is_case_insensitive_and_ordered():
    mapping = resolve_columns(SCHEMA, ["ID", "RAdeg", "DEdeg", "junk"])
    assert {k: b.source for k, b in mapping.items()} == {
        "source_id": "ID", "ra": "RAdeg", "dec": "DEdeg"
    }


def test_missing_required_role_raises_and_lists_the_real_columns():
    with pytest.raises(SchemaError) as excinfo:
        resolve_columns(SCHEMA, ["alpha", "delta"])
    message = str(excinfo.value)
    assert "ra" in message and "alpha" in message and "do not guess" in message


def test_ambiguous_role_raises_rather_than_picking_one():
    with pytest.raises(SchemaError, match="ambiguous"):
        resolve_columns(SCHEMA, ["ID", "ra", "RAdeg", "dec"])


def test_override_wins_over_candidate_matching():
    mapping = resolve_columns(
        SCHEMA, ["ID", "ra", "RAdeg", "dec"], overrides={"ra": "ra"}
    )
    assert mapping["ra"].source == "ra"


def test_override_pointing_at_a_missing_column_raises():
    with pytest.raises(SchemaError, match="not in the table"):
        resolve_columns(SCHEMA, ["ID", "ra", "dec"], overrides={"ra": "nope"})


def test_optional_role_is_silently_absent():
    mapping = resolve_columns(SCHEMA, ["ID", "ra", "dec"])
    assert "vlos" not in mapping


def _valid_table() -> Table:
    table = Table()
    table["source_id"] = np.arange(5, dtype=np.int64)
    table["ra"] = np.linspace(200, 201, 5) * u.deg
    table["dec"] = np.linspace(-48, -47, 5) * u.deg
    return table


def test_validate_accepts_a_conforming_table():
    report = validate_table(_valid_table(), SCHEMA)
    assert report["n_rows"] == 5
    assert report["unique_id_ok"] is True
    assert report["units"]["ra"] == "deg"
    assert report["n_nonfinite"]["ra"] == 0


def test_validate_rejects_a_wrong_unit():
    table = _valid_table()
    table["ra"] = table["ra"].value * u.km
    with pytest.raises(SchemaError, match="unit is km"):
        validate_table(table, SCHEMA)


def test_validate_rejects_a_merely_equivalent_unit():
    """A standardized product must be IN the schema unit, not merely convertible.

    standardize() converts, so radians reaching validation means the table was
    not built by this pipeline. Accepting it would leave every downstream
    consumer responsible for unit handling.
    """
    table = _valid_table()
    table["ra"] = np.radians(table["ra"].value) * u.rad
    with pytest.raises(SchemaError, match="expected deg"):
        validate_table(table, SCHEMA)


def test_validate_rejects_a_missing_unit():
    table = _valid_table()
    table["ra"].unit = None
    with pytest.raises(SchemaError, match="column has none"):
        validate_table(table, SCHEMA)


def test_validate_rejects_duplicate_ids():
    table = _valid_table()
    table["source_id"][-1] = table["source_id"][0]
    with pytest.raises(SchemaError, match="duplicates"):
        validate_table(table, SCHEMA)


def test_validate_counts_nonfinite_values():
    """NaNs are counted where the role tolerates missing data."""
    schema = TableSchema(
        name="tolerant",
        roles=(ColumnRole("ra", "deg", ("ra",), allow_missing=True),),
    )
    table = Table()
    table["ra"] = np.array([200.0, np.nan, 202.0]) * u.deg
    assert validate_table(table, schema)["n_nonfinite"]["ra"] == 1


def test_validate_rejects_a_missing_required_role():
    table = _valid_table()
    table.remove_column("dec")
    with pytest.raises(SchemaError, match="missing"):
        validate_table(table, SCHEMA)


def test_describe_table_reports_dtype_and_unit():
    rows = describe_table(_valid_table())
    entry = {r["name"]: r for r in rows}["ra"]
    assert entry["unit"] == "deg" and entry["dtype"].startswith("float")


def test_load_column_map_reads_yaml(tmp_path):
    path = tmp_path / "maps.yaml"
    path.write_text("demo:\n  ra: RAdeg\n  dec: null\n", encoding="utf-8")
    assert load_column_map(path, "demo") == {"ra": {"column": "RAdeg", "unit": None}}
    assert load_column_map(path, "absent") == {}
    assert load_column_map(tmp_path / "nofile.yaml", "demo") == {}
