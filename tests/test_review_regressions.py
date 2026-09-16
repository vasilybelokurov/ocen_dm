"""Regression tests for the twelve defects found in the 2026-09-16 review.

Each test is named for the finding it locks down. Every one of them fails
against the pre-review implementation.
"""

from __future__ import annotations

import hashlib
import io
import json
import urllib.error

import numpy as np
import pytest
from astropy import units as u
from astropy.table import MaskedColumn, Table

from ocen_dm.data import download as dl
from ocen_dm.data.base import standardize, write_processed
from ocen_dm.data.schema import ColumnRole, SchemaError, TableSchema, resolve_columns, validate_table
from ocen_dm.provenance import DatasetRecord, FileRecord, verify_file, write_manifest

PAYLOAD = b"synthetic catalogue bytes" * 100
PAYLOAD_MD5 = hashlib.md5(PAYLOAD).hexdigest()


class _FakeResponse(io.BytesIO):
    """Stand-in for an HTTP response with a controllable Content-Length."""

    def __init__(self, payload: bytes, claimed: int | None = None):
        super().__init__(payload)
        length = len(payload) if claimed is None else claimed
        self.headers = {"Content-Length": str(length)}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


# ---------------------------------------------------------------- finding 1 --
def test_f1_invalid_table_is_not_written(isolated_root, synthetic_periphery):
    """A violation caught by validate_table (not by standardize) must abort too.

    Duplicate identifiers pass column resolution and unit conversion and fail
    only in validation, so this exercises the fixed code path rather than an
    earlier one.
    """
    from ocen_dm.data import kuzma2025

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    bad = synthetic_periphery.copy()
    bad["source_id"][-1] = bad["source_id"][0]     # duplicate Gaia identifier
    bad.write(directory / "wCen_table.fits")

    with pytest.raises(SchemaError, match="duplicates"):
        kuzma2025.preprocess()
    assert not list((isolated_root / "data" / "processed").rglob("*.ecsv"))


def test_f1_failure_does_not_replace_a_good_product(isolated_root, synthetic_periphery):
    from ocen_dm.data import kuzma2025

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    synthetic_periphery.write(directory / "wCen_table.fits")
    good = kuzma2025.preprocess()
    first = Table.read(good["path"])

    bad = synthetic_periphery.copy()
    bad["dec"] = np.full(len(bad), 120.0) * u.deg   # impossible declination
    bad.write(directory / "wCen_table.fits", overwrite=True)
    with pytest.raises(SchemaError, match="physical bound"):
        kuzma2025.preprocess()

    reread = Table.read(good["path"])
    assert reread.colnames == first.colnames
    for name in first.colnames:
        assert np.array_equal(np.asarray(reread[name]), np.asarray(first[name]))


def test_f1_cli_preprocess_reports_failure(isolated_root, synthetic_periphery, capsys):
    from ocen_dm.cli import main

    directory = isolated_root / "data" / "raw" / "kuzma2025_pristine"
    directory.mkdir(parents=True)
    bad = synthetic_periphery.copy()
    bad["ra"] = bad["ra"].value * u.km
    bad.write(directory / "wCen_table.fits")
    assert main(["preprocess", "--dataset", "kuzma2025"]) == 1
    assert "FAILED" in capsys.readouterr().out


# ---------------------------------------------------------------- finding 2 --
def test_f2_discovery_keeps_publisher_checksum_and_size(isolated_root, monkeypatch):
    record = {
        "files": [
            {
                "key": "p.fits",
                "size": len(PAYLOAD),
                "checksum": f"md5:{PAYLOAD_MD5}",
                "links": {"self": "http://example/p.fits"},
            }
        ]
    }
    responses = [_FakeResponse(json.dumps(record).encode()), _FakeResponse(PAYLOAD)]
    monkeypatch.setattr(dl, "_open", lambda url, timeout: responses.pop(0))

    dataset = DatasetRecord(key="disc", title="t", role="r", kind="zenodo",
                            optional=False, data={"doi": "10.5281/zenodo.1"}, files=())
    result = dl.fetch_dataset(dataset, discover=True, log=lambda _: None)[0]
    assert result.record.published_md5 == PAYLOAD_MD5
    assert result.record.expected_bytes == len(PAYLOAD)
    assert result.record.checksum_source and "Zenodo API" in result.record.checksum_source
    assert result.verified is True


def test_f2_discovered_checksum_is_enforced(isolated_root, monkeypatch):
    """Serving different bytes than the discovered checksum must fail."""
    record = {
        "files": [
            {"key": "p.fits", "size": len(PAYLOAD), "checksum": "md5:" + "0" * 32,
             "links": {"self": "http://example/p.fits"}}
        ]
    }
    responses = [_FakeResponse(json.dumps(record).encode())] + [
        _FakeResponse(PAYLOAD) for _ in range(5)
    ]
    monkeypatch.setattr(dl, "_open", lambda url, timeout: responses.pop(0))
    dataset = DatasetRecord(key="disc", title="t", role="r", kind="zenodo",
                            optional=False, data={"doi": "10.5281/zenodo.1"}, files=())
    result = dl.fetch_dataset(dataset, discover=True, retries=1, log=lambda _: None)[0]
    assert result.status == "failed"
    assert not result.record.path.exists()


def test_f2b_empty_file_is_never_ok(tmp_path):
    path = tmp_path / "empty.bin"
    path.write_bytes(b"")
    result = verify_file(FileRecord(dataset="d", name="empty.bin", url=None), path)
    assert result["ok"] is False
    assert result["verification"] == "none"


# ---------------------------------------------------------------- finding 3 --
def test_f3_unknown_dataset_key_is_a_usage_error(isolated_root, registry_copy):
    from ocen_dm.cli import main

    assert main(["fetch-data", "--dataset", "kuzma2025_typo"]) == 2


def test_f3_failed_discovery_is_a_failed_result(isolated_root, monkeypatch):
    def boom(url, timeout):
        raise urllib.error.HTTPError(url, 504, "Gateway Time-out", {}, None)

    monkeypatch.setattr(dl, "_open", boom)
    dataset = DatasetRecord(key="empty", title="t", role="r", kind="zenodo",
                            optional=False, data={"doi": "10.5281/zenodo.14978551"}, files=())
    results = dl.fetch_dataset(dataset, discover=True, log=lambda _: None)
    assert [r.status for r in results] == ["failed"]
    assert results[0].record.required is True


def test_f3_missing_required_file_gives_nonzero_exit(isolated_root, registry_copy, monkeypatch):
    from ocen_dm.cli import main

    def boom(url, timeout):
        raise urllib.error.HTTPError(url, 504, "Gateway Time-out", {}, None)

    monkeypatch.setattr(dl, "_open", boom)
    assert main(["fetch-data", "--retries", "1", "--timeout", "1"]) == 1


def test_f3_unconfigured_profiles_are_not_success(isolated_root, registry_copy):
    from ocen_dm.cli import main

    assert main(["preprocess", "--dataset", "omegacat"]) == 1


# ---------------------------------------------------------------- finding 4 --
def test_f4_size_mismatch_is_a_failure(isolated_root, monkeypatch):
    monkeypatch.setattr(dl, "_open", lambda url, timeout: _FakeResponse(PAYLOAD))
    record = FileRecord(dataset="d", name="f.fits", url="http://example/f.fits",
                        expected_bytes=10 * len(PAYLOAD), expected_bytes_tolerance=0.0)
    dataset = DatasetRecord(key="d", title="t", role="r", kind="zenodo",
                            optional=False, data={}, files=(record,))
    result = dl.fetch_dataset(dataset, retries=1, log=lambda _: None)[0]
    assert result.status == "failed" and result.ok is False


def test_f4_manual_file_failing_its_checksum_is_not_cached(isolated_root):
    directory = isolated_root / "data" / "raw" / "manualset"
    directory.mkdir(parents=True)
    (directory / "table.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    record = FileRecord(dataset="manualset", name="table.*", url=None,
                        published_md5="0" * 32, checksum_source="test")
    dataset = DatasetRecord(key="manualset", title="t", role="r", kind="manual",
                            optional=False, data={}, files=(record,))
    result = dl.fetch_dataset(dataset, log=lambda _: None)[0]
    assert result.status == "failed"
    assert result.ok is False


def test_f4_ambiguous_manual_match_is_a_failure(isolated_root):
    directory = isolated_root / "data" / "raw" / "manualset"
    directory.mkdir(parents=True)
    (directory / "table.csv").write_text("a\n1\n", encoding="utf-8")
    (directory / "table.fits").write_bytes(b"x")
    record = FileRecord(dataset="manualset", name="table.*", url=None)
    dataset = DatasetRecord(key="manualset", title="t", role="r", kind="manual",
                            optional=False, data={}, files=(record,))
    result = dl.fetch_dataset(dataset, log=lambda _: None)[0]
    assert result.status == "failed" and "match" in result.message


def test_f4_manual_result_records_the_real_filename(isolated_root):
    directory = isolated_root / "data" / "raw" / "manualset"
    directory.mkdir(parents=True)
    (directory / "table.csv").write_text("a\n1\n", encoding="utf-8")
    record = FileRecord(dataset="manualset", name="table.*", url=None)
    dataset = DatasetRecord(key="manualset", title="t", role="r", kind="manual",
                            optional=False, data={}, files=(record,))
    result = dl.fetch_dataset(dataset, log=lambda _: None)[0]
    assert result.verification["name"] == "table.csv"   # not the wildcard


# ---------------------------------------------------------------- finding 5 --
SCHEMA_R = TableSchema("unit_probe", roles=(ColumnRole("r", "arcsec", ("r",)),))


def test_f5_unitless_column_is_rejected_not_assumed():
    table = Table()
    table["r"] = np.array([1.0, 2.0])          # arcmin in the paper, unitless in CSV
    with pytest.raises(SchemaError, match="carries no unit"):
        standardize(table, SCHEMA_R, overrides={"r": "r"})


def test_f5_declared_unit_is_converted_not_asserted():
    table = Table()
    table["r"] = np.array([1.0, 2.0])
    out = standardize(table, SCHEMA_R, overrides={"r": {"column": "r", "unit": "arcmin"}})
    assert out["r"].unit == u.arcsec
    assert np.allclose(out["r"].value, [60.0, 120.0])


def test_f5_source_unit_is_converted_to_the_schema_unit():
    table = Table()
    table["r"] = np.array([1.0, 2.0]) * u.arcmin
    out = standardize(table, SCHEMA_R, overrides={"r": "r"})
    assert out["r"].unit == u.arcsec and np.allclose(out["r"].value, [60.0, 120.0])
    assert out.meta["ocen_source_units"]["r"] == "arcmin"


def test_f5_velocity_and_angle_conversions():
    schema = TableSchema("v", roles=(ColumnRole("v", "km / s", ("v",)),
                                     ColumnRole("a", "deg", ("a",))))
    table = Table()
    table["v"] = np.array([1000.0]) * (u.m / u.s)
    table["a"] = np.array([np.pi]) * u.rad
    out = standardize(table, schema, overrides={"v": "v", "a": "a"})
    assert np.allclose(out["v"].value, [1.0]) and np.allclose(out["a"].value, [180.0])


def test_f5_inconvertible_unit_is_rejected():
    table = Table()
    table["r"] = np.array([1.0]) * u.km
    with pytest.raises(SchemaError, match="not convertible"):
        standardize(table, SCHEMA_R, overrides={"r": "r"})


# ---------------------------------------------------------------- finding 6 --
ID_SCHEMA = TableSchema(
    "id_probe", unique_id="source_id",
    roles=(ColumnRole("source_id", None, ("source_id",), kind="integer", is_identifier=True),),
)


def test_f6_float_identifier_is_rejected():
    table = Table()
    table["source_id"] = np.array([6.0e18, 6.1e18])
    with pytest.raises(SchemaError, match="floating point|not valid for a 'integer' role"):
        validate_table(table, ID_SCHEMA)


def test_f6_int64_identifiers_survive_the_round_trip(isolated_root):
    ids = np.array([6000000000000000001, 6000000000000000002], dtype=np.int64)
    table = Table()
    table["source_id"] = ids
    validate_table(table, ID_SCHEMA)
    path = write_processed(table, "id_probe")
    back = Table.read(path)
    assert np.array_equal(np.asarray(back["source_id"]), ids)


def test_f6_generic_id_column_cannot_pose_as_a_gaia_source_id():
    from ocen_dm.data import kuzma2026

    assert "ID" not in kuzma2026.SCHEMA.role("source_id").candidates
    assert "id" not in kuzma2026.SCHEMA.role("source_id").candidates


# ---------------------------------------------------------------- finding 7 --
def test_f7_product_records_raw_lineage(isolated_root, synthetic_periphery):
    from ocen_dm.data import kuzma2025
    from ocen_dm.provenance import file_digest

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    raw = directory / "wCen_table.fits"
    synthetic_periphery.write(raw)

    report = kuzma2025.preprocess()
    meta = Table.read(report["path"]).meta
    assert meta["ocen_raw"]["sha256"] == file_digest(raw, "sha256")
    assert meta["ocen_raw"]["bytes"] == raw.stat().st_size
    assert meta["ocen_package_version"]


def test_f7_altered_raw_file_is_detected(isolated_root, synthetic_periphery):
    from ocen_dm.data import kuzma2025
    from ocen_dm.provenance import file_digest

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    raw = directory / "wCen_table.fits"
    synthetic_periphery.write(raw)

    write_manifest([{  # pretend fetch-data verified this file
        "dataset": kuzma2025.DATASET, "name": "wCen_table.fits", "path": str(raw),
        "exists": True, "bytes": raw.stat().st_size, "sha256": file_digest(raw, "sha256"),
        "ok": True,
    }])
    edited = synthetic_periphery.copy()
    edited["ra"][0] = edited["ra"][0] + 0.001
    edited.write(raw, overwrite=True)

    with pytest.raises(SchemaError, match="changed since it was verified"):
        kuzma2025.preprocess()


# ---------------------------------------------------------------- finding 8 --
def test_f8_profile_schema_carries_bin_edges_and_asymmetric_errors():
    from ocen_dm.data import omegacat

    roles = omegacat.PROFILE_SCHEMAS["pm_radial"].role_names
    for expected in ("r_lower", "r_median", "r_upper",
                     "sigma_pmr", "sigma_pmr_err_lo", "sigma_pmr_err_hi"):
        assert expected in roles


def test_f8_unclaimed_columns_are_preserved(isolated_root, synthetic_periphery):
    from ocen_dm.data import kuzma2025

    table = synthetic_periphery.copy()
    table["pristine_ca_hk"] = np.linspace(0, 1, len(table))
    out = standardize(table, kuzma2025.SCHEMA, keep_extra_columns=True)
    assert "src_pristine_ca_hk" in out.colnames


def test_f8_omitted_columns_are_reported_when_not_kept(isolated_root, synthetic_periphery):
    from ocen_dm.data import kuzma2025

    table = synthetic_periphery.copy()
    table["extra_thing"] = np.zeros(len(table))
    out = standardize(table, kuzma2025.SCHEMA, keep_extra_columns=False)
    assert "extra_thing" in out.meta["ocen_omitted_columns"]


# ---------------------------------------------------------------- finding 9 --
def test_f9_truncated_transfer_is_rejected(isolated_root, monkeypatch):
    monkeypatch.setattr(dl, "_open", lambda url, t: _FakeResponse(b"short", claimed=100000))
    dest = isolated_root / "trunc.fits"
    with pytest.raises(dl.DownloadError):
        dl.download_file("http://example/x", dest, retries=1, backoff=0.0)
    assert not dest.exists()


def test_f9_verification_happens_before_publication(isolated_root, monkeypatch):
    monkeypatch.setattr(dl, "_open", lambda url, t: _FakeResponse(PAYLOAD))
    dest = isolated_root / "verified.fits"
    with pytest.raises(dl.DownloadError, match="verification failed"):
        dl.download_file("http://example/x", dest, retries=1, backoff=0.0,
                         verify=lambda p: (False, "bad checksum"))
    assert not dest.exists()


def test_f9_corrupt_refetch_does_not_destroy_a_good_copy(isolated_root, monkeypatch):
    record = FileRecord(dataset="d", name="f.fits", url="http://example/f.fits",
                        published_md5=PAYLOAD_MD5, checksum_source="test")
    dataset = DatasetRecord(key="d", title="t", role="r", kind="zenodo",
                            optional=False, data={}, files=(record,))
    monkeypatch.setattr(dl, "_open", lambda url, t: _FakeResponse(PAYLOAD))
    assert dl.fetch_dataset(dataset, log=lambda _: None)[0].ok

    monkeypatch.setattr(dl, "_open", lambda url, t: _FakeResponse(b"corrupted bytes"))
    result = dl.fetch_dataset(dataset, force=True, retries=1, log=lambda _: None)[0]
    assert result.status == "failed"
    assert record.path.read_bytes() == PAYLOAD          # the good copy survives


def test_f9_temporary_files_are_unique(isolated_root, monkeypatch):
    """Two concurrent transfers must not share one temporary path."""
    seen: list[str] = []

    def opener(url, timeout):
        # names of the temp files that exist while a transfer is in flight
        seen.extend(str(p.name) for p in isolated_root.glob(".*.part"))
        return _FakeResponse(PAYLOAD)

    monkeypatch.setattr(dl, "_open", opener)
    dl.download_file("http://example/a", isolated_root / "a.fits", retries=1)
    dl.download_file("http://example/a", isolated_root / "a.fits", retries=1)

    assert len(seen) == 2, f"expected one in-flight temp file per transfer, saw {seen}"
    assert seen[0] != seen[1], f"temporary name was reused: {seen}"
    assert not list(isolated_root.glob(".*part"))       # nothing left behind


def test_f9_no_partial_file_survives_failure(isolated_root, monkeypatch):
    def boom(url, timeout):
        raise urllib.error.URLError("nope")

    monkeypatch.setattr(dl, "_open", boom)
    dest = isolated_root / "f.fits"
    with pytest.raises(dl.DownloadError):
        dl.download_file("http://example/x", dest, retries=2, backoff=0.0)
    assert not dest.exists()
    assert not list(isolated_root.glob(".*"))


# --------------------------------------------------------------- finding 10 --
MASK_SCHEMA = TableSchema("mask_probe", roles=(ColumnRole("ra", "deg", ("ra",)),))


def test_f10_masked_required_value_is_rejected():
    table = Table()
    table["ra"] = MaskedColumn([200.0, 201.0], mask=[False, True], unit=u.deg)
    with pytest.raises(SchemaError, match="missing value"):
        validate_table(table, MASK_SCHEMA)


def test_f10_masked_optional_value_is_allowed_and_counted():
    schema = TableSchema("m2", roles=(ColumnRole("v", "km / s", ("v",),
                                                 required=False, allow_missing=True),))
    table = Table()
    table["v"] = MaskedColumn([1.0, 2.0], mask=[False, True], unit=u.km / u.s)
    assert validate_table(table, schema)["n_masked"]["v"] == 1


def test_f10_nan_in_a_required_role_is_rejected():
    table = Table()
    table["ra"] = np.array([200.0, np.nan]) * u.deg
    with pytest.raises(SchemaError, match="missing value"):
        validate_table(table, MASK_SCHEMA)


def test_f10_empty_table_is_rejected():
    table = Table()
    table["ra"] = np.array([], dtype=float) * u.deg
    with pytest.raises(SchemaError, match="no rows"):
        validate_table(table, MASK_SCHEMA)


def test_f10_out_of_range_values_are_rejected():
    schema = TableSchema("r", roles=(ColumnRole("p", None, ("p",), valid_range=(0.0, 1.0)),))
    table = Table()
    table["p"] = np.array([0.5, 1.5])
    with pytest.raises(SchemaError, match="physical bound"):
        validate_table(table, schema)


# --------------------------------------------------------------- finding 11 --
TWO_ROLES = TableSchema("two", roles=(ColumnRole("ra", "deg", ("ra",)),
                                      ColumnRole("dec", "deg", ("dec",))))


def test_f11_unknown_override_key_is_rejected():
    with pytest.raises(SchemaError, match="unknown role"):
        resolve_columns(TWO_ROLES, ["ra", "dec"], {"raa": "ra"})


def test_f11_two_roles_cannot_share_one_column():
    with pytest.raises(SchemaError, match="several roles"):
        resolve_columns(TWO_ROLES, ["ra", "dec"], {"dec": "ra"})


def test_f11_column_map_path_follows_the_project_root(isolated_root):
    from ocen_dm.data.base import column_map_path

    assert column_map_path() == isolated_root / "configs" / "column_maps.yaml"


def test_f11_column_map_is_read_from_the_active_root(isolated_root):
    from ocen_dm.data.base import standardize as std

    (isolated_root / "configs" / "column_maps.yaml").write_text(
        "unit_probe:\n  r: {column: R_arcmin, unit: arcmin}\n", encoding="utf-8"
    )
    table = Table()
    table["R_arcmin"] = np.array([1.0])
    out = std(table, SCHEMA_R)
    assert np.allclose(out["r"].value, [60.0])


# --------------------------------------------------------------- finding 12 --
def test_f12_manifest_merges_across_datasets(isolated_root):
    entry_a = {"dataset": "A", "name": "a.fits", "sha256": "a" * 64, "bytes": 1, "ok": True}
    entry_b = {"dataset": "B", "name": "b.fits", "sha256": "b" * 64, "bytes": 2, "ok": True}
    write_manifest([entry_a])
    manifest = write_manifest([entry_b])
    datasets = {f["dataset"] for f in json.loads(manifest.read_text())["files"]}
    assert datasets == {"A", "B"}


def test_f12_manifest_records_failures(isolated_root, monkeypatch, registry_copy):
    def boom(url, timeout):
        raise urllib.error.HTTPError(url, 504, "Gateway Time-out", {}, None)

    monkeypatch.setattr(dl, "_open", boom)
    dl.fetch_all(only=["kuzma2025_pristine"], retries=1, log=lambda _: None)
    manifest = json.loads(
        (isolated_root / "provenance" / "manifest.json").read_text()
    )
    entry = manifest["files"][0]
    assert entry["status"] == "failed" and entry["ok"] is False
    assert entry["url"]


def test_f12_processed_write_is_atomic(isolated_root, synthetic_periphery, monkeypatch):
    """A write that fails midway must leave the previous product byte-identical."""
    path = write_processed(synthetic_periphery, "atomic_probe")
    original = path.read_bytes()

    bigger = synthetic_periphery.copy()
    bigger["ra"] = bigger["ra"] * 0 + 12.3 * u.deg

    def explode(*args, **kwargs):
        raise OSError("disk full, halfway through")

    monkeypatch.setattr(type(bigger), "write", explode)
    with pytest.raises(OSError):
        write_processed(bigger, "atomic_probe")

    assert path.read_bytes() == original
    assert not list(path.parent.glob(".*tmp*"))


# ----------------------------------------------- defects found in self-audit --
def test_selfaudit_unit_conversion_preserves_the_mask():
    """Column.to() drops the mask, turning a missing value into a real number."""
    schema = TableSchema("mask_convert",
                         roles=(ColumnRole("v", "km / s", ("v",),
                                           required=False, allow_missing=True),))
    table = Table()
    table["v"] = MaskedColumn([1000.0, 2000.0], mask=[False, True], unit=u.m / u.s)
    out = standardize(table, schema, overrides={"v": "v"})
    assert np.array_equal(np.asarray(out["v"].mask), [False, True])
    assert out["v"].unit == u.km / u.s
    assert validate_table(out, schema)["n_masked"]["v"] == 1


def test_selfaudit_masked_and_nan_are_counted_as_a_union():
    """One masked row plus one NaN row is two missing values, not one."""
    schema = TableSchema("count", roles=(ColumnRole("ra", "deg", ("ra",)),))
    table = Table()
    table["ra"] = MaskedColumn([200.0, np.nan, 202.0], mask=[True, False, False], unit=u.deg)
    with pytest.raises(SchemaError, match="2 missing value"):
        validate_table(table, schema)


def test_selfaudit_range_check_ignores_masked_rows():
    schema = TableSchema("rng", roles=(ColumnRole("p", None, ("p",), valid_range=(0.0, 1.0),
                                                  required=False, allow_missing=True),))
    table = Table()
    table["p"] = MaskedColumn([0.5, 99.0], mask=[False, True])
    assert validate_table(table, schema)["n_masked"]["p"] == 1   # masked 99 is not a violation


# ------------------------------------------- second-round review regressions --
def test_r2_failed_checksum_is_not_trusted_at_consumption(isolated_root, synthetic_periphery):
    """A file recorded as checksum-FAILED must not become valid input later."""
    from ocen_dm.data import kuzma2025
    from ocen_dm.data.base import check_against_manifest
    from ocen_dm.provenance import file_digest

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    raw = directory / "wCen_table.fits"
    synthetic_periphery.write(raw)
    sha = file_digest(raw, "sha256")

    write_manifest([{
        "dataset": kuzma2025.DATASET, "name": "wCen_table.fits", "path": str(raw),
        "exists": True, "bytes": raw.stat().st_size, "sha256": sha,
        "md5": "0" * 32, "published_checksum": "f" * 32,
        "checksum_ok": False, "ok": False, "status": "failed",
    }])
    with pytest.raises(SchemaError, match="FAILED checksum"):
        check_against_manifest(kuzma2025.DATASET, raw, sha)
    with pytest.raises(SchemaError):
        kuzma2025.preprocess()


def test_r2_failed_refetch_keeps_the_verified_anchor(isolated_root):
    verified = {
        "dataset": "d", "name": "f.fits", "path": "/x/f.fits", "exists": True,
        "bytes": 10, "sha256": "a" * 64, "md5": "b" * 32, "ok": True, "status": "ok",
    }
    failed = {
        "dataset": "d", "name": "f.fits", "path": "/x/f.fits", "exists": False,
        "bytes": None, "sha256": None, "ok": False, "status": "failed",
    }
    write_manifest([verified])
    manifest = write_manifest([failed])
    entry = json.loads(manifest.read_text())["files"][0]
    assert entry["status"] == "failed"
    assert entry["last_verified"]["sha256"] == "a" * 64   # anchor survives


def test_r2_registry_honours_a_declared_checksum_algorithm(tmp_path):
    from ocen_dm.provenance import load_registry

    path = tmp_path / "reg.yaml"
    path.write_text(
        "schema_version: 1\n"
        "datasets:\n"
        "  d:\n"
        "    title: t\n"
        "    role: r\n"
        "    data: {kind: zenodo}\n"
        "    files:\n"
        "      - name: f.fits\n"
        "        url: http://e/f\n"
        "        published_md5: abc\n"
        "        published_algorithm: sha256\n"
        "        checksum_source: test\n",
        encoding="utf-8",
    )
    assert load_registry(path)["d"].files[0].published_algorithm == "sha256"


def test_r2_percentage_is_converted_not_taken_at_face_value():
    """0.5 percent is 0.005, not 0.5, in a probability role."""
    schema = TableSchema("p", roles=(ColumnRole("prob", None, ("p",),
                                                valid_range=(0.0, 1.0), dimensionless=True),))
    table = Table()
    table["p"] = np.array([0.5, 20.0]) * u.percent
    out = standardize(table, schema, overrides={"prob": "p"})
    assert np.allclose(np.asarray(out["prob"]), [0.005, 0.2])
    validate_table(out, schema)


def test_r2_dimensional_value_in_a_dimensionless_role_is_rejected():
    schema = TableSchema("p", roles=(ColumnRole("prob", None, ("p",), dimensionless=True),))
    table = Table()
    table["p"] = np.array([0.5]) * u.km
    with pytest.raises(SchemaError, match="pure number"):
        standardize(table, schema, overrides={"prob": "p"})


def test_r2_object_dtype_identifier_of_floats_is_rejected():
    schema = TableSchema("o", unique_id="name",
                         roles=(ColumnRole("name", None, ("name",),
                                           kind="string", is_identifier=True),))
    table = Table()
    table["name"] = np.array([1.5, 2.5], dtype=object)
    with pytest.raises(SchemaError, match="non-identifier values"):
        validate_table(table, schema)


def test_r2_blank_identifier_is_rejected():
    schema = TableSchema("o", unique_id="name",
                         roles=(ColumnRole("name", None, ("name",),
                                           kind="string", is_identifier=True),))
    table = Table()
    table["name"] = np.array(["Gaia DR3 1", "  "])
    with pytest.raises(SchemaError, match="blank identifier"):
        validate_table(table, schema)


def test_r2_inconsistent_bin_geometry_is_rejected():
    from ocen_dm.data import omegacat

    table = Table()
    table["r_lower"] = [100.0] * u.arcsec
    table["r_median"] = [10.0] * u.arcsec
    table["r_upper"] = [1.0] * u.arcsec
    table["sigma_pmr"] = [0.8] * u.mas / u.yr
    table["sigma_pmr_err_lo"] = [0.1] * u.mas / u.yr
    table["sigma_pmr_err_hi"] = [0.1] * u.mas / u.yr
    with pytest.raises(SchemaError, match="bin geometry"):
        validate_table(table, omegacat.PROFILE_SCHEMAS["pm_radial"])


def test_r2_negative_dispersion_is_rejected():
    from ocen_dm.data import omegacat

    table = Table()
    table["r_lower"] = [1.0] * u.arcsec
    table["r_median"] = [2.0] * u.arcsec
    table["r_upper"] = [3.0] * u.arcsec
    table["sigma_pmr"] = [-0.8] * u.mas / u.yr
    table["sigma_pmr_err_lo"] = [0.1] * u.mas / u.yr
    table["sigma_pmr_err_hi"] = [0.1] * u.mas / u.yr
    with pytest.raises(SchemaError, match="physical bound"):
        validate_table(table, omegacat.PROFILE_SCHEMAS["pm_radial"])


def test_r2_overlapping_bins_are_rejected():
    from ocen_dm.data import omegacat

    table = Table()
    table["r_lower"] = [1.0, 2.0] * u.arcsec
    table["r_median"] = [2.0, 3.0] * u.arcsec
    table["r_upper"] = [3.0, 4.0] * u.arcsec     # bin 1 ends at 3, bin 2 starts at 2
    table["sigma_pmr"] = [0.8, 0.7] * u.mas / u.yr
    table["sigma_pmr_err_lo"] = [0.1, 0.1] * u.mas / u.yr
    table["sigma_pmr_err_hi"] = [0.1, 0.1] * u.mas / u.yr
    with pytest.raises(SchemaError, match="overlap"):
        validate_table(table, omegacat.PROFILE_SCHEMAS["pm_radial"])


def test_r2_load_validates_by_default(isolated_root, synthetic_periphery):
    """load() must not hand back an invalid table to a downstream caller."""
    from ocen_dm.data import kuzma2025

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    bad = synthetic_periphery.copy()
    bad["dec"] = np.full(len(bad), 120.0) * u.deg
    bad.write(directory / "wCen_table.fits")

    with pytest.raises(SchemaError, match="physical bound"):
        kuzma2025.load()
    assert kuzma2025.load(validate=False) is not None    # inspection escape hatch


def test_r2_dataset_with_no_files_and_no_discovery_is_incomplete(isolated_root, registry_copy):
    from ocen_dm.cli import main

    assert main(["fetch-data", "--dataset", "omegacat_vi_kinematics", "--no-discover"]) == 1


def test_r2_empty_repository_record_is_a_failure(isolated_root, monkeypatch):
    monkeypatch.setattr(dl, "_open", lambda url, t: _FakeResponse(b'{"files": []}'))
    dataset = DatasetRecord(key="empty", title="t", role="r", kind="zenodo",
                            optional=False, data={"doi": "10.5281/zenodo.1"}, files=())
    results = dl.fetch_dataset(dataset, discover=True, log=lambda _: None)
    assert [r.status for r in results] == ["failed"]


def test_r2_download_does_not_close_a_descriptor_it_handed_over(isolated_root, monkeypatch):
    """The retry path must not close a descriptor os.fdopen() already owns."""
    import os

    def opener(url, timeout):
        return _FakeResponse(PAYLOAD)

    monkeypatch.setattr(dl, "_open", opener)
    probe = os.open(os.devnull, os.O_RDONLY)       # a descriptor we own
    try:
        with pytest.raises(dl.DownloadError):
            dl.download_file("http://example/x", isolated_root / "f.fits",
                             retries=1, backoff=0.0,
                             verify=lambda p: (False, "nope"))
        os.fstat(probe)                            # still open => not stolen
    finally:
        os.close(probe)
