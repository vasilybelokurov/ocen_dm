"""Tests for the dataset registry, checksums, and the manifest."""

from __future__ import annotations

import hashlib
import json

import pytest

from ocen_dm.paths import project_root, provenance_dir
from ocen_dm.provenance import (
    FileRecord,
    file_digest,
    load_registry,
    verify_file,
    write_manifest,
)

PRIMARY = [
    "kuzma2025_pristine",
    "kuzma2026_spectroscopy",
    "fimbulthul_ibata2019",
    "omegacat_vi_kinematics",
]


def test_registry_lists_the_four_primary_datasets():
    registry = load_registry()
    for key in PRIMARY:
        assert key in registry.datasets
        assert not registry[key].optional


def test_every_dataset_has_a_resolvable_reference():
    registry = load_registry()
    for dataset in registry:
        paper = dataset.paper
        assert paper.get("arxiv") or paper.get("doi") or paper.get("journal"), (
            f"{dataset.key} has no citable reference"
        )
        assert "verified" in paper, f"{dataset.key} does not state its verification status"


def test_checksums_declare_their_source():
    """No checksum may appear without saying where it came from."""
    registry = load_registry()
    for dataset in registry:
        for record in dataset.files:
            if record.published_md5 is not None:
                assert record.checksum_source, f"{dataset.key}/{record.name}"


def test_kuzma2025_published_md5_matches_the_specification():
    registry = load_registry()
    record = registry["kuzma2025_pristine"].files[0]
    assert record.name == "wCen_table.fits"
    assert record.published_md5 == "3dcd58e69767901fbbdf69385a14c768"


def test_registry_rejects_a_checksum_without_provenance(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        "schema_version: 1\n"
        "datasets:\n"
        "  x:\n"
        "    title: t\n"
        "    role: r\n"
        "    data: {kind: zenodo}\n"
        "    files:\n"
        "      - {name: f.fits, url: 'http://e/f', published_md5: 'deadbeef'}\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="checksum_source"):
        load_registry(bad)


def test_file_digest_matches_hashlib(tmp_path):
    path = tmp_path / "blob.bin"
    payload = b"omega centauri" * 1000
    path.write_bytes(payload)
    assert file_digest(path, "md5") == hashlib.md5(payload).hexdigest()
    assert file_digest(path, "sha256") == hashlib.sha256(payload).hexdigest()


def test_verify_file_detects_a_corrupted_download(tmp_path):
    payload = b"correct"
    path = tmp_path / "f.bin"
    path.write_bytes(payload)
    good = FileRecord(
        dataset="d", name="f.bin", url=None,
        published_md5=hashlib.md5(payload).hexdigest(), checksum_source="test",
    )
    assert verify_file(good, path)["ok"] is True

    path.write_bytes(b"corrupted")
    result = verify_file(good, path)
    assert result["checksum_ok"] is False
    assert result["ok"] is False


def test_verify_file_on_missing_file_is_not_ok(tmp_path):
    record = FileRecord(dataset="d", name="absent.bin", url=None)
    result = verify_file(record, tmp_path / "absent.bin")
    assert result["exists"] is False
    assert result["ok"] is False


def test_write_manifest_records_environment_and_checksums(isolated_root):
    results = [
        {"dataset": "d", "name": "f.bin", "path": "/x/f.bin", "exists": True,
         "bytes": 10, "md5": "a" * 32, "sha256": "b" * 64, "md5_ok": True,
         "size_ok": True, "ok": True}
    ]
    manifest = write_manifest(results)
    payload = json.loads(manifest.read_text())
    assert payload["environment"]["python"]
    assert payload["files"][0]["sha256"] == "b" * 64
    checksums = (manifest.parent / "checksums.txt").read_text()
    assert "b" * 64 in checksums and "d/f.bin" in checksums


def test_raw_data_is_not_committed():
    """data/raw must be gitignored (spec Milestone 1 acceptance criterion)."""
    ignore = (project_root() / ".gitignore").read_text()
    assert "data/raw/" in ignore


def test_provenance_directory_exists():
    assert (provenance_dir() / "datasets.yaml").is_file()
