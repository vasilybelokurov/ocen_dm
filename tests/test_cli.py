"""CLI tests: the inventory command must run with no data and no network."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

import ocen_dm
from ocen_dm.cli import main

#: Real repository root, resolved from the package location so that it stays
#: correct after ``OCEN_DM_ROOT`` has been redirected to a temporary directory.
REPO_ROOT = Path(ocen_dm.__file__).resolve().parents[2]


@pytest.fixture
def populated_root(isolated_root):
    """A temporary project root carrying a copy of the real registry."""
    shutil.copy(
        REPO_ROOT / "provenance" / "datasets.yaml",
        isolated_root / "provenance" / "datasets.yaml",
    )
    (isolated_root / ".gitignore").write_text("data/raw/\n", encoding="utf-8")
    return isolated_root


def test_inventory_runs_with_no_data(populated_root, capsys):
    assert main(["inventory"]) == 0
    out = capsys.readouterr().out
    for key in ("kuzma2025_pristine", "kuzma2026_spectroscopy",
                "fimbulthul_ibata2019", "omegacat_vi_kinematics"):
        assert key in out
    assert "MISSING" in out
    assert "processed products: 0" in out


def test_inventory_marks_unverified_references(populated_root, capsys):
    main(["inventory"])
    assert "UNVERIFIED" in capsys.readouterr().out


def test_inventory_performs_no_inference(populated_root, capsys):
    """The inventory command must not create any results product."""
    main(["inventory"])
    assert not list((populated_root / "results").rglob("*")) if (
        populated_root / "results"
    ).is_dir() else True


def test_preprocess_reports_missing_inputs_without_crashing(populated_root, capsys):
    status = main(["preprocess"])
    out = capsys.readouterr().out
    assert status == 1                      # inputs absent -> non-zero, but graceful
    assert "FAILED" in out or "skipped" in out
    report = json.loads(
        (populated_root / "data" / "processed" / "preprocess_report.json").read_text()
    )
    assert isinstance(report, list) and report


def test_preprocess_rejects_an_unknown_dataset(populated_root):
    assert main(["preprocess", "--dataset", "not_a_dataset"]) == 2


def test_inspect_omegacat_without_data_returns_nonzero(populated_root, capsys):
    assert main(["inspect-omegacat"]) == 1
    assert "fetch-data" in capsys.readouterr().out
