"""Downloader tests. The network is mocked; no test contacts a real server."""

from __future__ import annotations

import hashlib
import io
import json
import urllib.error

import pytest

from ocen_dm.data import download as dl
from ocen_dm.provenance import DatasetRecord, FileRecord

PAYLOAD = b"synthetic catalogue bytes" * 100
PAYLOAD_MD5 = hashlib.md5(PAYLOAD).hexdigest()


class _FakeResponse(io.BytesIO):
    """Minimal stand-in for an ``http.client.HTTPResponse``."""

    def __init__(self, payload: bytes):
        super().__init__(payload)
        self.headers = {"Content-Length": str(len(payload))}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _serve(payload: bytes, fail_times: int = 0):
    """Return an ``_open`` replacement that fails ``fail_times`` times first."""
    state = {"calls": 0}

    def opener(url, timeout):
        state["calls"] += 1
        if state["calls"] <= fail_times:
            raise urllib.error.URLError("simulated transient failure")
        return _FakeResponse(payload)

    opener.state = state
    return opener


def test_download_file_is_atomic_and_complete(tmp_path, monkeypatch):
    monkeypatch.setattr(dl, "_open", _serve(PAYLOAD))
    dest = tmp_path / "sub" / "f.fits"
    dl.download_file("http://example/f.fits", dest)
    assert dest.read_bytes() == PAYLOAD
    assert not (tmp_path / "sub" / "f.fits.part").exists()


def test_download_file_retries_transient_failures(tmp_path, monkeypatch):
    opener = _serve(PAYLOAD, fail_times=2)
    monkeypatch.setattr(dl, "_open", opener)
    dest = tmp_path / "f.fits"
    dl.download_file("http://example/f.fits", dest, retries=3, backoff=0.0)
    assert dest.read_bytes() == PAYLOAD
    assert opener.state["calls"] == 3


def test_download_file_gives_up_and_leaves_no_partial(tmp_path, monkeypatch):
    monkeypatch.setattr(dl, "_open", _serve(PAYLOAD, fail_times=99))
    dest = tmp_path / "f.fits"
    with pytest.raises(dl.DownloadError, match="after 3 attempts"):
        dl.download_file("http://example/f.fits", dest, retries=3, backoff=0.0)
    assert not dest.exists()
    assert not dest.with_suffix(".fits.part").exists()


def _dataset(md5: str | None, *, expected_bytes: int | None = None) -> DatasetRecord:
    return DatasetRecord(
        key="demo",
        title="demo",
        role="test",
        kind="zenodo",
        optional=False,
        data={"doi": "10.5281/zenodo.1"},
        files=(
            FileRecord(
                dataset="demo",
                name="f.fits",
                url="http://example/f.fits",
                published_md5=md5,
                checksum_source="test" if md5 else None,
                expected_bytes=expected_bytes,
            ),
        ),
    )


def test_fetch_dataset_verifies_the_published_checksum(isolated_root, monkeypatch):
    monkeypatch.setattr(dl, "_open", _serve(PAYLOAD))
    results = dl.fetch_dataset(_dataset(PAYLOAD_MD5), log=lambda _: None)
    assert [r.status for r in results] == ["ok"]
    assert results[0].verification["checksum_ok"] is True
    assert results[0].verified is True


def test_fetch_dataset_flags_a_checksum_mismatch_as_failed(isolated_root, monkeypatch):
    monkeypatch.setattr(dl, "_open", _serve(PAYLOAD))
    results = dl.fetch_dataset(_dataset("0" * 32), retries=1, log=lambda _: None)
    assert results[0].status == "failed"
    assert "verification failed" in results[0].message
    # nothing is published when the checksum does not match
    assert not results[0].record.path.exists()


def test_fetch_dataset_reuses_a_verified_local_copy(isolated_root, monkeypatch):
    opener = _serve(PAYLOAD)
    monkeypatch.setattr(dl, "_open", opener)
    dataset = _dataset(PAYLOAD_MD5)

    dl.fetch_dataset(dataset, log=lambda _: None)
    calls_after_first = opener.state["calls"]
    results = dl.fetch_dataset(dataset, log=lambda _: None)

    assert results[0].status == "cached"
    assert opener.state["calls"] == calls_after_first  # no second request


def test_force_refetches_even_a_verified_copy(isolated_root, monkeypatch):
    opener = _serve(PAYLOAD)
    monkeypatch.setattr(dl, "_open", opener)
    dataset = _dataset(PAYLOAD_MD5)
    dl.fetch_dataset(dataset, log=lambda _: None)
    dl.fetch_dataset(dataset, force=True, log=lambda _: None)
    assert opener.state["calls"] == 2


def test_manual_dataset_is_reported_not_faked(isolated_root):
    dataset = DatasetRecord(
        key="manualset", title="t", role="r", kind="manual", optional=False,
        data={"kind": "manual"},
        files=(FileRecord(dataset="manualset", name="table.*", url=None),),
    )
    results = dl.fetch_dataset(dataset, log=lambda _: None)
    assert results[0].status == "manual"
    assert "MANUAL_DOWNLOADS" in results[0].message
    assert not list((isolated_root / "data" / "raw").rglob("table.*"))


def test_zenodo_discovery_lists_files(monkeypatch):
    record = {
        "files": [
            {
                "key": "profiles.fits",
                "size": 1234,
                "checksum": "md5:" + "a" * 32,
                "links": {"self": "http://example/profiles.fits"},
            }
        ]
    }

    def opener(url, timeout):
        return _FakeResponse(json.dumps(record).encode())

    monkeypatch.setattr(dl, "_open", opener)
    listing = dl.zenodo_files(14978551)
    assert listing == [
        {
            "name": "profiles.fits",
            "url": "http://example/profiles.fits",
            "size": 1234,
            "checksum": "md5:" + "a" * 32,
        }
    ]


def test_zenodo_outage_is_surfaced_not_swallowed(monkeypatch):
    def opener(url, timeout):
        raise urllib.error.HTTPError(url, 504, "Gateway Time-out", {}, None)

    monkeypatch.setattr(dl, "_open", opener)
    with pytest.raises(RuntimeError, match="Zenodo API unreachable"):
        dl.zenodo_record(14791430)


def test_discovery_failure_does_not_abort_the_dataset(isolated_root, monkeypatch):
    def opener(url, timeout):
        raise urllib.error.HTTPError(url, 504, "Gateway Time-out", {}, None)

    monkeypatch.setattr(dl, "_open", opener)
    dataset = DatasetRecord(
        key="empty", title="t", role="r", kind="zenodo", optional=False,
        data={"doi": "10.5281/zenodo.14978551"}, files=(),
    )
    messages: list[str] = []
    results = dl.fetch_dataset(dataset, discover=True, log=messages.append)
    # An unreachable repository is a failure, never an empty dataset.
    assert [r.status for r in results] == ["failed"]
    assert any("discovery failed" in m for m in messages)
