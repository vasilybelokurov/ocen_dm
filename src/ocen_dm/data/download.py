"""Deterministic retrieval of external datasets.

Integrity rules, tightened after the 2026-09-16 review:

* a transfer shorter than the advertised ``Content-Length`` is a failure;
* the temporary file has a unique name, so concurrent runs cannot corrupt
  each other;
* the checksum is verified **on the temporary file**, and only a file that
  passes is moved into its final path, so a valid cached copy is never
  destroyed by a corrupt refetch;
* publisher checksums and sizes discovered from the Zenodo API are carried
  into the file record and enforced, not discarded;
* a status is never more optimistic than its verification.

Files that cannot be fetched programmatically (journal supplementary material
behind publisher sessions) are reported as manual actions, never faked.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from ..provenance import (
    DatasetRecord,
    FileRecord,
    Registry,
    load_registry,
    verify_file,
    write_manifest,
)

__all__ = [
    "DownloadResult",
    "DownloadError",
    "download_file",
    "fetch_dataset",
    "fetch_all",
    "zenodo_record",
    "zenodo_files",
]

USER_AGENT = "ocen-dm/0.1 (research pipeline; https://arxiv.org/abs/2502.01135)"
DEFAULT_TIMEOUT = 60.0
DEFAULT_RETRIES = 3
RETRY_BACKOFF = 5.0

#: statuses that mean "the file is present locally and passed its checks"
_OK_STATUSES = {"ok", "cached"}


class DownloadError(RuntimeError):
    """Raised when a file cannot be retrieved intact."""


@dataclass
class DownloadResult:
    """Outcome of attempting to obtain one file."""

    record: FileRecord
    status: str  # 'ok' | 'cached' | 'failed' | 'manual' | 'skipped'
    message: str = ""
    verification: dict[str, Any] | None = None

    @property
    def ok(self) -> bool:
        """True only when the file is present *and* its verification passed."""
        if self.status not in _OK_STATUSES:
            return False
        return bool(self.verification and self.verification.get("ok"))

    @property
    def verified(self) -> bool:
        """True when a publisher checksum was matched."""
        return bool(self.verification and self.verification.get("verified"))


def _open(url: str, timeout: float):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=timeout)


def zenodo_record(record_id: str | int, timeout: float = DEFAULT_TIMEOUT) -> dict[str, Any]:
    """Fetch a Zenodo record's JSON metadata.

    Parameters
    ----------
    record_id : str or int
        Numeric Zenodo record identifier, e.g. ``14791430``.
    timeout : float, optional
        Socket timeout in seconds.

    Returns
    -------
    dict
        The parsed Zenodo API record.

    Raises
    ------
    RuntimeError
        If Zenodo is unreachable or returns a non-JSON body. As of 2026-09-16
        the Zenodo API returns HTTP 504; this is surfaced rather than swallowed
        so that no dataset is silently skipped.
    """
    url = f"https://zenodo.org/api/records/{record_id}"
    try:
        with _open(url, timeout) as response:
            return json.load(response)
    except (urllib.error.URLError, json.JSONDecodeError, TimeoutError) as exc:
        raise RuntimeError(f"Zenodo API unreachable for record {record_id}: {exc}") from exc


def zenodo_files(record_id: str | int, timeout: float = DEFAULT_TIMEOUT) -> list[dict[str, Any]]:
    """List a Zenodo record's files as ``{name, url, size, checksum}`` dicts.

    ``checksum`` keeps Zenodo's ``'<algorithm>:<hexdigest>'`` form.
    """
    record = zenodo_record(record_id, timeout=timeout)
    out = []
    for item in record.get("files", []):
        out.append(
            {
                "name": item.get("key"),
                "url": (item.get("links") or {}).get("self"),
                "size": item.get("size"),
                "checksum": item.get("checksum"),
            }
        )
    return out


def _split_checksum(value: str | None) -> tuple[str | None, str | None]:
    """Split Zenodo's ``'md5:abc…'`` into ``('md5', 'abc…')``."""
    if not value:
        return None, None
    if ":" in value:
        algorithm, digest = value.split(":", 1)
        return algorithm.strip().lower(), digest.strip().lower()
    return "md5", value.strip().lower()


def download_file(
    url: str,
    dest: Path,
    *,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    backoff: float = RETRY_BACKOFF,
    verify: Callable[[Path], tuple[bool, str]] | None = None,
    progress: Callable[[int, int | None], None] | None = None,
) -> Path:
    """Download ``url`` to ``dest``, verifying before publishing it.

    Parameters
    ----------
    url : str
        Source URL.
    dest : pathlib.Path
        Final destination path; parent directories are created.
    timeout : float, optional
        Socket timeout in seconds.
    retries : int, optional
        Number of attempts before giving up.
    backoff : float, optional
        Seconds to wait after the first failure; doubled each further attempt.
    verify : callable, optional
        ``verify(temp_path) -> (ok, message)``, run on the completed temporary
        file. A failure is retried, and ``dest`` is left untouched.
    progress : callable, optional
        Called as ``progress(bytes_so_far, total_or_None)``.

    Returns
    -------
    pathlib.Path
        ``dest``.

    Raises
    ------
    DownloadError
        If every attempt fails, the transfer is short, or verification fails.
        Any existing file at ``dest`` is left in place.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        handle, tmp_name = tempfile.mkstemp(
            dir=dest.parent, prefix=f".{dest.name}.", suffix=".part"
        )
        tmp = Path(tmp_name)
        handle_owned = True          # ours until os.fdopen() takes it over
        try:
            response = _open(url, timeout)
            fh = os.fdopen(handle, "wb")
            handle_owned = False     # the file object owns the descriptor now
            with response, fh:
                total_header = response.headers.get("Content-Length")
                total = int(total_header) if total_header else None
                seen = 0
                while chunk := response.read(1 << 20):
                    fh.write(chunk)
                    seen += len(chunk)
                    if progress is not None:
                        progress(seen, total)

            if total is not None and seen != total:
                raise DownloadError(
                    f"truncated transfer: received {seen} of {total} advertised bytes"
                )
            if seen == 0:
                raise DownloadError("server returned an empty body")

            if verify is not None:
                passed, message = verify(tmp)
                if not passed:
                    raise DownloadError(f"verification failed: {message}")

            os.replace(tmp, dest)
            return dest
        except BaseException as exc:  # includes KeyboardInterrupt: always clean up
            if handle_owned:
                # Only close a descriptor we still own; closing one that
                # os.fdopen() already took could close an unrelated resource
                # that reused the number.
                try:
                    os.close(handle)
                except OSError:
                    pass
            tmp.unlink(missing_ok=True)
            if not isinstance(exc, Exception):
                raise
            last_error = exc
            if attempt < retries:
                time.sleep(backoff * 2 ** (attempt - 1))

    raise DownloadError(f"failed to download {url} after {retries} attempts: {last_error}")


def _discovered_records(
    dataset: DatasetRecord, timeout: float, log: Callable[[str], None]
) -> tuple[list[FileRecord], str | None]:
    """Query Zenodo for a dataset's files, keeping published checksums and sizes."""
    record_id = str(dataset.data.get("doi", "")).rsplit(".", 1)[-1]
    try:
        listing = zenodo_files(record_id, timeout=timeout)
    except RuntimeError as exc:
        log(f"  ! discovery failed: {exc}")
        return [], str(exc)

    files: list[FileRecord] = []
    for item in listing:
        algorithm, digest = _split_checksum(item.get("checksum"))
        log(f"  discovered: {item['name']}  {item['size']} B  {item.get('checksum')}")
        files.append(
            FileRecord(
                dataset=dataset.key,
                name=item["name"],
                url=item["url"],
                published_md5=digest,
                published_algorithm=algorithm or "md5",
                checksum_source=(
                    f"Zenodo API record {record_id}, files[].checksum, "
                    "read at download time"
                )
                if digest
                else None,
                expected_bytes=item.get("size"),
                expected_bytes_tolerance=0.0,
                required=not dataset.optional,
            )
        )
    return files, None


def fetch_dataset(
    dataset: DatasetRecord,
    *,
    force: bool = False,
    discover: bool = False,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    log: Callable[[str], None] = print,
) -> list[DownloadResult]:
    """Obtain every file of one dataset and verify it.

    Parameters
    ----------
    dataset : DatasetRecord
        Registry entry.
    force : bool, optional
        Re-download even when a verified local copy exists. The existing file is
        replaced only after the new one passes verification.
    discover : bool, optional
        For Zenodo datasets with no declared files, query the Zenodo API and
        download everything the record contains, enforcing its checksums.
    timeout, retries : optional
        Passed to :func:`download_file`.
    log : callable, optional
        Line logger.

    Returns
    -------
    list of DownloadResult
        Includes a ``failed`` entry when discovery itself failed, so that an
        unreachable repository can never be mistaken for an empty dataset.
    """
    results: list[DownloadResult] = []
    files: list[FileRecord] = list(dataset.files)

    if not files and not (discover and dataset.kind == "zenodo"):
        # A required dataset with nothing to fetch is incomplete, not a success.
        return [
            DownloadResult(
                FileRecord(dataset=dataset.key, name="<no files declared>", url=None,
                           required=not dataset.optional),
                "failed",
                "no files are declared for this dataset and discovery is disabled",
            )
        ]

    if discover and dataset.kind == "zenodo" and not files:
        discovered, error = _discovered_records(dataset, timeout, log)
        if error is not None:
            placeholder = FileRecord(
                dataset=dataset.key,
                name="<discovery>",
                url=None,
                required=not dataset.optional,
            )
            results.append(DownloadResult(placeholder, "failed", error))
            return results
        files = discovered
        if not files:
            log("  ! record lists no files")
            results.append(
                DownloadResult(
                    FileRecord(dataset=dataset.key, name="<no files>", url=None,
                               required=not dataset.optional),
                    "failed",
                    "the repository record lists no files",
                )
            )
            return results

    for record in files:
        if record.is_glob:
            present = sorted(dataset.raw_dir.glob(record.name)) if dataset.raw_dir.is_dir() else []
            if len(present) > 1:
                results.append(
                    DownloadResult(
                        record,
                        "failed",
                        f"{len(present)} files match {record.name!r} "
                        f"({[p.name for p in present]}); leave exactly one",
                    )
                )
                continue
            if present:
                verification = verify_file(record, present[0])
                verification["name"] = present[0].name  # record the real name, not the pattern
                status = "cached" if verification["ok"] else "failed"
                message = (
                    f"found {present[0].name}"
                    if verification["ok"]
                    else f"{present[0].name} failed verification: {verification}"
                )
                results.append(DownloadResult(record, status, message, verification))
            else:
                results.append(
                    DownloadResult(
                        record,
                        "manual",
                        f"retrieve by hand into {dataset.raw_dir} "
                        "(see docs/MANUAL_DOWNLOADS.md)",
                    )
                )
            continue

        if record.url is None:
            results.append(DownloadResult(record, "skipped", "no URL in registry"))
            continue

        if record.path.is_file() and not force:
            verification = verify_file(record)
            if verification["ok"]:
                results.append(DownloadResult(record, "cached", "verified local copy", verification))
                continue
            log(f"  local copy of {record.name} failed verification; refetching")

        log(f"  downloading {record.name} ...")

        def _verify(temp_path: Path, record: FileRecord = record) -> tuple[bool, str]:
            check = verify_file(record, temp_path)
            if check["ok"]:
                return True, "ok"
            return False, (
                f"checksum_ok={check['checksum_ok']} size_ok={check['size_ok']} "
                f"bytes={check['bytes']} md5={check['md5']}"
            )

        try:
            download_file(
                record.url, record.path, timeout=timeout, retries=retries, verify=_verify
            )
        except DownloadError as exc:
            results.append(DownloadResult(record, "failed", str(exc)))
            continue

        verification = verify_file(record)
        status = "ok" if verification["ok"] else "failed"
        results.append(
            DownloadResult(
                record,
                status,
                "downloaded" if status == "ok" else f"post-write verification failed: {verification}",
                verification,
            )
        )

    return results


def fetch_all(
    registry: Registry | None = None,
    *,
    only: list[str] | None = None,
    force: bool = False,
    discover: bool = True,
    include_optional: bool = False,
    timeout: float = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    log: Callable[[str], None] = print,
) -> list[DownloadResult]:
    """Fetch every registered dataset and update the provenance manifest.

    Parameters
    ----------
    registry : Registry, optional
        Defaults to the on-disk registry.
    only : list of str, optional
        Restrict to these dataset keys. Unknown keys raise.
    force : bool, optional
        Ignore cached copies.
    discover : bool, optional
        Query Zenodo for datasets with no declared files.
    include_optional : bool, optional
        Also fetch datasets flagged ``optional: true``.
    timeout, retries : optional
        Passed through to :func:`fetch_dataset`.
    log : callable, optional
        Line logger.

    Returns
    -------
    list of DownloadResult

    Raises
    ------
    KeyError
        If ``only`` names a dataset that is not in the registry: a misspelled
        key must not look like a successful no-op.
    """
    registry = registry or load_registry()
    if only:
        unknown = sorted(set(only) - set(registry.datasets))
        if unknown:
            raise KeyError(
                f"unknown dataset key(s) {unknown}; registry has {sorted(registry.datasets)}"
            )

    results: list[DownloadResult] = []
    for dataset in registry:
        if only and dataset.key not in only:
            continue
        if dataset.optional and not include_optional and not only:
            continue
        if dataset.kind == "none":
            continue

        log(f"[{dataset.key}] {dataset.title}")
        dataset.raw_dir.mkdir(parents=True, exist_ok=True)
        results.extend(
            fetch_dataset(
                dataset,
                force=force,
                discover=discover,
                timeout=timeout,
                retries=retries,
                log=log,
            )
        )

    recorded = []
    for result in results:
        entry = dict(result.verification) if result.verification else {
            "dataset": result.record.dataset,
            "name": result.record.name,
            "path": str(result.record.path),
            "exists": False,
            "bytes": None,
            "sha256": None,
            "ok": False,
        }
        entry["status"] = result.status
        entry["message"] = result.message
        entry["url"] = result.record.url
        recorded.append(entry)

    manifest = write_manifest(recorded)
    log(f"manifest -> {manifest}")
    return results
