"""Command-line interface for the oCen_dm pipeline.

Milestone 1 commands only::

    ocen fetch-data          # download raw data, verify checksums, write manifest
    ocen preprocess          # standardize + validate -> data/processed
    ocen inventory           # print the dataset inventory (no inference)
    ocen inspect-omegacat    # dump real column structure of oMEGACat raw files
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, Sequence

from .paths import ensure_data_tree, processed_dir, project_root
from .provenance import load_registry, verify_file

__all__ = ["main", "build_parser"]


def _fmt_bytes(n: int | None) -> str:
    """Format a byte count as a short human-readable string."""
    if n is None:
        return "-"
    value = float(n)
    for unit in ("B", "kB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return f"{value:.1f} GB"


def cmd_fetch_data(args: argparse.Namespace) -> int:
    """Download every registered dataset and update the provenance manifest.

    Returns 0 only when every required file is present and passed verification.
    """
    from .data.download import fetch_all

    ensure_data_tree()
    try:
        results = fetch_all(
            only=args.dataset or None,
            force=args.force,
            discover=not args.no_discover,
            include_optional=args.include_optional,
            timeout=args.timeout,
            retries=args.retries,
        )
    except KeyError as exc:
        print(f"usage error: {exc}", file=sys.stderr)
        return 2

    print("\nsummary")
    incomplete = 0
    for result in results:
        level = (result.verification or {}).get("verification", "none")
        mark = "verified" if result.verified else (
            f"present, {level}-checked only" if result.ok else "")
        print(
            f"  {result.status:8s} {result.record.dataset}/{result.record.name}"
            f"{'  [' + mark + ']' if mark else ''}"
            f"{'  -- ' + result.message if result.message else ''}"
        )
        if result.record.required and not result.ok:
            incomplete += 1

    manual = [r for r in results if r.status == "manual"]
    if manual:
        print(f"\n{len(manual)} file(s) need manual retrieval: see docs/MANUAL_DOWNLOADS.md")
    if incomplete:
        print(
            f"\n{incomplete} required file(s) are missing or unverified; "
            "the dataset is incomplete."
        )
        return 1
    return 0


def cmd_preprocess(args: argparse.Namespace) -> int:
    """Standardize the raw catalogues into validated processed products.

    A schema violation aborts that dataset's product: nothing is written, and
    the command exits non-zero. An invalid table is never published.
    """
    from .data import baumgardt_catalogue, fimbulthul, kuzma2025, kuzma2026, omegacat, vasiliev2021

    ensure_data_tree()
    reports: list[dict[str, Any]] = []
    modules = {
        "kuzma2025": kuzma2025,
        "kuzma2026": kuzma2026,
        "fimbulthul": fimbulthul,
        "omegacat": omegacat,
        "vasiliev2021": vasiliev2021,
        "baumgardt": baumgardt_catalogue,
    }
    selected = args.dataset or list(modules)

    status = 0
    for key in selected:
        module = modules.get(key)
        if module is None:
            print(f"unknown dataset {key!r}; known: {sorted(modules)}", file=sys.stderr)
            status = 2
            continue
        try:
            result = module.preprocess()
        except Exception as exc:  # noqa: BLE001 - reported per dataset, never silent
            print(f"[{key}] FAILED: {type(exc).__name__}: {exc}")
            reports.append({"dataset": key, "error": f"{type(exc).__name__}: {exc}"})
            status = max(status, 1)
            continue
        for item in result if isinstance(result, list) else [result]:
            reports.append(item)
            if "error" in item:
                print(f"[{key}] FAILED {item.get('profile', '')}: {item['error'].splitlines()[0]}")
                status = max(status, 1)
            elif "skipped" in item:
                print(f"[{key}] skipped {item.get('profile', '')}: {item['skipped'].splitlines()[0]}")
                status = max(status, 1)
            else:
                print(f"[{key}] {item.get('n_rows', '?')} rows -> {item.get('path')}")

    report_path = processed_dir() / "preprocess_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(reports, fh, indent=2, default=str)
        fh.write("\n")
    print(f"report -> {report_path}")
    if status:
        print("one or more products could not be built; nothing invalid was written.")
    return status


def cmd_inventory(args: argparse.Namespace) -> int:
    """Print a concise dataset inventory. Performs no scientific inference."""
    registry = load_registry()
    print(f"oCen_dm dataset inventory   root={project_root()}")
    print(f"registry: {registry.path}  (schema v{registry.schema_version})\n")

    for dataset in registry:
        flag = " [optional]" if dataset.optional else ""
        print(f"{dataset.key}{flag}")
        print(f"  {dataset.title}")
        print(f"  role: {dataset.role.strip()}")
        paper = dataset.paper
        link = paper.get("arxiv") or paper.get("journal") or (
            f"https://doi.org/{paper['doi']}" if paper.get("doi") else None
        )
        verified = paper.get("verified")
        if link:
            mark = {True: "verified", False: "UNVERIFIED", None: "unchecked"}[verified]
            print(f"  paper: {link}  [{mark}]")
        if dataset.data.get("doi"):
            print(f"  data:  https://doi.org/{dataset.data['doi']}  (kind={dataset.kind})")
        else:
            print(f"  data:  kind={dataset.kind}")

        if not dataset.files:
            if dataset.kind == "zenodo":
                print("  files: none declared -- run `ocen fetch-data` to discover")
            print()
            continue

        for record in dataset.files:
            # A wildcard record names a pattern, not a file: resolve it the way
            # fetch-data does, or the inventory reports a present file as MISSING.
            target = None
            if record.is_glob:
                matches = sorted(dataset.raw_dir.glob(record.name)) if dataset.raw_dir.is_dir() else []
                target = matches[0] if len(matches) == 1 else None
            check = verify_file(record, target)
            if check["exists"]:
                digest = {
                    True: "checksum ok",
                    False: "CHECKSUM MISMATCH",
                    None: "no published checksum",
                }[check["checksum_ok"]]
                state = f"present  {_fmt_bytes(check['bytes'])}  {digest}"
                if not check["ok"]:
                    state += "  [NOT USABLE]"
            else:
                state = "MISSING"
            shown = target.name if target is not None else record.name
            print(f"  file:  {shown:32s} {state}")
        print()

    processed = sorted(processed_dir().rglob("*.ecsv")) if processed_dir().is_dir() else []
    print(f"processed products: {len(processed)}")
    for path in processed:
        print(f"  {path.relative_to(project_root())}")
    return 0


def cmd_inspect_omegacat(args: argparse.Namespace) -> int:
    """List downloaded oMEGACat files and dump their real column structure."""
    from .data import omegacat

    files = omegacat.inventory()
    if not files:
        print("no raw oMEGACat VI files found; run `ocen fetch-data` first")
        return 1

    for item in files:
        print(f"{item['name']}  {_fmt_bytes(item['bytes'])}")
        if args.columns:
            try:
                report = omegacat.describe(item["name"])
            except Exception as exc:  # noqa: BLE001 - some files are not tables
                print(f"    (not readable as a table: {type(exc).__name__}: {exc})")
                continue
            print(f"    rows: {report['n_rows']}")
            for column in report["columns"]:
                unit = f" [{column['unit']}]" if column["unit"] else ""
                print(f"    - {column['name']}  {column['dtype']}{unit}")
    return 0


def cmd_plot_data(args: argparse.Namespace) -> int:
    """Write overview PNGs of every ingested dataset into plots/."""
    import matplotlib

    matplotlib.use("Agg")
    from .plotting.data_overview import plot_all

    for path in plot_all(only=args.only or None):
        print(f"  wrote {path.relative_to(project_root())}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the ``ocen`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="ocen", description="Omega Cen dark-matter / tidal-tail pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch = subparsers.add_parser("fetch-data", help="download and verify raw data")
    fetch.add_argument("--dataset", action="append", help="restrict to this dataset key")
    fetch.add_argument("--force", action="store_true", help="ignore cached copies")
    fetch.add_argument("--no-discover", action="store_true", help="skip Zenodo file discovery")
    fetch.add_argument("--include-optional", action="store_true")
    fetch.add_argument("--timeout", type=float, default=60.0, help="socket timeout (s)")
    fetch.add_argument("--retries", type=int, default=3, help="attempts per file")
    fetch.set_defaults(func=cmd_fetch_data)

    pre = subparsers.add_parser("preprocess", help="build validated processed products")
    pre.add_argument("--dataset", action="append", help="kuzma2025|kuzma2026|fimbulthul|omegacat|vasiliev2021|baumgardt")
    pre.set_defaults(func=cmd_preprocess)

    inv = subparsers.add_parser("inventory", help="print the dataset inventory")
    inv.set_defaults(func=cmd_inventory)

    plot = subparsers.add_parser("plot-data", help="write overview PNGs of the ingested data")
    plot.add_argument("--only", action="append", help="function name, e.g. plot_sky_overview")
    plot.set_defaults(func=cmd_plot_data)

    insp = subparsers.add_parser("inspect-omegacat", help="dump oMEGACat raw file structure")
    insp.add_argument("--columns", action="store_true", help="also list columns")
    insp.set_defaults(func=cmd_inspect_omegacat)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point. Returns a process exit status."""
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
