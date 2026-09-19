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
    from .data import (baumgardt2019, baumgardt_catalogue, fimbulthul, kuzma2025, kuzma2026,
                       omegacat, trager1995, vasiliev2021)

    ensure_data_tree()
    reports: list[dict[str, Any]] = []
    modules = {
        "kuzma2025": kuzma2025,
        "kuzma2026": kuzma2026,
        "fimbulthul": fimbulthul,
        "omegacat": omegacat,
        "vasiliev2021": vasiliev2021,
        "baumgardt": baumgardt_catalogue,
        "trager1995": trager1995,
        "baumgardt2019": baumgardt2019,
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


def cmd_fit(args: argparse.Namespace) -> int:
    """Run one nested-sampling fit (K1 or K2) and write it under results/fits/<label>."""
    import numpy as np

    from .kinematics import DarkMatterModel, FitProblem, KinematicData, NoDarkMatterModel, run_nested
    from .paths import results_dir

    preset = None
    if args.preset:
        from .kinematics.presets import build_preset
        preset, family, datasets = build_preset(args.preset)
        print(f"preset {preset.name}: {preset.reference}\n  {preset.url}\n  {preset.notes}")
    else:
        datasets = args.datasets.split(",")
    data = KinematicData.load(datasets, gaia_edr3_pm={"r_min_arcsec": args.gaia_r_min})
    kw = dict(tracer=args.tracer, backend=args.backend)
    if getattr(args, "no_scales", False):
        kw["instruments"] = ()
    if preset is not None:
        pass
    elif args.family == "K1":
        family = NoDarkMatterModel(**kw)
    elif args.family in ("K2-cored", "K2-nfw"):
        family = DarkMatterModel(gamma=0.0 if args.family == "K2-cored" else 1.0, **kw)
    else:
        print(f"unknown family {args.family!r}", file=sys.stderr)
        return 2
    problem = FitProblem(family, data)
    data_provenance: dict[str, object] = {"kind": "real"}
    if args.mock_from:
        # the mock is generated by the family that produced the parameter vector
        # (default K1), then fitted with --family; this is the false-positive test
        x_true = np.load(args.mock_from)
        if args.mock_family == "K1":
            gen = NoDarkMatterModel(tracer=args.tracer)
        else:
            gen = DarkMatterModel(gamma=0.0 if args.mock_family == "K2-cored" else 1.0, tracer=args.tracer)
        if len(x_true) != len(gen.names):
            print(f"--mock-from vector has {len(x_true)} entries; {gen.label} has {len(gen.names)} parameters", file=sys.stderr)
            return 2
        rng = np.random.default_rng(args.seed)
        problem = FitProblem(family, FitProblem(gen, data).mock_data(x_true, rng))
        data_provenance = {"kind": "mock", "generating_family": gen.label, "seed": args.seed,
                           "source_file": str(args.mock_from), "tracer": args.tracer,
                           "truth": {n: float(v) for n, v in zip(gen.names, x_true)}}
        print(f"mock data generated from {gen.label} with seed {args.seed}: " + ", ".join(f"{n}={v:.4g}" for n, v in zip(gen.names, x_true)))
    label = args.label or (f"preset_{preset.name}" if preset else family.label + ("_mock" if args.mock_from else ""))
    out = results_dir() / "fits" / label
    print(f"fit {family.label}: {len(family.names)} parameters, {problem.data.n_points} points -> {out}")
    summary = run_nested(problem, out, n_live=args.n_live, dlogz=args.dlogz, seed=args.seed,
                         max_ncalls=args.max_ncalls, verbose=args.verbose, step_sampler=args.step_sampler,
                         data_provenance=data_provenance)
    print(f"logZ = {summary['logz']:.2f} +- {summary['logzerr']:.2f}; chi2_ml = {summary['chi2_ml_total']:.1f} "
          f"/ {summary['n_points']} points; {summary['n_calls']} calls in {summary['elapsed_s']:.0f} s")
    for name, q in summary["parameters"].items():
        print(f"  {name:10s} {q['p50']:12.4g}  [{q['p16']:.4g}, {q['p84']:.4g}]  {q['unit']}")
    if preset is not None:
        from astropy.table import Table
        from .kinematics.presets import derived_quantities
        post = Table.read(out / "posterior.ecsv")
        rows = [derived_quantities(preset, {n: float(r[n]) for n in family.names}, family) for r in post]
        red = summary["chi2_ml_total"] / max(summary["n_points"] - len(family.names), 1)
        print(f"\nlike-for-like comparison (chi2/dof of the preset model on our data = {red:.1f}"
              + ("; the pulls below use only the published error and are NOT meaningful for a model this poor)" if red > 2 else ")"))
        for key, (val, err, note) in preset.published.items():
            ours = np.array([r[key] for r in rows if key in r])
            if len(ours) == 0:
                continue
            lo, med, hi = np.percentile(ours, [16, 50, 84])
            if np.isfinite(err):
                pull = (med - val) / np.hypot(err, 0.5 * (hi - lo))
                print(f"  {key:9s} ours {med:.4g} [{lo:.4g}, {hi:.4g}]  published {val:.4g} +- {err:.2g}  ({pull:+.1f} sigma)  -- {note}")
            else:
                print(f"  {key:9s} ours 95% upper limit {np.percentile(ours, 95):.4g}  published limit {val:.4g}  -- {note}")
    return 0


#: Gaia EDR3 now enters as our own measurement beyond 460 arcsec rather than the published
#: spline, whose inner points continue inward over a region with no usable Gaia star
#: (JOURNAL 2026-09-19). The old key "gaia_edr3_pm" remains available via --datasets.
DEFAULT_DATASETS = ("hst_pm_radial_ours,hst_pm_tangential_ours,muse_los_dispersion,"
                    "gaia_dr2_pm,gaia_edr3_ours_radial,gaia_edr3_ours_tangential")


def cmd_plot_constraints(args: argparse.Namespace) -> int:
    """Write the constraint-map and outer-tracer-audit figures."""
    from .plotting.constraints import plot_constraint_map, plot_contamination_model, plot_outer_tracer_audit

    from .plotting.constraints import (fit_quality_table, plot_annulus_fits, plot_dataset_step,
                                       plot_extended_profile, plot_hst_gaia_star_by_star,
                                       plot_master_datasets, plot_method_comparison,
                                       plot_offset_explained, plot_periphery,
                                       plot_periphery_density, plot_pm_datasets,
                                       plot_residual_significance)

    t = fit_quality_table()
    print("  per-annulus fit quality (chi2 per bin of the projected histogram):")
    print("    r [arcsec]     N    f_field   sigma_R        sigma_T        chi2/bin R (all, peak)   T (all, peak)")
    for row in t:
        print("    %4.0f-%4.0f %7d   %.3f   %.3f±%.3f  %.3f±%.3f     %.2f, %.2f          %.2f, %.2f" % (
            row["r_lower"], row["r_upper"], row["n_stars"], row["f_field"], row["sigma_pmr"], row["sigma_pmr_err"],
            row["sigma_pmt"], row["sigma_pmt_err"], row["chi2_r_wide"], row["chi2_r_peak"],
            row["chi2_t_wide"], row["chi2_t_peak"]))
    for path in (plot_master_datasets(),
                 plot_pm_datasets(),
                 plot_extended_profile(),
                 plot_periphery(),
                 plot_periphery_density(),
                 plot_hst_gaia_star_by_star(),
                 plot_offset_explained(),
                 plot_dataset_step(),
                 plot_residual_significance(),
                 plot_method_comparison(),
                 plot_annulus_fits("plots/outer_fit_annuli_radial.png", component="r"),
                 plot_annulus_fits("plots/outer_fit_annuli_tangential.png", component="t"),
                 # four constraint maps: {field modelled, P > 0.9} x {rotation from our own
                 # fit, rotation from the published curve}
                 plot_constraint_map("plots/constraint_map_selfconsistent.png", k1=args.k1, k2=args.k2,
                                     contamination_modelled=True, streaming="self"),
                 plot_constraint_map("plots/constraint_map.png", k1=args.k1, k2=args.k2,
                                     contamination_modelled=True, streaming="published"),
                 plot_constraint_map("plots/constraint_map_pcut_selfconsistent.png", k1=args.k1, k2=args.k2,
                                     contamination_modelled=False, streaming="self"),
                 plot_constraint_map("plots/constraint_map_pcut.png", k1=args.k1, k2=args.k2,
                                     contamination_modelled=False, streaming="published"),
                 plot_outer_tracer_audit(reference="pcut"),
                 plot_outer_tracer_audit("plots/outer_tracer_audit_mixture.png", reference="mixture"),
                 plot_contamination_model()):
        print(f"  wrote {path}")
    return 0


def cmd_fetch_field_template(args: argparse.Namespace) -> int:
    """Fetch the independent Gaia DR3 field annulus used as the contamination template."""
    from .selection.field_template import build_product

    path = build_product(g_max=args.g_max)
    from astropy.table import Table
    t = Table.read(path)
    print(f"wrote {path}: {t.meta['n_kept']:,} field stars from {t.meta['n_fetched']:,} fetched, "
          f"{t.meta['surface_density_per_arcmin2']:.2f} per arcmin^2")
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
    pre.add_argument("--dataset", action="append", help="kuzma2025|kuzma2026|fimbulthul|omegacat|vasiliev2021|baumgardt|trager1995|baumgardt2019")
    pre.set_defaults(func=cmd_preprocess)

    inv = subparsers.add_parser("inventory", help="print the dataset inventory")
    inv.set_defaults(func=cmd_inventory)

    plot = subparsers.add_parser("plot-data", help="write overview PNGs of the ingested data")
    plot.add_argument("--only", action="append", help="function name, e.g. plot_sky_overview")
    plot.set_defaults(func=cmd_plot_data)

    fit = subparsers.add_parser("fit", help="run a nested-sampling Jeans fit (K1 / K2-cored / K2-nfw)")
    fit.add_argument("--family", default="K1", choices=["K1", "K2-cored", "K2-nfw"])
    fit.add_argument("--preset", default=None, choices=["watkins2013", "omegacat6", "baumgardt2018", "imbh_limit"],
                     help="run under a published analysis's assumptions and compare with its numbers")
    fit.add_argument("--datasets", default=DEFAULT_DATASETS, help="comma-separated dataset keys")
    fit.add_argument("--no-scales", action="store_true",
                     help="remove the per-instrument multiplicative nuisances (every dataset compared with the same model)")
    fit.add_argument("--tracer", default="composite", choices=["composite", "trager"],
                     help="tracer density: HST star counts inside 25 arcsec + Trager light (default), or Trager only")
    fit.add_argument("--gaia-r-min", type=float, default=300.0, help="inner cut for the Gaia EDR3 profile, arcsec")
    fit.add_argument("--backend", default="jeans", choices=["jeans", "jam", "agama"],
                     help="dynamical engine: our spherical Jeans solver (default), JamPy, or an AGAMA Cuddeford-Osipkov-Merritt DF")
    fit.add_argument("--n-live", type=int, default=400)
    fit.add_argument("--dlogz", type=float, default=0.5)
    fit.add_argument("--max-ncalls", type=int, default=None)
    fit.add_argument("--seed", type=int, default=42)
    fit.add_argument("--label", default=None, help="output folder name under results/fits (default: family label)")
    fit.add_argument("--mock-from", default=None, help=".npy parameter vector: fit mock data generated from it")
    fit.add_argument("--mock-family", default="K1", choices=["K1", "K2-cored", "K2-nfw"],
                     help="family that generated the --mock-from vector (default K1)")
    fit.add_argument("--verbose", action="store_true")
    fit.add_argument("--step-sampler", action="store_true", help="use a slice step sampler instead of MLFriends rejection")
    fit.set_defaults(func=cmd_fit)

    fft = subparsers.add_parser("fetch-field-template", help="WSDB: Gaia DR3 field annulus for the contamination model")
    fft.add_argument("--g-max", type=float, default=20.5)
    fft.set_defaults(func=cmd_fetch_field_template)

    pcon = subparsers.add_parser("plot-constraints", help="figures: what constrains the mass where, and the outer-tracer audit")
    pcon.add_argument("--k1", default="K1_noDM_composite", help="finished no-DM run to draw")
    pcon.add_argument("--k2", default="K2_cored_composite", help="finished dark-matter run to draw")
    pcon.set_defaults(func=cmd_plot_constraints)

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
