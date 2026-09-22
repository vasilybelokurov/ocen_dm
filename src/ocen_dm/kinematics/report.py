"""Summaries and figures for saved nested-sampling runs.

New runs replay saved likelihood arrays and resolved models. Legacy runs use
their recorded options and current products, with an explicit warning. Reports
check that the reconstructed best-sample likelihood matches the saved result.
"""

from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.table import Table

from ..paths import results_dir

__all__ = ["load_run", "comparison_table", "posterior_upper_limit", "write_report", "engine_crosscheck"]

_LOG_PARAMS = ("M_star", "M_rem", "a_rem", "M_bh", "r_beta", "M_dm_100", "r_s")


def load_run_metadata(d: Path) -> dict:
    """Read a run's configuration without requiring posterior/profile products."""
    d = Path(d)
    if not (d / "summary.json").exists():
        raise FileNotFoundError(f"run has no summary.json under {d}")
    out = {"label": d.name, "dir": d, "summary": json.loads((d / "summary.json").read_text())}
    if (d / "run.yaml").exists():          # input hashes and mock provenance, for comparability
        import yaml
        out["run"] = yaml.safe_load((d / "run.yaml").read_text()) or {}
    for key in ("model", "data_snapshot", "dataset_options", "data"):
        if key in out.get("run", {}):
            if key in out["summary"] and out["summary"][key] != out["run"][key]:
                raise ValueError(f"{d}: summary.json and run.yaml disagree on {key}")
            out["summary"][key] = out["run"][key]
    return out


def load_run(label: str) -> dict:
    """Read metadata, posterior samples and the profile stack of one run."""
    d = results_dir() / "fits" / label
    out = load_run_metadata(d)
    out["label"] = label
    out["posterior"] = Table.read(d / "posterior.ecsv")
    if (d / "profiles.npz").exists():
        out["profiles"] = dict(np.load(d / "profiles.npz"))
    return out


def _fmt(name: str, q: dict) -> str:
    if name in _LOG_PARAMS:
        return f"{q['p50']:.3g} [{q['p16']:.3g}, {q['p84']:.3g}]"
    return f"{q['p50']:.3f} [{q['p16']:.3f}, {q['p84']:.3f}]"


def posterior_upper_limit(run: dict, name: str, level: float = 0.95) -> float:
    """One-sided posterior upper limit on a parameter (e.g. ``M_dm_100`` at 95 per cent)."""
    return float(np.percentile(np.asarray(run["posterior"][name]), 100 * level))


def comparison_table(runs: Iterable[dict]) -> str:
    """Markdown table: evidence, chi2 at the best sample, and every parameter's median and 16-84 interval."""
    runs = list(runs)
    if not runs:
        raise ValueError("at least one run is required")
    names: list[str] = []
    for r in runs:
        for n in r["summary"]["parameters"]:
            if n not in names:
                names.append(n)
    # New runs compare actual arrays. For legacy runs retain conservative file hashes,
    # observational switches and the writer's real mock schema (seed, truth, generator).
    def _fingerprint(r):
        run = r.get("run", {}) or {}
        summary = r["summary"]
        snapshot = summary.get("data_snapshot") or run.get("data_snapshot")
        if snapshot:
            from .run_io import data_fingerprint, read_data_snapshot
            # Validate the file too: a damaged/missing snapshot is not evidence of identity.
            return ("snapshot", data_fingerprint(read_data_snapshot(r["dir"], snapshot)))
        data = run.get("data") or summary.get("data") or {}
        inputs = run.get("inputs")
        if not inputs or data.get("kind") not in ("real", "mock"):
            return None
        mock = None
        if data["kind"] == "mock":
            if any(k not in data for k in ("source_file", "generating_family", "truth", "seed")):
                return None
            mock = {k: data.get(k) for k in ("source_file", "generating_family", "truth", "seed",
                                            "tracer", "model", "noise_model")}
        opts = summary.get("dataset_options", {}) | run.get("dataset_options", {})
        observing = {"gaia_errors": opts.get("gaia_errors", "raw"),
                     "gaia_rotation": opts.get("gaia_rotation", "published"),
                     "gaia_r_min": opts.get("gaia_r_min", 300.0)}
        return (tuple(sorted(r["summary"]["datasets"])),
                int(r["summary"].get("n_points", -1)),
                tuple(sorted(inputs.items())), json.dumps(observing, sort_keys=True),
                data["kind"], json.dumps(mock, sort_keys=True))
    prints = [_fingerprint(r) for r in runs]
    datasets = [tuple(sorted(r["summary"]["datasets"])) for r in runs]
    common = datasets[0] if prints[0] is not None and all(p == prints[0] for p in prints) else None
    if common is None and all(d == datasets[0] for d in datasets):
        lines_note = ("same dataset keys but observations differ or identity is unverified "
                      "(input hashes, point count or mock provenance): evidences are NOT comparable")
    else:
        lines_note = None
    ref = max(r["summary"]["logz"] for r in runs) if common else None
    lines = ["| quantity | " + " | ".join(r["label"] for r in runs) + " |",
             "|---|" + "---|" * len(runs)]
    lines.append("| ln Z | " + " | ".join(f"{r['summary']['logz']:.2f} ± {r['summary']['logzerr']:.2f}" for r in runs) + " |")
    lines.append("| Δ ln Z vs best | " + " | ".join((f"{r['summary']['logz'] - ref:+.2f}" if ref is not None else "n/a (different data)") for r in runs) + " |")
    lines.append("| max ln L | " + " | ".join(f"{r['summary']['lnL_max']:.1f}" for r in runs) + " |")
    lines.append("| χ² at max L / N | " + " | ".join(f"{r['summary']['chi2_ml_total']:.0f} / {r['summary']['n_points']}" for r in runs) + " |")
    for ds in runs[0]["summary"]["chi2_ml"]:
        lines.append(f"| χ² {ds} | " + " | ".join(f"{r['summary']['chi2_ml'].get(ds, {}).get('chi2', float('nan')):.0f} / {r['summary']['chi2_ml'].get(ds, {}).get('n', 0)}" for r in runs) + " |")
    if lines_note:
        lines.append("| **warning** | " + lines_note + " |" + " |" * (len(runs) - 1))
    lines.append("| likelihood calls / time | " + " | ".join(f"{r['summary']['n_calls']:,} / {r['summary']['elapsed_s'] / 60:.0f} min" for r in runs) + " |")
    for n in names:
        cells = []
        for r in runs:
            q = r["summary"]["parameters"].get(n)
            cells.append(_fmt(n, q) if q else "—")
        unit = next((r["summary"]["parameters"][n]["unit"] for r in runs if n in r["summary"]["parameters"]), "")
        lines.append(f"| {n} {('[' + unit + ']') if unit else ''} | " + " | ".join(cells) + " |")
    for r in runs:
        if "M_dm_100" in r["posterior"].colnames:
            lines.append(f"| M_DM(<100 pc) 95 % upper limit, {r['label']} | {posterior_upper_limit(r, 'M_dm_100'):.3g} M☉ |" + " |" * (len(runs) - 1))
    return "\n".join(lines)


def write_report(labels: Iterable[str], path: Path | None = None, plots_dir: Path | None = None) -> Path:
    """Comparison table plus the posterior-profile and data-vs-model figures for finished runs."""
    from ..plotting.fits import plot_posterior_profiles, plot_profile_fit

    runs = [load_run(l) for l in labels]
    problems = [problem_for(r) for r in runs]  # validate replay before writing report products
    path = path or results_dir() / "fits" / "comparison.md"
    plots_dir = plots_dir or (results_dir().parent / "plots")
    md = ["# Nested-sampling runs: comparison\n", comparison_table(runs), ""]
    # independent-engine check of every best sample (JamPy shares the parametrisation)
    try:
        import jampy  # noqa: F401
        md.append("## Engine cross-check of the best samples (our Jeans solver vs JamPy)\n")
        for r, problem in zip(runs, problems):
            if problem.family.backend == "agama":
                md.append(f"{r['label']}: AGAMA has a different parametrisation; engine cross-check skipped.\n")
                continue
            md += [f"**{r['label']}**\n", engine_crosscheck(r["label"], ("jeans", "jam")), ""]
    except ImportError:
        md.append("(JamPy not installed: engine cross-check skipped)\n")
    fig = plot_posterior_profiles({r["label"]: r["dir"] for r in runs if "profiles" in r},
                                  plots_dir / "fit_posterior_profiles.png", title="Enclosed mass, dark fraction, circular speed (16-84 %)")
    md.append(f"![Posterior mass profiles]({os.path.relpath(fig, path.parent)})\n")
    for r, P in zip(runs, problems):
        s = r["summary"]
        prov = s.get("data", {"kind": "real"})
        fam = P.family
        x_ml = np.array([s["parameters"][n]["ml"] for n in fam.names])
        samples = np.array([np.asarray(r["posterior"][n]) for n in fam.names]).T
        suffix = "" if prov.get("kind") != "mock" else (
            f"  [MOCK data from {prov['generating_family']}, seed {prov['seed']}]")
        out = plot_profile_fit(P, x_ml, plots_dir / f"fit_{r['label']}_posterior_profiles.png",
                               title=f"{r['label']}: best sample and 16-84 % posterior band{suffix}", samples=samples)
        md.append(f"![{r['label']} fitted profiles]({os.path.relpath(out, path.parent)})")
    path.write_text("\n".join(md))
    return path


def data_for(summary: dict, *, run_dir: Path | None = None):
    """Replay a snapshot, or load legacy profiles with the recorded Gaia options."""
    from .fit import FitProblem
    from .likelihood import KinematicData

    prov = summary.get("data", {"kind": "real"})
    if "data_snapshot" in summary:
        from .run_io import read_data_snapshot
        if run_dir is None:
            raise ValueError("run_dir is required to read this run's data snapshot")
        return read_data_snapshot(Path(run_dir), summary["data_snapshot"]), prov
    warnings.warn("Legacy run has no data snapshot; replay uses recorded options and current "
                  "products. Verify its saved likelihood before interpreting the report.",
                  UserWarning, stacklevel=2)
    opts = summary.get("dataset_options", {})
    gaia = {"error_model": opts.get("gaia_errors", "raw"),
            "rotation": opts.get("gaia_rotation", "published")}
    data = KinematicData.load(list(summary["datasets"]),
                             gaia_edr3_pm={"r_min_arcsec": opts.get("gaia_r_min", 300.0)},
                             gaia_edr3_ours_radial=gaia, gaia_edr3_ours_tangential=gaia)
    if prov.get("kind") != "mock":
        return data, prov
    gen = _family_for({"family": prov["generating_family"],
                       "dataset_options": {"tracer": prov["tracer"]} if "tracer" in prov else {},
                       **({"model": prov["model"]} if "model" in prov else {})})
    if set(gen.names) != set(prov["truth"]):
        raise ValueError("legacy mock provenance does not reconstruct its generating model")
    x_true = np.array([prov["truth"][n] for n in gen.names])
    rng = np.random.default_rng(prov["seed"])
    return FitProblem(gen, data).mock_data(x_true, rng), prov


def _family_for(summary: dict, backend: str | None = None, options: dict | None = None):
    """Restore resolved models; use recorded ladder switches for legacy runs."""
    from .fit import DarkMatterModel, NoDarkMatterModel
    if "model" in summary:
        from .run_io import family_from_config
        return family_from_config(summary["model"], backend=backend)
    fam_label = summary["family"]
    options = summary.get("dataset_options", {}) | (options or {})
    tracer = options.get("tracer", "composite" if "_composite" in fam_label else "trager")
    inferred_backend = next((b for b in ("jam", "agama") if fam_label.endswith("_" + b)), "jeans")
    kw = dict(tracer=tracer, backend=backend or options.get("backend", inferred_backend))
    if options.get("no_scales"):
        kw["instruments"] = ()
    if options.get("isotropic"):
        kw["constant_beta"] = True
        kw["fixed"] = {"beta_0": 0.0}
    elif options.get("constant_beta"):
        kw["constant_beta"] = True
    family = (NoDarkMatterModel(**kw) if fam_label.startswith("K1") else
              DarkMatterModel(gamma=0.0 if "cored" in fam_label else 1.0, **kw))
    if "parameters" in summary and set(family.names) != set(summary["parameters"]):
        raise ValueError("legacy run lacks the model configuration needed for its saved parameters; "
                         "do not reconstruct it from default settings")
    return family


def problem_for(run: dict):
    """Restore a saved problem and check its maximum-likelihood sample."""
    from .fit import FitProblem
    s = run["summary"]
    data, _ = data_for(s, run_dir=run["dir"])
    family = _family_for(s, options=run.get("run", {}).get("dataset_options"))
    problem = FitProblem(family, data)
    x = np.array([s["parameters"][n]["ml"] for n in family.names])
    actual = problem.loglike_vector(x)
    if not np.isclose(actual, s["lnL_max"], rtol=1e-9, atol=1e-7):
        raise ValueError(f"{run['label']}: replay ln L {actual:.12g} differs from saved "
                         f"{s['lnL_max']:.12g}; check model, inputs and code version")
    return problem


def engine_crosscheck(label: str, backends: Iterable[str] = ("jeans", "jam")) -> str:
    """Evaluate a run's best sample with other engines: chi2 per dataset and total ln L.

    Only engines sharing the parametrisation can be compared on the same vector
    (``jeans`` and ``jam``); the AGAMA DF family has its own anisotropy parameters
    and is compared through its own fit instead.
    """
    from .fit import FitProblem

    run = load_run(label)
    s = run["summary"]
    data = problem_for(run).data
    lines = ["| engine | " + " | ".join(f"χ² {d}" for d in s["datasets"]) + " | ln L |", "|---|" + "---|" * (len(s["datasets"]) + 1)]
    for b in backends:
        fam = _family_for(s, b, options=run.get("run", {}).get("dataset_options"))
        x = np.array([s["parameters"][n]["ml"] for n in fam.names])
        P = FitProblem(fam, data)
        chi = P.chi2(x)
        lines.append(f"| {b} | " + " | ".join(f"{chi[d][0]:.1f}" for d in s["datasets"]) + f" | {P.loglike_vector(x):.2f} |")
    return "\n".join(lines)
