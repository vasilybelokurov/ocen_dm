"""Summaries and figures for finished nested-sampling runs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.table import Table

from ..paths import results_dir

__all__ = ["load_run", "comparison_table", "posterior_upper_limit", "write_report"]

_LOG_PARAMS = ("M_star", "M_rem", "a_rem", "M_bh", "r_beta", "M_dm_100", "r_s")


def load_run(label: str) -> dict:
    """Read ``summary.json``, the posterior table and the profile stack of one run."""
    d = results_dir() / "fits" / label
    if not (d / "summary.json").exists():
        raise FileNotFoundError(f"run {label!r} has no summary.json under {d}")
    out = {"label": label, "dir": d, "summary": json.loads((d / "summary.json").read_text()),
           "posterior": Table.read(d / "posterior.ecsv")}
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
    names: list[str] = []
    for r in runs:
        for n in r["summary"]["parameters"]:
            if n not in names:
                names.append(n)
    ref = max(r["summary"]["logz"] for r in runs)
    lines = ["| quantity | " + " | ".join(r["label"] for r in runs) + " |",
             "|---|" + "---|" * len(runs)]
    lines.append("| ln Z | " + " | ".join(f"{r['summary']['logz']:.2f} ± {r['summary']['logzerr']:.2f}" for r in runs) + " |")
    lines.append("| Δ ln Z vs best | " + " | ".join(f"{r['summary']['logz'] - ref:+.2f}" for r in runs) + " |")
    lines.append("| max ln L | " + " | ".join(f"{r['summary']['lnL_max']:.1f}" for r in runs) + " |")
    lines.append("| χ² at max L / N | " + " | ".join(f"{r['summary']['chi2_ml_total']:.0f} / {r['summary']['n_points']}" for r in runs) + " |")
    for ds in runs[0]["summary"]["chi2_ml"]:
        lines.append(f"| χ² {ds} | " + " | ".join(f"{r['summary']['chi2_ml'].get(ds, {}).get('chi2', float('nan')):.0f} / {r['summary']['chi2_ml'].get(ds, {}).get('n', 0)}" for r in runs) + " |")
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
    from .fit import DarkMatterModel, FitProblem, NoDarkMatterModel
    from .likelihood import KinematicData

    runs = [load_run(l) for l in labels]
    path = path or results_dir() / "fits" / "comparison.md"
    plots_dir = plots_dir or (results_dir().parent / "plots")
    md = ["# Nested-sampling runs: comparison\n", comparison_table(runs), ""]
    fig = plot_posterior_profiles({r["label"]: r["dir"] for r in runs if "profiles" in r},
                                  plots_dir / "fit_posterior_profiles.png", title="Enclosed mass, dark fraction, circular speed (16-84 %)")
    md.append(f"Figure: `{fig.relative_to(plots_dir.parent) if fig.is_relative_to(plots_dir.parent) else fig}`\n")
    for r in runs:
        s = r["summary"]
        datasets = list(s["datasets"])
        data = KinematicData.load(datasets, gaia_edr3_pm={"r_min_arcsec": 300.0})
        tracer = "composite" if s["family"].endswith("_composite") else "trager"
        if s["family"].startswith("K1"):
            fam = NoDarkMatterModel(tracer=tracer)
        else:
            fam = DarkMatterModel(gamma=0.0 if "cored" in s["family"] else 1.0, tracer=tracer)
        P = FitProblem(fam, data)
        x_ml = np.array([s["parameters"][n]["ml"] for n in fam.names])
        samples = np.array([np.asarray(r["posterior"][n]) for n in fam.names]).T
        out = plot_profile_fit(P, x_ml, plots_dir / f"fit_{r['label']}_posterior_profiles.png",
                               title=f"{r['label']}: best sample and 16-84 % posterior band", samples=samples)
        md.append(f"Figure: `{out.name}`")
    path.write_text("\n".join(md))
    return path
