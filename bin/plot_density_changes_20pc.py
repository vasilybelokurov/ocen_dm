"""Summarise 20-pc DM shell-density changes in the audited pilot batch.

Read saved diagnostics only. Show tidal runs relative to matched isolation,
and isolation runs relative to their own sampled initial conditions. Bars are
temporal percentiles over the final 10 Myr, not confidence intervals.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from astropy.table import Table
import yaml

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "results/dynamics/batch_analysis_20260921_v2"
OUTPUT = ROOT / "plots/density_changes_20pc"
OUTPUT.parent.mkdir(exist_ok=True)
DATA_OUTPUT = ROOT / "results/plot_data/density_changes_20pc.json"
DATA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
FIELD = "dm_20pc_rho_msun_pc3"
COLOURS = {"cusp": "#2869a4", "core": "#c25e24"}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def statistics(values, times):
    values, times = np.asarray(values, float), np.asarray(times, float)
    assert np.all(np.isfinite(values)) and np.all(np.diff(times) > 0)
    late = times >= times[-1] - 10
    p16, median, p84 = np.percentile(values[late], [16, 50, 84])
    return dict(endpoint_pct=float(values[-1]), last_10_myr_median_pct=float(median),
                last_10_myr_p16_pct=float(p16), last_10_myr_p84_pct=float(p84),
                last_10_myr_n_snapshots=int(late.sum()))


def draw(ax, rows, title, reference, limits):
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=28)
    ax.text(0, 1.025, reference, transform=ax.transAxes, fontsize=9.5, color=".3")
    for y, (label, row) in enumerate(rows):
        colour = COLOURS[row["shape"]]
        stats = row["plotted_change"]
        lo, med, hi = [stats[f"last_10_myr_{q}_pct"] for q in ["p16", "median", "p84"]]
        if row.get("unresolved"):
            ax.axhspan(y - .44, y + .44, color="#fff1dd", zorder=0)
        ax.errorbar(med, y, xerr=[[med-lo], [hi-med]], fmt="o", color=colour,
                    capsize=3, markersize=6, elinewidth=2, zorder=3)
        ax.plot(stats["endpoint_pct"], y, "x", color=colour, ms=7, mew=1.4, zorder=4)
    ax.axvline(0, color=".35", ls="--", lw=1)
    ax.set(yticks=range(len(rows)), yticklabels=[label for label, _ in rows],
           ylim=(len(rows)-.5, -.5), xlim=limits, xlabel="DM shell-density change [%]")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.grid(axis="x", color=".87", linewidth=.7)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(".7")
    ax.set_axisbelow(True)


def main():
    checks = json.loads((ANALYSIS / "artifact_checksums.json").read_text())
    summary_path = ANALYSIS / "summary.json"
    assert digest(summary_path) == checks["summary.json"]
    summary = json.loads(summary_path.read_text())
    inputs = {str(summary_path.relative_to(ROOT)): digest(summary_path)}
    tables, records = {}, {}
    for inventory in summary["inventory"]:
        key = inventory["run"]
        candidates = [p for p in summary["inputs_sha256"] if p.endswith(key + "/diagnostics.ecsv")]
        assert len(candidates) == 1, key
        path = ROOT / candidates[0]
        for source in [path, path.with_name("run.yaml")]:
            relative = str(source.relative_to(ROOT))
            inputs[relative] = digest(source)
            assert inputs[relative] == summary["inputs_sha256"][relative], relative
        run = yaml.safe_load(path.with_name("run.yaml").read_text())
        assert run["status"] == "complete"
        table = Table.read(path)
        times = np.asarray(table["elapsed_myr"], float)
        density = np.asarray(table[FIELD], float)
        assert len(table) == run["n_snapshots"] and density[0] > 0
        assert abs(times[0]) < .001
        shape = "cusp" if "cusp" in key else "core"
        own_change = statistics(100*(density/density[0]-1), times)
        records[key] = dict(run=key, shape=shape, engine=run["engine"],
                            programme=inventory["programme"], orbit=inventory["orbit"],
                            duration_myr=float(times[-1]),
                            initial_density_msun_pc3=float(density[0]),
                            initial_shell_neff=float(table["dm_20pc_shell_neff"][0]),
                            relative_to_initial=own_change, plotted_change=own_change,
                            plotted_reference="own initial shell density",
                            unresolved=key.startswith("progenitor_core/"))
        tables[key] = table
    assert len(records) == 24

    tidal, progenitor = [], []
    for shape in COLOURS:
        for engine in ["live", "frozen"]:
            key = f"remnant_{shape}_passage/evolution_{engine}"
            reference = f"remnant_{shape}_isolation_reference/evolution_{engine}"
            table, control = tables[key], tables[reference]
            time, control_time = np.asarray(table["elapsed_myr"]), np.asarray(control["elapsed_myr"])
            assert time[0] >= control_time[0] - 1e-7 and time[-1] <= control_time[-1] + 1e-7
            assert np.isclose(table[FIELD][0], control[FIELD][0], rtol=1e-12)
            ratio = np.asarray(table[FIELD])/np.interp(time, control_time, control[FIELD])
            stats = statistics(100*(ratio-1), time)
            saved = summary["remnant_comparisons"][f"{shape}_{engine}"]["20"]["rho_msun_pc3"]
            assert all(np.isclose(value, saved[name], atol=1e-10, rtol=1e-10)
                       for name, value in stats.items())
            row = records[key]
            row.update(plotted_change=stats, plotted_reference=reference,
                       relative_to_matched_isolation=stats)
            tidal.append((f"{shape.title()}, {engine}", row))
            row = records[f"progenitor_{shape}/evolution_{engine}"]
            count = round(row["initial_shell_neff"])
            progenitor.append((f"{shape.title()}, {engine}\n" + r"$N_{\rm eff,0}=$" + str(count), row))

    controls = {}
    for shape in COLOURS:
        cases = [
            ("_isolation_reference", "live", "88 Myr: live"),
            ("_isolation_reference", "frozen", "88 Myr: frozen"),
            ("", "live", "20 Myr: baseline live"),
            ("", "frozen", "20 Myr: baseline frozen"),
            ("_dt_half", "live", "20 Myr: half timestep"),
            ("_eps_half", "live_p1_fixed" if shape == "cusp" else "live", "20 Myr: half softening"),
            ("_n2", "live", "20 Myr: double N, live"),
            ("_n2", "frozen", "20 Myr: double N, frozen"),
        ]
        controls[shape] = [(label, records[f"remnant_{shape}{suffix}/evolution_{engine}"])
                           for suffix, engine, label in cases]
    displayed = [row["run"] for _, row in tidal+progenitor+controls["cusp"]+controls["core"]]
    assert len(displayed) == len(set(displayed)) == 24 and set(displayed) == set(records)

    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10,
                         "xtick.labelsize": 9, "ytick.labelsize": 10, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.5),
                             gridspec_kw={"height_ratios": [1, 1.5]})
    fig.subplots_adjust(left=.185, right=.97, top=.84, bottom=.16, wspace=.69, hspace=.62)
    draw(axes[0, 0], tidal, "A  Tidal remnants (88 Myr)",
         "Relative to matched isolation at the same time", (-7, 1))
    draw(axes[0, 1], progenitor, "B  Progenitors in isolation (100 Myr)",
         "Relative to own initial density; core unresolved", (-80, 25))
    for ax, shape, letter in zip(axes[1], COLOURS, ["C", "D"]):
        draw(ax, controls[shape], f"{letter}  Remnants in isolation: {shape}",
             "Relative to own initial density", (-3, 2))
        for y in [1.5, 3.5, 5.5]:
            ax.axhline(y, color=".88", lw=.8)
    fig.suptitle("Dark-matter density changes at 20 pc", x=.52, y=.975,
                 fontsize=19, fontweight="bold")
    fig.text(.52, .939, "All 24 completed pilot evolutions | Density measured in the 18.1-22.1 pc shell",
             ha="center", fontsize=11, color=".25")
    legend = [Line2D([], [], marker="o", color=".3", lw=0, label="Final 10 Myr median"),
              Line2D([], [], color=".3", lw=2, marker="|", markersize=8, label="16-84% temporal range"),
              Line2D([], [], marker="x", color=".3", lw=0, label="Final snapshot")]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(.52, .914),
               ncol=3, frameon=False, fontsize=10)
    fig.text(.025, .095,
             "Negative values indicate depletion; positive values indicate increased shell density. "
             "Each panel states its reference; horizontal scales differ.", fontsize=10)
    fig.text(.025, .070,
             "Bars describe time variation in a single realisation, not statistical confidence. "
             "Doubled-N runs use a different initial particle sample.", fontsize=10)
    fig.text(.025, .045,
             "Shaded progenitor-core rows start with only 15 effective shell particles: "
             "large changes are unresolved. This batch samples one tidal orbit.", fontsize=10)
    paths = [OUTPUT.with_suffix(".png")]
    fig.savefig(paths[0], dpi=200)
    plt.close(fig)
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                  source_sha256=digest(Path(__file__)), inputs_sha256=inputs,
                  field=FIELD, shell_edges_pc=[float(20*np.exp(-.1)), float(20*np.exp(.1))],
                  n_completed_evolutions=len(records),
                  interpretation="Temporal ranges are not confidence intervals. Only one tidal orbit; progenitor core unresolved at 20 pc.",
                  runs=list(records.values()),
                  figures_sha256={str(p.relative_to(ROOT)): digest(p) for p in paths})
    DATA_OUTPUT.write_text(json.dumps(record, indent=2)+"\n")
    for label, row in tidal:
        print(label, row["plotted_change"])
    print("Verified and plotted all", len(records), "completed evolutions:", OUTPUT)


if __name__ == "__main__":
    main()
