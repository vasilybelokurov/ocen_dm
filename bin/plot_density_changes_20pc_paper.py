"""Lay out the verified 20-pc batch summary at the report's text width."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/plot_data/density_changes_20pc.json"
OUTPUT = ROOT / "plots/density_changes_20pc_paper.png"
PROVENANCE = ROOT / "results/plot_data/density_changes_20pc_paper.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUTPUT.parent.mkdir(exist_ok=True)
    PROVENANCE.parent.mkdir(parents=True, exist_ok=True)
    data = json.loads(SOURCE.read_text())
    for name, checksum in data["inputs_sha256"].items():
        assert digest(ROOT / name) == checksum, name
    records = {r["run"]: r for r in data["runs"]}
    panels = [[], [], [], []]
    for column, shape in enumerate(["cusp", "core"]):
        for engine in ["live", "frozen"]:
            panels[0].append((f"{shape.title()}, {engine}", records[f"remnant_{shape}_passage/evolution_{engine}"]))
            row = records[f"progenitor_{shape}/evolution_{engine}"]
            panels[1].append((f"{shape.title()}, {engine}\n" + r"$N_{\rm eff,0}=$" + str(round(row["initial_shell_neff"])), row))
        for suffix, engine, label in [
            ("_isolation_reference", "live", "88 Myr: live"),
            ("_isolation_reference", "frozen", "88 Myr: frozen"),
            ("", "live", "20 Myr: baseline live"),
            ("", "frozen", "20 Myr: baseline frozen"),
            ("_dt_half", "live", "20 Myr: half timestep"),
            ("_eps_half", "live_p1_fixed" if shape == "cusp" else "live", "20 Myr: half softening"),
            ("_n2", "live", "20 Myr: double N, live"),
            ("_n2", "frozen", "20 Myr: double N, frozen"),
        ]:
            panels[column+2].append((label, records[f"remnant_{shape}{suffix}/evolution_{engine}"]))
    assert {r["run"] for panel in panels for _, r in panel} == set(records)
    assert sum(map(len, panels)) == len(records) == 24

    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(6.65, 5.3),
                             gridspec_kw={"height_ratios": [1, 1.6]})
    fig.subplots_adjust(left=.22, right=.985, top=.83, bottom=.09, wspace=.91, hspace=.85)
    titles = ["A  Tidal remnants", "B  Progenitor isolation",
              "C  Remnant isolation: cusp", "D  Remnant isolation: core"]
    references = ["88 Myr; relative to matched isolation", "100 Myr; relative to initial density",
                  "Relative to initial density", "Relative to initial density"]
    colours = {"cusp": "#2869a4", "core": "#c25e24"}
    for index, (ax, panel) in enumerate(zip(axes.flat, panels)):
        ax.set_title(titles[index], loc="left", fontsize=8.5, fontweight="bold", pad=22)
        ax.text(0, 1.025, references[index], transform=ax.transAxes, fontsize=7.2, color=".3")
        for y, (_, row) in enumerate(panel):
            stats = row["plotted_change"]
            lo, med, hi = [stats[f"last_10_myr_{q}_pct"] for q in ["p16", "median", "p84"]]
            colour = colours[row["shape"]]
            if row["unresolved"]:
                ax.axhspan(y-.44, y+.44, color="#fff1dd", zorder=0)
            ax.errorbar(med, y, xerr=[[med-lo], [hi-med]], fmt="o", color=colour,
                        capsize=2, markersize=3.8, elinewidth=1.2, zorder=3)
            ax.plot(stats["endpoint_pct"], y, "x", color=colour, ms=4.7, mew=1, zorder=4)
        limits = (-7, 1) if index == 0 else (-80, 25) if index == 1 else (-3, 2)
        ticks = [-6, -4, -2, 0] if index == 0 else [-80, -40, 0, 20] if index == 1 else [-3, -2, -1, 0, 1, 2]
        ax.set(yticks=range(len(panel)), yticklabels=[label for label, _ in panel],
               ylim=(len(panel)-.5, -.5), xlim=limits, xticks=ticks,
               xlabel="Shell-density change [%]")
        ax.tick_params(axis="y", length=0, pad=4, labelsize=7.4)
        ax.tick_params(axis="x", labelsize=7.5)
        ax.axvline(0, color=".4", ls="--", lw=.8)
        ax.grid(axis="x", color=".87", linewidth=.6)
        if index >= 2:
            for y in [1.5, 3.5, 5.5]:
                ax.axhline(y, color=".88", lw=.6)
        for spine in ["top", "right", "left"]:
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(".7")
        ax.set_axisbelow(True)
    fig.text(.52, .972, "DM shell density at 20 pc: all 24 completed evolutions", ha="center", fontsize=10)
    legend = [Line2D([], [], marker="o", color=".3", lw=0, ms=4, label="Final 10 Myr median"),
              Line2D([], [], color=".3", lw=1.2, label="16-84% temporal range"),
              Line2D([], [], marker="x", color=".3", lw=0, ms=5, label="Final snapshot")]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(.52, .952),
               ncol=3, frameon=False, fontsize=7.5, handlelength=1.4, columnspacing=1.2)
    fig.savefig(OUTPUT, dpi=300)
    plt.close(fig)
    provenance = dict(source_sha256=digest(Path(__file__)),
                      summary_sha256=digest(SOURCE), figure_sha256=digest(OUTPUT),
                      n_evolutions=24, summary_path=str(SOURCE.relative_to(ROOT)))
    PROVENANCE.write_text(json.dumps(provenance, indent=2)+"\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
