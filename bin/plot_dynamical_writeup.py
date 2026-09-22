"""Typeset the audited pilot diagnostics as legible, page-width PNG figures."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT/"results/dynamics/batch_analysis_20260921_v2"
PASSAGE = ROOT/"results/dynamics/remnant_passage_20260921T114939Z"
ISOLATION = ROOT/"results/dynamics/isolation_20260921T113007Z"
OUTPUT = ROOT/"plots"
OUTPUT.mkdir(exist_ok=True)
DATA_OUTPUT = ROOT/"results/plot_data"
DATA_OUTPUT.mkdir(parents=True, exist_ok=True)
summary = json.loads((ANALYSIS/"summary.json").read_text())
profiles = Table.read(ANALYSIS/"endpoint_profiles.ecsv")
inputs = {}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(batch, model, engine="live"):
    path = batch/model/("evolution_"+engine)/"diagnostics.ecsv"
    checksum = digest(path)
    assert checksum == summary["inputs_sha256"][str(path.relative_to(ROOT))]
    inputs[str(path.relative_to(ROOT))] = checksum
    return Table.read(path)


def save(fig, name):
    fig.savefig(OUTPUT/("dynamical_batch_"+name+".png"), dpi=300)
    plt.close(fig)


def events(ax):
    ax.axvline(summary["orbit"]["pericentre_elapsed_myr"], color=".5", ls=":", lw=.8)
    for event in summary["disc_crossings"]:
        ax.axvline(event["elapsed_myr"], color=".7", ls="-.", lw=.7)


plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 9.5, "axes.labelsize": 8.5,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
                     "axes.grid": True, "grid.alpha": .2, "lines.linewidth": 1.2,
                     "pdf.fonttype": 42, "savefig.bbox": "tight"})
colours = {"cusp": "#2869a4", "core": "#c25e24"}
radii = [3, 10, 20, 35, 50, 70]

fig, axes = plt.subplots(2, 2, figsize=(6.65, 4.65), layout="constrained")
for shape, colour in colours.items():
    for engine, style in [("live", "-"), ("frozen", "--")]:
        t = read(PASSAGE, f"remnant_{shape}_passage", engine)
        z = read(PASSAGE, f"remnant_{shape}_isolation_reference", engine)
        for ax, field in zip(axes.flat[:3], ["dm_20pc_mass_msun", "dm_70pc_mass_msun", "dm_bound_mass_msun"]):
            control = np.interp(t["elapsed_myr"], z["elapsed_myr"], z[field])
            ax.plot(t["elapsed_myr"], 100*(t[field]/control-1), style, color=colour, label=f"{shape.title()}, {engine}")
        if engine == "live":
            axes[1, 1].plot(t["elapsed_myr"], t["tidal_scale_pc"], color=colour)
for ax, title in zip(axes.flat, ["Mass inside 20 pc", "Mass inside 70 pc", "Spherical bound-mass proxy", "Instantaneous tidal scale"]):
    ax.set(title=title, xlabel="Elapsed time [Myr]", ylabel="Relative to isolation [%]")
    events(ax)
axes[1, 1].set_ylabel("Tidal scale [pc]")
axes[0, 0].legend(loc="lower left", framealpha=.9)
save(fig, "remnant_passage")

fig = plt.figure(figsize=(6.65, 4.55), layout="constrained")
grid = fig.add_gridspec(2, 2)
axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, :])]
for shape, colour in colours.items():
    for engine, style in [("live", "-"), ("frozen", "--")]:
        result = summary["remnant_comparisons"][shape+"_"+engine]
        for ax, field in zip(axes[:2], ["mass_msun", "rho_msun_pc3"]):
            ax.plot(radii, [result[str(r)][field]["last_10_myr_median_pct"] for r in radii], style, marker=".", color=colour)
            if engine == "live":
                ax.fill_between(radii, [result[str(r)][field]["last_10_myr_p16_pct"] for r in radii], [result[str(r)][field]["last_10_myr_p84_pct"] for r in radii], color=colour, alpha=.12)
        axes[2].plot(radii, [summary["late_anisotropy"][shape+"_"+engine][str(r)] for r in radii], style, marker=".", color=colour, label=f"{shape.title()}, {engine}")
for ax, title in zip(axes, ["Enclosed mass", "Shell density", "DM velocity anisotropy"]):
    ax.set(xscale="log", xlabel="Radius [pc]", title=title)
    ax.axhline(0, color=".6", lw=.6)
for ax in axes[:2]:
    ax.set_ylabel("Relative to isolation [%]")
axes[2].set_ylabel(r"$\beta$")
axes[2].legend(ncol=2, loc="lower left")
save(fig, "remnant_radial_response")

orbit = Table.read(ANALYSIS/"orbit_context.ecsv")
fig, axes = plt.subplots(3, 1, figsize=(6.65, 5.8), sharex=True, layout="constrained")
axes[0].plot(orbit["elapsed_myr"], orbit["radius_kpc"], label="Galactocentric radius")
axes[0].plot(orbit["elapsed_myr"], orbit["z_kpc"], label="Height above disc")
axes[0].set_ylabel("Distance [kpc]")
axes[1].plot(orbit["elapsed_myr"], -orbit["tidal_lambda_min_myr2"], label="Strongest compression")
axes[1].plot(orbit["elapsed_myr"], orbit["tidal_lambda_max_myr2"], label="Strongest stretching")
axes[1].set_ylabel("Host tidal gradient\n"+r"[Myr$^{-2}$]")
t = read(PASSAGE, "remnant_cusp_passage", "frozen")
for r in [20, 70]:
    axes[2].plot(t["elapsed_myr"], t[f"dm_initial_{r}pc_median_delta_reference_energy_kms2"], label=f"Initially in {r} pc shell")
axes[2].set(xlabel="Elapsed time [Myr]", ylabel="Median reference energy change\n"+r"[km$^2$ s$^{-2}$]")
for ax in axes:
    events(ax)
    ax.legend(loc="upper left" if ax == axes[2] else "upper right")
save(fig, "orbit_and_heating")

fig, axes = plt.subplots(2, 2, figsize=(6.65, 4.6), layout="constrained")
for column, shape in enumerate(colours):
    for suffix, engine, label in [("", "live", "Baseline"), ("_dt_half", "live", "Half timestep"), ("_eps_half", "live_p1_fixed" if shape == "cusp" else "live", "Half softening"), ("_n2", "live", "Double N")]:
        batch = ISOLATION if not suffix else PASSAGE
        model = f"remnant_{shape}{suffix}"
        t = read(batch, model, engine)
        axes[0, column].plot(t["elapsed_myr"], 100*(t["dm_20pc_mass_msun"]/t["dm_20pc_mass_msun"][0]-1), label=label)
        path = batch/model/("evolution_"+engine)/"energy.log"
        energy = np.loadtxt(path, usecols=(0, 1))
        assert digest(path) == summary["inputs_sha256"][str(path.relative_to(ROOT))]
        inputs[str(path.relative_to(ROOT))] = digest(path)
        axes[1, column].plot(energy[:, 0]*977.7922216807892, 100*(energy[:, 1]/energy[0, 1]-1))
    axes[0, column].set(title=shape.title(), ylabel="Mass inside 20 pc\nchange from own IC [%]")
    axes[1, column].set(xlabel="Elapsed time [Myr]", ylabel=r"$(E-E_0)/E_0$ [%]")
axes[0, 0].legend(loc="lower left", ncol=2, fontsize=7)
save(fig, "numerical_controls")

fig, axes = plt.subplots(2, 3, figsize=(6.65, 4.05), layout="constrained")
for row, shape in enumerate(colours):
    for engine, style in [("live", "-"), ("frozen", "--")]:
        t = read(ISOLATION, f"progenitor_{shape}", engine)
        for component, colour in [("dm", "#2869a4"), ("stars", "#a55075")]:
            for column, r in enumerate([20, 70]):
                field = f"{component}_{r}pc_mass_msun"
                axes[row, column].plot(t["elapsed_myr"], 100*(t[field]/t[field][0]-1), style, color=colour, label=f"{component.upper() if component == 'dm' else 'Stars'}, {engine}")
    t = read(ISOLATION, f"progenitor_{shape}")
    for component, colour in [("dm", "#2869a4"), ("stars", "#a55075")]:
        axes[row, 2].plot(radii, [np.median(t[f"{component}_{r}pc_shell_neff"]) for r in radii], marker=".", color=colour)
    for column in [0, 1]:
        axes[row, column].set(title=f"{shape.title()}: <{[20, 70][column]} pc", xlabel="Time [Myr]")
    axes[row, 0].set_ylabel("Mass change from IC [%]")
    axes[row, 2].set(title=f"{shape.title()}: shell counts", xlabel="Radius [pc]", ylabel=r"$N_{\rm eff}$", xscale="log", yscale="log")
    axes[row, 2].axhline(100, color=".5", ls=":")
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="outside upper center", ncol=4, fontsize=7.5)
save(fig, "progenitor_isolation")

fig, axes = plt.subplots(2, 2, figsize=(6.65, 4.65), layout="constrained")
for row, shape in enumerate(colours):
    for column, component in enumerate(["dm", "stars"]):
        ax = axes[row, column]
        for engine, style in [("initial", ":"), ("live", "-"), ("frozen", "--")]:
            p = profiles[(profiles["model"] == f"progenitor_{shape}") & (profiles["component"] == component) & (profiles["engine"] == engine) & (profiles["radius_pc"] >= 50)]
            ax.plot(p["radius_pc"], 100*(p["mass_msun"]/p["analytic_mass_msun"]-1), style, marker=".", label=engine.title())
        ax.set(xscale="log", xlabel="Radius [pc]", title=f"{shape.title()}: {'DM' if component == 'dm' else 'stars'}")
        ax.axhline(0, color=".6", lw=.7)
    axes[row, 0].set_ylabel("Mass difference from\nanalytic IC [%]")
axes[0, 0].legend(loc="lower right")
save(fig, "progenitor_global_profiles")

for name in ["summary.json", "endpoint_profiles.ecsv", "orbit_context.ecsv"]:
    path = ANALYSIS/name
    inputs[str(path.relative_to(ROOT))] = digest(path)
record = dict(source_sha256=digest(Path(__file__)), inputs_sha256=inputs,
              figures_sha256={str(p.relative_to(ROOT)): digest(p) for p in sorted(OUTPUT.glob("dynamical_batch_*.png"))})
(DATA_OUTPUT/"dynamical_batch_provenance.json").write_text(json.dumps(record, indent=2)+"\n")
print("Created six paper figures in", OUTPUT)
