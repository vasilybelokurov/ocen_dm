"""Audit and summarise the 2026-09-21 pilot batch without changing saved runs.

Run from the project root with PYTHONPATH=src and the project Python environment.
The output directory must be new. Temporal percentiles describe variation within
one realisation, not confidence intervals or independent repeat experiments.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from scipy.integrate import quad
from scipy.optimize import brentq
import yaml

from ocen_dm.dynamics import TIME_UNIT_MYR, nemo
from ocen_dm.dynamics.diagnostics import effective_number, velocity_moments, tidal_tensor
from ocen_dm.dynamics.experiment import load_prepared
from ocen_dm.dynamics.models import agama_module

ROOT = Path(__file__).resolve().parents[1]
ISOLATION = ROOT / "results/dynamics/isolation_20260921T113007Z"
PASSAGE = ROOT / "results/dynamics/remnant_passage_20260921T114939Z"
RADII = [3, 10, 20, 35, 50, 70]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def describe(values):
    a = np.asarray(values, float)
    finite = a[np.isfinite(a)]
    return dict(initial=a[0], final=a[-1], minimum=np.min(finite) if len(finite) else np.nan,
                maximum=np.max(finite) if len(finite) else np.nan,
                median=np.median(finite) if len(finite) else np.nan,
                n_nonfinite=int(np.sum(~np.isfinite(a))))


def at(table, field, times):
    x = np.asarray(table["elapsed_myr"])
    if min(times) < x[0]-1e-7 or max(times) > x[-1]+1e-7:
        raise ValueError("Comparison would extrapolate beyond a saved evolution")
    return np.interp(times, x, table[field])


def comparison(table, reference, field):
    time = np.asarray(table["elapsed_myr"])
    ratio = np.asarray(table[field])/at(reference, field, time)
    late = time >= time[-1]-10
    pct = 100*(ratio-1)
    return dict(endpoint_pct=pct[-1], last_10_myr_median_pct=np.median(pct[late]),
                last_10_myr_p16_pct=np.percentile(pct[late], 16),
                last_10_myr_p84_pct=np.percentile(pct[late], 84),
                last_10_myr_n_snapshots=int(late.sum()))


def final_snapshot(directory, record, particles, track):
    if record["engine"] == "frozen":
        source = directory/"snapshots.h5"
        with h5py.File(source) as f:
            if f.attrs["snapshots_written"] != len(f["time"]):
                raise ValueError("Incomplete frozen snapshot file")
            time, ids, mass, xv = f["time"][-1], f["particle_id"][:], f["mass"][:], f["xv"][-1]
    else:
        source = directory/"snapshots.nemo"
        # snapprint's time selector has a tolerance. Select near the endpoint,
        # then check the full-precision time carried by the particle rows.
        command = nemo.command_args([
            "snapprint", f"in={source}", f"times={track.time[-1]:.17g}",
            "options=t,key,m,x,y,z,vx,vy,vz", "header=nbody", "newline=t", "format=%.17g"])
        matches = []
        with tempfile.TemporaryFile(mode="w+t") as error:
            with subprocess.Popen(command, cwd=directory, env=nemo.runtime_environment(directory),
                                  stdout=subprocess.PIPE, stderr=error, text=True) as process:
                for item in nemo.parse_snapshots(process.stdout):
                    if abs(item[0]-track.time[-1])*TIME_UNIT_MYR < 1e-5:
                        matches.append(item)
                result = process.wait()
                if result:
                    error.seek(0)
                    raise RuntimeError(error.read())
        if not matches:
            raise ValueError(f"Final snapshot missing: {directory}")
        time, ids, mass, xv = matches[-1]
    assert np.array_equal(ids, particles["particle_id"])
    assert np.allclose(mass, particles["mass"], rtol=2e-6, atol=0)
    assert xv.shape == particles["xv"].shape and np.all(np.isfinite(xv))
    assert abs(time-track.time[-1])*TIME_UNIT_MYR < 1e-5
    return xv-track.at(time), source


def profile(xv, particles, component, radii):
    member = particles["component"] == component
    r = np.linalg.norm(xv[:, :3], axis=1)*1000
    mass = particles["mass"]
    rows = []
    for radius in radii:
        inside = member & (r <= radius)
        shell = member & (r >= radius*np.exp(-.1)) & (r < radius*np.exp(.1))
        sigma, beta = velocity_moments(xv[shell], mass[shell])
        rows.append(dict(radius_pc=radius, mass_msun=mass[inside].sum(),
                         aperture_neff=effective_number(mass[inside]),
                         shell_neff=effective_number(mass[shell]),
                         rho_msun_pc3=mass[shell].sum()/(4*np.pi/3*radius**3*(np.exp(.3)-np.exp(-.3))),
                         sigma_r_kms=sigma, beta=beta))
    return rows


def analytic_mass(component, equilibrium, radius_pc):
    a, cut = component.scale_pc/1000, component.cutoff_pc/1000
    norm = equilibrium["density_norm_msun_kpc3"]
    def integrand(r):
        return r*r*(r/a)**(-component.gamma)*(1+r/a)**(component.gamma-component.outer_slope)*np.exp(-(r/cut)**2)
    return 4*np.pi*norm*quad(integrand, 0, radius_pc/1000, epsabs=1e-15, epsrel=1e-9)[0]


def save_figure(fig, output, name):
    plots = ROOT/"plots"
    plots.mkdir(exist_ok=True)
    path = plots/(output.name+"_"+name+".png")
    if path.exists():
        raise FileExistsError(f"Figure already exists: {path}")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    tables, records, preparations, inventory, checks, excluded = {}, {}, {}, [], {}, []
    for batch in [ISOLATION, PASSAGE]:
        for prepared in sorted(p for p in batch.iterdir() if (p/"config.yaml").exists()):
            config, particles, track, manifest = load_prepared(prepared)
            eq = json.loads((prepared/"equilibrium.json").read_text())
            preparations[prepared.name] = (config, particles, track, eq)
            for name, checksum in manifest["files"].items():
                checks[str((prepared/name).relative_to(ROOT))] = checksum
            checks[str((prepared/"run.yaml").relative_to(ROOT))] = digest(prepared/"run.yaml")
            for directory in sorted(prepared.glob("evolution*")):
                record = yaml.safe_load((directory/"run.yaml").read_text())
                if record["status"] != "complete":
                    excluded.append(dict(path=str(directory.relative_to(ROOT)), status=record["status"], error=record.get("error")))
                    continue
                key = prepared.name+"/"+directory.name
                source = directory/"diagnostics.ecsv"
                checksum = digest(source)
                assert checksum == record["diagnostics_sha256"], source
                assert record["prepared_manifest_sha256"] == digest(prepared/"run.yaml")
                table = Table.read(source)
                time = np.asarray(table["elapsed_myr"])
                assert len(table) == record["n_snapshots"] and np.all(np.diff(time) > 0)
                assert abs(time[0]) < .001 and abs(time[-1]-(track.time[-1]-track.time[0])*TIME_UNIT_MYR) < .001
                for c in config.components:
                    total = np.asarray(table[f"{c.name}_total_mass_msun"])
                    assert np.all(total == total[0])
                    masses = np.array([table[f"{c.name}_{r}pc_mass_msun"] for r in RADII])
                    assert np.all(np.isfinite(masses)) and np.all(np.diff(masses, axis=0) >= 0)
                    assert np.all(masses >= 0) and np.all(masses <= total)
                checks[str(source.relative_to(ROOT))] = checksum
                checks[str((directory/"run.yaml").relative_to(ROOT))] = digest(directory/"run.yaml")
                tables[key], records[key] = table, (directory, record)
                row = dict(run=key, programme=config.programme, engine=record["engine"],
                           orbit=config.orbit.kind, duration_myr=time[-1], snapshots=len(table),
                           columns={name: describe(table[name]) for name in table.colnames})
                energy = directory/"energy.log"
                if energy.exists():
                    raw = np.loadtxt(energy, usecols=(0, 1))
                    checks[str(energy.relative_to(ROOT))] = digest(energy)
                    # The prescribed moving nucleus does work in tidal runs;
                    # only static isolation energy is a conservation diagnostic.
                    if config.orbit.kind == "isolation":
                        row["max_logged_fractional_energy_drift"] = np.max(abs(raw[:, 1]/raw[0, 1]-1))
                inventory.append(row)
                print("Verified", key, flush=True)
    assert len(inventory) == 24 and len(preparations) == 14

    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), n_evolutions=len(inventory),
                   n_preparations=len(preparations), excluded_attempts=excluded, inventory=inventory,
                   orbit=json.loads((PASSAGE/"orbit_selection.json").read_text()),
                   interpretation="Temporal percentiles are descriptive, not confidence intervals. Only one seed; sensitivity runs are in isolation.",
                   remnant_comparisons={}, numerical_controls={}, progenitor_comparisons={})
    checks[str((PASSAGE/"orbit_selection.json").relative_to(ROOT))] = digest(PASSAGE/"orbit_selection.json")
    for shape in ["cusp", "core"]:
        for engine in ["live", "frozen"]:
            t = tables[f"remnant_{shape}_passage/evolution_{engine}"]
            z = tables[f"remnant_{shape}_isolation_reference/evolution_{engine}"]
            summary["remnant_comparisons"][shape+"_"+engine] = {
                str(r): {field: comparison(t, z, f"dm_{r}pc_{field}") for field in ["mass_msun", "rho_msun_pc3"]}
                for r in RADII}
        base = tables[f"remnant_{shape}/evolution_live"]
        for suffix, evolution in [("dt_half", "live"), ("eps_half", "live_p1_fixed" if shape == "cusp" else "live"), ("n2", "live")]:
            t = tables[f"remnant_{shape}_{suffix}/evolution_{evolution}"]
            result = {str(r): comparison(t, base, f"dm_{r}pc_mass_msun") for r in RADII}
            if suffix == "n2":
                # The doubled-N sample is a different realisation: compare it
                # with its own frozen evolution as well as the original run.
                result["own_frozen"] = {str(r): comparison(t, tables[f"remnant_{shape}_n2/evolution_frozen"], f"dm_{r}pc_mass_msun") for r in RADII}
            summary["numerical_controls"][shape+"_"+suffix] = result
        live = tables[f"progenitor_{shape}/evolution_live"]
        frozen = tables[f"progenitor_{shape}/evolution_frozen"]
        summary["progenitor_comparisons"][shape] = {
            c: {str(r): comparison(live, frozen, f"{c}_{r}pc_mass_msun") for r in RADII if np.all(np.asarray(frozen[f"{c}_{r}pc_mass_msun"]) > 0)}
            for c in ["dm", "stars"]}

    # Independently remeasure the final snapshots, extending progenitor profiles
    # to the scales containing most of their mass (absent from the pilot tables).
    profile_rows = []
    for shape in ["cusp", "core"]:
        for model in [f"progenitor_{shape}", f"remnant_{shape}_passage", f"remnant_{shape}_isolation_reference"]:
            config, particles, track, eq = preparations[model]
            radii = RADII + ([100, 200, 500, 1000, 2000, 5000, 10000] if config.programme == "progenitor" else [100, 150, 200])
            for engine in ["initial", "live", "frozen"]:
                if engine == "initial":
                    xv = particles["xv"]
                else:
                    key = model+"/evolution_"+engine
                    directory, record = records[key]
                    xv, source = final_snapshot(directory, record, particles, track)
                    checks[str(source.relative_to(ROOT))] = digest(source)
                for n, component in enumerate(config.components):
                    for row in profile(xv, particles, n, radii):
                        row.update(model=model, engine=engine, component=component.name,
                                   analytic_mass_msun=analytic_mass(component, eq["components"][component.name], row["radius_pc"]))
                        if engine != "initial" and row["radius_pc"] in RADII:
                            expected = tables[key][f'{component.name}_{row["radius_pc"]}pc_mass_msun'][-1]
                            assert np.isclose(row["mass_msun"], expected, rtol=1e-12, atol=1e-6)
                        profile_rows.append(row)
                print("Remeasured", model, engine, flush=True)
    profile_table = Table(rows=profile_rows)
    profile_table.meta["note"] = "Final snapshot (or initial IC); exact aperture sums; shell width +/-0.1 in ln r."
    profile_table.write(output/"endpoint_profiles.ecsv")

    # Verify that all nominally matched controls use identical phase space and masses.
    matched = []
    for shape in ["cusp", "core"]:
        base = preparations[f"remnant_{shape}"][1]
        for suffix in ["dt_half", "eps_half", "passage", "isolation_reference"]:
            other = preparations[f"remnant_{shape}_{suffix}"][1]
            for field in ["xv", "mass", "particle_id", "component"]:
                assert np.array_equal(base[field], other[field]), (shape, suffix, field)
            matched.append(f"remnant_{shape}_{suffix}")
    for asset in ["orbit.npz", "host.ini"]:
        assert digest(PASSAGE/"remnant_cusp_passage"/asset) == digest(PASSAGE/"remnant_core_passage"/asset)
    track = preparations["remnant_cusp_passage"][2]
    host = agama_module().Potential(str(PASSAGE/"remnant_cusp_passage/host.ini"))
    crossings = []
    for i in np.where(track.xv[:-1, 2]*track.xv[1:, 2] < 0)[0]:
        time = brentq(lambda tm: track.at(tm)[2], track.time[i], track.time[i+1], xtol=1e-14)
        xv = track.at(time)
        crossings.append(dict(elapsed_myr=(time-track.time[0])*TIME_UNIT_MYR,
                              cylindrical_radius_kpc=np.linalg.norm(xv[:2]), vz_kms=xv[5],
                              tidal_tensor_eigenvalues_kms2_kpc2=np.linalg.eigvalsh(tidal_tensor(host, xv[:3], time))))
    summary["disc_crossings"] = crossings
    summary["late_anisotropy"] = {}
    for shape in ["cusp", "core"]:
        for engine in ["live", "frozen"]:
            table = tables[f"remnant_{shape}_passage/evolution_{engine}"]
            late = np.asarray(table["elapsed_myr"]) >= table["elapsed_myr"][-1]-10
            summary["late_anisotropy"][shape+"_"+engine] = {
                str(r): np.nanmedian(table[f"dm_{r}pc_beta"][late]) for r in RADII}
    summary["matched_initial_conditions_verified"] = matched
    summary["source_sha256"] = digest(Path(__file__))
    summary["inputs_sha256"] = checks
    (output/"summary.json").write_text(json.dumps(clean(summary), indent=2, allow_nan=False)+"\n")
    (output/"analysis_source.py").write_bytes(Path(__file__).read_bytes())

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .2})
    colours = {"cusp": "#2869a4", "core": "#c25e24"}
    peri = summary["orbit"]["pericentre_elapsed_myr"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for shape, color in colours.items():
        for engine, style in [("live", "-"), ("frozen", "--")]:
            t = tables[f"remnant_{shape}_passage/evolution_{engine}"]
            z = tables[f"remnant_{shape}_isolation_reference/evolution_{engine}"]
            time = np.asarray(t["elapsed_myr"])
            for ax, field in zip(axes.flat[:3], ["dm_20pc_mass_msun", "dm_70pc_mass_msun", "dm_bound_mass_msun"]):
                ax.plot(time, 100*(np.asarray(t[field])/at(z, field, time)-1), style, color=color, label=f"{shape.title()}, {engine}")
            if engine == "live":
                axes[1, 1].plot(time, t["tidal_scale_pc"], color=color, label=shape.title())
    for ax, title in zip(axes.flat, ["DM mass inside 20 pc", "DM mass inside 70 pc", "Spherical bound-mass proxy", "Instantaneous tidal scale"]):
        ax.axvline(peri, color=".5", ls=":", lw=1)
        for crossing in crossings:
            ax.axvline(crossing["elapsed_myr"], color=".75", ls="-.", lw=.8)
        ax.set(title=title, xlabel="Elapsed time [Myr]", ylabel="Change relative to isolation [%]")
    axes[1, 1].set_ylabel("Tidal scale [pc]")
    axes[0, 0].legend(fontsize=9)
    fig.suptitle("One Galactic passage: compact remnants\nDotted line: pericentre; dot-dashed lines: disc crossings", fontsize=12)
    save_figure(fig, output, "remnant_passage")

    sample_time = track.time[::10]
    phase = track.at(sample_time)
    elapsed = (sample_time-track.time[0])*TIME_UNIT_MYR
    eigenvalues = np.array([np.linalg.eigvalsh(tidal_tensor(host, xv[:3], tm)) for tm, xv in zip(sample_time, phase)]) / TIME_UNIT_MYR**2
    orbit_table = Table(dict(elapsed_myr=elapsed, radius_kpc=np.linalg.norm(phase[:, :3], axis=1),
                            z_kpc=phase[:, 2], tidal_lambda_min_myr2=eigenvalues[:, 0],
                            tidal_lambda_max_myr2=eigenvalues[:, -1]))
    orbit_table.meta["tensor"] = "Inertial host force gradient; no centrifugal term."
    orbit_table.write(output/"orbit_context.ecsv")
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True, layout="constrained")
    axes[0].plot(elapsed, orbit_table["radius_kpc"], label="Galactocentric radius")
    axes[0].plot(elapsed, phase[:, 2], label="Height above disc")
    axes[0].set_ylabel("Distance [kpc]")
    axes[1].plot(elapsed, -eigenvalues[:, 0], label="Strongest compression")
    axes[1].plot(elapsed, eigenvalues[:, -1], label="Strongest stretching")
    axes[1].set_ylabel(r"Host tidal gradient [Myr$^{-2}$]")
    for r in [20, 70]:
        table = tables["remnant_cusp_passage/evolution_frozen"]
        axes[2].plot(table["elapsed_myr"], table[f"dm_initial_{r}pc_median_delta_reference_energy_kms2"], label=f"Initially in {r} pc shell")
    axes[2].set(xlabel="Elapsed time [Myr]", ylabel=r"Median $\Delta E_{\rm initial\ potential}$ [km$^2$ s$^{-2}$]")
    for ax in axes:
        ax.axvline(peri, color=".5", ls=":", lw=1)
        for crossing in crossings:
            ax.axvline(crossing["elapsed_myr"], color=".7", ls="-.", lw=1)
        ax.legend(fontsize=9)
    fig.suptitle("Orbit and tidal forcing; energy response shown for the frozen cusp\nDotted: pericentre; dot-dashed: disc crossings", fontsize=12)
    save_figure(fig, output, "orbit_and_heating")

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), layout="constrained")
    for shape, color in colours.items():
        for engine, style in [("live", "-"), ("frozen", "--")]:
            result = summary["remnant_comparisons"][shape+"_"+engine]
            for ax, field in zip(axes[:2], ["mass_msun", "rho_msun_pc3"]):
                med = [result[str(r)][field]["last_10_myr_median_pct"] for r in RADII]
                ax.plot(RADII, med, style, marker="o", color=color, label=f"{shape.title()}, {engine}")
                if engine == "live":
                    ax.fill_between(RADII, [result[str(r)][field]["last_10_myr_p16_pct"] for r in RADII], [result[str(r)][field]["last_10_myr_p84_pct"] for r in RADII], color=color, alpha=.12)
            t = tables[f"remnant_{shape}_passage/evolution_{engine}"]
            late = np.asarray(t["elapsed_myr"]) >= t["elapsed_myr"][-1]-10
            axes[2].plot(RADII, [np.nanmedian(t[f"dm_{r}pc_beta"][late]) for r in RADII], style, marker="o", color=color)
    for ax, title, ylabel in zip(axes, ["Enclosed mass", "Shell density", "DM velocity anisotropy"], ["Change relative to isolation [%]", "Change relative to isolation [%]", r"$\beta$"]):
        ax.set(xscale="log", xlabel="Radius [pc]", title=title, ylabel=ylabel)
        ax.axhline(0, color=".6", lw=.8)
    axes[0].legend(fontsize=8)
    fig.suptitle("Final 10 Myr medians; shading is temporal 16–84% variation, not an uncertainty", fontsize=11)
    save_figure(fig, output, "remnant_radial_response")

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    for column, shape in enumerate(colours):
        for suffix, evolution, label in [("", "live", "Baseline: 100k, 0.2 pc"), ("_dt_half", "live", "Half timestep"), ("_eps_half", "live_p1_fixed" if shape == "cusp" else "live", "Half softening"), ("_n2", "live", "200k particles")]:
            key = f"remnant_{shape}{suffix}/evolution_{evolution}"
            t = tables[key]
            axes[0, column].plot(t["elapsed_myr"], 100*(t["dm_20pc_mass_msun"]/t["dm_20pc_mass_msun"][0]-1), label=label)
            directory, _ = records[key]
            energy = np.loadtxt(directory/"energy.log", usecols=(0, 1))
            axes[1, column].plot(energy[:, 0]*TIME_UNIT_MYR, 100*(energy[:, 1]/energy[0, 1]-1))
        axes[0, column].set(title=shape.title(), ylabel=r"$M_{\rm DM}(<20\,pc)$ change from own IC [%]")
        axes[1, column].set(ylabel=r"$(E-E_0)/E_0$ [%]", xlabel="Elapsed time [Myr]")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("20 Myr isolation sensitivity checks; doubled N uses a different particle sample", fontsize=12)
    save_figure(fig, output, "numerical_controls")

    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained")
    for row, shape in enumerate(colours):
        for engine, style in [("live", "-"), ("frozen", "--")]:
            t = tables[f"progenitor_{shape}/evolution_{engine}"]
            for component, color in [("dm", "#2869a4"), ("stars", "#a55075")]:
                for column, r in enumerate([20, 70]):
                    field = f"{component}_{r}pc_mass_msun"
                    axes[row, column].plot(t["elapsed_myr"], 100*(t[field]/t[field][0]-1), style, color=color, label=f"{component}, {engine}")
        for component, color in [("dm", "#2869a4"), ("stars", "#a55075")]:
            t = tables[f"progenitor_{shape}/evolution_live"]
            axes[row, 2].plot(RADII, [np.median(t[f"{component}_{r}pc_shell_neff"]) for r in RADII], marker="o", color=color)
        for column in [0, 1]:
            axes[row, column].set(title=f"{shape.title()}: mass inside {[20, 70][column]} pc", xlabel="Elapsed time [Myr]", ylabel="Change from initial sample [%]")
        axes[row, 2].set(title=f"{shape.title()}: median shell counts", xlabel="Radius [pc]", ylabel=r"$N_{\rm eff}$", xscale="log", yscale="log")
        axes[row, 2].axhline(100, color=".5", ls=":")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("100 Myr progenitor isolation: fluctuations must be read with particle counts", fontsize=12)
    save_figure(fig, output, "progenitor_isolation")

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
    for row, shape in enumerate(colours):
        for column, component in enumerate(["dm", "stars"]):
            ax = axes[row, column]
            for engine, style in [("initial", ":"), ("live", "-"), ("frozen", "--")]:
                p = profile_table[(profile_table["model"] == f"progenitor_{shape}") & (profile_table["component"] == component) & (profile_table["engine"] == engine) & (profile_table["radius_pc"] >= 50)]
                ax.plot(p["radius_pc"], 100*(p["mass_msun"]/p["analytic_mass_msun"]-1), style, marker=".", label=engine)
            ax.set(xscale="log", xlabel="Radius [pc]", ylabel="Mass difference from analytic IC [%]", title=f"{shape.title()}: {component}")
            ax.axhline(0, color=".6", lw=.8)
    axes[0, 0].legend(fontsize=9)
    fig.suptitle("Progenitor mass profiles: initial sample and 100 Myr endpoints", fontsize=12)
    save_figure(fig, output, "progenitor_global_profiles")
    (output/"artifact_checksums.json").write_text(json.dumps({p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file()}, indent=2)+"\n")
    print("Completed analysis:", output, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output.resolve())
