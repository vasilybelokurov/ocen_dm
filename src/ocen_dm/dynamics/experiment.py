"""Prepare once, evolve with live gravity or a matched frozen potential."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tarfile

import h5py
import numpy as np
import yaml

from ..kinematics.run_io import code_state
from . import G, TIME_UNIT_MYR
from .config import Experiment, safe_label
from .diagnostics import snapshot_diagnostics, weighted_median
from .models import agama_module, build_equilibrium
from .orbits import Track, build_track, write_external_potential, timestep
from . import nemo


def _now():
    return datetime.now(timezone.utc).isoformat()


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_manifest(path, record):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(yaml.safe_dump(record, sort_keys=False))
    temporary.replace(path)


def _source_archive(directory, name="source.tar.gz"):
    root = Path(__file__).resolve().parents[3]
    with tarfile.open(directory/name, "w:gz") as archive:
        for source in sorted((root/"src").rglob("*.py")):
            archive.add(source, arcname=str(source.relative_to(root)))
        archive.add(root/"pyproject.toml", arcname="pyproject.toml")
        for name in ("nemo_env.sh", "build_nbody_keys.sh", "build_nbody_kernel.py", "ocen_keep_keys.cc",
                     "nbody_compat/mm_malloc.h", "nbody_compat/p1_coefficients.h"):
            archive.add(root/"bin"/name, arcname="bin/"+name)


def prepare(config, directory, stellar_data=None):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    record = {"schema_version": 1, "status": "preparing", "started_utc": _now(),
              "config": config.to_dict(), "code": code_state(),
              "units": {"length": "kpc", "velocity": "km/s", "mass": "Msun",
                        "time_unit_myr": TIME_UNIT_MYR, "G": G}}
    path = directory/"run.yaml"
    _write_manifest(path, record)
    try:
        _source_archive(directory)
        if stellar_data is not None:
            from ..kinematics.run_io import write_data_snapshot
            record['stellar_observations'] = write_data_snapshot(stellar_data, directory)
        (directory/"config.yaml").write_text(yaml.safe_dump(config.to_dict(), sort_keys=False))
        particles, equilibrium = build_equilibrium(config, directory)
        (directory/"equilibrium.json").write_text(json.dumps(equilibrium, indent=2)+"\n")
        track, orbit = build_track(config, directory)
        for frozen in (False, True):
            write_external_potential(directory, frozen, config=config)
        record.update(status="prepared", prepared_utc=_now(), orbit=orbit,
                      n_particles=len(particles["mass"]),
                      files={str(p.relative_to(directory)): _hash(p)
                             for p in sorted(directory.rglob("*")) if p.is_file() and p != path})
        _write_manifest(path, record)
    except BaseException as exc:
        record.update(status="failed", failed_utc=_now(), error=f"{type(exc).__name__}: {exc}")
        _write_manifest(path, record)
        raise
    return directory


def load_prepared(directory):
    directory = Path(directory).resolve()
    record = yaml.safe_load((directory/"run.yaml").read_text())
    if record.get("schema_version") != 1 or record.get("status") != "prepared":
        raise ValueError("Initial conditions have not completed preparation")
    for name, digest in record["files"].items():
        path = directory/name
        if not path.is_file() or _hash(path) != digest:
            raise ValueError(f"Prepared input missing or changed: {name}")
    config = Experiment.from_dict(record["config"])
    if Experiment.load(directory/"config.yaml") != config:
        raise ValueError("Prepared manifest and saved configuration disagree")
    with np.load(directory/"initial_local.npz", allow_pickle=False) as data:
        particles = {name: data[name] for name in data.files}
    return config, particles, Track.load(directory), record


def _run_frozen(config, particles, track, prepared, evolution):
    agama = agama_module()
    potential = agama.Potential(str(prepared/"external_frozen.ini"))
    times = np.arange(track.time[0], track.time[-1], config.integrator.output_step_myr/TIME_UNIT_MYR)
    times = np.append(times, track.time[-1])
    initial = particles["xv"]+track.xv[0]
    with h5py.File(evolution/"snapshots.h5", "x") as output:
        output.attrs.update(time_unit_myr=TIME_UNIT_MYR, engine="frozen",
                            mass_role="diagnostic weights; no particle self-gravity")
        output.create_dataset("time", data=times)
        output.create_dataset("particle_id", data=particles["particle_id"])
        output.create_dataset("mass", data=particles["mass"])
        positions = output.create_dataset("xv", shape=(len(times), len(initial), 6), dtype="f8",
                                           chunks=(1, min(4096, len(initial)), 6), compression="gzip")
        positions[0] = initial
        output.attrs["snapshots_written"] = 1
        output.flush()
        current = initial
        for i in range(1, len(times)):
            with agama.setNumThreads(1):
                results = agama.orbit(ic=current, potential=potential, timestart=times[i-1],
                                      time=times[i]-times[i-1], trajsize=1,
                                      accuracy=config.integrator.accuracy, dtype="float64", verbose=False)
            current = np.vstack(results[:, 1])
            if not np.all(np.isfinite(current)):
                raise ValueError("Frozen integration produced non-finite phase-space coordinates")
            positions[i] = current
            output.attrs["snapshots_written"] = i+1
            output.flush()
            print(f"frozen: snapshot {i+1}/{len(times)} at {times[i]*TIME_UNIT_MYR:.6g} Myr", flush=True)


def run(prepared, engine="live", label=None, env_script=None):
    if engine not in ("live", "frozen"):
        raise ValueError("engine must be live or frozen")
    prepared = Path(prepared).resolve()
    config, particles, track, inputs = load_prepared(prepared)
    label = safe_label(label or engine)
    directory = prepared/("evolution_"+label)
    directory.mkdir(exist_ok=False)
    kmax, dt = timestep(config)
    record = {"schema_version": 1, "status": "running", "started_utc": _now(),
              "engine": engine, "prepared_manifest_sha256": _hash(prepared/"run.yaml"),
              "code": code_state(), "frame": "host inertial Cartesian",
              "friction": False,
              "nucleus_motion": "self-gravitating particles" if config.nucleus.live and engine == "live" else "prescribed",
              "nucleus_representation": "live Plummer particles" if config.nucleus.live else "rigid Plummer",
              "physical_stellar_relaxation": False}
    if engine == "live":
        record.update(integrator="gyrfalcON block steps", kmax=kmax,
                      actual_max_step_myr=dt*TIME_UNIT_MYR,
                      actual_min_step_myr=dt*TIME_UNIT_MYR/2**(config.integrator.levels-1))
    else:
        record.update(integrator="AGAMA adaptive DOP853", accuracy=config.integrator.accuracy)
    path = directory/"run.yaml"
    _write_manifest(path, record)
    try:
        _source_archive(directory)
        if engine == "live":
            record["kernel_runtime"] = nemo.build_safe_kernel(directory, env_script)
            record["executables"] = nemo.executable_info(directory, env_script)
            plugin = nemo.build_key_plugin(directory, env_script)
            record["key_plugin_sha256"] = _hash(plugin)
            nemo.write_input(particles, track, directory, env_script)
            arguments = nemo.gyrfalcon_args(config, prepared)
            record["command"] = arguments
            _write_manifest(path, record)
            nemo.call(arguments, directory, "gyrfalcon.log", env_script)
        else:
            _run_frozen(config, particles, track, prepared, directory)
        record.update(status="evolved", evolved_utc=_now())
        _write_manifest(path, record)
        analyse(directory, env_script=env_script)
    except BaseException as exc:
        status = "analysis_failed" if record["status"] == "evolved" else "failed"
        record.update(status=status, failed_utc=_now(), error=f"{type(exc).__name__}: {exc}")
        _write_manifest(path, record)
        raise
    return directory


def _frozen_snapshots(path):
    with h5py.File(path) as data:
        count = data.attrs["snapshots_written"]
        if count != len(data["time"]):
            raise ValueError("Frozen evolution was interrupted before all snapshots were written")
        for i in range(count):
            yield data["time"][i], data["particle_id"][:], data["mass"][:], data["xv"][i]


def analyse(directory, env_script=None):
    from astropy.table import Table
    directory = Path(directory).resolve()
    prepared = directory.parent
    config, particles, track, inputs = load_prepared(prepared)
    path = directory/"run.yaml"
    record = yaml.safe_load(path.read_text())
    if record.get("status") not in ("evolved", "analysis_failed", "complete"):
        raise ValueError("Only a finished evolution can be analysed")
    if record["prepared_manifest_sha256"] != _hash(prepared/"run.yaml"):
        raise ValueError("Prepared manifest changed since evolution")
    agama = agama_module()
    host = agama.Potential(str(prepared/"host.ini")) if (prepared/"host.ini").exists() else None
    reference = agama.Potential(str(prepared/"satellite_initial.ini"))
    initial_radius = np.linalg.norm(particles["xv"][:, :3], axis=1)
    initial_energy = reference.potential(particles["xv"][:, :3])+.5*np.sum(particles["xv"][:, 3:]**2, axis=1)
    iterator = (_frozen_snapshots(directory/"snapshots.h5") if record["engine"] == "frozen" else
                nemo.snapshots(directory/"snapshots.nemo", directory, env_script))
    rows = []
    previous_time = -np.inf
    for time, ids, mass, xv in iterator:
        if not np.array_equal(ids, particles["particle_id"]):
            raise ValueError("Evolution lost, duplicated or changed particle identities")
        if not np.allclose(mass, particles["mass"], rtol=2e-6, atol=0):
            raise ValueError("Evolution changed particle masses")
        if xv.shape != particles["xv"].shape or not np.all(np.isfinite(xv)):
            raise ValueError("Invalid snapshot phase space")
        if time < previous_time-1e-10:
            raise ValueError("Snapshot times decreased")
        if time <= previous_time+1e-12:
            continue  # Some NEMO versions repeat the final output.
        centre = track.at(time)
        row = snapshot_diagnostics(xv, particles, config, centre, time, host)
        row.update(time_myr=float(time*TIME_UNIT_MYR), elapsed_myr=float((time-track.time[0])*TIME_UNIT_MYR))
        measured_centre = np.array([row[f"centre_{name}"] for name in
                                   ["x_kpc", "y_kpc", "z_kpc", "vx_kms", "vy_kms", "vz_kms"]])
        local = xv-measured_centre
        energy = reference.potential(local[:, :3])+.5*np.sum(local[:, 3:]**2, axis=1)
        row["median_abs_fractional_reference_energy_change"] = weighted_median(
            np.abs((energy-initial_energy)/initial_energy), particles["mass"])
        for number, component in enumerate(config.particle_components):
            for radius in config.radii_pc:
                shell = ((initial_radius >= radius/1000*np.exp(-.1)) &
                         (initial_radius < radius/1000*np.exp(.1)) & (particles["component"] == number))
                # Fixed initial-potential energy, separated from changing live binding energy.
                row[f"{component.name}_initial_{radius:g}pc_median_delta_reference_energy_kms2"] = (
                    weighted_median(energy[shell]-initial_energy[shell], particles["mass"][shell]))
        rows.append(row)
        previous_time = time
    if len(rows) < 2 or abs(rows[0]["time_myr"]-track.time[0]*TIME_UNIT_MYR) > .001:
        raise ValueError("Evolution does not contain both initial and later snapshots")
    if abs(rows[-1]["time_myr"]-track.time[-1]*TIME_UNIT_MYR) > .001:
        raise ValueError("Evolution did not reach its requested end time")
    table = Table(rows=rows)
    table.meta.update(programme=config.programme, engine=record["engine"],
                      binding="iterative self-excluded spherical satellite energy about measured stellar centre" if config.nucleus.live else "iterative self-excluded spherical satellite energy about rigid nucleus",
                      tidal_scale="largest eigenvalue of host force gradient plus instantaneous centrifugal matrix",
                      density="finite shell r*exp(-0.1) to r*exp(0.1), not a pointwise estimate",
                      energy_change="mass-weighted specific energy in the fixed initial satellite potential",
                      interpretation="collisionless live nucleus in imposed host; no physical stellar relaxation or host wake" if config.nucleus.live else "controlled prescribed-orbit experiment; no dynamical friction",
                      projected_radii="half-mass radii along three Cartesian axes; half-light only for constant M/L")
    table.write(directory/"diagnostics.ecsv", overwrite=True)
    if config.nucleus.live and 'stellar_observations' in inputs:
        from ..kinematics.run_io import read_data_snapshot
        from .lifetime import projected_stellar_profiles
        data = read_data_snapshot(prepared, inputs['stellar_observations'])
        observable_rows = []
        for epoch, coordinates, reference_centre in [
                ('initial', particles['xv']+track.xv[0], track.xv[0]),
                ('final', xv, track.xv[-1])]:
            for row in projected_stellar_profiles(coordinates, particles, config, reference_centre, data):
                row['epoch'] = epoch
                observable_rows.append(row)
        observables = Table(rows=observable_rows)
        observables.meta.update(distance_kpc=5.34, minimum_neff=30,
            interpretation='Diagnostic Cartesian views, complete annuli, constant stellar M/L; not an observational survival classification',
            missing_physics='No rotation, physical stellar relaxation, mass evolution or measured spatial selection')
        observables.write(directory/'stellar_observables.ecsv', overwrite=True)
        record['stellar_observables_sha256'] = _hash(directory/'stellar_observables.ecsv')
    _source_archive(directory, "analysis_source.tar.gz")
    record.update(status="complete", completed_utc=_now(), n_snapshots=len(rows), analysis_code=code_state(),
                  diagnostics_sha256=_hash(directory/"diagnostics.ecsv"))
    _write_manifest(path, record)
    return table
