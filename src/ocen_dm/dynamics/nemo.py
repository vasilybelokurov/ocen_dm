"""NEMO processes and identity-preserving snapshot streams."""
from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from ..paths import project_root
from . import G, TIME_UNIT_MYR
from .orbits import end_time, timestep


def command_args(arguments, env_script=None):
    script = Path(env_script or os.environ.get("OCEN_NEMO_ENV", project_root()/"bin/nemo_env.sh"))
    if not script.is_file():
        raise FileNotFoundError(f"NEMO environment script not found: {script}")
    # Paths and arguments are positional shell parameters, never shell source.
    shell = ('source "$1" && shift && '
             'if [ -n "${OCEN_NBODY_FALCON_LIBRARY:-}" ]; then '
             '_ocen_libdir="${OCEN_NBODY_FALCON_LIBRARY%/*}"; '
             'export DYLD_LIBRARY_PATH="$_ocen_libdir:${DYLD_LIBRARY_PATH:-}"; '
             'export LD_LIBRARY_PATH="$_ocen_libdir:${LD_LIBRARY_PATH:-}"; fi && exec "$@"')
    return ["/bin/zsh", "-c", shell, "ocen-nbody",
            str(script.resolve()), *map(str, arguments)]


def runtime_environment(directory):
    environment = os.environ.copy()
    scratch = directory/"runtime"
    scratch.mkdir(exist_ok=True)
    environment.update(OCEN_NBODY_SCRATCH=str(scratch.resolve()), OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    library = scratch/"libfalcON.so"
    manifest = scratch/"kernel_build.json"
    environment.pop("OCEN_NBODY_FALCON_LIBRARY", None)
    if manifest.exists() or library.exists():
        if not manifest.is_file() or not library.is_file():
            raise ValueError("Run-local falcON library or build manifest missing")
        record = json.loads(manifest.read_text())
        if hashlib.sha256(library.read_bytes()).hexdigest() != record["library_sha256"]:
            raise ValueError("Run-local falcON library missing or checksum mismatch")
    if library.is_file():
        environment["OCEN_NBODY_FALCON_LIBRARY"] = str(library.resolve())
    return environment


def call(arguments, directory, log_name, env_script=None):
    command = command_args(arguments, env_script)
    (directory/(log_name+".command.json")).write_text(json.dumps(command, indent=2)+"\n")
    with (directory/log_name).open("xb") as log:
        subprocess.run(command, cwd=directory, env=runtime_environment(directory),
                       stdout=log, stderr=subprocess.STDOUT, check=True)


def write_input(particles, track, directory, env_script=None):
    xv = particles["xv"] + track.xv[0]
    # NEMO keys preserve component labels even when the tree reorders particles.
    table = np.column_stack((particles["mass"], xv, particles["softening"], particles["particle_id"]))
    np.savetxt(directory/"initial.txt", table, fmt=["%.17g"]*8+["%d"])
    call(["a2s", "in=initial.txt", "out=initial.nemo", f"N={len(xv)}", "read=mxvek",
          f"time={track.time[0]:.17g}"], directory, "a2s.log", env_script)


def build_safe_kernel(directory, env_script=None):
    """Replace only the unsafe scalar P1 arithmetic in a run-local library."""
    runtime_environment(directory)
    manifest = directory/"runtime/kernel_build.json"
    if not manifest.exists():
        root = Path(__file__).resolve().parents[3]
        call([sys.executable, root/"bin/build_nbody_kernel.py", directory/"runtime"],
             directory, "build_kernel.log", env_script)
    record = json.loads(manifest.read_text())
    library = directory/"runtime/libfalcON.so"
    if hashlib.sha256(library.read_bytes()).hexdigest() != record["library_sha256"]:
        raise ValueError("Run-local falcON library checksum mismatch")
    return record


def build_key_plugin(directory, env_script=None):
    root = Path(__file__).resolve().parents[3]
    build_safe_kernel(directory, env_script)
    plugin = directory/"runtime/ocen_keep_keys.so"
    call(["/bin/zsh", root/"bin/build_nbody_keys.sh", root/"bin/ocen_keep_keys.cc", plugin],
         directory, "build_keys.log", env_script)
    return plugin


def executable_info(directory, env_script=None):
    names = ("gyrfalcON", "a2s", "snapprint")
    result = subprocess.run(command_args(["which", *names], env_script), cwd=directory,
                            env=runtime_environment(directory), capture_output=True, text=True, check=True)
    paths = result.stdout.strip().splitlines()
    if len(paths) != len(names):
        raise RuntimeError("Could not identify the NEMO executables")
    return {name: {"path": path, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()}
            for name, path in zip(names, paths)}


def gyrfalcon_args(config, prepared):
    kmax, dt = timestep(config)
    arguments = ["gyrfalcON", "in=initial.nemo", "out=snapshots.nemo", "eps=-1", "kernel=1",
            f"kmax={kmax}", f"Nlev={config.integrator.levels}", "fea=0.2",
            f"theta={config.integrator.theta:.17g}", f"Grav={G:.17g}",
            f"tstop={end_time(config):.17g}",
            f"step={config.integrator.output_step_myr/TIME_UNIT_MYR:.17g}",
            "give=mxvek", "startout=t", "lastout=t", "logfile=energy.log", "logstep=32",
            "manipname=ocen_keep_keys", "manippath=runtime", "manipinit=t"]
    if not config.nucleus.live or config.orbit.kind != "isolation":
        arguments += ["accname=agama", f"accpars=0,{G:.17g}", f"accfile={prepared/'external_live.ini'}"]
    return arguments


def parse_snapshots(lines):
    """Parse nbody headers and t,id,mass,x,y,z,vx,vy,vz rows lazily.

    Time is read from the explicitly formatted rows: snapprint's standalone
    time header uses only six significant digits, regardless of format=.
    """
    numeric = (line.strip() for line in lines if line.strip() and not line.lstrip().startswith("#"))
    while True:
        try:
            header = next(numeric)
        except StopIteration:
            return
        try:
            n_float = float(header)
            n = int(n_float)
            if n < 1 or n != n_float:
                raise ValueError("Invalid particle count")
            data = np.array([[float(x) for x in next(numeric).split()] for _ in range(n)])
        except (StopIteration, ValueError) as exc:
            raise ValueError("Malformed or truncated NEMO snapshot stream") from exc
        if data.shape != (n, 9) or not np.all(np.isfinite(data)):
            raise ValueError("NEMO snapshot has wrong dimensions or non-finite values")
        time = data[0, 0]
        if not np.all(data[:, 0] == time):
            raise ValueError("NEMO snapshot contains inconsistent times")
        ids = data[:, 1].astype(np.int64)
        if not np.array_equal(ids, data[:, 1]) or len(np.unique(ids)) != n:
            raise ValueError("NEMO particle keys are not unique integers")
        if np.any(data[:, 2] <= 0):
            raise ValueError("NEMO particle masses must be positive")
        order = np.argsort(ids)
        yield time, ids[order], data[order, 2], data[order, 3:]


def snapshots(path, directory, env_script=None):
    command = command_args(["snapprint", f"in={Path(path).resolve()}",
                            "options=t,key,m,x,y,z,vx,vy,vz", "header=nbody",
                            "newline=t", "format=%.17g"], env_script)
    with tempfile.TemporaryFile(mode="w+t") as error:
        process = subprocess.Popen(command, cwd=directory, env=runtime_environment(directory),
                                   stdout=subprocess.PIPE, stderr=error, text=True)
        try:
            yield from parse_snapshots(process.stdout)
            result = process.wait()
            if result:
                error.seek(0)
                raise RuntimeError(f"snapprint failed ({result}): {error.read()[-4000:]}")
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()
            process.stdout.close()
