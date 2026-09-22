"""Regression checks for the observed single-precision P1 overflow."""
import os
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from ocen_dm.dynamics import G
from ocen_dm.dynamics import nemo
from ocen_dm.dynamics.orbits import Track


@pytest.mark.parametrize("damage", ["library_changed", "library_missing", "manifest_missing"])
def test_runtime_refuses_damaged_library_instead_of_falling_back(tmp_path, damage):
    runtime = tmp_path/"runtime"
    runtime.mkdir()
    library = runtime/"libfalcON.so"
    manifest = runtime/"kernel_build.json"
    library.write_bytes(b"test runtime library")
    manifest.write_text(json.dumps({"library_sha256": hashlib.sha256(library.read_bytes()).hexdigest()}))
    assert nemo.runtime_environment(tmp_path)["OCEN_NBODY_FALCON_LIBRARY"] == str(library.resolve())
    if damage == "library_changed":
        library.write_bytes(b"changed runtime library")
    elif damage == "library_missing":
        library.unlink()
    else:
        manifest.unlink()
    with pytest.raises(ValueError, match="Run-local falcON"):
        nemo.runtime_environment(tmp_path)


def test_p1_coefficients_match_extended_precision_without_unused_overflow(tmp_path):
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        pytest.skip("C++ compiler unavailable")
    source = tmp_path/"check.cc"
    source.write_text(r'''
#include "p1_coefficients.h"
#include <limits>
#include <iostream>
int main() {
    // Values recorded immediately before the real 100k-particle failure.
    volatile float x=5.43901e7f, previous=1.3316e30f;
    const float unsafe = 7*x*previous;
    if (std::isfinite(unsafe)) return 1;
    for (double mass: {1e-6, 1., 75., 1e4})
    for (double inverse: {1., 1e4, 1e6, 5.43901e7, 1e8})
    for (double fraction: {0., .1, .5}) {
        float xx=inverse, h=fraction/inverse, d[4]={float(mass),0,0,0};
        if (!ocen_p1_coefficients<3>(xx,h,d)) return 2;
        // Independent reference: form the extra raw derivative in long double,
        // then apply the original softened-potential derivative formula.
        long double raw[5]; raw[0]=static_cast<long double>(float(mass))*std::sqrt(static_cast<long double>(xx));
        for (int n=1;n<5;++n) raw[n]=(2*n-1)*static_cast<long double>(xx)*raw[n-1];
        for (int n=0;n<4;++n) {
            const long double reference=raw[n]+static_cast<long double>(h)*raw[n+1];
            if (std::abs(d[n]/reference-1)>2e-6) return 3;
        }
    }
    float too_large[4]={1e30f,0,0,0};
    if (ocen_p1_coefficients<3>(1e8f, 5e-9f, too_large)) return 4;
    return 0;
}
''')
    include = Path(__file__).resolve().parents[1]/"bin/nbody_compat"
    subprocess.run([compiler, "-std=c++11", "-O2", "-I"+str(include), str(source),
                    "-o", str(tmp_path/"check")], check=True, capture_output=True, text=True)
    subprocess.run([str(tmp_path/"check")], check=True, capture_output=True, text=True)


@pytest.mark.skipif(os.environ.get("OCEN_RUN_NEMO_TESTS") != "1", reason="local compiled NEMO required")
def test_close_cell_p1_forces_are_finite_and_match_direct_sum(tmp_path):
    # Two compact cells at a separation/softening that overflow the old float
    # Taylor calculation. Verify actual tree forces against direct P1 forces.
    rng = np.random.default_rng(12)
    position = rng.normal(size=(64, 3))*1e-6
    position[:32, 0] -= 4e-5
    position[32:, 0] += 4e-5
    particles = dict(xv=np.column_stack((position, np.zeros_like(position))),
                     mass=np.full(64, 4.), softening=np.full(64, 1e-4),
                     particle_id=np.arange(64))
    track = Track(np.array([0., 1e-8]), np.zeros((2, 6)))
    nemo.build_key_plugin(tmp_path)
    nemo.write_input(particles, track, tmp_path)
    nemo.call(["gyrfalcON", "in=initial.nemo", "out=snapshots.nemo", "eps=-1", "kernel=1",
               "kmax=30", "Nlev=1", "theta=0.5", f"Grav={G:.17g}", "tstop=1e-9", "step=1e-9",
               "give=mxveka", "startout=t", "lastout=t", "manipname=ocen_keep_keys",
               "manippath=runtime", "manipinit=t"], tmp_path, "gyrfalcon.log")
    nemo.call(["snapprint", "in=snapshots.nemo", "times=0", "options=t,key,ax,ay,az",
               "header=", "format=%.17g"], tmp_path, "forces.txt")
    forces = np.loadtxt(tmp_path/"forces.txt")
    # NEMO's time selector has a tolerance that includes this tiny final step.
    forces = forces[forces[:, 0] == 0]
    forces = forces[np.argsort(forces[:, 1]), 2:]
    assert forces.shape == (64, 3) and np.all(np.isfinite(forces))
    delta = position[:, None, :] - position[None, :, :]
    q = np.sum(delta**2, axis=2)+1e-8
    factor = G*4*(q**(-1.5)+1.5e-8*q**(-2.5))
    exact = -np.sum(factor[:, :, None]*delta, axis=1)
    assert np.max(np.linalg.norm(forces-exact, axis=1)/np.linalg.norm(exact, axis=1)) < .005
