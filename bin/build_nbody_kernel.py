#!/usr/bin/env python3
"""Build a run-local falcON library with stable P1 scalar Taylor coefficients.

Uses the installed objects/ABI and replaces only kernel.o. No installed file is
modified. Run inside the selected NEMO environment, as for build_nbody_keys.sh.
An unknown upstream implementation is rejected instead of patched heuristically.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import sys


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output / "kernel_build.json").exists():
        raise FileExistsError("Kernel build already exists; use a fresh evolution directory")
    falcon = Path(os.environ["FALCON"])
    nemo = Path(os.environ["NEMO"])
    compat = Path(__file__).resolve().parent / "nbody_compat"
    upstream = falcon / "src/public/lib/kernel.cc"
    source = upstream.read_text()
    start = source.index("  template<> struct _block<p1,1>")
    stop = source.index("  //////////////////////////////////////////////////////////////////////////////\n  // kern_type = p2", start)
    original = source[start:stop]
    expected = "19ce6a8807b50e791b402ad5d5f6c18b1028f9a3062d944c8e4c8a21601128ad"
    if hashlib.sha256(original.encode()).hexdigest() != expected:
        raise ValueError("Unrecognized falcON P1 kernel; review the source before applying this patch")
    precision = os.environ.get("OCEN_NEMO_PRECISION", "SINGLE")
    if precision not in ("SINGLE", "DOUBLE"):
        raise ValueError("OCEN_NEMO_PRECISION must match the installed SINGLE or DOUBLE ABI")
    replacement = '''  template<int K> struct _block<p1,K> : public _setE<p1> {
    enum { ND=K+1 };
    sv b(real&X, real D[ND], real, real HQ, real) {
      if(!ocen_p1_coefficients<K>(X,HQ,D))
        falcON_THROW("P1 Taylor coefficients exceed native precision; rescale units or use a double-precision build");
    } };
'''
    source = source[:start] + replacement + source[stop:]
    source = source.replace("#include <public/kernel.h>",
                            '#include <public/kernel.h>\n#include "p1_coefficients.h"')
    (output / "kernel.cc").write_text(source)
    compiler = shlex.split(os.environ.get("OCEN_NEMO_LIBRARY_CXX", "clang++" if platform.system() == "Darwin" else "g++"))
    compile_command = compiler + ["-std=c++11", "-O2", "-fPIC", "-D_FILE_OFFSET_BITS=64",
        "-DfalcON_NEMO", "-DfalcON_" + precision, "-I" + str(compat), "-I" + str(falcon / "inc"),
        "-I" + str(falcon / "inc/utils"), "-I" + str(nemo / "inc"),
        "-c", str(output / "kernel.cc"), "-o", str(output / "kernel.o")]
    names = ("basic", "nemo++", "body", "tree", "gravity", "kernel", "partner", "nbody",
             "forcesC", "tools", "sample", "manip", "profile", "bodyfunc", "neighbours", "PotExp")
    objects = [output / "kernel.o" if n == "kernel" else falcon / "lib" / (n + ".o") for n in names]
    library = output / "libfalcON.so"
    temporary = output / "libfalcON.building.so"
    link_command = compiler + (["-dynamiclib"] if platform.system() == "Darwin" else ["-shared"])
    link_command += list(map(str, objects)) + ["-L" + str(falcon.parent / "utils/lib"), "-lWDutils",
                                             "-L" + str(nemo / "lib"), "-lnemo"]
    if platform.system() == "Darwin":
        link_command += ["-Wl,-undefined,dynamic_lookup", "-Wl,-install_name," + str(library)]
    link_command += ["-o", str(temporary)]
    for command in (compile_command, link_command):
        print(shlex.join(command), flush=True)
        subprocess.run(command, check=True)
    temporary.replace(library)
    record = {"patch": "stable P1 Taylor coefficients; no overflowing unused derivative",
        "precision": precision, "upstream_kernel_sha256": digest(upstream),
        "original_p1_block_sha256": expected, "patched_kernel_sha256": digest(output / "kernel.cc"),
        "helper_sha256": digest(compat / "p1_coefficients.h"), "builder_sha256": digest(__file__),
        "input_object_sha256": {str(p): digest(p) for p in objects},
        "original_library_sha256": digest(falcon / "lib/libfalcON.so"),
        "wdutils_library_sha256": digest(falcon.parent / "utils/lib/libWDutils.so"),
        "library": str(library), "library_sha256": digest(library),
        "compile_command": compile_command, "link_command": link_command,
        "compiler_version": subprocess.check_output(compiler + ["--version"], text=True)}
    (output / "kernel_build.json").write_text(json.dumps(record, indent=2) + "\n")
    print("Built " + str(library), flush=True)


if __name__ == "__main__":
    main()
