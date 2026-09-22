#!/usr/bin/env python3
"""Replay the preserved spherical-mapper input in polar/equatorial planes.

The exported potential slightly changes rounding; the original nonfinite
polar output is recorded separately in the fixture, not asserted on replay.
"""
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ocen-df-mpl")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

import numpy as np
from ocen_dm.kinematics.positive_df import agama_pc


def main():
    fixture = ROOT/"tests/fixtures/compact_frequency"
    record = json.loads((fixture/"orbit.json").read_text())
    ag = agama_pc()
    pot = ag.Potential(str(fixture/"potential.ini"))
    mapper = ag.ActionMapper(pot)
    polar = np.array([record["action_angles"]])
    equatorial = polar.copy()
    equatorial[:, 2], equatorial[:, 1] = polar[:, 1], 0.
    print("Originally recorded:", record["original_polar_xv"])
    print("Replayed with exported potential (rounding differs):")
    for label, aa in (("polar", polar), ("equatorial", equatorial)):
        xv = mapper(aa)
        print(label, "phase space:", xv.tolist())
        if np.all(np.isfinite(xv)):
            actions, omega = ag.actions(pot, xv, frequencies=True)
            print(label, "recovered Jr, L:", actions[:, 0].tolist(),
                  (actions[:, 1]+abs(actions[:, 2])).tolist())
            print(label, "Omega_r/Omega_t:", (omega[:, 0]/omega[:, 1]).tolist())


if __name__ == "__main__":
    main()
