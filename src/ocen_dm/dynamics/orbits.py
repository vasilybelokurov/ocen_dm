"""Full Cartesian forward tracks with explicit host clocks and units."""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import shutil

import numpy as np
from scipy.interpolate import CubicHermiteSpline

from . import TIME_UNIT_MYR
from .models import agama_module


def timestep(config):
    kmax = math.ceil(math.log2(TIME_UNIT_MYR/config.integrator.max_step_myr))
    return kmax, 2.**(-kmax)


def end_time(config):
    _, dt = timestep(config)
    # End on a full NEMO step; retain the host's absolute clock (especially bar phase).
    duration = math.ceil(config.orbit.duration_myr/TIME_UNIT_MYR/dt)*dt
    return config.orbit.start_time_myr/TIME_UNIT_MYR + duration


@dataclass
class Track:
    time: np.ndarray  # kpc/(km/s), increasing
    xv: np.ndarray    # kpc, km/s in the host's inertial Cartesian frame

    def __post_init__(self):
        self.time = np.asarray(self.time, float)
        self.xv = np.asarray(self.xv, float)
        if self.time.ndim != 1 or self.xv.shape != (len(self.time), 6) or len(self.time) < 2:
            raise ValueError("Track needs N increasing times and an N x 6 Cartesian array")
        if not np.all(np.isfinite(self.xv)) or not np.all(np.isfinite(self.time)):
            raise ValueError("Track contains non-finite values")
        if np.any(np.diff(self.time) <= 0):
            raise ValueError("Track times must increase; reorder rows without reversing velocities")
        self._spline = CubicHermiteSpline(self.time, self.xv[:, :3], self.xv[:, 3:])

    def at(self, time):
        if np.any(np.asarray(time) < self.time[0]-1e-10) or np.any(np.asarray(time) > self.time[-1]+1e-10):
            raise ValueError("Track extrapolation is forbidden")
        return np.concatenate((self._spline(time), self._spline(time, 1)), axis=-1)

    def save(self, directory):
        np.savez_compressed(directory / "orbit.npz", time=self.time, xv=self.xv,
                            time_unit_myr=TIME_UNIT_MYR)
        np.savetxt(directory / "center.txt", np.column_stack((self.time, self.xv)),
                   fmt="%.17g", header="time[kpc/(km/s)] x y z[kpc] vx vy vz[km/s]")

    @classmethod
    def load(cls, directory):
        with np.load(Path(directory)/"orbit.npz", allow_pickle=False) as data:
            if not np.isclose(data["time_unit_myr"], TIME_UNIT_MYR, rtol=1e-12):
                raise ValueError("Unrecognized orbit time unit")
            return cls(data["time"], data["xv"])


def _orbit(agama, potential, initial, start, stop, count, accuracy):
    time, xv = agama.orbit(potential=potential, ic=initial, timestart=start,
                           time=stop-start, trajsize=count, accuracy=accuracy,
                           dtype="float64", verbose=False)
    order = np.argsort(time)
    return np.asarray(time)[order], np.asarray(xv)[order]


def build_track(config, directory):
    agama = agama_module()
    spec = config.orbit
    start, stop = spec.start_time_myr/TIME_UNIT_MYR, end_time(config)
    count = max(3, int(np.ceil((stop-start)*TIME_UNIT_MYR/spec.sample_step_myr))+1)
    notes = {"kind": spec.kind, "friction": False, "nucleus_motion": "prescribed",
             "time_unit_myr": TIME_UNIT_MYR, "start_time_myr": start*TIME_UNIT_MYR,
             "end_time_myr": stop*TIME_UNIT_MYR,
             "requested_duration_myr": spec.duration_myr,
             "effective_duration_myr": (stop-start)*TIME_UNIT_MYR}
    if spec.kind == "isolation":
        track = Track(np.linspace(start, stop, count), np.zeros((count, 6)))
        host = None
    else:
        initial = spec.initial_xv
        if spec.kind == "class1":
            source = Path(agama.__file__).parent/"data/McMillan17.ini"
            host = agama.Potential(str(source))
            shutil.copyfile(source, directory/"host.ini")
            if initial is None:
                from ..tails.progenitor_orbits import present_day_orbit
                orbit = present_day_orbit()
                present = np.array([orbit.x(), orbit.y(), orbit.z(), orbit.vx(), orbit.vy(), orbit.vz()])
                if spec.end_at_present:
                    _, xb = _orbit(agama, host, present, stop, start, 2, config.integrator.accuracy)
                    initial = xb[0]
                    notes["initial_state_source"] = "observed present orbit rewound by effective duration"
                    notes["target_present_xv"] = present.tolist()
                else:
                    tb, xb = _orbit(agama, host, present, 0., -.3, 12001, config.integrator.accuracy)
                    radius = np.linalg.norm(xb[:, :3], axis=1)
                    maxima = np.where((radius[1:-1] > radius[:-2]) & (radius[1:-1] > radius[2:]))[0]+1
                    if not len(maxima):
                        raise ValueError("Could not find a past apocentre for class1")
                    initial = xb[maxima[-1]]
                    notes["initial_state_source"] = "nearest past apocentre of nominal present orbit"
                    notes["rewind_myr"] = float(-tb[maxima[-1]]*TIME_UNIT_MYR)
        elif spec.kind == "class3":
            from ..tails import bar_migration as bm
            history = bm.BarHistory(spec.omega_final)
            if stop > history.tf:
                raise ValueError("class3 end time exceeds the bar history; preserve its absolute clock")
            if spec.end_at_present and abs(stop-history.tf) > 1e-10:
                raise ValueError("class3 end_at_present requires the exact end of the bar history")
            host, _ = bm.slowing_bar_potential(history)
            for suffix, target in (("axi", "host_axis.ini"), ("baryon_full", "host_bar.ini"),
                                   ("baryon_axi", "host_bar_axis.ini")):
                shutil.copyfile(bm.POT_DIR/f"MWPotentialHunter24_{suffix}.ini", directory/target)
            size = history.omega_f/history.omega
            amplitude = history.frac*size
            np.savetxt(directory/"bar_scale.txt", np.column_stack((history.t, amplitude, size)), fmt="%.17g")
            np.savetxt(directory/"bar_negative_scale.txt", np.column_stack((history.t, -amplitude, size)), fmt="%.17g")
            np.savetxt(directory/"bar_rotation.txt", np.column_stack((history.t, history.phi)), fmt="%.17g")
            (directory/"host.ini").write_text(
                "[Potential axis]\nfile=host_axis.ini\n\n"
                "[Potential bar]\nfile=host_bar.ini\nscale=bar_scale.txt\nrotation=bar_rotation.txt\n\n"
                "[Potential subtract]\nfile=host_bar_axis.ini\nscale=bar_negative_scale.txt\n")
            if initial is None:
                _, backwards = _orbit(agama, host, bm.present_day_samples(1)[0],
                                      history.tf, start, 2, config.integrator.accuracy)
                initial = backwards[0]
                notes["initial_state_source"] = "nominal present orbit rewound in this bar potential"
                if spec.end_at_present:
                    notes["target_present_xv"] = bm.present_day_samples(1)[0].tolist()
            notes["bar_omega_final_kms_kpc"] = spec.omega_final
            notes["bar_present_time_myr"] = history.tf*TIME_UNIT_MYR
        else:
            # Analytic, portable host for integration tests and controlled demonstrations.
            host = agama.Potential(type="Plummer", mass=1e11, scaleRadius=2.)
            (directory/"host.ini").write_text("[Potential]\ntype=Plummer\nmass=1e11\nscaleRadius=2\n")
            if initial is None:
                vc = np.sqrt(-5*host.force([5., 0., 0.])[0])
                initial = [5., 0., 0., 0., vc, 0.]
            notes["host"] = {"type": "Plummer", "mass_msun": 1e11, "scale_kpc": 2.}
        # Reload the exported potential before using it: the simulation uses this exact asset.
        restored = agama.Potential(str(directory / "host.ini"))
        for t in (start, (start+stop)/2, stop):
            probes = np.array([[1., .2, .1], [5., -2., .3], [15., 1., 2.]])
            if not np.allclose(restored.force(probes, t=t), host.force(probes, t=t), rtol=1e-6, atol=1e-7):
                raise ValueError("Exported host does not reproduce the original time-dependent force")
        host = restored
        time, xv = _orbit(agama, host, initial, start, stop, count, config.integrator.accuracy)
        track = Track(time, xv)
        if "target_present_xv" in notes:
            offset = track.xv[-1]-notes["target_present_xv"]
            notes["present_endpoint_position_error_pc"] = float(np.linalg.norm(offset[:3])*1000)
            notes["present_endpoint_velocity_error_kms"] = float(np.linalg.norm(offset[3:]))
            if np.linalg.norm(offset[:3]) > .01 or np.linalg.norm(offset[3:]) > .1:
                raise ValueError("Historical track failed present-day endpoint replay")
        midpoint = (track.time[1:]+track.time[:-1])/2
        # A coarse trajectory can give a moving nucleus spurious accelerations.
        selection = np.unique(np.linspace(0, len(midpoint)-1, min(len(midpoint), 3000)).astype(int))
        errors = []
        for i in selection:
            t = midpoint[i]
            actual = host.force(track.at(t)[:3], t=t)
            errors.append(np.linalg.norm(track._spline(t, 2)-actual)/max(np.linalg.norm(actual), 1e-20))
        notes["max_interpolated_acceleration_relative_error"] = float(max(errors))
        if max(errors) > .01:
            raise ValueError("Orbit interpolation acceleration error exceeds 1%; reduce sample_step_myr")
    track.save(directory)
    return track, notes


def write_external_potential(directory, frozen=False, config=None):
    """Inertial host + moving nucleus, optionally the frozen initial DM/stars.

    Frozen controls retain the same initial satellite potential as live runs.
    They measure the effect of suppressing its gravitational response.
    """
    parts = []
    if (directory/"host.ini").exists():
        parts.append("[Potential host]\nfile=host.ini\n")
    if frozen or config is None or not config.nucleus.live:
        asset = "satellite_initial.ini" if frozen else "nucleus.ini"
        parts.append(f"[Potential satellite]\nfile={asset}\ncenter=center.txt\n")
    path = directory/("external_frozen.ini" if frozen else "external_live.ini")
    path.write_text("\n".join(parts) if parts else "# Isolated live system: no external force.\n")
    return path
