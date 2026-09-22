"""Strict, portable configuration for the two dynamical programmes."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
import math
import re

import yaml


def _positive(name, value):
    if isinstance(value, bool) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _construct(cls, values):
    if not isinstance(values, dict):
        raise ValueError(f"{cls.__name__} must be a mapping")
    unknown = set(values) - {f.name for f in fields(cls)}
    if unknown:
        raise ValueError(f"Unknown {cls.__name__} settings: {sorted(unknown)}")
    return cls(**values)


@dataclass(frozen=True)
class Nucleus:
    mass_msun: float = 3.55e6
    scale_pc: float = 5.4
    live: bool = False
    n_particles: int = 100000
    softening_pc: float = .1
    beta: float = 0.0

    def __post_init__(self):
        _positive("nucleus mass", self.mass_msun)
        _positive("nucleus scale", self.scale_pc)
        if type(self.live) is not bool or type(self.n_particles) is not int or self.n_particles < 16:
            raise ValueError("nucleus live must be boolean and n_particles >= 16")
        _positive("nucleus softening", self.softening_pc)
        if not -.49 <= self.beta <= 0:
            raise ValueError("Plummer nucleus requires -0.49 <= beta <= 0")

    def component(self):
        return Component(name="nucleus", scale_pc=self.scale_pc, cutoff_pc=100*self.scale_pc,
                         n_particles=self.n_particles, softening_pc=self.softening_pc,
                         gamma=0., outer_slope=5., beta=self.beta, mass_msun=self.mass_msun,
                         profile="plummer")


@dataclass(frozen=True)
class Component:
    """Spherical double power law with exp[-(r/r_cut)^2] density taper.

    A core means gamma=0, alpha=1 (finite central density, not a Plummer core).
    Choose either a total mass or an aperture mass, never both.
    """
    name: str
    scale_pc: float
    cutoff_pc: float
    n_particles: int
    softening_pc: float
    gamma: float = 1.0
    outer_slope: float = 3.0
    beta: float = 0.0
    mass_msun: float | None = None
    mass_within_msun: float | None = None
    radius_pc: float | None = None
    refinement_radius_pc: float | None = None
    refinement_floor: float = .001
    profile: str = "spheroid"

    def __post_init__(self):
        if self.name not in ("dm", "stars", "nucleus"):
            raise ValueError("Component name must be dm, stars or nucleus")
        if self.profile not in ("spheroid", "plummer"):
            raise ValueError("Unknown density profile")
        if self.profile == "plummer" and (self.gamma != 0 or self.outer_slope != 5 or self.mass_msun is None):
            raise ValueError("Plummer component requires gamma=0, outer_slope=5 and total mass")
        for name in ("scale_pc", "cutoff_pc", "softening_pc"):
            _positive(name, getattr(self, name))
        if self.refinement_radius_pc is not None:
            _positive("refinement_radius_pc", self.refinement_radius_pc)
        if not 0 < self.refinement_floor <= 1:
            raise ValueError("refinement_floor must lie in (0, 1]")
        if type(self.n_particles) is not int or self.n_particles < 16:
            raise ValueError("n_particles must be an integer >= 16")
        if not 0 <= self.gamma < 2 or not 2 < self.outer_slope <= 8:
            raise ValueError("Require 0 <= gamma < 2 and 2 < outer_slope <= 8")
        if not -.49 <= self.beta <= .49:
            raise ValueError("The DF validator supports -0.49 <= beta <= 0.49")
        if self.gamma < 2 * self.beta:
            raise ValueError("Central slope must satisfy gamma >= 2 beta")
        if (self.mass_msun is None) == (self.mass_within_msun is None):
            raise ValueError("Choose exactly one of mass_msun and mass_within_msun")
        if self.mass_msun is not None:
            _positive("mass_msun", self.mass_msun)
            if self.radius_pc is not None:
                raise ValueError("radius_pc belongs to mass_within_msun")
        else:
            _positive("mass_within_msun", self.mass_within_msun)
            if self.radius_pc is None:
                raise ValueError("mass_within_msun requires radius_pc")
            _positive("radius_pc", self.radius_pc)


@dataclass(frozen=True)
class Orbit:
    kind: str = "isolation"
    duration_myr: float = 10.0
    sample_step_myr: float = .05
    start_time_myr: float = 0.0
    initial_xv: list[float] | None = None
    omega_final: float = 24.0
    end_at_present: bool = False

    def __post_init__(self):
        if self.kind not in ("isolation", "class1", "class3", "plummer_host"):
            raise ValueError("Orbit kind must be isolation, class1, class3 or plummer_host; "
                             "mass-dependent friction is not implemented")
        _positive("duration_myr", self.duration_myr)
        _positive("sample_step_myr", self.sample_step_myr)
        _positive("omega_final", self.omega_final)
        if type(self.end_at_present) is not bool:
            raise ValueError("end_at_present must be boolean")
        if self.end_at_present and (self.kind not in ("class1", "class3") or self.initial_xv is not None):
            raise ValueError("end_at_present requires class1/class3 and no explicit initial_xv")
        if not math.isfinite(self.start_time_myr) or self.start_time_myr < 0:
            raise ValueError("start_time_myr must be finite and nonnegative")
        if self.initial_xv is not None:
            if len(self.initial_xv) != 6 or not all(map(math.isfinite, self.initial_xv)):
                raise ValueError("initial_xv must contain six finite numbers (kpc, km/s)")
            if self.kind == "isolation":
                raise ValueError("An isolated nucleus is stationary at the origin")
        if self.sample_step_myr > self.duration_myr:
            raise ValueError("Orbit sampling step exceeds duration")


@dataclass(frozen=True)
class Integrator:
    max_step_myr: float = .02
    levels: int = 4
    output_step_myr: float = 1.0
    theta: float = .5
    accuracy: float = 1e-8

    def __post_init__(self):
        for name in ("max_step_myr", "output_step_myr", "theta", "accuracy"):
            _positive(name, getattr(self, name))
        if type(self.levels) is not int or not 1 <= self.levels <= 20:
            raise ValueError("levels must be an integer from 1 to 20")
        if self.output_step_myr < self.max_step_myr:
            raise ValueError("output_step_myr must be >= max_step_myr")
        if self.theta > 1 or self.accuracy >= .01:
            raise ValueError("Require theta <= 1 and accuracy < 0.01")


@dataclass(frozen=True)
class Experiment:
    programme: str
    components: tuple[Component, ...]
    seed: int = 42
    nucleus: Nucleus = field(default_factory=Nucleus)
    orbit: Orbit = field(default_factory=Orbit)
    integrator: Integrator = field(default_factory=Integrator)
    radii_pc: tuple[float, ...] = (3., 10., 20., 35., 50., 70.)

    def __post_init__(self):
        names = [c.name for c in self.components]
        expected = {"remnant": ["dm"], "progenitor": ["dm", "stars"]}
        if self.programme not in expected or sorted(names) != expected[self.programme]:
            raise ValueError("remnant requires dm; progenitor requires dm and stars, once each")
        if type(self.seed) is not int or not 0 <= self.seed < 2**31:
            raise ValueError("seed must be a nonnegative 31-bit integer")
        if not self.radii_pc or any(not math.isfinite(r) or r <= 0 for r in self.radii_pc):
            raise ValueError("radii_pc must be positive finite radii")
        if sorted(set(self.radii_pc)) != list(self.radii_pc):
            raise ValueError("radii_pc must be unique and increasing")

    def to_dict(self):
        # Round trip through YAML uses only portable list/dict/scalar types.
        out = asdict(self)
        out["components"] = list(out["components"])
        out["radii_pc"] = list(out["radii_pc"])
        return out

    @property
    def particle_components(self):
        return self.components + ((self.nucleus.component(),) if self.nucleus.live else ())

    @classmethod
    def from_dict(cls, data):
        data = dict(data)
        data["components"] = tuple(_construct(Component, c) for c in data["components"])
        for key, kind in (("nucleus", Nucleus), ("orbit", Orbit), ("integrator", Integrator)):
            data[key] = _construct(kind, data.get(key, {}))
        if "radii_pc" in data:
            data["radii_pc"] = tuple(data["radii_pc"])
        return _construct(cls, data)

    @classmethod
    def load(cls, path):
        return cls.from_dict(yaml.safe_load(Path(path).read_text()))


def safe_label(label):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label):
        raise ValueError("Run labels must be a single name using letters, digits, _, - or .")
    return label
