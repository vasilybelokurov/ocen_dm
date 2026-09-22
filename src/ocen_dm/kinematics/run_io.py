"""Versioned snapshots of the resolved inputs to a kinematic fit.

Store the arrays and MGE actually in memory, rather than infer them later from
filenames or CLI defaults. Legacy runs remain readable through report.py.
"""

from __future__ import annotations

from dataclasses import asdict
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import subprocess

import numpy as np

from .likelihood import BinnedProfile, KinematicData

_ARRAY_FIELDS = ("r", "r_lower", "r_upper", "value", "err_lo", "err_hi",
                 "streaming2", "r_nodes")


def _json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def data_payload(data: KinematicData) -> dict:
    profiles = []
    for p in data.profiles:
        row = asdict(p)
        for key in _ARRAY_FIELDS:
            row[key] = None if row[key] is None else np.asarray(row[key], float).tolist()
        profiles.append(row)
    return {"schema_version": 1, "profiles": profiles}


def data_fingerprint(data: KinematicData) -> str:
    """Hash likelihood inputs; prose notes and profile order do not change identity."""
    rows = data_payload(data)["profiles"]
    for row in rows:
        row.pop("note")
        row.pop("shares_stars_with")
    return hashlib.sha256(_json(sorted(rows, key=lambda p: p["name"])).encode()).hexdigest()


def write_data_snapshot(data: KinematicData, run_dir: Path) -> dict:
    path = run_dir / "data_snapshot.json"
    path.write_text(_json(data_payload(data)) + "\n")
    return {"file": path.name, "sha256": sha256(path), "fingerprint": data_fingerprint(data)}


def read_data_snapshot(run_dir: Path, record: dict) -> KinematicData:
    path = run_dir / record["file"]
    if sha256(path) != record["sha256"]:
        raise ValueError(f"data snapshot checksum mismatch: {path}")
    payload = json.loads(path.read_text())
    if payload["schema_version"] != 1:
        raise ValueError("unsupported data snapshot version")
    profiles = []
    for row in payload["profiles"]:
        for key in _ARRAY_FIELDS:
            row[key] = None if row[key] is None else np.asarray(row[key], float)
        row["shares_stars_with"] = tuple(row["shares_stars_with"])
        profiles.append(BinnedProfile(**row))
    data = KinematicData(tuple(profiles))
    if data_fingerprint(data) != record["fingerprint"]:
        raise ValueError(f"data snapshot fingerprint mismatch: {path}")
    return data


def _prior_config(prior) -> dict | None:
    if prior is None:
        return None
    keys = {"uniform": ("lo", "hi"), "loguniform": ("lo", "hi"),
            "normal": ("mu", "sigma"), "truncnormal": ("lo", "hi", "mu", "sigma")}
    return {"kind": prior.kind, **{k: float(getattr(prior, k)) for k in keys[prior.kind]}}


def family_config(family) -> dict:
    """Capture all resolved constructor settings, MGE coefficients and priors."""
    from .fit import DarkMatterModel, NoDarkMatterModel
    from .density_conditioning import DensityConditionedModel

    if type(family) is DensityConditionedModel:
        return {"schema_version": 1, "class": "DensityConditionedModel",
                "parent_model": family_config(family.parent),
                "density": family.density, "radius_pc": family.radius_pc,
                "conditioning": "parent_loguniform_mass_full_radius_support"}

    if type(family) not in (NoDarkMatterModel, DarkMatterModel):
        raise TypeError(f"no snapshot implementation for {type(family).__name__}")
    mge = asdict(family.mge_fit)
    for key in ("sigmas_arcsec", "fractions"):
        mge[key] = np.asarray(mge[key], float).tolist()
    config = {
        "schema_version": 1, "class": type(family).__name__, "label": family.label,
        "tracer": family.tracer, "backend": family.backend,
        "instruments": list(family.nuisance_instruments),
        "fixed": {k: float(v) for k, v in family.fixed.items()},
        "constant_beta": bool(family.constant_beta), "beta0_max": family.beta0_max,
        "anisotropy_profile": family.anisotropy_profile,
        "distance_kpc": family.fixed_distance_kpc,
        "distance_prior": _prior_config(family.distance_prior), "mge_fit": mge,
        "parameters": [{"name": p.name, "prior": _prior_config(p.prior),
                        "unit": p.unit, "latex": p.latex} for p in family.parameters],
    }
    if isinstance(family, DarkMatterModel):
        config.update(gamma=family.gamma, r_t=family.r_t)
    return config


def family_from_config(config: dict, backend: str | None = None):
    from ..light_model import MGEFit
    from .fit import DarkMatterModel, NoDarkMatterModel, Parameter, Prior

    if config["schema_version"] != 1:
        raise ValueError("unsupported model snapshot version")
    if config["class"] == "DensityConditionedModel":
        from .density_conditioning import DensityConditionedModel
        if config["conditioning"] != "parent_loguniform_mass_full_radius_support":
            raise ValueError("unknown density-conditioning convention")
        return DensityConditionedModel(family_from_config(config["parent_model"], backend),
                                       config["density"], config["radius_pc"])
    classes = {"NoDarkMatterModel": NoDarkMatterModel, "DarkMatterModel": DarkMatterModel}
    cls = classes[config["class"]]
    mge = dict(config["mge_fit"])
    for key in ("sigmas_arcsec", "fractions"):
        mge[key] = np.asarray(mge[key], float)
    kw = {k: config[k] for k in ("tracer", "instruments", "fixed", "constant_beta",
                                "beta0_max", "distance_kpc")}
    kw.update(mge_fit=MGEFit(**mge), backend=backend or config["backend"],
              anisotropy_profile=config.get("anisotropy_profile", "monotonic"),
              distance_prior=Prior(**config["distance_prior"]) if config["distance_prior"] else None)
    if cls is DarkMatterModel:
        kw.update(gamma=config["gamma"], r_t=config["r_t"])
    family = cls(**kw)
    parameters = tuple(Parameter(p["name"], Prior(**p["prior"]), p["unit"], p["latex"])
                       for p in config["parameters"])
    if set(family.names) != {p.name for p in parameters}:
        raise ValueError("saved parameter names do not match the reconstructed model/backend")
    family.parameters = parameters
    family.names = tuple(p.name for p in parameters)
    return family


def code_state() -> dict:
    """Record source identity and installed numerical packages at driver entry."""
    root = Path(__file__).resolve().parents[3]
    def git(*args):
        try:
            return subprocess.run(["git", *args], cwd=root, capture_output=True,
                                  text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return None
    status = git("status", "--porcelain")
    packages = {}
    for name in ("numpy", "scipy", "astropy", "ultranest", "jampy", "agama"):
        try:
            packages[name] = version(name)
        except (PackageNotFoundError, TypeError, KeyError):
            # Some optional local installations have incomplete distribution metadata.
            packages[name] = None
    paths = sorted((root / "src").rglob("*.py")) + [root / "pyproject.toml"]
    return {"git_commit": git("rev-parse", "HEAD"),
            "dirty": None if status is None else bool(status),
            "source_hashes": {str(p.relative_to(root)): sha256(p) for p in paths if p.is_file()},
            "packages": packages}
