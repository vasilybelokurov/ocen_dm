"""Omega Centauri dark-matter / tidal-tail inference pipeline.

Milestone 1 scope: data provenance, deterministic ingestion, validated
loaders. No scientific inference lives in this package yet.
"""

__version__ = "0.1.0"

from .paths import load_local_env

#: machine-specific locations of external inputs, from the gitignored .env (see kb/getting_started)
load_local_env()
