"""Vasiliev & Baumgardt (2021) Gaia EDR3 globular-cluster members: omega Cen.

Source: Zenodo 10.5281/zenodo.4891252 (``clusters.zip``). The working file is the
user's FITS compilation of the zip's per-cluster tables in the canonical store,
``~/data/catalogues/gc_members_gaia_vasiliev.fits``; its omega Cen block was
checked identical to the zip's ``catalogues/NGC_5139_oCen.txt`` on 2026-09-16.

Column meanings come from the zip's ``!readme.txt``. The FITS carries no units,
so every input unit is declared in ``configs/column_maps.yaml`` from that readme.

Two products:

* ``vasiliev2021_ocen_members`` -- 228,055 stars to G = 21 within 0.67 deg, with
  per-star PM covariance and membership probability;
* ``vasiliev2021_ocen_pm_profiles`` -- the authors' PM dispersion and rotation
  profiles (percentiles vs radius) from ``profiles/NGC_5139_oCen.txt``, read from
  the zip at its canonical location. That table has 12 columns where the readme
  documents 11; the twelfth is carried as ``undocumented_col12`` and must not be
  used until its meaning is established.

Likelihood rule: membership probabilities are astrometric only (no photometric
or chemical information) and are model-dependent; use them for sample selection
with a stated threshold, and treat that threshold as part of the selection
function, not as a per-star weight in the science likelihood.
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.table import Table

from ..provenance import load_registry
from ._product import build_product
from .base import raw_lineage, schema_report, write_processed
from .schema import ColumnRole, SchemaError, TableSchema

__all__ = ["DATASET", "CLUSTER_NAME", "SCHEMA", "PROFILE_SCHEMA", "members_path",
           "zip_path", "load", "load_profiles", "preprocess"]

DATASET = "vasiliev2021_gc_members"
CLUSTER_NAME = "NGC_5139_oCen"

LIKELIHOOD_RULE = (
    "astrometric membership probability: use as a stated selection threshold that "
    "belongs to the selection function, not as a per-star likelihood weight"
)

SCHEMA = TableSchema(
    name="vasiliev2021_ocen_members",
    unique_id="source_id",
    notes="Gaia EDR3 stars within 0.67 deg of omega Cen with astrometric membership.",
    roles=(
        ColumnRole("source_id", None, ("SOURCE_ID", "source_id"), kind="integer",
                   is_identifier=True, description="Gaia EDR3 source_id"),
        ColumnRole("ra", "deg", ("RA", "ra"), valid_range=(0.0, 360.0)),
        ColumnRole("dec", "deg", ("DEC", "dec"), valid_range=(-90.0, 90.0)),
        ColumnRole("x", "deg", ("X", "x"), description="cluster-centred tangent-plane x"),
        ColumnRole("y", "deg", ("Y", "y"), description="cluster-centred tangent-plane y"),
        ColumnRole("parallax", "mas", ("PLX", "plx"), allow_missing=True,
                   description="parallax with Lindegren et al. zero-point correction"),
        ColumnRole("pmra", "mas / yr", ("PMRA", "pmra")),
        ColumnRole("pmdec", "mas / yr", ("PMDEC", "pmdec")),
        ColumnRole("parallax_error", "mas", ("PLXE", "plxe"), allow_missing=True,
                   valid_range=(0.0, None)),
        ColumnRole("pmra_error", "mas / yr", ("PMRAE", "pmrae"), valid_range=(0.0, None)),
        ColumnRole("pmdec_error", "mas / yr", ("PMDECE", "pmdece"), valid_range=(0.0, None)),
        ColumnRole("pmra_pmdec_corr", None, ("PMCORR", "pmcorr"), dimensionless=True,
                   valid_range=(-1.0, 1.0)),
        ColumnRole("g_mag", "mag", ("G_MAG", "g_mag"), allow_missing=True),
        ColumnRole("bp_rp", "mag", ("BP_RP", "bp_rp"), allow_missing=True),
        ColumnRole("source_density", "1 / arcmin2", ("SIGMA", "Sigma"), allow_missing=True,
                   description="local source density used for the error scaling"),
        ColumnRole("quality_flag", None, ("GFLAG", "qflag"), kind="integer",
                   description="bit0: 1 = 5-parameter solution; bit1: 2 = passed all quality filters"),
        ColumnRole("membership_prob", None, ("PROB", "memberprob"), dimensionless=True,
                   valid_range=(0.0, 1.0)),
    ),
)

PROFILE_SCHEMA = TableSchema(
    name="vasiliev2021_ocen_pm_profiles",
    notes="Authors' PM dispersion and rotation percentiles versus radius.",
    roles=(
        ColumnRole("r", "arcsec", ("r",), valid_range=(0.0, None)),
        ColumnRole("sigma_pm_p02", "mas / yr", ("s02",), valid_range=(0.0, None)),
        ColumnRole("sigma_pm_p16", "mas / yr", ("s16",), valid_range=(0.0, None)),
        ColumnRole("sigma_pm", "mas / yr", ("s50",), valid_range=(0.0, None),
                   description="median PM dispersion"),
        ColumnRole("sigma_pm_p84", "mas / yr", ("s84",), valid_range=(0.0, None)),
        ColumnRole("sigma_pm_p98", "mas / yr", ("s98",), valid_range=(0.0, None)),
        ColumnRole("vrot_pm_p02", "mas / yr", ("v02",)),
        ColumnRole("vrot_pm_p16", "mas / yr", ("v16",)),
        ColumnRole("vrot_pm", "mas / yr", ("v50",), description="median PM rotation amplitude"),
        ColumnRole("vrot_pm_p84", "mas / yr", ("v84",)),
        ColumnRole("vrot_pm_p98", "mas / yr", ("v98",)),
        ColumnRole("undocumented_col12", None, ("c12",), dimensionless=True, required=False,
                   description="12th column of profiles/*.txt, absent from the readme; meaning unknown"),
    ),
)


def _record(name_fragment: str):
    dataset = load_registry()[DATASET]
    for record in dataset.files:
        if name_fragment in record.name:
            return record
    raise KeyError(f"{DATASET}: no file matching {name_fragment!r} in the registry")


def members_path() -> Path:
    """Canonical-store path of the FITS compilation."""
    return _record("gc_members_gaia_vasiliev").path


def zip_path() -> Path:
    """Canonical-store path of the Zenodo zip (for the profiles)."""
    return _record("clusters.zip").path


def _is_ocen(table: Table) -> np.ndarray:
    names = np.char.strip(np.asarray(table["NAME"]).astype(str))
    return names == CLUSTER_NAME


def load(validate: bool = True) -> Any:
    """Standardized omega Cen member table (all 228,055 stars, no probability cut)."""
    from .base import read_table, standardize
    from .schema import validate_table

    table = read_table(members_path())
    table = table[_is_ocen(table)]
    out = standardize(table, SCHEMA, keep_extra_columns=False)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_likelihood_rule"] = LIKELIHOOD_RULE
    if validate:
        validate_table(out, SCHEMA)
    return out


def load_profiles() -> Table:
    """Read ``profiles/NGC_5139_oCen.txt`` from the zip into a raw table.

    Columns are positional (the file has no header); names follow the readme's
    order, with the undocumented twelfth column kept but flagged.
    """
    path = zip_path()
    with zipfile.ZipFile(path) as archive:
        raw = archive.read(f"profiles/{CLUSTER_NAME}.txt").decode("ascii")
    data = np.loadtxt(io.StringIO(raw))
    if data.ndim != 2 or data.shape[1] not in (11, 12):
        raise SchemaError(
            f"{DATASET}: profiles/{CLUSTER_NAME}.txt has shape {data.shape}; the readme "
            "documents 11 columns (12 observed on 2026-09-16)"
        )
    names = ["r", "s02", "s16", "s50", "s84", "s98", "v02", "v16", "v50", "v84", "v98"]
    if data.shape[1] == 12:
        names.append("c12")
    table = Table(data, names=names)
    table["r"].unit = u.deg
    for n in names[1:11]:
        table[n].unit = u.mas / u.yr
    return table


def preprocess() -> list[dict[str, Any]]:
    """Write both products: the member catalogue and the PM profiles."""
    reports: list[dict[str, Any]] = []

    reports.append(build_product(
        dataset=DATASET,
        schema=SCHEMA,
        path=members_path(),
        subdir="tails",
        likelihood_rule=LIKELIHOOD_RULE,
        row_filter=_is_ocen,
        keep_extra_columns=False,
        extra_meta={"ocen_cluster_name": CLUSTER_NAME,
                    "ocen_source_doi": "10.5281/zenodo.4891252"},
    ))

    from .base import standardize

    table = load_profiles()
    out = standardize(table, PROFILE_SCHEMA, keep_extra_columns=False)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_likelihood_rule"] = (
        "authors' binned PM profiles from the same stars as the member catalogue: "
        "not independent of a likelihood built on those stars"
    )
    lineage = raw_lineage(zip_path())
    lineage["member"] = f"profiles/{CLUSTER_NAME}.txt"
    out.meta["ocen_raw"] = lineage
    report = schema_report(out, PROFILE_SCHEMA)
    report["path"] = str(write_processed(out, PROFILE_SCHEMA.name, subdir="kinematics"))
    report["dataset"] = DATASET
    report["lineage"] = lineage
    reports.append(report)
    return reports
