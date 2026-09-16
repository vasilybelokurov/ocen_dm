"""Loaders for the oMEGACat internal-kinematics data products.

oMEGACat VI (https://doi.org/10.3847/1538-4357/adbe67, data
https://doi.org/10.5281/zenodo.14978551) supplies the binned proper-motion and
line-of-sight profiles that are the v1 internal-kinematics likelihood data.

**File names are verified, column names are not.** Appendix A of the paper
(arXiv:2503.04903v2, fetched 2026-09-16) names the deposited products; those
names are listed in :data:`PRODUCT_FILES`. The literal FITS column names could
not be inspected because Zenodo returned HTTP 504, so every profile schema
requires an explicit map in ``configs/column_maps.yaml`` before it will load.

The schema structure follows the published Tables 2 and 3, which give
``r_lower, r_median, r_upper``, a star count, and **asymmetric** errors on each
quantity. A single symmetric error column cannot represent these data, so the
schemas carry ``*_err_lo`` / ``*_err_hi`` and the bin edges: the likelihood
must integrate the model over the same bin, and a symmetric approximation is a
modelling choice to be made explicitly later, not silently at ingestion.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from ..paths import configs_dir, raw_dir
from ._product import build_product
from .base import read_table, schema_report
from .schema import ColumnRole, SchemaError, TableSchema

__all__ = [
    "DATASET",
    "PRODUCT_FILES",
    "PROFILE_SCHEMAS",
    "inventory",
    "describe",
    "product_paths",
    "load_profile",
    "preprocess",
]

DATASET = "omegacat_vi_kinematics"

LIKELIHOOD_RULE = (
    "binned profile: compare a model integrated over the same radial bin and "
    "tracer selection; bins share stars with other binnings of the same data, "
    "so alternative binnings are not independent likelihood factors"
)

#: File names documented in Appendix A of arXiv:2503.04903v2 (verified
#: 2026-09-16 from the paper's own HTML; NOT from the Zenodo record, which was
#: unreachable). Column names inside these files remain unknown.
PRODUCT_FILES: dict[str, str] = {
    "pm_log_bins": "proper_motion_dispersion_log_bins.fits",
    "pm_lin_bins": "proper_motion_dispersion_lin_bins.fits",
    "pm_equaln_bins": "proper_motion_dispersion_equaln_bins.fits",
    "pm_voronoi": "proper_motion_dispersion_voronoi_bins.fits",
    "los_profile": "los_profile.fits",
    "los_voronoi": "los_dispersion_and_rotation_voronoi_bins.fits",
    "equipartition_profiles": "energy_equipartition_profiles.fits",
    "equipartition_massbins": "energy_equipartition_massbins.fits",
    "catalog_and_selections": "catalog_and_selections.fits",
}

_BIN_ROLES = (
    ColumnRole("r_lower", "arcsec", (), valid_range=(0.0, None),
               description="inner edge of the radial bin"),
    ColumnRole("r_median", "arcsec", (), valid_range=(0.0, None),
               description="median radius of the stars in the bin"),
    ColumnRole("r_upper", "arcsec", (), valid_range=(0.0, None),
               description="outer edge of the radial bin"),
    ColumnRole("n_stars", None, (), kind="integer", required=False,
               valid_range=(1, None), description="number of tracers in the bin"),
)


def _asymmetric(value: str, unit: str | None, description: str) -> tuple[ColumnRole, ...]:
    """Return the value role plus its lower and upper error roles."""
    return (
        ColumnRole(
            value, unit, (),
            # A dispersion is non-negative; a rotation amplitude or angle is not,
            # so only dispersions carry the bound.
            valid_range=(0.0, None) if value.startswith("sigma") else None,
            description=description,
        ),
        ColumnRole(f"{value}_err_lo", unit, (), valid_range=(0.0, None),
                   description=f"lower (minus) uncertainty on {value}"),
        ColumnRole(f"{value}_err_hi", unit, (), valid_range=(0.0, None),
                   description=f"upper (plus) uncertainty on {value}"),
    )


def _profile_schema(name: str, value_roles: tuple[ColumnRole, ...], notes: str) -> TableSchema:
    """Build a binned-profile schema: bin geometry, counts, asymmetric errors.

    The bin edges must satisfy ``r_lower <= r_median <= r_upper`` and must not
    overlap: a profile whose geometry is inconsistent cannot be compared with a
    model integrated over the same bins.
    """
    return TableSchema(
        name=name,
        notes=notes,
        roles=_BIN_ROLES + value_roles,
        ordered_triples=(("r_lower", "r_median", "r_upper"),),
    )


#: Profile products the kinematic likelihood expects. Candidate name lists are
#: empty on purpose: every role must be mapped explicitly once the real column
#: names are known.
PROFILE_SCHEMAS: dict[str, TableSchema] = {
    "pm_radial": _profile_schema(
        "omegacat_vi_pm_radial",
        _asymmetric("sigma_pmr", "mas / yr", "radial proper-motion dispersion"),
        "Radial PM dispersion profile (oMEGACat VI, Table 2 structure).",
    ),
    "pm_tangential": _profile_schema(
        "omegacat_vi_pm_tangential",
        _asymmetric("sigma_pmt", "mas / yr", "tangential proper-motion dispersion"),
        "Tangential PM dispersion profile (oMEGACat VI, Table 2 structure).",
    ),
    "pm_combined": _profile_schema(
        "omegacat_vi_pm_combined",
        _asymmetric("sigma_pmc", "mas / yr", "combined proper-motion dispersion"),
        "Combined PM dispersion profile (oMEGACat VI, Table 2 structure).",
    ),
    "los_dispersion": _profile_schema(
        "omegacat_vi_los_dispersion",
        _asymmetric("sigma_los", "km / s", "line-of-sight velocity dispersion"),
        "LOS dispersion profile (oMEGACat VI, Table 3 structure).",
    ),
    "los_rotation": _profile_schema(
        "omegacat_vi_los_rotation",
        _asymmetric("v_rot", "km / s", "line-of-sight rotation amplitude")
        + _asymmetric("theta_0", "deg", "position angle of the rotation axis"),
        "LOS rotation profile with rotation-axis angle (oMEGACat VI, Table 3).",
    ),
}


def inventory() -> list[dict[str, Any]]:
    """List the raw oMEGACat VI files present locally.

    Returns
    -------
    list of dict
        ``name``, ``bytes``, ``suffix`` and ``documented`` (whether the name
        appears in :data:`PRODUCT_FILES`) for each downloaded file.
    """
    directory = raw_dir() / DATASET
    if not directory.is_dir():
        return []
    documented = set(PRODUCT_FILES.values())
    return [
        {
            "name": p.name,
            "bytes": p.stat().st_size,
            "suffix": p.suffix.lower(),
            "documented": p.name in documented,
        }
        for p in sorted(directory.iterdir())
        if p.is_file()
    ]


def describe(name: str, hdu: int | str | None = None) -> dict[str, Any]:
    """Dump the real column structure of one raw oMEGACat file.

    Parameters
    ----------
    name : str
        File name under ``data/raw/omegacat_vi_kinematics/``.
    hdu : int or str, optional
        FITS extension.

    Returns
    -------
    dict
        Schema report with no schema validation, since the mapping is unknown.
    """
    table = read_table(raw_dir() / DATASET / name, hdu=hdu)
    report = schema_report(table)
    report["file"] = name
    return report


def product_paths() -> dict[str, dict[str, Any]]:
    """Return the configured ``profile -> {file, hdu}`` mapping.

    Read from ``configs/data.yaml`` under ``omegacat_vi.products``. Empty when
    the config does not yet name the files.
    """
    path = configs_dir() / "data.yaml"
    if not path.is_file():
        return {}
    with open(path, "r", encoding="utf-8") as fh:
        doc = yaml.safe_load(fh) or {}
    products = ((doc.get("omegacat_vi") or {}).get("products")) or {}
    return {str(k): dict(v) for k, v in products.items() if v}


def load_profile(kind: str, *, validate: bool = True) -> Any:
    """Load one standardized oMEGACat VI profile.

    Parameters
    ----------
    kind : str
        One of :data:`PROFILE_SCHEMAS`.

    Returns
    -------
    astropy.table.Table

    Raises
    ------
    KeyError
        If ``kind`` is not a known profile.
    SchemaError
        If ``configs/data.yaml`` does not yet say which file holds the profile.
    """
    if kind not in PROFILE_SCHEMAS:
        raise KeyError(f"unknown profile {kind!r}; known: {sorted(PROFILE_SCHEMAS)}")

    products = product_paths()
    if kind not in products:
        raise SchemaError(
            f"configs/data.yaml does not name the file for profile {kind!r}. "
            f"Appendix A of arXiv:2503.04903v2 documents the deposited names "
            f"(see PRODUCT_FILES), but not the columns inside them. Run "
            "`ocen inspect-omegacat --columns` and record both the file and its "
            "column map. Do not guess."
        )

    entry = products[kind]
    from .base import standardize

    table = read_table(raw_dir() / DATASET / entry["file"], hdu=entry.get("hdu"))
    out = standardize(table, PROFILE_SCHEMAS[kind], keep_extra_columns=True)
    out.meta["ocen_dataset"] = DATASET
    out.meta["ocen_profile"] = kind
    out.meta["ocen_likelihood_rule"] = LIKELIHOOD_RULE
    if validate:
        from .schema import validate_table

        validate_table(out, PROFILE_SCHEMAS[kind])
    return out


def preprocess() -> list[dict[str, Any]]:
    """Standardize every configured profile and write the processed products.

    Returns
    -------
    list of dict
        One report per profile. Profiles that are not yet configured are
        reported with a ``'skipped'`` reason; the caller treats a skipped
        required product as a failure, not a success.
    """
    products = product_paths()
    reports: list[dict[str, Any]] = []
    for kind, schema in PROFILE_SCHEMAS.items():
        if kind not in products:
            reports.append(
                {
                    "dataset": DATASET,
                    "profile": kind,
                    "skipped": (
                        f"configs/data.yaml does not name the file for profile "
                        f"{kind!r}; run `ocen inspect-omegacat --columns` first"
                    ),
                }
            )
            continue
        entry = products[kind]
        try:
            report = build_product(
                dataset=DATASET,
                schema=schema,
                path=raw_dir() / DATASET / entry["file"],
                hdu=entry.get("hdu"),
                subdir="kinematics",
                likelihood_rule=LIKELIHOOD_RULE,
                extra_meta={"ocen_profile": kind},
            )
        except (SchemaError, FileNotFoundError) as exc:
            reports.append({"dataset": DATASET, "profile": kind, "error": str(exc)})
            continue
        report["profile"] = kind
        reports.append(report)
    return reports
