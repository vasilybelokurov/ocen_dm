"""Loader tests using synthetic files written to a temporary project root.

No test touches the network or the real catalogues.
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table

from ocen_dm.data import fimbulthul, kuzma2025, kuzma2026, omegacat
from ocen_dm.data.base import find_raw, read_table, standardize, write_processed
from ocen_dm.data.schema import SchemaError, validate_table


def test_read_table_round_trips_fits_and_csv(tmp_path, synthetic_periphery):
    for suffix in (".fits", ".csv", ".ecsv"):
        path = tmp_path / f"t{suffix}"
        synthetic_periphery.write(path, overwrite=True)
        back = read_table(path)
        assert len(back) == len(synthetic_periphery)
        assert np.allclose(np.asarray(back["ra"]), np.asarray(synthetic_periphery["ra"]))


def test_read_table_missing_file_points_at_fetch_data(tmp_path):
    with pytest.raises(FileNotFoundError, match="fetch-data"):
        read_table(tmp_path / "nothing.fits")


def test_find_raw_rejects_ambiguous_matches(isolated_root, synthetic_periphery):
    directory = isolated_root / "data" / "raw" / "demo"
    directory.mkdir(parents=True)
    synthetic_periphery.write(directory / "kuzma2026_members.csv")
    synthetic_periphery.write(directory / "kuzma2026_members.fits")
    with pytest.raises(SchemaError, match="exactly one raw file"):
        find_raw("demo", "kuzma2026_members.*")


def test_kuzma2025_standardization_and_validation(isolated_root, synthetic_periphery):
    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    synthetic_periphery.write(directory / "wCen_table.fits")

    table = kuzma2025.load()
    assert set(["source_id", "ra", "dec", "pmra", "pmdec"]) <= set(table.colnames)
    assert table["membership_prob"] is not None  # resolved from the 'prob' column
    assert table.meta["ocen_source_columns"]["membership_prob"] == "prob"
    assert "diagnostic only" in table.meta["ocen_likelihood_rule"]

    report = validate_table(table, kuzma2025.SCHEMA)
    assert report["n_rows"] == len(synthetic_periphery)
    assert report["unique_id_ok"] is True


def test_kuzma2025_preprocess_writes_a_readable_product(isolated_root, synthetic_periphery):
    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    synthetic_periphery.write(directory / "wCen_table.fits")

    report = kuzma2025.preprocess()
    written = read_table(report["path"])
    assert len(written) == len(synthetic_periphery)
    assert written["ra"].unit == u.deg
    assert written.meta["ocen_schema"] == "kuzma2025_periphery"


def test_kuzma2026_requires_velocity_columns(isolated_root, synthetic_periphery):
    """A spectroscopic table without v_los must fail, not silently drop it."""
    directory = isolated_root / "data" / "raw" / kuzma2026.DATASET
    directory.mkdir(parents=True)
    synthetic_periphery.write(directory / "kuzma2026_members.ecsv")
    with pytest.raises(SchemaError, match="vlos"):
        kuzma2026.load()


def test_kuzma2026_loads_a_conforming_table(isolated_root, synthetic_periphery, rng):
    table = synthetic_periphery.copy()
    n = len(table)
    table["HRV"] = (232.0 + 10.0 * rng.standard_normal(n)) * u.km / u.s
    table["e_RV"] = np.full(n, 1.5) * u.km / u.s
    directory = isolated_root / "data" / "raw" / kuzma2026.DATASET
    directory.mkdir(parents=True)
    # ECSV, not CSV: plain CSV drops units, which the loader now refuses to
    # guess back (see test_f5_* in test_review_regressions.py).
    table.write(directory / "kuzma2026_members.ecsv")

    loaded = kuzma2026.load()
    assert loaded["vlos"].unit.is_equivalent(u.km / u.s)
    assert "condition on target positions" in loaded.meta["ocen_likelihood_rule"]
    validate_table(loaded, kuzma2026.SCHEMA)


def test_fimbulthul_track_only_rule_is_recorded(isolated_root, synthetic_vizier_fimbulthul):
    """Column names follow the real VizieR catalogue J/other/NatAs/3.667."""
    directory = isolated_root / "data" / "raw" / fimbulthul.DATASET
    directory.mkdir(parents=True)
    synthetic_vizier_fimbulthul.write(directory / "fimbulthul_supp_table1.ecsv")
    loaded = fimbulthul.load()
    assert "track constraint only" in loaded.meta["ocen_likelihood_rule"]
    assert loaded["name"][0].startswith("Gaia DR3 ")
    validate_table(loaded, fimbulthul.SCHEMA)


def test_fimbulthul_identifier_is_a_string_not_an_integer():
    """306 of 309 real rows carry a Gaia name; three carry CRTS/UCAC4 names."""
    role = fimbulthul.SCHEMA.role("name")
    assert role.kind == "string" and role.is_identifier


def test_omegacat_profile_requires_an_explicit_config(isolated_root, synthetic_profile):
    directory = isolated_root / "data" / "raw" / omegacat.DATASET
    directory.mkdir(parents=True)
    synthetic_profile.write(directory / "profiles.fits")

    assert omegacat.inventory()[0]["name"] == "profiles.fits"
    with pytest.raises(SchemaError, match="does not name the file"):
        omegacat.load_profile("pm_radial")


def test_omegacat_profile_loads_once_configured(isolated_root, synthetic_profile):
    directory = isolated_root / "data" / "raw" / omegacat.DATASET
    directory.mkdir(parents=True)
    synthetic_profile.write(directory / "profiles.ecsv")

    (isolated_root / "configs" / "data.yaml").write_text(
        "omegacat_vi:\n  products:\n    pm_radial: {file: profiles.ecsv}\n",
        encoding="utf-8",
    )
    (isolated_root / "configs" / "column_maps.yaml").write_text(
        "omegacat_vi_pm_radial:\n"
        "  r_lower: R_low\n"
        "  r_median: R_med\n"
        "  r_upper: R_up\n"
        "  sigma_pmr: sigma_R\n"
        "  sigma_pmr_err_lo: sigma_R_lo\n"
        "  sigma_pmr_err_hi: sigma_R_hi\n"
        "  n_stars: N\n",
        encoding="utf-8",
    )
    table = omegacat.load_profile("pm_radial")
    assert table["r_median"].unit == u.arcsec
    assert table["sigma_pmr"].unit == u.mas / u.yr
    assert table["sigma_pmr_err_lo"].unit == u.mas / u.yr
    validate_table(table, omegacat.PROFILE_SCHEMAS["pm_radial"])


def test_omegacat_preprocess_skips_unconfigured_profiles(isolated_root):
    reports = omegacat.preprocess()
    assert len(reports) == len(omegacat.PROFILE_SCHEMAS)
    assert all("skipped" in r for r in reports)


def test_omegacat_product_filenames_come_from_the_paper():
    """Appendix A of arXiv:2503.04903v2, verified live on 2026-09-16."""
    assert omegacat.PRODUCT_FILES["los_profile"] == "los_profile.fits"
    assert (omegacat.PRODUCT_FILES["pm_log_bins"]
            == "proper_motion_dispersion_log_bins.fits")


def test_write_processed_preserves_units_and_metadata(isolated_root, synthetic_periphery):
    table = standardize(synthetic_periphery, kuzma2025.SCHEMA)
    path = write_processed(table, "demo_product")
    back = read_table(path)
    assert back["ra"].unit == u.deg
    assert back.meta["ocen_schema"] == "kuzma2025_periphery"
