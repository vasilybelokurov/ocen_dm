"""Loaders for datasets kept in the canonical catalogue store (kind: local).

Synthetic fixtures only: a multi-cluster FITS with the real Vasiliev column names,
and a zip with a 12-column profile table.
"""

from __future__ import annotations

import io
import zipfile

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table

from ocen_dm.data import baumgardt_catalogue, vasiliev2021
from ocen_dm.data._product import build_product
from ocen_dm.data.schema import SchemaError, validate_table
from ocen_dm.provenance import DatasetRecord, FileRecord, load_registry


def _multi_cluster_fits(path, n_ocen=40, rng=np.random.default_rng(1)):
    """Two clusters in one file, real column names, no units (as in the real FITS)."""
    n = 2 * n_ocen
    t = Table()
    t["SOURCE_ID"] = np.arange(n, dtype=np.int64) + 6_080_000_000_000_000_000
    t["RA"] = 201.7 + rng.normal(0, 0.2, n)
    t["DEC"] = -47.5 + rng.normal(0, 0.2, n)
    t["X"] = rng.normal(0, 0.2, n); t["Y"] = rng.normal(0, 0.2, n)
    t["PLX"] = rng.normal(0.18, 0.1, n)
    t["PMRA"] = rng.normal(-3.2, 0.5, n); t["PMDEC"] = rng.normal(-6.7, 0.5, n)
    t["PLXE"] = np.full(n, 0.1); t["PMRAE"] = np.full(n, 0.1); t["PMDECE"] = np.full(n, 0.1)
    t["PMCORR"] = rng.uniform(-0.5, 0.5, n)
    t["G_MAG"] = rng.uniform(12, 21, n); t["BP_RP"] = rng.uniform(0.5, 1.5, n)
    t["SIGMA"] = np.full(n, 13.0); t["GFLAG"] = np.full(n, 2, dtype=np.int64)
    t["PROB"] = rng.uniform(0, 1, n)
    t["NAME"] = np.array(["NGC_5139_oCen  "] * n_ocen + ["NGC_104_47Tuc  "] * n_ocen)
    t.write(path)
    return t


def _profile_zip(path, ncol=12):
    r = np.linspace(0, 0.66, 100)
    cols = [r] + [0.5 - 0.3 * r + 0.02 * k for k in range(5)] + [0.05 * r * (1 + 0.1 * k) for k in range(5)]
    if ncol == 12:
        cols.append(0.07 * r)
    data = np.column_stack(cols)
    buf = io.StringIO()
    np.savetxt(buf, data, fmt="%.6g")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("!readme.txt", "synthetic")
        z.writestr("profiles/NGC_5139_oCen.txt", buf.getvalue())
        z.writestr("catalogues/NGC_5139_oCen.txt", "# x\n")


@pytest.fixture
def store(isolated_root, registry_copy, tmp_path, monkeypatch):
    """Point the registry's local paths at synthetic files."""
    fits = tmp_path / "gc_members.fits"
    zipf = tmp_path / "clusters.zip"
    cat = tmp_path / "gc_catalog.fits"
    _multi_cluster_fits(fits)
    _profile_zip(zipf)
    c = Table()
    c["NAME"] = np.array(["NGC_5139_oCen  ", "NGC_104_47Tuc  "])
    c["RA"] = [201.697, 6.02]; c["DEC"] = [-47.48, -72.08]
    c["PMRA"] = [-3.236, 5.25]; c["PMDEC"] = [-6.731, -2.53]
    c["PMRA_ERR"] = [0.011, 0.01]; c["PMDEC_ERR"] = [0.011, 0.01]
    c["DIST"] = [5.43, 4.52]; c["DIST_ERR"] = [0.05, 0.03]
    c["RV"] = [232.78, -17.2]; c["RV_ERR"] = [0.21, 0.2]
    c["MASS"] = [3.94e6, 8.9e5]; c["RH"] = [7.56, 6.3]; c["RC"] = [4.3, 0.5]
    c["RPERI"] = [1.3, 5.5]; c["RAPO"] = [6.8, 7.4]; c["FEH"] = [-1.53, -0.72]
    c.write(cat)

    reg = isolated_root / "provenance" / "datasets.yaml"
    text = reg.read_text()
    text = text.replace("~/data/catalogues/gc_members_gaia_vasiliev.fits", str(fits))
    text = text.replace("~/Work/Code/apogee_halo_rotation/data/external/vasiliev_baumgardt_2021_clusters.zip", str(zipf))
    text = text.replace("~/data/catalogues/gc_catalog_updated.fits", str(cat))
    # The registry's expected sizes and the zip's published md5 describe the REAL
    # files; the synthetic stand-ins must not be held to them here. Size and
    # checksum enforcement have their own tests.
    import re
    text = re.sub(r"expected_bytes: \d+", "expected_bytes: null", text)
    text = text.replace('published_md5: "3a07093be9a1fcbfb0e3c545ca175c99"', "published_md5: null")
    text = text.replace('checksum_source: "Zenodo API record 4891252, files[].checksum, read 2026-09-16"',
                        "checksum_source: null")
    reg.write_text(text)
    # the real configs are needed for the declared units
    import shutil
    from pathlib import Path
    import ocen_dm
    repo = Path(ocen_dm.__file__).resolve().parents[2]
    shutil.copy(repo / "configs" / "column_maps.yaml", isolated_root / "configs" / "column_maps.yaml")
    return {"fits": fits, "zip": zipf, "cat": cat}


def test_local_kind_files_are_verified_in_place(store):
    from ocen_dm.data import download as dl

    results = dl.fetch_all(only=["vasiliev2021_gc_members", "baumgardt_gc_catalogue"],
                           log=lambda _: None)
    assert [r.status for r in results] == ["cached", "cached", "cached"]
    assert all("canonical store" in r.message for r in results)
    assert all(r.verification["sha256"] for r in results)     # identity recorded


def test_local_file_with_wrong_size_fails_verification(store, isolated_root):
    from ocen_dm.data import download as dl

    reg = isolated_root / "provenance" / "datasets.yaml"
    # The fixture nulled every expected size; re-arm only the gc_catalog file's,
    # located by its name rather than by position in the registry.
    import re
    text = reg.read_text()
    text = re.sub(r'(name: "gc_catalog_updated.fits"(?:(?!name:).)*?expected_bytes: )null',
                  r"\g<1>12345", text, count=1, flags=re.S)
    assert "expected_bytes: 12345" in text
    reg.write_text(text)
    results = dl.fetch_all(only=["baumgardt_gc_catalogue"], log=lambda _: None)
    assert results[0].status == "failed"


def test_local_file_missing_is_a_failure(store):
    from ocen_dm.data import download as dl

    store["cat"].unlink()
    results = dl.fetch_all(only=["baumgardt_gc_catalogue"], log=lambda _: None)
    assert results[0].status == "failed" and "canonical location" in results[0].message


def test_local_files_are_never_copied_into_the_repo(store, isolated_root):
    from ocen_dm.data import download as dl

    dl.fetch_all(only=["vasiliev2021_gc_members"], log=lambda _: None)
    assert not (isolated_root / "data" / "raw" / "vasiliev2021_gc_members").exists()


def test_vasiliev_loader_selects_only_omega_cen(store):
    table = vasiliev2021.load()
    assert len(table) == 40
    assert table["source_id"].dtype.kind == "i"
    assert table["pmra"].unit == u.mas / u.yr
    assert np.all((table["pmra_pmdec_corr"] >= -1) & (table["pmra_pmdec_corr"] <= 1))
    validate_table(table, vasiliev2021.SCHEMA)


def test_vasiliev_preprocess_records_the_row_selection(store):
    reports = vasiliev2021.preprocess()
    members = reports[0]
    assert members["n_rows"] == 40
    assert members["lineage"]["rows_in_file"] == 80
    assert members["lineage"]["rows_selected"] == 40


def test_vasiliev_profiles_keep_the_undocumented_column_flagged(store):
    raw = vasiliev2021.load_profiles()
    assert raw.colnames[-1] == "c12"
    reports = vasiliev2021.preprocess()
    prof = Table.read(reports[1]["path"])
    assert "undocumented_col12" in prof.colnames
    assert prof["r"].unit == u.arcsec            # converted from the file's degrees
    assert np.isclose(prof["r"][-1], 0.66 * 3600, rtol=1e-6)
    assert prof["sigma_pm"].unit == u.mas / u.yr


def test_vasiliev_profiles_with_eleven_columns_also_load(store, tmp_path):
    _profile_zip(store["zip"], ncol=11)
    raw = vasiliev2021.load_profiles()
    assert "c12" not in raw.colnames


def test_vasiliev_profiles_reject_an_unexpected_shape(store):
    with zipfile.ZipFile(store["zip"], "w") as z:
        z.writestr("profiles/NGC_5139_oCen.txt", "0 1 2\n1 2 3\n")
    with pytest.raises(SchemaError, match="documents 11 columns"):
        vasiliev2021.load_profiles()


def test_baumgardt_loader_returns_one_row_with_units(store):
    table = baumgardt_catalogue.load()
    assert len(table) == 1
    assert table["distance"].unit == u.kpc
    assert float(table["distance"][0]) == pytest.approx(5.43)
    assert table["mass"].unit == u.Msun
    assert "5139" in str(table["name"][0])


def test_baumgardt_product_is_flagged_unverified(store):
    report = baumgardt_catalogue.preprocess()
    meta = Table.read(report["path"]).meta
    assert "UNVERIFIED" in meta["ocen_likelihood_rule"]


def test_registry_carries_the_zenodo_checksum_for_the_zip():
    record = [f for f in load_registry()["vasiliev2021_gc_members"].files if "zip" in f.name][0]
    assert record.published_md5 == "3a07093be9a1fcbfb0e3c545ca175c99"
    assert record.local_path is not None


def test_row_filter_lineage_in_build_product(isolated_root, synthetic_periphery):
    from ocen_dm.data import kuzma2025

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    synthetic_periphery.write(directory / "wCen_table.fits")
    report = build_product(
        dataset=kuzma2025.DATASET, schema=kuzma2025.SCHEMA,
        path=directory / "wCen_table.fits", hdu=1, subdir="tails", likelihood_rule="t",
        row_filter=lambda t: np.asarray(t["ra"]) > np.median(np.asarray(t["ra"])),
    )
    assert report["lineage"]["rows_in_file"] == len(synthetic_periphery)
    assert report["lineage"]["rows_selected"] == report["n_rows"]


def test_transform_is_recorded_in_lineage(isolated_root, synthetic_periphery):
    from ocen_dm.data import kuzma2025

    directory = isolated_root / "data" / "raw" / kuzma2025.DATASET
    directory.mkdir(parents=True)
    synthetic_periphery.write(directory / "wCen_table.fits")

    def shift(table):
        table = table.copy(); table["ra"] = table["ra"] + 0 * u.deg; return table

    report = build_product(
        dataset=kuzma2025.DATASET, schema=kuzma2025.SCHEMA, path=directory / "wCen_table.fits",
        hdu=1, subdir="tails", likelihood_rule="t", transform=shift, transform_note="identity shift",
    )
    assert report["lineage"]["transform"] == "identity shift"
    assert Table.read(report["path"]).meta["ocen_raw"]["transform"] == "identity shift"
