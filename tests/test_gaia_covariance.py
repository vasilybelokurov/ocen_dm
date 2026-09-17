"""Gaia covariance product: verification logic tested with a mocked WSDB join."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table

from ocen_dm.selection import gaia_covariance as gc


def _periphery(isolated_root, n=30, rng=np.random.default_rng(5)):
    t = Table()
    t["source_id"] = rng.permutation(np.arange(n, dtype=np.int64) + 6_000_000_000_000_000_000)
    t["pmra"] = rng.normal(-3.2, 0.5, n) * u.mas / u.yr
    t["pmdec"] = rng.normal(-6.7, 0.5, n) * u.mas / u.yr
    path = isolated_root / "data" / "processed" / "tails" / "kuzma2025_periphery.ecsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    t.write(path)
    return t, path


def _fake_join(table, corr=0.1, pm_offset=0.0, drop_last=False):
    ids = np.sort(np.asarray(table["source_id"], dtype=np.int64))
    order = np.argsort(np.asarray(table["source_id"]))
    n = len(ids)
    pmra = np.asarray(table["pmra"])[order] + pm_offset
    if drop_last:
        pmra = pmra.copy(); pmra[-1] = np.nan
    return {
        "source_id": ids, "pmra": pmra, "pmdec": np.asarray(table["pmdec"])[order],
        "pmra_error": np.full(n, 0.05), "pmdec_error": np.full(n, 0.05),
        "pmra_pmdec_corr": np.full(n, corr), "parallax": np.full(n, 0.18),
        "parallax_error": np.full(n, 0.03), "parallax_pmra_corr": np.zeros(n),
        "parallax_pmdec_corr": np.zeros(n), "ruwe": np.full(n, 1.0),
        "phot_g_mean_mag": np.full(n, 15.0), "astrometric_params_solved": np.full(n, 31),
    }


def test_product_is_written_and_checks_recorded(isolated_root, monkeypatch):
    table, path = _periphery(isolated_root)
    monkeypatch.setattr(gc, "fetch_gaia_covariance", lambda ids: _fake_join(table))
    report = gc.build_product(path)
    out = Table.read(report["path"])
    assert len(out) == len(table)
    assert out["pmra"].unit == u.mas / u.yr
    assert out.meta["ocen_checks"]["matched_fraction"] == 1.0
    assert "pmra_pmdec_corr" in out.colnames
    assert np.array_equal(np.asarray(out["source_id"]), np.sort(np.asarray(table["source_id"])))


def test_unmatched_star_is_refused(isolated_root, monkeypatch):
    table, path = _periphery(isolated_root)
    monkeypatch.setattr(gc, "fetch_gaia_covariance", lambda ids: _fake_join(table, drop_last=True))
    with pytest.raises(RuntimeError, match="no Gaia DR3 row"):
        gc.build_product(path)


def test_different_release_is_refused(isolated_root, monkeypatch):
    """Catalogue PMs that differ from Gaia's beyond rounding mean a different release."""
    table, path = _periphery(isolated_root)
    monkeypatch.setattr(gc, "fetch_gaia_covariance", lambda ids: _fake_join(table, pm_offset=0.05))
    with pytest.raises(RuntimeError, match="not the same release"):
        gc.build_product(path)


def test_out_of_range_correlation_is_refused(isolated_root, monkeypatch):
    table, path = _periphery(isolated_root)
    monkeypatch.setattr(gc, "fetch_gaia_covariance", lambda ids: _fake_join(table, corr=1.5))
    with pytest.raises(RuntimeError, match="outside"):
        gc.build_product(path)


def test_query_is_a_left_join_on_source_id():
    assert "LEFT JOIN gaia_dr3.gaia_source" in gc.QUERY
    assert "pmra_pmdec_corr" in gc.QUERY
