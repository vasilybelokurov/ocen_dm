"""End-to-end smoke test of ``ocen fit`` on mock data with a deliberately tiny sampler budget."""

from __future__ import annotations

import json

import numpy as np
import pytest

from ocen_dm.cli import main
from ocen_dm.paths import processed_dir, results_dir

_HAS_DATA = (processed_dir() / "kinematics" / "omegacat_vi_pm_radial.ecsv").exists()


@pytest.mark.skipif(not _HAS_DATA, reason="processed kinematics not present")
def test_fit_cli_writes_all_products(tmp_path, monkeypatch):
    from ocen_dm.kinematics import NoDarkMatterModel
    fam = NoDarkMatterModel(tracer="composite")      # the CLI's default tracer, not the class's
    x = fam.transform(np.full(len(fam.names), 0.5))
    np.save(tmp_path / "x.npy", x)
    monkeypatch.setattr("ocen_dm.paths.results_dir", lambda: tmp_path / "results")   # results -> tmp, data untouched
    rc = main(["fit", "--family", "K1", "--datasets", "hst_pm_radial,muse_los_dispersion",
               "--mock-from", str(tmp_path / "x.npy"), "--n-live", "40", "--max-ncalls", "1200",
               "--label", "smoke", "--seed", "3"])
    assert rc == 0
    out = tmp_path / "results" / "fits" / "smoke"
    for f in ("posterior.ecsv", "summary.json", "profiles.npz", "run.yaml", "ml_x.npy"):
        assert (out / f).exists(), f
    summary = json.loads((out / "summary.json").read_text())
    assert summary["family"] == "K1_noDM_composite" and set(summary["parameters"]) == set(fam.names)   # default tracer
    # the run must record that it was fitted to mock data, and the report must rebuild
    # that exact realisation rather than drawing the model against the real profiles
    prov = summary["data"]
    assert prov["kind"] == "mock" and prov["seed"] == 3 and prov["generating_family"] == "K1_noDM_composite"
    assert prov["truth"] == pytest.approx({n: float(v) for n, v in zip(fam.names, x)})
    from ocen_dm.kinematics import FitProblem
    from ocen_dm.kinematics.report import data_for
    mock, prov2 = data_for(summary)
    real, _ = data_for({**summary, "data": {"kind": "real"}})
    assert prov2["kind"] == "mock"
    by_mock = {p.name: p for p in mock.profiles}; by_real = {p.name: p for p in real.profiles}
    for name in by_mock:
        assert not np.allclose(by_mock[name].value, by_real[name].value)      # a different realisation
        np.testing.assert_allclose(by_mock[name].err_lo, by_real[name].err_lo)  # same errors and bins
        np.testing.assert_allclose(by_mock[name].r, by_real[name].r)
    # rebuilt exactly: the maximum-likelihood chi2 matches what the run itself recorded
    x_ml = np.array([summary["parameters"][n]["ml"] for n in fam.names])
    chi2 = sum(c for c, _ in FitProblem(fam, mock).chi2(x_ml).values())
    assert chi2 == pytest.approx(summary["chi2_ml_total"], abs=0.01)   # summary stores it rounded to 2 dp
    assert summary["n_points"] == 69                       # 40 + 29 bins
    prof = np.load(out / "profiles.npz")
    assert prof["M_total"].shape[1] == 60 and np.all(np.isfinite(prof["v_circ"]))
    # the report machinery runs on it
    from ocen_dm.kinematics.report import comparison_table, load_run, write_report
    run = load_run("smoke")
    table = comparison_table([run, run])
    assert "ln Z" in table and "M_star" in table and "χ² hst_pm_radial" in table
    report = write_report(["smoke"], path=tmp_path / "cmp.md", plots_dir=tmp_path / "plots")
    assert report.exists() and (tmp_path / "plots" / "fit_smoke_posterior_profiles.png").exists()


def test_class_and_cli_tracer_defaults_differ_deliberately():
    """Guard against silent drift: the class default is the light profile, the CLI's is the
    star-count composite. Reports must take the tracer from the run label, never from a
    freshly constructed family (this mismatch produced a wrong mock comparison on 2026-09-18)."""
    from ocen_dm.kinematics import NoDarkMatterModel
    from ocen_dm.kinematics.report import _family_for
    assert NoDarkMatterModel().tracer == "trager"
    assert _family_for({"family": "K1_noDM_composite"}).tracer == "composite"
    assert _family_for({"family": "K1_noDM"}).tracer == "trager"
    assert len(NoDarkMatterModel(tracer="trager").mge_fit.sigmas_arcsec) != len(
        NoDarkMatterModel(tracer="composite").mge_fit.sigmas_arcsec)


def test_fit_cli_rejects_unknown_family():
    with pytest.raises(SystemExit):
        main(["fit", "--family", "K9"])


@pytest.mark.skipif(not _HAS_DATA, reason="processed kinematics not present")
def test_fit_cli_preset_runs_and_prints_comparison(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("ocen_dm.paths.results_dir", lambda: tmp_path / "results")
    rc = main(["fit", "--preset", "watkins2013", "--n-live", "30", "--max-ncalls", "600", "--seed", "1"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "preset watkins2013" in out and "like-for-like comparison" in out
    assert "ML_V" in out and "published 2.71" in out and "beta_0" in out
    assert (tmp_path / "results" / "fits" / "preset_watkins2013" / "summary.json").exists()
