"""Regression tests for saved observations, mock identity and preservation of runs."""

from copy import deepcopy
from dataclasses import replace
import json
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.table import Table
import yaml

from ocen_dm.cli import main
from ocen_dm.kinematics.fit import DarkMatterModel, FitProblem, NoDarkMatterModel, Prior, run_nested
from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData
from ocen_dm.kinematics.report import (_family_for, comparison_table, data_for, load_run_metadata,
                                       problem_for)
from ocen_dm.kinematics.run_io import (data_fingerprint, family_config, family_from_config,
                                      read_data_snapshot, write_data_snapshot)
from ocen_dm.light_model import MGEFit

MGE = MGEFit(np.array([50., 250.]), np.array([0.7, 0.3]), 0., 0., 2)


@pytest.fixture(autouse=True)
def no_tracer_files(monkeypatch):
    monkeypatch.setattr("ocen_dm.kinematics.fit.load_tracer_profile",
                        lambda *a: SimpleNamespace(r_arcsec=np.array([1., 300.])))
    monkeypatch.setattr("ocen_dm.kinematics.fit.fit_mge_projected", lambda *a, **k: MGE)


@pytest.fixture
def data():
    r = np.array([50., 100.])
    p = BinnedProfile("gaia_edr3_ours_tangential", "pmt", r, r-10, r+10,
                      np.array([0.4, 0.3]), np.array([0.02, 0.03]), np.array([0.04, 0.05]),
                      "GaiaEDR3", streaming2=np.array([0.002, 0.003]),
                      r_nodes=np.column_stack((r-5, r+5)))
    q = BinnedProfile("muse_los_dispersion", "los", r, None, None,
                      np.array([18., 16.]), np.ones(2), np.ones(2)*2, "MUSE")
    return KinematicData((p, q))


def baseline(**kw):
    return NoDarkMatterModel(mge_fit=MGE, instruments=(), constant_beta=True,
                             fixed={"beta_0": 0.0}, tracer="composite", **kw)


def summary_for(family, data):
    x = family.transform(np.full(len(family.names), 0.5))
    problem = FitProblem(family, data)
    return {"family": family.label, "model": family_config(family),
            "parameters": {n: {"ml": float(v)} for n, v in zip(family.names, x)},
            "datasets": {p.name: {} for p in data.profiles}, "data": {"kind": "real"},
            "lnL_max": problem.loglike_vector(x)}


def test_snapshot_replays_every_array_without_current_products(data, tmp_path, monkeypatch):
    record = write_data_snapshot(data, tmp_path)
    monkeypatch.setattr(KinematicData, "load", lambda *a, **k: pytest.fail("read live products"))
    restored, _ = data_for({"data_snapshot": record}, run_dir=tmp_path)
    for before, after in zip(data.profiles, restored.profiles):
        for name in ("r", "r_lower", "r_upper", "value", "err_lo", "err_hi", "streaming2", "r_nodes"):
            np.testing.assert_array_equal(getattr(before, name), getattr(after, name))
    path = tmp_path / record["file"]
    path.write_text(path.read_text() + " ")
    with pytest.raises(ValueError, match="checksum"):
        read_data_snapshot(tmp_path, record)


@pytest.mark.parametrize("field", ["r", "r_lower", "r_upper", "value", "err_lo", "err_hi", "streaming2", "r_nodes"])
def test_identity_includes_every_likelihood_array(data, field):
    p, q = data.profiles
    changed = replace(p, **{field: getattr(p, field) * 1.01})
    assert data_fingerprint(data) != data_fingerprint(KinematicData((changed, q)))
    assert data_fingerprint(data) == data_fingerprint(KinematicData((q, replace(p, note="new prose"))))


@pytest.mark.parametrize("backend", ["jeans", "jam", "agama"])
@pytest.mark.parametrize("halo", [False, True])
def test_resolved_family_roundtrip_retains_configuration(backend, halo):
    cls = DarkMatterModel if halo else NoDarkMatterModel
    kw = dict(gamma=0.7, r_t=75.) if halo else {}
    fam = cls(mge_fit=MGE, instruments=("MUSE",), backend=backend, tracer="composite",
              fixed={"M_bh": 1234.}, distance_prior=Prior("uniform", 4.1, 6.2),
              distance_kpc=5.1, beta0_max=-0.6, **kw)
    config = json.loads(json.dumps(family_config(fam)))
    restored = family_from_config(config)
    assert family_config(restored) == config
    np.testing.assert_array_equal(fam.transform(np.full(len(fam.names), 0.3)),
                                  restored.transform(np.full(len(fam.names), 0.3)))


@pytest.mark.parametrize("fixed_distance", [False, True])
def test_snapshot_model_replays_identical_likelihood(data, fixed_distance):
    fam = baseline(fix_distance=fixed_distance, distance_kpc=4.59)
    x = fam.transform(np.full(len(fam.names), 0.5))
    restored = _family_for({"model": family_config(fam)})
    assert FitProblem(restored, data).loglike_vector(x) == FitProblem(fam, data).loglike_vector(x)
    assert restored.fixed_distance_kpc == 4.59
    assert _family_for({"model": family_config(fam)}, backend="jam").backend == "jam"


def test_legacy_loader_honours_gaia_variants_and_radius_cut(data, monkeypatch):
    seen = {}
    def loader(names, **opts):
        seen.update(opts)
        return data
    monkeypatch.setattr(KinematicData, "load", loader)
    s = {"datasets": [p.name for p in data.profiles],
         "dataset_options": {"gaia_errors": "eta", "gaia_rotation": "ours", "gaia_r_min": 450.}}
    with pytest.warns(UserWarning, match="Legacy run"):
        data_for(s)
    assert seen["gaia_edr3_ours_radial"] == {"error_model": "eta", "rotation": "ours"}
    assert seen["gaia_edr3_ours_tangential"] == seen["gaia_edr3_ours_radial"]
    assert seen["gaia_edr3_pm"] == {"r_min_arcsec": 450.}


def test_legacy_family_reads_summary_ladder_switches():
    s = {"family": "K1_noDM_composite", "dataset_options": {"isotropic": True, "no_scales": True}}
    assert _family_for(s).names == baseline().names


def comparison_run():
    return {"label": "test", "posterior": Table(),
            "summary": {"datasets": {"x": {}}, "n_points": 1, "parameters": {},
                        "logz": 1., "logzerr": 0.1, "lnL_max": 2., "chi2_ml_total": 0.,
                        "chi2_ml": {"x": {"chi2": 0., "n": 1}}, "n_calls": 5, "elapsed_s": 1.},
            "run": {"inputs": {"x.ecsv": "abc"}, "dataset_options": {},
                    "data": {"kind": "mock", "source_file": "ml_x.npy", "generating_family": "K1",
                             "truth": {"M_star": 2e6}, "seed": 42}}}


@pytest.mark.parametrize("key,value", [("generating_family", "K2"), ("truth", {"M_star": 3e6}),
                                       ("seed", 43), ("source_file", "other.npy")])
def test_mock_comparison_rejects_different_generators_or_realisations(key, value):
    a = comparison_run()
    b = deepcopy(a)
    b["run"]["data"][key] = value
    assert "NOT comparable" in comparison_table([a, b])
    assert "NOT comparable" not in comparison_table([a, deepcopy(a)])


def test_comparison_rejects_missing_provenance_and_ignores_model_switches():
    a = comparison_run()
    b = deepcopy(a)
    b["run"]["dataset_options"] = {"isotropic": True, "no_scales": True, "backend": "jam"}
    assert "NOT comparable" not in comparison_table([a, b])
    b["run"]["dataset_options"]["gaia_rotation"] = "ours"
    assert "NOT comparable" in comparison_table([a, b])
    a.pop("run")
    assert "NOT comparable" in comparison_table([a, deepcopy(a)])


def test_new_comparison_uses_observations_not_unrelated_file_metadata(data, tmp_path):
    a = comparison_run()
    a["dir"] = tmp_path
    a["summary"]["data_snapshot"] = write_data_snapshot(data, tmp_path)
    b = deepcopy(a)
    b["run"]["inputs"] = {"unused.ecsv": "different"}
    b["run"]["dataset_options"] = {"constant_beta": True}
    assert "NOT comparable" not in comparison_table([a, b])
    b["dir"] = tmp_path / "different_data"
    b["dir"].mkdir()
    p, q = data.profiles
    changed = KinematicData((replace(p, value=p.value * 1.01), q))
    b["summary"]["data_snapshot"] = write_data_snapshot(changed, b["dir"])
    assert "NOT comparable" in comparison_table([a, b])


def test_reports_detect_replay_mismatch(data, tmp_path):
    s = summary_for(baseline(), data)
    s["data_snapshot"] = write_data_snapshot(data, tmp_path)
    run = {"label": "test", "summary": s, "dir": tmp_path}
    problem_for(run)
    s["lnL_max"] += 1
    with pytest.raises(ValueError, match="differs from saved"):
        problem_for(run)


@pytest.mark.parametrize("rung", [0, 1, 2])
@pytest.mark.parametrize("halo", [False, True])
def test_source_run_defines_mock_generator_independently_of_target(data, tmp_path, monkeypatch, rung, halo):
    source = tmp_path / "source"
    source.mkdir()
    cls = DarkMatterModel if halo else NoDarkMatterModel
    fam = cls(mge_fit=MGE, instruments=(), constant_beta=rung < 2,
              fixed={"beta_0": 0.0} if rung == 0 else {}, tracer="composite",
              **({"r_t": 87.} if halo else {}))
    s = summary_for(fam, data)
    s["data_snapshot"] = write_data_snapshot(data, source)
    (source / "summary.json").write_text(json.dumps(s))
    x = np.array([s["parameters"][n]["ml"] for n in fam.names])
    np.save(source / "ml_x.npy", x)
    monkeypatch.setattr(KinematicData, "load", lambda *a, **k: data)
    monkeypatch.setattr("ocen_dm.paths.results_dir", lambda: tmp_path / "results")
    captured = {}
    def driver(problem, out, **kw):
        captured.update(problem=problem, **kw)
        return {"logz": 0., "logzerr": 0., "chi2_ml_total": 0., "n_points": data.n_points,
                "n_calls": 1, "elapsed_s": 0., "parameters": {}}
    monkeypatch.setattr("ocen_dm.kinematics.run_nested", driver)
    target_rung = "--isotropic" if rung == 1 else "--constant-beta"
    assert main(["fit", "--family", "K2-cored", target_rung, "--no-scales",
                 "--mock-from", str(source), "--seed", "7"]) == 0
    assert captured["problem"].family.names != fam.names
    assert captured["data_provenance"]["model"] == family_config(fam)
    expected = FitProblem(fam, data).mock_data(x, np.random.default_rng(7))
    assert data_fingerprint(captured["problem"].data) == data_fingerprint(expected)


def test_cli_requires_metadata_or_explicit_bare_vector_model(data, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(KinematicData, "load", lambda *a, **k: data)
    np.save(tmp_path / "x.npy", np.ones(5))
    assert main(["fit", "--mock-from", str(tmp_path / "x.npy")]) == 2
    assert "need --mock-family" in capsys.readouterr().err


def test_existing_run_is_untouched_by_driver_and_cli(data, tmp_path, monkeypatch, capsys):
    out = tmp_path / "fits" / "keep"
    out.mkdir(parents=True)
    sentinel = out / "summary.json"
    sentinel.write_bytes(b"original posterior summary")
    with pytest.raises(FileExistsError, match="choose a new label"):
        run_nested(FitProblem(baseline(), data), out)
    monkeypatch.setattr(KinematicData, "load", lambda *a, **k: data)
    monkeypatch.setattr("ocen_dm.paths.results_dir", lambda: tmp_path)
    assert main(["fit", "--label", "keep"]) == 2
    assert "choose a new --label" in capsys.readouterr().err
    assert sentinel.read_bytes() == b"original posterior summary"
    assert list(out.iterdir()) == [sentinel]


@pytest.mark.parametrize("nsteps", [None, 48])
def test_driver_snapshots_before_sampling_and_retains_launch_state(data, tmp_path, monkeypatch, nsteps):
    out = tmp_path / "new"
    problem = FitProblem(baseline(fix_distance=True, distance_kpc=4.59), data)
    monkeypatch.setattr("ocen_dm.kinematics.run_io.code_state",
                        lambda: {"git_commit": "launch", "dirty": True})
    class Sampler:
        def __init__(self, names, loglike, transform, **kw):
            self.x = transform(np.full(len(names), 0.5))
            self.loglike = loglike
            launch = yaml.safe_load((out / "run.yaml").read_text())
            assert launch["status"] == "started"
            assert launch["distance_kpc"] == 4.59
            assert data_fingerprint(read_data_snapshot(out, launch["data_snapshot"])) == data_fingerprint(data)
            monkeypatch.setattr("ocen_dm.kinematics.run_io.code_state", lambda: pytest.fail("late code capture"))
        def run(self, **kw):
            assert self.stepsampler.nsteps == (nsteps or 2*len(problem.family.names))
            return {"samples": np.array([self.x, self.x]),
                    "weighted_samples": {"points": np.array([self.x]), "logl": [self.loglike(self.x)]},
                    "logz": 0., "logzerr": 0.1, "ncall": 1}
    monkeypatch.setattr("ultranest.ReactiveNestedSampler", Sampler)
    s = run_nested(problem, out, n_live=10, n_profile_samples=1, step_sampler=True,
                   slice_nsteps=nsteps,
                   dataset_options={"gaia_errors": "eta", "gaia_rotation": "ours", "gaia_r_min": 400.})
    meta = yaml.safe_load((out / "run.yaml").read_text())
    assert meta["status"] == "complete" and meta["git_commit"] == "launch"
    assert meta["sampler"]["step_sampler"] is True
    assert meta["sampler"]["slice_nsteps"] == (nsteps or 2*len(problem.family.names))
    assert meta["dataset_options"]["gaia_r_min"] == 400.
    run = load_run_metadata(out)
    assert problem_for(run).loglike_vector(np.load(out / "ml_x.npy")) == s["lnL_max"]


def test_interrupted_sampler_leaves_launch_snapshot(data, tmp_path, monkeypatch):
    out = tmp_path / "interrupted"
    def fail(*a, **k):
        raise RuntimeError("simulated interruption")
    monkeypatch.setattr("ultranest.ReactiveNestedSampler", fail)
    with pytest.raises(RuntimeError, match="simulated interruption"):
        run_nested(FitProblem(baseline(), data), out)
    assert yaml.safe_load((out / "run.yaml").read_text())["status"] == "started"
    assert (out / "data_snapshot.json").is_file()
    with pytest.raises(FileExistsError):
        run_nested(FitProblem(baseline(), data), out)
