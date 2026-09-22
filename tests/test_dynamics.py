"""Physical reference checks and end-to-end tests of controlled dynamics."""
from dataclasses import replace
import io
import os
from pathlib import Path

import numpy as np
import pytest
import yaml

from ocen_dm.dynamics import G, TIME_UNIT_MYR
from ocen_dm.dynamics.config import Component, Experiment, Integrator, Nucleus, Orbit, safe_label
from ocen_dm.dynamics.diagnostics import (bound_particles, effective_number,
    snapshot_diagnostics, spherical_particle_potential, tidal_scale, tidal_tensor, weighted_median)
from ocen_dm.dynamics.nemo import gyrfalcon_args, parse_snapshots
from ocen_dm.dynamics.orbits import Track, end_time, timestep


def config(**kwargs):
    component = Component(name="dm", scale_pc=20., cutoff_pc=100., n_particles=512,
                          softening_pc=.2, mass_within_msun=1e5, radius_pc=20.)
    return Experiment(programme="remnant", components=(component,),
                      orbit=Orbit(duration_myr=.1, sample_step_myr=.01),
                      integrator=Integrator(max_step_myr=.001, levels=2, output_step_myr=.05),
                      **kwargs)


def test_config_round_trip_and_strict_unknown_options():
    original = config()
    assert Experiment.from_dict(yaml.safe_load(yaml.safe_dump(original.to_dict()))) == original
    broken = original.to_dict()
    broken["orbit"]["friction"] = True
    with pytest.raises(ValueError, match="Unknown"):
        Experiment.from_dict(broken)


def test_cli_exposes_separate_preparation_and_evolution():
    from ocen_dm.cli import build_parser
    parser = build_parser()
    prepared = parser.parse_args(["nbody", "prepare", "case.yaml", "--label", "case"])
    assert prepared.nbody_action == "prepare" and prepared.label == "case"
    evolved = parser.parse_args(["nbody", "run", "case", "--engine", "frozen"])
    assert evolved.engine == "frozen"


@pytest.mark.parametrize("change", [dict(gamma=0, beta=.3), dict(mass_msun=1e6),
    dict(scale_pc=-1), dict(n_particles=1), dict(beta=.7), dict(cutoff_pc=np.inf),
    dict(softening_pc=np.nan), dict(radius_pc=None)])
def test_invalid_component_rejected(change):
    with pytest.raises(ValueError):
        replace(config().components[0], **change)


@pytest.mark.parametrize("name", ["../escape", "/absolute", "a/b", "a b", "", ".."])
def test_run_labels_are_single_names(name):
    with pytest.raises(ValueError):
        safe_label(name)


def test_programmes_require_the_correct_components():
    with pytest.raises(ValueError, match="progenitor"):
        replace(config(), programme="progenitor")
    with pytest.raises(ValueError, match="friction"):
        Orbit(kind="class2")


def test_time_units_and_binary_timestep_preserve_absolute_epoch():
    c = replace(config(), orbit=Orbit(duration_myr=.1, sample_step_myr=.01, start_time_myr=2000))
    kmax, dt = timestep(c)
    assert dt == 2**(-kmax)
    assert dt*TIME_UNIT_MYR <= c.integrator.max_step_myr
    assert 0 <= end_time(c)*TIME_UNIT_MYR-2000-.1 < dt*TIME_UNIT_MYR
    command = gyrfalcon_args(c, Path("/space in path/inputs"))
    assert f"Grav={G:.17g}" in command and f"accpars=0,{G:.17g}" in command
    assert "eps=-1" in command and "manipname=ocen_keep_keys" in command
    assert "accfile=/space in path/inputs/external_live.ini" in command


def test_track_uses_velocity_derivatives_and_forbids_extrapolation(tmp_path):
    t = np.array([2., 3., 4.])
    v = np.array([1., 2., -3.])
    track = Track(t, np.column_stack((t[:, None]*v, np.tile(v, (3, 1)))))
    np.testing.assert_allclose(track.at(2.5), np.r_[2.5*v, v])
    track.save(tmp_path)
    np.testing.assert_allclose(Track.load(tmp_path).at(2.5), track.at(2.5))
    with pytest.raises(ValueError, match="extrapolation"):
        track.at(4.1)
    with pytest.raises(ValueError, match="increase"):
        Track(t[::-1], track.xv[::-1])


class PointHost:
    def __init__(self, mass):
        self.mass = mass

    def force(self, x, t=0):
        x = np.asarray(x)
        return -G*self.mass*x/np.linalg.norm(x)**3


class FlatHost:
    def __init__(self, speed):
        self.speed = speed

    def force(self, x, t=0):
        x = np.asarray(x)
        return -self.speed**2*x/np.linalg.norm(x)**2


def test_tidal_tensor_point_mass_and_circular_jacobi_radius():
    host, radius, satellite = PointHost(1e11), 8., 3.55e6
    tensor = tidal_tensor(host, [radius, 0., 0.])
    np.testing.assert_allclose(tensor, np.diag([2., -1., -1.])*G*host.mass/radius**3, rtol=1e-8, atol=1e-9)
    centre = np.array([radius, 0., 0., 0., np.sqrt(G*host.mass/radius), 0.])
    rj, lam = tidal_scale(host, centre, 0, lambda r: satellite)
    assert lam == pytest.approx(3*G*host.mass/radius**3, rel=1e-8)
    assert rj == pytest.approx(radius*(satellite/(3*host.mass))**(1/3), rel=1e-8)


def test_circular_flat_rotation_uses_coefficient_two():
    host = FlatHost(220.)
    rj, lam = tidal_scale(host, np.array([8., 0., 0., 0., 220., 0.]), 0, lambda r: 3.55e6)
    assert lam == pytest.approx(2*220**2/8**2, rel=1e-8)
    assert rj == pytest.approx((G*3.55e6/lam)**(1/3))


def test_particle_monopole_excludes_self_and_unbinding_keeps_arrays():
    radius, mass = np.array([1., 2., 3.]), np.array([2., 3., 4.])
    potential = spherical_particle_potential(radius, mass, np.ones(3, bool))
    np.testing.assert_allclose(potential, -G*np.array([3/2+4/3, 2/2+4/3, 2/3+3/3]))
    local = np.array([[.01, 0, 0, 0, 0, 0], [.02, 0, 0, 10000, 0, 0]], float)
    original = local.copy()
    keep, energy = bound_particles(local, np.ones(2), Nucleus())
    assert list(keep) == [True, False]
    np.testing.assert_array_equal(local, original)
    assert energy[1] > 0


def test_shell_density_and_aperture_mass_have_distinct_meanings():
    c = config()
    xv = np.zeros((3, 6))
    xv[:, 0] = [.005, .02, .1]
    particles = {"mass": np.array([1., 2., 7.]), "component": np.zeros(3, int)}
    row = snapshot_diagnostics(xv, particles, c, np.zeros(6), 0)
    assert row["dm_20pc_mass_msun"] == 3
    volume = 4*np.pi/3*20**3*(np.exp(.3)-np.exp(-.3))
    assert row["dm_20pc_rho_msun_pc3"] == pytest.approx(2/volume)
    assert row["dm_20pc_shell_neff"] == 1
    assert np.isnan(row["dm_20pc_beta"])
    assert effective_number([1, 3]) == 1.6
    assert weighted_median([0., 10., 20.], [1., 100., 1.]) == 10.
    assert np.isnan(weighted_median([], []))


def test_snapshot_parser_sorts_keys_and_rejects_truncation():
    data = "# header\n2\n0.25 1 3 1 2 3 4 5 6\n0.25 0 2 6 5 4 3 2 1\n"
    time, ids, mass, xv = list(parse_snapshots(io.StringIO(data)))[0]
    assert time == .25 and list(ids) == [0, 1] and list(mass) == [2, 3]
    assert list(xv[0]) == [6, 5, 4, 3, 2, 1]
    with pytest.raises(ValueError, match="truncated"):
        list(parse_snapshots(io.StringIO(data.rsplit("0.25 0 2", 1)[0])))
    with pytest.raises(ValueError, match="keys"):
        list(parse_snapshots(io.StringIO(data.replace("1 3 1 2", "0 3 1 2"))))


def test_signed_eddington_inversion_matches_analytic_plummer():
    pytest.importorskip("agama")
    from ocen_dm.dynamics.models import agama_module, signed_df_grid
    agama = agama_module()
    mass, scale = 2e6, .02
    potential = agama.Potential(type="Plummer", mass=mass, scaleRadius=scale)

    def derivatives(radius):
        s = (radius/scale)**2
        return -5*s/(1+s), -10*s/(1+s)**2

    result = signed_df_grid(potential, potential, derivatives, 0., scale*.01, scale*30)
    analytic = 24*np.sqrt(2)/(7*np.pi**3)*scale**2/(G**5*mass**4)*result["binding_energy_kms2"]**3.5
    np.testing.assert_allclose(result["f_energy"], analytic, rtol=3e-4)
    assert np.all(result["signed_over_absolute"] > 0)


@pytest.fixture
def prepared(tmp_path):
    pytest.importorskip("agama")
    from ocen_dm.dynamics.experiment import prepare
    return prepare(config(), tmp_path/"prepared")


def test_prepared_potential_density_recovery_and_preservation(prepared):
    from ocen_dm.dynamics.experiment import load_prepared, prepare
    from ocen_dm.dynamics.models import agama_module
    c, particles, track, manifest = load_prepared(prepared)
    assert manifest["status"] == "prepared"
    assert len(particles["mass"]) == 512
    p = agama_module().Potential(str(prepared/"nucleus.ini"))
    expected = -G*c.nucleus.mass_msun/np.sqrt(.02**2+(.001*c.nucleus.scale_pc)**2)
    assert p.potential([.02, 0., 0.]) == pytest.approx(expected)
    with pytest.raises(FileExistsError):
        prepare(c, prepared)
    manifest["config"]["nucleus"]["mass_msun"] *= 2
    (prepared/"run.yaml").write_text(yaml.safe_dump(manifest))
    with pytest.raises(ValueError, match="disagree"):
        load_prepared(prepared)
    manifest["config"]["nucleus"]["mass_msun"] /= 2
    (prepared/"run.yaml").write_text(yaml.safe_dump(manifest))
    (prepared/"nucleus.ini").write_text("tampered")
    with pytest.raises(ValueError, match="changed"):
        load_prepared(prepared)


def test_frozen_isolation_conserves_specific_energy(prepared):
    from ocen_dm.dynamics.experiment import run, analyse
    directory = run(prepared, "frozen")
    table = analyse(directory)
    assert len(table) >= 3
    assert table["median_abs_fractional_reference_energy_change"][-1] < 1e-6
    np.testing.assert_allclose(table["dm_total_mass_msun"], table["dm_total_mass_msun"][0])
    with pytest.raises(FileExistsError):
        run(prepared, "frozen")


def test_repeated_sampling_is_reproducible_and_core_normalization_is_matched(tmp_path):
    pytest.importorskip("agama")
    from ocen_dm.dynamics.models import build_equilibrium
    c = config()
    samples, reports = [], []
    for label, gamma in (("first", 1.), ("repeat", 1.), ("core", 0.)):
        path = tmp_path/label
        path.mkdir()
        particles, report = build_equilibrium(replace(c, components=(replace(c.components[0], gamma=gamma),)), path)
        samples.append(particles["xv"])
        reports.append(report)
    np.testing.assert_array_equal(samples[0], samples[1])
    for report in reports:
        assert report["components"]["dm"]["mass_inside_pc"]["20.0"] == pytest.approx(1e5, rel=1e-8)


def test_impossible_isotropic_progenitor_core_fails_without_clipping(tmp_path):
    pytest.importorskip("agama")
    from ocen_dm.dynamics.models import build_equilibrium
    root = Path(__file__).resolve().parents[1]
    c = Experiment.load(root/"configs/dynamics/progenitor_core.yaml")
    c = replace(c, components=(replace(c.components[0], beta=0., n_particles=128), c.components[1]))
    with pytest.raises(ValueError, match="negative distribution function"):
        build_equilibrium(c, tmp_path)


def test_pericentre_refinement_recovers_weighted_density_and_protects_centre(tmp_path):
    pytest.importorskip("agama")
    from ocen_dm.dynamics.models import build_equilibrium, agama_module
    c = config()
    component = replace(c.components[0], n_particles=8000, refinement_radius_pc=10.,
                        refinement_floor=.01)
    particles, report = build_equilibrium(replace(c, components=(component,)), tmp_path)
    radius = np.linalg.norm(particles["xv"][:, :3], axis=1)
    weight = particles["mass"]
    assert weight.max() > 3*weight.min()
    np.testing.assert_allclose(weight[radius < .01], weight.min())
    density = agama_module().Potential(str(tmp_path/"dm_potential.ini"))
    for r in (.01, .02, .05, .1):
        fraction = density.enclosedMass(r)/density.totalMass()
        indicator = radius < r
        measured = weight[indicator].sum()/weight.sum()
        error = np.sqrt(np.sum(weight**2*(indicator-fraction)**2))/weight.sum()
        assert abs(measured-fraction) < 6*error
    assert report["components"]["dm"]["sampling"].startswith("pericentre")


@pytest.mark.parametrize("filename", ["progenitor_cusp.yaml", "progenitor_core.yaml"])
def test_joint_progenitor_templates_are_admissible(tmp_path, filename):
    pytest.importorskip("agama")
    from ocen_dm.dynamics.models import build_equilibrium
    root = Path(__file__).resolve().parents[1]
    c = Experiment.load(root/"configs/dynamics"/filename)
    c = replace(c, components=tuple(replace(p, n_particles=1000) for p in c.components))
    particles, report = build_equilibrium(c, tmp_path)
    assert len(particles["mass"]) == 2000
    for component in c.components:
        result = report["components"][component.name]
        assert result["max_density_fractional_error"] < .001
        assert result["max_beta_error"] < .002
        assert result["minimum_signed_df_ratio"] >= 0


@pytest.mark.parametrize("kind, epoch", [("plummer_host", 0.), ("class1", 0.), ("class3", 2000.)])
def test_host_assets_reproduce_forward_orbit_and_clock(tmp_path, kind, epoch):
    pytest.importorskip("agama")
    from ocen_dm.dynamics.orbits import build_track
    if kind == "class3" and not (Path.home()/"Work/Code/oCen_bar/agama_potentials").exists():
        pytest.skip("Hunter bar assets unavailable")
    c = replace(config(), orbit=Orbit(kind=kind, duration_myr=.1, sample_step_myr=.01, start_time_myr=epoch))
    track, record = build_track(c, tmp_path)
    assert track.time[0]*TIME_UNIT_MYR == pytest.approx(epoch)
    assert record["max_interpolated_acceleration_relative_error"] < .001
    assert np.all(np.diff(track.time) > 0)


@pytest.mark.skipif(os.environ.get("OCEN_RUN_NEMO_TESTS") != "1", reason="set OCEN_RUN_NEMO_TESTS=1 for compiled NEMO tests")
def test_live_integrator_units_particle_keys_and_analytic_nucleus(tmp_path):
    pytest.importorskip("agama")
    from ocen_dm.dynamics.experiment import prepare, run, analyse
    from ocen_dm.dynamics import nemo
    c = config()
    c = replace(c, components=(replace(c.components[0], mass_within_msun=1e-3),))
    path = prepare(c, tmp_path/"tracers")
    evolution = run(path, "live")
    table = analyse(evolution)
    assert table["median_abs_fractional_reference_energy_change"][-1] < 1e-3
    times = []
    for time, ids, mass, xv in nemo.snapshots(evolution/"snapshots.nemo", evolution):
        np.testing.assert_array_equal(ids, np.arange(512))
        times.append(time)
    assert times[-1]*TIME_UNIT_MYR == pytest.approx(end_time(c)*TIME_UNIT_MYR, abs=1e-5)


@pytest.mark.skipif(os.environ.get("OCEN_RUN_NEMO_TESTS") != "1", reason="set OCEN_RUN_NEMO_TESTS=1 for compiled NEMO tests")
def test_live_self_gravity_reproduces_a_circular_binary(tmp_path):
    from ocen_dm.dynamics import nemo
    radius, mass = .01, 1e6
    omega = np.sqrt(G*mass/(4*radius**3))
    quarter_period_myr = np.pi/(2*omega)*TIME_UNIT_MYR
    c = replace(config(), orbit=Orbit(duration_myr=quarter_period_myr, sample_step_myr=.01),
                integrator=Integrator(max_step_myr=.002, levels=1, output_step_myr=.2))
    track = Track(np.array([0., end_time(c)]), np.zeros((2, 6)))
    particles = dict(xv=np.array([[radius, 0, 0, 0, omega*radius, 0],
                                  [-radius, 0, 0, 0, -omega*radius, 0]]),
                     mass=np.array([mass, mass]), softening=np.full(2, 1e-6),
                     particle_id=np.array([101, 202]))
    # A negligible external mass leaves a direct, independent self-gravity reference.
    (tmp_path/"external_live.ini").write_text("[Potential]\ntype=Plummer\nmass=1e-8\nscaleRadius=1\n")
    nemo.build_key_plugin(tmp_path)
    nemo.write_input(particles, track, tmp_path)
    nemo.call(nemo.gyrfalcon_args(c, tmp_path), tmp_path, "gyrfalcon.log")
    snapshots = list(nemo.snapshots(tmp_path/"snapshots.nemo", tmp_path))
    time, keys, actual_mass, xv = snapshots[-1]
    np.testing.assert_array_equal(keys, [101, 202])
    np.testing.assert_allclose(actual_mass, mass)
    expected = radius*np.array([np.cos(omega*time), np.sin(omega*time), 0.])
    np.testing.assert_allclose(xv[0, :3], expected, atol=radius*2e-4)
    np.testing.assert_allclose(xv[1, :3], -expected, atol=radius*2e-4)
