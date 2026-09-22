"""Reference tests for responsive nuclei and historical endpoints."""
from dataclasses import replace
from pathlib import Path
import os

import numpy as np
import pytest
import yaml

from ocen_dm.dynamics import G, TIME_UNIT_MYR
from ocen_dm.dynamics.config import Experiment, Orbit
from ocen_dm.dynamics.diagnostics import snapshot_diagnostics, nucleus_potential
from ocen_dm.dynamics.experiment import prepare, load_prepared, run
from ocen_dm.dynamics.models import agama_module
from ocen_dm.dynamics.nemo import gyrfalcon_args
from ocen_dm.dynamics.orbits import build_track


def model(gamma=1.):
    base = Experiment.load(Path(__file__).resolve().parents[1]/'configs/dynamics/remnant_cusp.yaml')
    return replace(base, nucleus=replace(base.nucleus, live=True, n_particles=4096, beta=-.1),
                   components=(replace(base.components[0], n_particles=2048, gamma=gamma),),
                   orbit=Orbit(duration_myr=.1, sample_step_myr=.01),
                   integrator=replace(base.integrator, max_step_myr=.001, output_step_myr=.05))


@pytest.fixture
def live_prepared(tmp_path):
    pytest.importorskip('agama')
    from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData
    datum = BinnedProfile(name='test_los', kind='los', instrument='test', r=np.array([200.]),
                         r_lower=np.array([0.]), r_upper=np.array([1000.]),
                         value=np.array([15.]), err_lo=np.array([1.]), err_hi=np.array([1.]))
    return prepare(model(), tmp_path/'live', stellar_data=KinematicData((datum,)))


def test_live_nucleus_has_gravity_exactly_once(live_prepared):
    c, particles, track, record = load_prepared(live_prepared)
    assert len(particles['mass']) == 6144
    agama = agama_module()
    combined = agama.Potential(str(live_prepared/'satellite_initial.ini'))
    dm = agama.Potential(str(live_prepared/'dm_potential.ini'))
    position = [.02, 0., 0.]
    expected = -G*c.nucleus.mass_msun/np.sqrt(.02**2+.0054**2)
    assert combined.potential(position)-dm.potential(position) == pytest.approx(expected)
    assert nucleus_potential(.02, c.nucleus) == 0
    assert not any(arg.startswith('accname=') for arg in gyrfalcon_args(c, live_prepared))
    assert 'nucleus.ini' not in (live_prepared/'external_live.ini').read_text()
    assert 'satellite_initial.ini' in (live_prepared/'external_frozen.ini').read_text()
    np.testing.assert_allclose(np.average(particles['xv'][:, 3:], weights=particles['mass'], axis=0), 0, atol=1e-12)


def test_stellar_diagnostics_are_translation_and_boost_invariant(live_prepared):
    c, particles, _, _ = load_prepared(live_prepared)
    a = snapshot_diagnostics(particles['xv'], particles, c, np.zeros(6), 0)
    shift = np.array([3., -4., .7, 100., -50., 12.])
    b = snapshot_diagnostics(particles['xv']+shift, particles, c, np.zeros(6), 0)
    for key in a:
        if key.startswith(('dm_', 'nucleus_')):
            assert b[key] == pytest.approx(a[key], rel=1e-9, abs=1e-9, nan_ok=True)
    assert a['nucleus_half_mass_radius_pc'] == pytest.approx(5.4/(2**(2/3)-1)**.5, rel=.05)
    for axis in 'xyz':
        assert a[f'nucleus_projected_half_mass_{axis}_pc'] == pytest.approx(5.4, rel=.06)
    assert a['nucleus_bound_fraction'] > .995
    hot = particles['xv'].copy()
    hot[:, 3:] *= 1000
    disrupted = snapshot_diagnostics(hot, particles, c, np.zeros(6), 0)
    assert disrupted['nucleus_bound_fraction'] < .01


def test_isotropic_plummer_in_cusp_is_rejected(tmp_path):
    from ocen_dm.dynamics.models import build_equilibrium
    c = model()
    with pytest.raises(ValueError, match='nucleus: negative distribution function'):
        build_equilibrium(replace(c, nucleus=replace(c.nucleus, beta=0)), tmp_path)


def test_live_nucleus_frozen_control_completes(live_prepared):
    from astropy.table import Table
    directory = run(live_prepared, 'frozen')
    table = Table.read(directory/'diagnostics.ecsv')
    assert table['nucleus_bound_fraction'][-1] > .99
    assert table['nucleus_half_mass_radius_pc'][-1]/table['nucleus_half_mass_radius_pc'][0] == pytest.approx(1, abs=.03)
    observables = Table.read(directory/'stellar_observables.ecsv')
    assert len(observables) == 6
    assert set(observables['epoch']) == {'initial', 'final'}
    assert np.all(np.isfinite(observables['prediction']))


@pytest.mark.parametrize('kind', ['class1', 'class3'])
def test_historical_orbit_returns_to_present(tmp_path, kind):
    duration = 50.
    c = model()
    from ocen_dm.dynamics.orbits import timestep
    dt = timestep(c)[1]*TIME_UNIT_MYR
    duration = np.ceil(duration/dt)*dt
    start = 8*TIME_UNIT_MYR-duration if kind == 'class3' else 0.
    c = replace(c, orbit=Orbit(kind=kind, start_time_myr=start, duration_myr=duration,
                               sample_step_myr=.02, end_at_present=True),
                integrator=replace(c.integrator, accuracy=1e-11))
    _, record = build_track(c, tmp_path)
    assert record['present_endpoint_position_error_pc'] < .01
    assert record['present_endpoint_velocity_error_kms'] < .001


@pytest.mark.skipif(os.environ.get('OCEN_RUN_NEMO_TESTS') != '1', reason='compiled NEMO integration')
def test_live_self_gravitating_nucleus_completes(live_prepared):
    from astropy.table import Table
    from ocen_dm.dynamics.nemo import snapshots
    directory = run(live_prepared, 'live')
    table = Table.read(directory/'diagnostics.ecsv')
    assert table['nucleus_bound_fraction'][-1] > .99
    assert table['nucleus_half_mass_radius_pc'][-1]/table['nucleus_half_mass_radius_pc'][0] == pytest.approx(1, abs=.04)
    # The shrinking-sphere centre has finite-sample noise; conservation of
    # total momentum is the physical reference for this isolated calculation.
    centres = [np.average(xv, weights=mass, axis=0) for _, _, mass, xv in
               snapshots(directory/'snapshots.nemo', directory)]
    assert np.linalg.norm(centres[-1][3:]-centres[0][3:]) < 1e-3
    assert np.linalg.norm(centres[-1][:3]-centres[0][:3])*1000 < .01


def test_heating_time_units_and_mass_scaling():
    from ocen_dm.dynamics.lifetime import heating_time_gyr
    # Independent pc/(km/s) conversion, rather than production kpc units.
    reference = .814*20**3/((G*1000)**2*.5*100*12)*.9777922216807892/1000
    assert heating_time_gyr(20., 100., .5, 12.) == pytest.approx(reference)
    assert heating_time_gyr(20., 100., 1., 12.) == pytest.approx(reference/2)
    with pytest.raises(ValueError):
        heating_time_gyr(20., 0., .5, 12.)


def test_projected_profiles_follow_units_and_flag_sparse_bins(live_prepared):
    from ocen_dm.dynamics.lifetime import projected_stellar_profiles
    from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData
    c, particles, _, _ = load_prepared(live_prepared)
    profiles = []
    for kind in ['los', 'pmr', 'pmt']:
        profiles.append(BinnedProfile(name=kind, instrument='test', kind=kind,
                         r=np.array([200., 2e6]), r_lower=np.array([0., 1e6]),
                         r_upper=np.array([1e6, 3e6]), value=np.array([1., 1.]),
                         err_lo=np.ones(2), err_hi=np.ones(2)))
    rows = projected_stellar_profiles(particles['xv'], particles, c, np.zeros(6), KinematicData(tuple(profiles)))
    resolved = [r for r in rows if r['radius_arcsec'] == 200.]
    assert all(r['neff'] > 4000 for r in resolved)
    assert all(10 < r['prediction'] < 25 for r in resolved if r['kind'] == 'los')
    assert all(.3 < r['prediction'] < 1 for r in resolved if r['kind'] != 'los')
    assert all(np.isnan(r['prediction']) for r in rows if r['radius_arcsec'] == 2e6)


@pytest.mark.skipif(os.environ.get('OCEN_RUN_NEMO_TESTS') != '1', reason='compiled NEMO integration')
def test_live_nucleus_in_host_has_no_duplicate_moving_nucleus(tmp_path):
    from ocen_dm.dynamics.nemo import snapshots
    c = replace(model(), orbit=Orbit(kind='plummer_host', duration_myr=.1, sample_step_myr=.01))
    prepared = prepare(c, tmp_path/'host')
    external = (prepared/'external_live.ini').read_text()
    assert 'host.ini' in external and 'nucleus' not in external and 'center=' not in external
    directory = run(prepared, 'live')
    states = [(time, np.average(xv, weights=mass, axis=0)) for time, _, mass, xv in
              snapshots(directory/'snapshots.nemo', directory)]
    host = agama_module().Potential(str(prepared/'host.ini'))
    elapsed = states[-1][0]-states[0][0]
    expected = host.force(states[0][1][:3])*elapsed
    measured = states[-1][1][3:]-states[0][1][3:]
    np.testing.assert_allclose(measured, expected, atol=.005)


@pytest.mark.parametrize('stage', ['passage', 'bridge'])
def test_queue_numerical_gate_rejects_drift_and_changed_inputs(tmp_path, stage):
    import hashlib
    import importlib.util
    from astropy.table import Table
    path = Path(__file__).resolve().parents[1]/'bin/run_lifetime_batch.py'
    spec = importlib.util.spec_from_file_location('lifetime_queue', path)
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)
    controls = ['base', 'dt_half', 'n2']+(['seed43'] if stage == 'bridge' else [])
    suffix = '_live_nucleus' if stage == 'bridge' else ''
    def write_case(label, density):
        p = tmp_path/'runs'/label/'evolution_live'
        p.mkdir(parents=True, exist_ok=True)
        table = Table(dict(elapsed_myr=np.arange(6.), dm_20pc_rho_msun_pc3=density,
            dm_20pc_shell_neff=np.full(6, 2000.), nucleus_half_mass_radius_pc=np.full(6, 7.),
            nucleus_bound_mass_msun=np.full(6, 3.55e6), nucleus_sigma_x_kms=np.full(6, 20.),
            nucleus_sigma_y_kms=np.full(6, 20.), nucleus_sigma_z_kms=np.full(6, 20.)))
        table.write(p/'diagnostics.ecsv', overwrite=True)
        (p/'run.yaml').write_text(yaml.safe_dump(dict(status='complete',
            diagnostics_sha256=hashlib.sha256((p/'diagnostics.ecsv').read_bytes()).hexdigest())))
        return p
    for shape in ['cusp', 'core']:
        for control in controls:
            for host in ['isolation', 'class1']:
                write_case(f'{shape}_{control}_{host}{suffix}', np.full(6, 100. if host == 'isolation' else 97.))
    assert queue.gate(tmp_path, stage)['passed']
    p = write_case(f'cusp_base_isolation{suffix}', np.array([100., 100., 100., 130., 130., 130.]))
    assert not queue.gate(tmp_path, stage)['passed']
    with (p/'diagnostics.ecsv').open('a') as stream:
        stream.write('# changed\n')
    with pytest.raises(ValueError, match='changed diagnostics'):
        queue.gate(tmp_path, stage)
