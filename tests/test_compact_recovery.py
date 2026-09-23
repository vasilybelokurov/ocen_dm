from dataclasses import replace

import numpy as np
import pytest

from ocen_dm.kinematics.compact_recovery import (
    mock_problem, recovery_gates, recovery_metrics, refined_config, residual_vector,
)
from ocen_dm.kinematics.df_fit import PhotometricData
from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData
from ocen_dm.kinematics.regularized_df import RegularizedDFConfig


class ToyModel:
    config = RegularizedDFConfig()

    def _check_inside_tracer(self, r):
        assert np.all(np.asarray(r) > 0)

    def projected_moments(self, r):
        r = np.asarray(r)
        sigma = np.exp(-r/20)
        return dict(Sigma=sigma, los=100*sigma, pmr=90*sigma, pmt=80*sigma)


def problem():
    profile = BinnedProfile("mock", "los", np.array([10., 30.]), None, None,
                            np.array([8., 9.]), np.array([.5, .6]), np.array([1., 1.2]),
                            "test", streaming2=np.ones(2))
    photo = PhotometricData(np.array([10., 30.]), np.array([20., 21.]),
                            np.array([.1, .2]), "test", "adopted")
    return mock_problem(KinematicData((profile,)), photo, ToyModel()), profile


def test_mock_has_zero_residual_and_preserves_template():
    p, original = problem()
    assert p.data.profiles[0].streaming2 is None
    np.testing.assert_array_equal(original.streaming2, 1)
    np.testing.assert_array_equal(p.data.profiles[0].err_lo, original.err_lo)
    ev = p.evaluate(ToyModel.config, ToyModel())
    np.testing.assert_allclose(residual_vector(p, ev), 0, atol=1e-13)


def test_residual_norm_matches_objective_with_asymmetric_errors():
    p, _ = problem()
    altered = replace(p.data.profiles[0], value=np.array([9., 11.]))
    p.data = KinematicData((altered,))
    ev = p.evaluate(ToyModel.config, ToyModel())
    r = residual_vector(p, ev)
    assert r @ r == pytest.approx(ev["objective"])
    np.testing.assert_allclose(r[:2], [1., -1/.6])


def test_refinement_includes_contour_map():
    c = RegularizedDFConfig()
    f = refined_config(c)
    assert f.contours.action_nodes == 2*c.contours.action_nodes
    assert f.contours.circularity_nodes == 2*c.contours.circularity_nodes-1
    assert f.numerics.velocity_nodes == 2*c.numerics.velocity_nodes
    assert f.M_star == c.M_star


def test_good_observables_do_not_certify_mass_recovery():
    truth = dict(r_pc=[1., 20.], total=[1e4, 1e6], M_star=3e6, rho20=0.)
    recovered = dict(truth, rho20=2.)
    metrics = recovery_metrics(truth, recovered)
    gates = recovery_gates(True, True, np.zeros(10), metrics)
    assert gates["observable_accuracy_passed"]
    assert not gates["passed"]
    assert metrics["rho20_relative_error"] is None


def test_capped_run_cannot_pass_and_injection_uses_relative_error():
    truth = dict(r_pc=[1., 20.], total=[1e4, 1e6], M_star=3e6, rho20=2.)
    metrics = recovery_metrics(truth, dict(truth, rho20=2.1))
    assert metrics["physical_accuracy_passed"]
    assert not recovery_gates(False, True, np.zeros(10), metrics)["passed"]
    assert recovery_gates(True, True, np.zeros(10), metrics)["passed"]


def test_absolute_stop_needs_full_window_of_small_changes():
    from ocen_dm.kinematics.compact_recovery import AbsoluteStop, Converged
    stop = AbsoluteStop(delta=.01, window=2)
    for value in (100., 10., 1., .995):  # last window change is 0.005 + big step before
        stop.update(value)
    with pytest.raises(Converged):
        stop.update(.994)  # 1.0 -> 0.994 over two steps < 0.01


def test_absolute_stop_continues_while_improving_and_validates():
    from ocen_dm.kinematics.compact_recovery import AbsoluteStop
    stop = AbsoluteStop(delta=.01, window=2)
    for value in (1., .5, .25, .12, .06, .03):
        stop.update(value)
    with pytest.raises(ValueError):
        AbsoluteStop(delta=0)


def _evaluation(kin_chi2, phot_residual, sigma):
    phot_residual, sigma = np.asarray(phot_residual, float), np.asarray(sigma, float)
    return dict(terms=dict(a=dict(chi2=kin_chi2, n=10)),
                photometry=dict(chi2=float(np.sum((phot_residual/sigma)**2)), prediction=phot_residual),
                photometry_mu=np.zeros_like(phot_residual), photometry_sigma=sigma)


GATES = dict(chi2_kin_per_point=1.3, chi2_per_point_any_dataset=2., chi2_phot_per_point=2.)


def test_data_gate_uses_weighted_photometry_not_raw_rms():
    from ocen_dm.kinematics.compact_recovery import data_fit_gates
    # A 1 mag outlier with 0.58 mag error is within 2 sigma: raw RMS is large, fit is fine.
    ev = _evaluation(10., [0., .05, -.05, 1.], [.1, .1, .1, .58])
    g = data_fit_gates(True, True, ev, GATES)
    assert g["chi2_phot_per_point"] < 2 and g["data_fit_passed"] and g["passed"]
    raw_rms = float(np.sqrt(np.mean(np.array([0., .05, .05, 1.])**2)))
    assert g["photometric_weighted_rms_mag"] < raw_rms


def test_data_gate_fails_on_any_dataset_or_optimizer():
    from ocen_dm.kinematics.compact_recovery import data_fit_gates
    assert not data_fit_gates(True, True, _evaluation(25., [0.], [.1]), GATES)["data_fit_passed"]
    g = data_fit_gates(False, True, _evaluation(10., [0.], [.1]), GATES)
    assert g["data_fit_passed"] and not g["passed"]


def test_data_radius_covers_bin_edges_and_photometry():
    from ocen_dm.kinematics.compact_recovery import data_radius_pc
    from ocen_dm.kinematics.likelihood import ARCSEC_PER_RAD
    p, _ = problem()
    edged = replace(p.data.profiles[0], r_lower=np.array([5., 20.]), r_upper=np.array([20., 40.]))
    p.data = KinematicData((edged,))
    assert data_radius_pc(p, 5.) == pytest.approx(40.*5000/ARCSEC_PER_RAD)


def test_difference_column_falls_back_and_never_uses_rejections():
    from ocen_dm.kinematics.compact_recovery import difference_column
    base = np.array([1., 2.])
    np.testing.assert_allclose(difference_column(base, .5, forward=[2., 2.]), [2., 0.])
    np.testing.assert_allclose(difference_column(base, .5, backward=[0., 2.]), [2., 0.])
    np.testing.assert_array_equal(difference_column(base, .5), [0., 0.])


def test_stop_near_rejection_window():
    from ocen_dm.kinematics.compact_recovery import stopped_near_rejection
    assert stopped_near_rejection([95], 100, 20)
    assert not stopped_near_rejection([10, 50], 100, 20)
    assert not stopped_near_rejection([], 100, 20)


def test_outer_ratio_coordinate_keeps_transitions_ordered():
    from ocen_dm.kinematics.df_fit import FitCoordinate, config_at
    base = RegularizedDFConfig.from_dict(dict(RegularizedDFConfig().to_dict(),
        stellar=dict(RegularizedDFConfig().to_dict()["stellar"], b_outer=-1., J_outer=900.)))
    coords = [FitCoordinate("stellar.J_a", 5., 500., log=True),
              FitCoordinate("stellar.log_J_outer_ratio", np.log(2.), np.log(100.))]
    assert coords[1].get(base) == pytest.approx(np.log(900./base.stellar.J_a))
    c = config_at(base, coords, [np.log(400.), np.log(3.)])
    assert c.stellar.J_a == pytest.approx(400.) and c.stellar.J_outer == pytest.approx(1200.)
    with pytest.raises(ValueError):
        FitCoordinate("stellar.log_J_outer_ratio", 0., 2.)       # must stay positive
    with pytest.raises(ValueError):
        config_at(RegularizedDFConfig(), coords[1:], [np.log(3.)])  # no second transition


def _count_profile():
    from ocen_dm.kinematics.counts import CountProfile
    lo = np.array([5., 10., 20.]); hi = np.array([10., 20., 40.])
    nodes = np.sqrt(lo*hi)[:, None]*np.array([.8, 1., 1.2])[None, :]
    return CountProfile("toy_counts", lo, hi, nodes, np.array([120, 80, 30]), np.pi*(hi**2-lo**2), "s", "sel", fit_field=False)


def test_joint_problem_with_counts_and_prior_residual_layout():
    from ocen_dm.kinematics.compact_recovery import residual_vector
    from ocen_dm.kinematics.df_fit import DFJointProblem
    p, profile = problem()
    joint = DFJointProblem(p.data, None, [_count_profile()])
    assert joint.n_residuals == 2+3
    ev = joint.evaluate(ToyModel.config, ToyModel())
    assert ev["photometry"] is None and set(ev["counts"]) == {"toy_counts"}
    assert ev["objective"] == pytest.approx(ev["chi2_kinematic"]+ev["deviance_counts"])
    r = residual_vector(joint, ev, priors={"distance_kpc": (5.0, 0.1)})
    assert r.size == 6
    assert r[-1] == pytest.approx((ToyModel.config.distance_kpc-5.0)/0.1)
    assert (r[2:5] @ r[2:5]) == pytest.approx(ev["deviance_counts"])
    with pytest.raises(ValueError):
        DFJointProblem(p.data, None, [])


def test_data_gates_with_counts_only():
    from ocen_dm.kinematics.compact_recovery import data_fit_gates
    ev = dict(terms=dict(a=dict(chi2=10., n=10)), photometry=None,
              counts=dict(c=dict(deviance=15., n=10, amplitude=1., field_density_per_arcmin2=0.)))
    g = data_fit_gates(True, True, ev, dict(chi2_kin_per_point=1.3, chi2_per_point_any_dataset=2., deviance_per_bin_counts=2.))
    assert g["deviance_per_bin_by_counts"]["c"] == pytest.approx(1.5) and g["passed"]
    ev["counts"]["c"]["deviance"] = 25.
    assert not data_fit_gates(True, True, ev, dict(chi2_kin_per_point=1.3, chi2_per_point_any_dataset=2., deviance_per_bin_counts=2.))["passed"]


def test_mock_counts_problem_noiseless_and_poisson():
    from ocen_dm.kinematics.compact_recovery import mock_counts_problem
    p, profile = problem()
    c = _count_profile()
    truth = ToyModel()
    noiseless = mock_counts_problem(p.data, [c], truth, {"toy_counts": 50.})
    assert noiseless.photometry is None and noiseless.counts[0].counts.sum() > 0
    ev = noiseless.evaluate(ToyModel.config, truth)
    assert ev["deviance_counts"] < .5          # rounding only
    assert ev["chi2_kinematic"] == pytest.approx(0., abs=1e-10)
    rng = np.random.default_rng(5)
    noisy = mock_counts_problem(p.data, [c], truth, {"toy_counts": 50.}, rng=rng)
    assert not np.array_equal(noisy.counts[0].counts, noiseless.counts[0].counts)
    assert noisy.data.profiles[0].value.shape == profile.value.shape
    np.testing.assert_array_equal(noisy.counts[0].r_nodes, c.r_nodes)


def test_noisy_recovery_gates_require_fit_and_physics():
    from ocen_dm.kinematics.compact_recovery import noisy_recovery_gates
    ev = dict(terms=dict(a=dict(chi2=10., n=10)), photometry=None,
              counts=dict(c=dict(deviance=12., n=10, amplitude=1., field_density_per_arcmin2=0.)))
    thr = dict(chi2_kin_per_point=1.3, chi2_per_point_any_dataset=2., deviance_per_bin_counts=2.)
    good = dict(physical_accuracy_passed=True); bad = dict(physical_accuracy_passed=False)
    assert noisy_recovery_gates(True, True, ev, thr, good)["passed"]
    assert not noisy_recovery_gates(True, True, ev, thr, bad)["passed"]
    assert not noisy_recovery_gates(False, True, ev, thr, good)["passed"]
