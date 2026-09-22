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
