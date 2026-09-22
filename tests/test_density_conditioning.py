"""Fixed-density sampling preserves the parent model and conditional prior."""
from dataclasses import replace
import json
from pathlib import Path
import runpy

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest

from ocen_dm.kinematics.density_conditioning import DensityConditionedModel
from ocen_dm.kinematics.fit import DarkMatterModel, FitProblem, Prior
from ocen_dm.kinematics.likelihood import BinnedProfile, KinematicData
from ocen_dm.kinematics.run_io import family_config, family_from_config
from ocen_dm.light_model import MGEFit


def parent(gamma=0.):
    return DarkMatterModel(gamma=gamma, instruments=(), anisotropy_profile="turnover",
        mge_fit=MGEFit(np.array([30.,300.,1500.]), np.array([.1,.8,.1]), 0., 0., 3))


@pytest.mark.parametrize("gamma", [0., 1.])
def test_conditioning_preserves_density_likelihood_and_snapshot(gamma):
    original = parent(gamma)
    conditional = DensityConditionedModel(original, 2.)
    assert len(conditional.names) == 11
    assert conditional.parameters == tuple(p for p in original.parameters if p.name!="M_dm_100")
    data = KinematicData((BinnedProfile("mock", "pmt", np.array([30.,100.,500.]),
        None, None, np.array([.5,.4,.3]), np.full(3,.02), np.full(3,.04), "mock"),))
    u = np.linspace(.15,.85,11)
    x = conditional.transform(u)
    theta = conditional.complete(conditional.to_dict(x))
    ref = np.array([theta[n] for n in original.names])
    assert_allclose(FitProblem(conditional,data).loglike_vector(x),
                    FitProblem(original,data).loglike_vector(ref), atol=1e-10)
    config = json.loads(json.dumps(family_config(conditional)))
    restored = family_from_config(config)
    assert family_config(restored)==config
    assert_array_equal(restored.transform(u),x)
    assert_allclose(FitProblem(restored,data).loglike_vector(x),
                    FitProblem(conditional,data).loglike_vector(x),atol=1e-10)
    for rs in np.geomspace(3,1000,50):
        theta["r_s"]=rs
        full=conditional.complete(theta)
        assert_allclose(original.halo(full).density(20.),2.,rtol=1e-12)
        assert conditional.mass_support_envelope[0] <= full["M_dm_100"] <= conditional.mass_support_envelope[1]


def test_loguniform_mass_jacobian_leaves_radius_prior_unchanged():
    c=DensityConditionedModel(parent(),2.)
    mass_range=np.log(c.mass_prior.hi/c.mass_prior.lo)
    # p(M | rs) dM/drho is independent of rs: verify across the full support.
    transformed=[]
    for rs in np.geomspace(3,1000,21):
        mass=c.complete({"r_s":rs})["M_dm_100"]
        jacobian=mass/c.density
        transformed.append(jacobian/(mass*mass_range))
    assert_allclose(transformed,c.prior_density_at_constraint,rtol=1e-14)


def test_unsupported_conditioning_fails_instead_of_silently_changing_prior():
    with pytest.raises(ValueError,match="truncate"):
        DensityConditionedModel(parent(),1e4)
    with pytest.raises(ValueError,match="positive finite"):
        DensityConditionedModel(parent(),0.)
    p=parent()
    p.parameters=tuple(replace(q,prior=Prior("uniform",1e3,3e7))
                       if q.name=="M_dm_100" else q for q in p.parameters)
    with pytest.raises(ValueError,match="log-uniform"):
        DensityConditionedModel(p,2.)


def test_launch_preserves_existing_batch_and_fit(tmp_path,monkeypatch):
    driver=runpy.run_path(str(Path(__file__).resolve().parents[1]/"bin/run_density_posterior_checks.py"))
    launch=driver["launch"]
    monkeypatch.setitem(launch.__globals__,"ROOT",tmp_path)
    with pytest.raises(FileExistsError,match="Batch already exists"):
        launch(tmp_path)
    batch=tmp_path/"new_batch"
    label=driver["labels_for"](batch)["rho2"]
    existing=tmp_path/"results/fits"/label
    existing.mkdir(parents=True)
    marker=existing/"preserve.txt"
    marker.write_text("original")
    with pytest.raises(FileExistsError,match="Fit already exists"):
        launch(batch)
    assert marker.read_text()=="original"
    assert not batch.exists()


def test_evidence_ratio_recovers_density_not_probability_atom(tmp_path,monkeypatch):
    driver=runpy.run_path(str(Path(__file__).resolve().parents[1]/"bin/run_density_posterior_checks.py"))
    compare=driver["compare_if_complete"]
    monkeypatch.setitem(compare.__globals__,"ROOT",tmp_path)
    batch=tmp_path/"batch"
    batch.mkdir()
    for kind,label in driver["labels_for"](batch).items():
        path=tmp_path/"results/fits"/label
        path.mkdir(parents=True)
        result=dict(logz=10.+(np.log(3) if kind=="rho2" else 0.),logzerr=.1,
                    parent_prior_pdf_at_fixed_density=.05)
        (path/"density_diagnostics.json").write_text(json.dumps(result))
    compare(batch)
    result=json.loads((batch/"density_comparison.json").read_text())
    assert_allclose(result["core_posterior_pdf_at_rho2"],.15)
