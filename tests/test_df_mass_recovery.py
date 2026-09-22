"""Physical controls for independent mass/light weighting and template gravity."""
from dataclasses import replace
import json

import numpy as np
import pytest
from scipy.integrate import simpson

from ocen_dm.kinematics.df_mass_recovery import (
    MassLightConfig, MassLightDFModel, TemplateMassConfig, DensityTemplate,
    RecoveryProblem, gaussian_mock_noise,
)
from ocen_dm.kinematics.positive_df import (
    ActionComponent, DFModelConfig, DFNumerics, PositiveDFModel, agama_pc,
    spherical_velocity_moments,
)
from ocen_dm.kinematics.df_mock import challenge_mock
from ocen_dm.kinematics.df_capacity import CapacityObservations
from ocen_dm.kinematics.df_fit import PhotometricData
from ocen_dm.kinematics.likelihood import BinnedProfile,KinematicData


@pytest.fixture(scope="module")
def config():
    pytest.importorskip("agama")
    return MassLightConfig(DFModelConfig(M_bh=0.,M_rem=0.,components=(
        ActionComponent(fraction=.65,J0=240.),ActionComponent(fraction=.35,J0=450.,g_r=1.6))),(.4,.6))


@pytest.fixture(scope="module")
def model(config):
    return MassLightDFModel(config)


def test_mass_light_config_roundtrip_and_rejection(config):
    assert MassLightConfig.from_dict(json.loads(json.dumps(config.to_dict())))==config
    for fractions in ((.5,),(.2,.2),(-.1,1.1),(np.nan,.5)):
        with pytest.raises(ValueError): replace(config,light_fractions=fractions)
    with pytest.raises(ValueError): replace(config,mass=replace(config.mass,M_rem=1.),dark=TemplateMassConfig())


def test_equal_weights_reproduce_original_model(config):
    same=replace(config,light_fractions=tuple(c.fraction for c in config.mass.components))
    original=PositiveDFModel(config.mass)
    new=MassLightDFModel(same)
    r=np.geomspace(.05,70,16)
    for key,v in original.projected_moments(r).items():
        np.testing.assert_allclose(v,new.projected_moments(r)[key],rtol=1e-13)


def test_light_changes_moments_without_changing_gravity(model,config):
    other=model.reweight_light(replace(config,light_fractions=(.8,.2)))
    r=np.geomspace(.1,70,20)
    np.testing.assert_array_equal(model.mass_profile(r),other.mass_profile(r))
    assert np.max(abs(model.projected_moments(r)['Sigma']/other.projected_moments(r)['Sigma']-1))>.1
    assert model.mass_light_config.light_fractions==(.4,.6)
    assert model.diagnostics['light_fractions']==[.4,.6]
    np.testing.assert_allclose(model.mass_df.totalMass(),config.mass.M_star,rtol=1e-7)
    # Native action-space normalization is itself a numerical integral.
    np.testing.assert_allclose(other.df.totalMass(),config.mass.M_star,rtol=3e-7)
    actions=10**np.random.default_rng(4).uniform(-6,6,(100,3))
    assert np.all(other.df(actions)>0) and np.all(other.mass_df(actions)>0)


def test_mass_df_still_sources_the_potential(model):
    r=np.geomspace(.05,70,30)
    density=spherical_velocity_moments(model.potential,model.mass_df,model.af,r,96)[:,0]
    xyz=np.column_stack((r,r*0,r*0))
    np.testing.assert_allclose(density,model.stellar_density.density(xyz),rtol=.003)


def test_warm_and_cold_solutions_agree(config,model):
    warm=MassLightDFModel(config,initial_stellar_potential=model.stellar_potential)
    r=np.geomspace(.1,60,16)
    for key,v in model.projected_moments(r).items():
        np.testing.assert_allclose(v,warm.projected_moments(r)[key],rtol=.001)


def test_density_template_scaling_and_mass():
    mock=challenge_mock('mixed_with_dm')
    template=DensityTemplate.from_mock(mock,'extended_dark_component')
    r=np.geomspace(1e-6,template.edge*1.7,5001)
    rho=template.density(r,4e5,1.7)
    assert simpson(4*np.pi*r**3*rho,x=np.log(r))==pytest.approx(4e5,rel=3e-6)
    assert template.density(np.array([2*template.edge]),4e5)[0]==0
    pot=template.potential(4e5,1.7)
    assert pot.totalMass()==pytest.approx(4e5,rel=2e-5)
    restored=DensityTemplate.from_dict(template.to_dict())
    np.testing.assert_array_equal(restored.density(r,4e5,1.7),rho)


def test_recovery_coordinate_noise_and_zero_dark_boundaries(config):
    mock=challenge_mock('mixed_with_dm')
    r=np.geomspace(4,2000,10)
    profile=BinnedProfile('test','los',r,None,None,r*0+10,r*0+1,r*0+1,'mock')
    photo=PhotometricData(r,r*0,r*0+.1,'mock','adopted')
    obs=CapacityObservations(KinematicData((profile,)),photo,mock)
    templates={name:DensityTemplate.from_mock(mock,pop) for name,pop in
               [('remnant','dark_remnants'),('halo','extended_dark_component')]}
    config=replace(config,dark=TemplateMassConfig(remnant_mass=0.,halo_mass=0.))
    problem=RecoveryProblem(obs,config,templates=templates,shapes='weights')
    decoded=problem.decode(problem.encode())
    np.testing.assert_allclose(problem.encode(decoded),problem.encode(),rtol=2e-15,atol=2e-15)
    assert decoded.dark==config.dark
    noise=gaussian_mock_noise(obs,17)
    np.testing.assert_array_equal(noise,gaussian_mock_noise(obs,17))
    noisy=RecoveryProblem(obs,config,templates=templates,shapes='weights',noise=noise)
    np.testing.assert_array_equal(noisy.target,obs.truth+noise)
    assert not np.array_equal(noise,gaussian_mock_noise(obs,18))
    shared=RecoveryProblem(obs,config,templates=templates,shapes='weights',shared_ml=True)
    decoded_shared=shared.decode(shared.encode())
    assert decoded_shared.light_fractions==tuple(c.fraction for c in decoded_shared.mass.components)
    assert len(shared.names)==len(problem.names)-1
    with pytest.raises(ValueError): RecoveryProblem(obs,config,noise=np.array([1.]))
