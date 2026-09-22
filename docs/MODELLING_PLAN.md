# Modelling plan: add freedom when the residuals justify it

*Adopted 2026-09-20; status and implementation checked 2026-09-21.*

## Current position

Rung 0 is complete for K1, K2-cored and K2-NFW on the same 89 measurements.
The cored-minus-K1 evidence difference is +1.08 ± 0.60; NFW-minus-K1 is
−0.08 ± 0.59. This completes the intended baseline stage: the three families provide
a common reference, and their outer residuals identify the next modelling questions.
The [rung-0 analysis](RUNG0_ANALYSIS.md) motivates anisotropy and rotation tests.
Developing those controls is the planned research sequence.

Keep the observations, tracer and error/rotation options fixed within each
family comparison. Report goodness of fit as well as evidence.

## Implemented models

| Stage | Anisotropy | K1 / K2 sampled parameters | Status |
|---|---|---:|---|
| Rung 0 | β = 0 | 5 / 7 | Complete for all three families |
| Rung 1 | constant β, uniform [−1, 0.5] | 6 / 8 | Complete for K1, NFW and cored; see the rung-1 comparison |
| Simple companion | single-transition β(r), three parameters; rung-1 endpoint range | 8 / 10 | Complete for K1, NFW and cored |
| Legacy varying model | monotonic β(r), three parameters | 8 / 10 | Preserved as an optional control |
| Flexible rung 2 | two-transition β(r), five parameters | 10 / 12 | Complete for K1, NFW and cored |

All comparisons exclude instrument scales. Counts include distance. K1 fits
`M_star`, `M_rem`, `a_rem`, `M_bh` and the selected anisotropy; K2 adds
`M_dm_100` and `r_s`. The cored family is γ = 0 gNFW, not Burkert.

The [data/model comparison across all rungs](RUNG_MODEL_COMPARISON.md) shows
the best-fit no-DM and lowest-χ² DM models, including both rung-2 variants,
with the same axes and observations in each column.

The existing form is

`beta(r) = beta_0 + (beta_inf - beta_0) r² / (r² + r_beta²)`.

It can cross zero once, but cannot rise from a tangential centre to a radial
intermediate region and then fall to tangential outskirts. Rung 2 therefore
uses the two-transition extension below. The legacy form remains available
without changing archived models or their replay.

The completed [rung-1 comparison](RUNG1_COMPARISON.md) confirms the tension:
for the cored model, Gaia tangential χ² rises from 23.56 to 58.50 while HST
tangential χ² falls from 483.75 to 88.59. At fixed rung-1 masses and distance,
setting β to zero reduces Gaia tangential χ² to 14.93 but worsens HST strongly.
This motivates radial anisotropy freedom. The existing central prior below
prevents a monotonic profile from reproducing a tangential centre, radial
intermediate region and an outer decline. Rung 2 allows that pattern, without
requiring a peak or fixing the sign of intermediate or outer anisotropy.

### Simple companion: one smooth transition

The authorised intermediate comparison uses

`beta(r) = beta_in + (beta_out-beta_in) r²/(r²+r_a²)`.

The code names are `beta_0 = beta_in`, `beta_inf = beta_out`, and `r_beta = r_a`.
Both endpoint values have independent uniform priors on [-1, 0.5]; the radius
has a log-uniform prior on [0.5, 100] pc. The transition power is fixed at two.
At `r_a` the anisotropy is halfway between the endpoints. The profile can
increase or decrease, and becomes the rung-1 constant profile exactly when
the endpoints coincide. It cannot have an intermediate peak. This is a
fixed-shape member of the generalised Osipkov–Merritt family
([Baes & Van Hese 2007](https://arxiv.org/abs/0705.4109)).

This gives a controlled extension from one to three anisotropy parameters,
with the same diagnostic endpoint range as rung 1. It allows mild radial
anisotropy in the HST region to decline towards isotropy in the Gaia region.
The endpoint range deliberately does not enforce the central DF condition:
it preserves rung 1's diagnostic assumptions while testing radial variation.
A central-prior-restricted control is a separate experiment and is not part
of this launch. The completed turnover fits retain their different central
prior; their evidence comparison with the simple companion therefore tests
both the shape family and the stated prior choices.

The companion uses the exact saved rung-1 observations, MGE, mass/distance
priors, fixed parameters, and error/rotation settings. No instrument scales
are fitted. The launcher records the resolved priors explicitly, checks
saved-model reconstruction, and verifies that each parent maximum likelihood
is reproduced in the constant limit at `r_a = 0.5, 10, 100` pc. Global model
defaults and existing run archives retain their original settings.

### Adopted two-transition profile

Set `f_i(r) = r²/(r²+r_i²)`, with `r_2 > r_1`, and use

`beta(r) = beta_0 (1-f_1) + beta_mid (f_1-f_2) + beta_inf f_2`.

The three weights are non-negative and sum to one. Therefore the curve stays
within its control levels and has `beta(r) < 1`. It admits a radial maximum,
a tangential minimum or a monotonic transition. `beta_mid` is a control level,
not necessarily the attained peak. Transition powers remain fixed at two for
this first fit, avoiding two additional weakly constrained shape parameters.
The Jeans integrating factor is analytic:

`ln g = 2 beta_0 ln r + (beta_mid-beta_0) ln(r²+r_1²) + (beta_inf-beta_mid) ln(r²+r_2²)`.

The independent priors are:

| Parameter | Prior |
|---|---|
| `beta_0` | uniform [−1, −0.5] |
| `beta_mid`, `beta_inf` | uniform [−1, 1), independently |
| `r_beta = r_1` | log-uniform [0.05, 30] pc |
| `delta_r_beta = r_2-r_1` | log-uniform [0.5, 150] pc |

The positive gap orders the scales without a hidden prior rejection. The
inner transition can lie below the innermost kinematic bin; the outer
transition can lie within or beyond the Gaia coverage. These are explicit
modelling priors, not quantities inferred from the pilot simulations. Their
effect should be checked if posterior mass accumulates at a boundary.

A radial intermediate region followed by more isotropic or tangential
outskirts is physically motivated by evolution and escape in tidally limited
clusters ([Tiongco, Vesperini & Varri 2016](https://academic.oup.com/mnras/article/455/4/3693/1268368)).
The fit is allowed to choose another shape. It does not import a DM-particle
anisotropy profile as the stellar profile or fit away rotation by definition.

## Physical meaning of the baselines

The [22 September audit](DF_CONSISTENCY_AUDIT.md) checks all saved runs and
reruns. Spatial density is positive throughout. The completed rung-0,
rung-1 and simple-beta rung-2 posteriors fail the central DF condition below.
Turnover rung 2 passes that condition but its three best samples fail a
finite-radius necessary condition **if a separable augmented density is
assumed**. This does not rule out every nonseparable DF. The audit records the
distinction and posterior failure fractions; no fit is certified DF-positive.

For a cored tracer in a point-mass-dominated central potential, the
[An & Evans necessary condition](https://arxiv.org/abs/astro-ph/0511686) is
`beta(0) <= -1/2`. It applies to a constant β as well as a radially varying one.
The broad constant-β prior and the isotropic baseline deliberately relax this
condition to diagnose the observed profiles. They are Jeans moment models;
they are not guaranteed to admit a non-negative distribution function.

The legacy varying model uses `beta_0 ~ U[-1, -0.5]`, `beta_inf ~ U[-1, 1]`, and
`r_beta ~ logU[0.5, 100] pc`. Meeting the central bound alone does not prove
global DF positivity. The baseline central-point-mass estimates are conditional
on these assumptions too.

Rung 2 retains this necessary central condition because the current MGE is
cored asymptotically and the fitted potential includes an unsoftened point
mass. This condition depends on that central extrapolation; a different
tracer cusp or a no-BH model would require its own treatment. The new smooth
profile and positive Jeans moments do not certify global DF positivity or
stability. The existing JAM logistic and AGAMA Cuddeford implementations
cannot represent this profile, so requesting them is rejected explicitly.
Because the central prior differs from the diagnostic rung-1 prior, this is
not a strictly nested comparison across rungs; compare residuals as well as
evidence and retain that prior distinction when interpreting improvements.

## Commands and output handling

Use `python -m ocen_dm.cli` with `PYTHONPATH=src`, or the installed `ocen`
entry point. `ocen-dm` is not declared in `pyproject.toml`. These commands are
recipes, not a request to launch fits.

```bash
# Reproduction recipes, with fresh labels to preserve archived runs.
python -m ocen_dm.cli fit --family K1 --isotropic --no-scales --label rung0_K1_repeat
python -m ocen_dm.cli fit --family K2-cored --isotropic --no-scales --label rung0_K2_cored_repeat
python -m ocen_dm.cli fit --family K2-nfw --isotropic --no-scales --step-sampler --label rung0_K2_nfw_repeat

# Rung-1 reproduction recipes, with fresh labels.
python -m ocen_dm.cli fit --family K1 --constant-beta --no-scales --step-sampler --label rung1_K1_repeat
python -m ocen_dm.cli fit --family K2-cored --constant-beta --no-scales --step-sampler --label rung1_K2_cored_repeat
python -m ocen_dm.cli fit --family K2-nfw --constant-beta --no-scales --step-sampler --label rung1_K2_nfw_repeat

# Rung-2 reproduction recipes, with fresh labels.
python -m ocen_dm.cli fit --family K1 --anisotropy-profile turnover --no-scales --step-sampler --label rung2_K1_repeat
python -m ocen_dm.cli fit --family K2-nfw --anisotropy-profile turnover --no-scales --step-sampler --label rung2_K2_nfw_repeat
python -m ocen_dm.cli fit --family K2-cored --anisotropy-profile turnover --no-scales --step-sampler --label rung2_K2_cored_repeat
```

The authorised batch uses `bin/run_rung2_fits.py`: it replays each completed
rung-1 maximum likelihood, retains its exact MGE and 89-point observation
snapshot, checks 258 prior points per new model, and launches three workers
using archived source. The run labels are `rung2_K1_turnover`, `rung2_K2_nfw_turnover`
and `rung2_K2_cored_turnover`. Sampler settings are 400 live points,
`dlogz=0.5`, seed 42 and the slice sampler with `2*ndim` steps. Each process
uses one numerical-library thread. Launch manifests and logs are under
`results/launches/rung2_turnover_20260921/`; existing runs are never reused.

The simple companion uses the same driver with an explicit profile option:

```bash
PYTHONPATH=src python bin/run_rung2_fits.py launch results/launches/rung2_simple_beta_20260921 --profile simple
```

This selects `rung2_K1_simple_beta`, `rung2_K2_nfw_simple_beta`, and
`rung2_K2_cored_simple_beta`, with 8/10/10 sampled parameters. The sampler
settings remain 400 live points, `dlogz=0.5`, seed 42, and `2*ndim` slice
steps. Each worker uses archived source and one numerical-library thread.
The command refuses an existing batch or any existing output run. A bare
`fit --anisotropy-profile monotonic` uses the legacy priors and does not
reproduce this companion's endpoint priors.

All three companion workers passed launch verification and began sampling.
The batch contains `verification.json` with prior checks, constant-limit
likelihood replay, data identity, and source/configuration checksums, plus
`tests.log` recording 115 passing focused tests. Each family passed 258 finite
prior likelihood evaluations before launch. The existing turnover source
archive and configurations were also verified unchanged.

The completed replacements are `rung1_K1_slice` and `rung1_K2_cored_slice`; the completed
NFW run is `rung1_K2_nfw`. The replacements reconstruct the exact model and likelihood
arrays from the original run snapshots. Only the sampling method changes. The original
K1/cored runs are marked interrupted and retain their checkpoints; the replacements
start fresh. Their launch records are in
`results/repairs/numerics_20260921T125417Z/slice_launch.json`.

**Bare `fit` does not select rung 0.** Without flags it retains the monotonic
β(r) family and three legacy instrument scales (including Gaia DR2 even though
DR2 is absent from the default data). Specify the stage and a fresh label.
The driver refuses an existing run directory, even if it is empty or contains an
interrupted run. Resume of an existing project run is not implemented; preserve
the old directory and use a fresh label.

Current defaults are the composite tracer, local Jeans solver, raw Gaia errors
and published rotation. `--gaia-errors eta` and `--gaia-rotation ours` change
the Gaia components only. There is no HST rotation switch. A valid HST
alternative needs an external absolute rotation estimate: locally corrected
HST mean motion does not restore removed streaming. `--gaia-r-min` applies to
the legacy published EDR3 profile, not our prebuilt bins.

## Reporting and mock experiments

New runs write `run.yaml` and `data_snapshot.json` before the sampler starts.
They capture the resolved model and parameter order, priors, fixed values, distance,
halo taper, backend, tracer MGE, likelihood arrays, sampler settings and code identity.
Reports use these snapshots and check the saved best-sample likelihood before plotting.
Missing or altered snapshots raise an error instead of falling back to current products.

Legacy runs have no snapshots. Their reports load current profiles using the recorded
Gaia options and check the replayed likelihood. This preserves the existing rung-0
results without claiming to recover missing historical inputs. A legacy configuration
that cannot be reconstructed from its saved parameters is rejected explicitly.

Use a source run directory or its `ml_x.npy` for an injection. Its saved model defines
the generator independently of the fitted family and stage; the current dataset/selection
flags define the mock's bins, errors and streaming terms. For example, when requested:

```bash
python -m ocen_dm.cli fit --family K2-cored --constant-beta --no-scales \
  --mock-from results/fits/rung0_K1 --seed 42 --label mock_rung0_K1_fit_rung1_K2_cored
```

A standalone `.npy` without adjacent run metadata requires `--mock-family`; in that
case the tracer, backend and ladder flags also define the generator. Prefer the saved
run form when source and target configurations differ. New runs preserve the generated
arrays and full generating-model metadata. The historical symmetric Gaussian noise
with averaged errors and a positive floor is unchanged and explicitly recorded as
`gaussian_average_errors_clipped_v1`; it is not an exact split-normal recovery experiment.

Evidence comparisons use the saved likelihood arrays for new runs. Legacy comparisons
use file hashes, observational options and full mock provenance, and withhold numerical
differences when input identity is missing or inconsistent. Model settings such as the
anisotropy stage do not themselves make the observations different.

## Decisions after each stage

1. Inspect per-dataset and per-radius residuals. Constant β is a useful control,
   but is unlikely to describe the full radial pattern.
2. Compare family evidences after checking the likelihood inputs, assumptions
   and priors. A matching fingerprint is not proof of comparability.
3. Require an adequate fit before interpreting the halo posterior. A small
   evidence advantage among misspecified models is not a detection.
4. Test rotation, error calibration, prior sensitivity and injection recovery
   before adopting a DM constraint. The turnover must not simply absorb an
   injected halo signal.

Rung-0 K1, cored and NFW took 4.0, 10.0 and 20.6 hours, with 3.2, 6.6 and
14.3 million likelihood calls. Earlier estimates of minutes were too optimistic.
Use the step sampler as a convergence-tested alternative for the NFW degeneracy;
do not narrow a physical prior merely to speed up sampling.

## Instrument scales and other limitations

Scales remain available for one diagnostic after the main sequence. Its evidence
is not part of the comparison. A displaced scale should lead back to calibration,
selection, tracer populations or model inadequacy.

Rotation uncertainty is propagated per bin, although the curve uncertainty is
coherent. A shared curve would need a defensible covariance model. Other open
issues are sphericity, outer equilibrium, tracer populations and independence
of the profiles.

AGAMA supplies a limited DF diagnostic: it clips negative values and its
Cuddeford family differs from the Jeans anisotropy. Successful construction is
not a positivity certificate, and its evidence is not compared with these fits.

## Standing rules

- Launch fits only when explicitly requested, one stage at a time.
- Keep instrument scales out of the compared sequence.
- Preserve run directories and record changed inputs and options.
- Record evidence, per-dataset χ² and residual figures in the journal before
  choosing the next stage.
