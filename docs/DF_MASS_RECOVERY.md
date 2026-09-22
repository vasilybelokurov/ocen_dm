# Free-potential DF mass recovery

The [fixed-potential challenge](DF_CAPACITY_CHALLENGE.md) established that three
positive action DFs can reproduce the light and velocity moments of the
independent equilibrium mocks. The next experiment fits the gravitational
potential as well. It separates stellar mass and light weights, then asks
whether the known total and component masses are recovered.

The [subsequent literature review](DF_MODEL_LITERATURE_REVIEW.md) considers
more compact stellar families in the same stars-plus-DM problem. The present
mixture remains a flexibility and recovery experiment; the capacity test does
not establish that the real data require three stellar basis components.

The current batch is
`results/df/mass_recovery_20260922T151408Z/`.
See the [generated status report](DF_MASS_RECOVERY_STATUS.md) for the latest
saved fits and plots. Results from running or unconverged fits are preliminary.

## Model and the meaning of its weights

For three positive unit-integrated action shapes `F_i(J)`, the model uses

\[
 f_{\rm mass}(\boldsymbol J)=M_\star\sum_{i=1}^3 w_iF_i(\boldsymbol J),\qquad
 f_{\rm light}(\boldsymbol J)\propto\sum_{i=1}^3 \ell_iF_i(\boldsymbol J).
\]

Both weight sets are positive and separately normalized. Their ratio `w_i/l_i`
is a relative mass-to-light ratio for a basis component. The basis components
are not assumed to identify individual physical stellar populations. Summing
different weights over the same positive orbital shapes permits a radial
change in the stellar mass-to-light ratio.

Each trial solves the stellar Poisson problem from **the mass DF**, including
the trial dark gravity in the potential. Photometry and all three projected
velocity moments are then computed from **the light DF** in that converged
potential. Absolute light normalization is profiled out. For numerical
convenience the light DF integrates to `M_star`, but its amplitude is a
luminosity proxy; it does not supply the gravitational density. The model
retains separate `mass_df` and light `df` objects.

The implementation is `src/ocen_dm/kinematics/df_mass_recovery.py`. Its
`MassLightDFModel` extends the existing solver without changing the default
behaviour of `PositiveDFModel`. The latter now optionally accepts a starting
stellar potential and additional static potentials. A warm starting potential
changes the initial guess only: the same convergence, mass and density
closure guards still apply.

## A controlled first recovery test

The truth is unchanged: the independently constructed positive energy-based
cluster with compact stars, extended stars and dark remnants, with and without
the extended dark component. All observations sample the same light-weighted
mixture. The distance remains 5.43 kpc, there is no black hole or rotation, and
the radial selection/binning remains that of the 89 kinematic and 82 photometric
points used in the capacity challenge.

The remnant and halo profile families in this first test contain the injected
shapes. A tabulated unit-mass density `u(r)` can change its total mass `M` and
radial scale `a` independently:

\[
 \rho(r;M,a)=\frac{M}{a^3}u(r/a),\qquad \int4\pi r^2u(r)\,dr=1.
\]

The saved tables are generated from the physical mock and interpolated on
1,200 radii. Their normalization is divided out before the free mass is applied.
Thus the optimizer receives the profile family but not a fixed dark mass or
radius. Both masses include zero in their allowed range. In the no-halo truth,
the fit still permits a halo, using the halo-case template as its shape family.

This is deliberately an optimistic, matched-family mass-recovery control.
It avoids making the first test simultaneously a test of an incorrect halo
profile. A successful result would not establish robustness to unknown
remnant/halo shapes. Those need subsequent shape-mismatch experiments.
The template mass and scale parameters have their own explicit configuration
fields; they do not reuse the legacy `M_dm_100` meaning.

The stellar mass and light DFs are positive and self-consistent. Fitted dark
templates have positive spatial density; their phase-space DFs have **not**
been certified in each new fitted potential. This retains the stated scope
of the existing stellar-DF branch. The generating mocks do have positive DFs
for all their gravitating components.

## Search and numerical controls

Each mock has three noiseless starting points with different stellar masses,
dark masses and relative mass/light weights. Their initial action shapes come
from the previous known-potential capacity fits; this is a favourable
initialization, not a blind search from an unknown tracer structure. All shape
coordinates become free in the final stage.

The optimizer works through three stages:

| Stage | Parameters allowed to vary | Number |
|---|---|---:|
| Weights | Stellar mass, two stellar mass ratios, two light ratios, remnant and halo masses/scales | 9 |
| Scales | The above plus each stellar component's action scale and outer orbital coefficient | 15 |
| Full | The above plus each component's inner/outer slopes, transition steepness and inner orbital coefficient | 27 |

Each stage has at most 60 main residual evaluations. Finite-difference derivative
probes are additional physical model calls and are recorded separately. Failed
equilibrium or accuracy checks reject a trial; negative DFs or moments are not
clipped into valid observations. Accepted points use cold, reproducible
equilibria. Derivative probes start from the same converged reference potential
and retain the same accuracy checks. Explicit finite steps resolve the numerical
accuracy of the equilibrium calculation instead of using machine-scale
differences through an iterative solver.

The base numerical settings are 160 potential nodes, 32 nodes per velocity
dimension, 180 radial moment nodes and 80 projection nodes, with equilibrium
convergence tolerance `1e-5`. Every final selected model is rebuilt with all
four resolutions doubled and a stricter convergence tolerance. It must change
the predicted data by less than 0.05 adopted errors. Native AGAMA projected
integrals must agree within 0.5%, and a cold replay of the saved model must
reproduce its residuals within `1e-8` adopted errors. Stellar mass and density
closure are checked during every trial.

The controller runs two workers with one numerical thread each. Source files,
tests, observational snapshots and capacity-fit seeds are copied into the batch
before jobs launch; workers import the frozen scientific source. Each job saves
its true model, normalized dark templates, noise realization, starting model,
best model, complete evaluation histories, stage results and numerical checks.
Existing run directories are refused. Progress JSON files are replaced atomically.

## Noisy experiments and what they can establish

The initial six noiseless fits run first. A mock proceeds to the next part only
if at least one full-stage solution meets an optimizer termination criterion,
passes numerical validation and meets the previously defined approximation
target at doubled resolution: RMS below 0.1 errors and every bin below 0.3.
A failed gate is recorded explicitly; it does not silently launch noisy fits
from an inadequate model.

For each case passing the gate, the controller adds a shared-M/L control and
eight independent Gaussian noise realizations, with two starts per realization.
The planned complete batch therefore contains 40 fits: six noiseless fits,
two shared-M/L controls and 32 noisy fits. Noise uses the explicitly adopted
kinematic error scales and `0.1 mag / sqrt(relative photometric weight)`.
It is calibrated only for that synthetic error model; the real photometric
weights are not measured uncertainties. The seeds and noise vectors are saved.

The main quantities are the total enclosed mass profile, stellar and dark
component masses, and halo density at 20 pc. Several comparably good fits with
different component masses would indicate ambiguity in the decomposition,
even if the total mass profile is recovered. The shared-M/L control measures
how much that restriction changes the result.

Eight noisy realizations per case provide an initial view of recovery bias and
scatter. They do **not** constitute a posterior-coverage calibration. No posterior
sampler, evidence calculation or new fit to the real cluster data is launched
by this batch. Formal interval coverage requires a larger, converged inference
ensemble after these controls establish an adequate model and search.

## Verification before launch

The 111 focused tests passed. These cover the new independent mass/light weights,
unchanged gravity under light reweighting, reproduction of the previous model
when weights agree, mass-DF density closure, warm/cold equilibria, exact dark
template mass/radius scaling, zero dark-mass boundaries and reproducible noise,
as well as the earlier DF, binning, consistency, anisotropy and replay controls.
A subsequent focused rerun of the seven new tests also passed after adding a
check of the shared-M/L coordinate path.

Both mock starting models passed end-to-end numerical validation. In a separate
four-evaluation optimizer smoke test, with action shapes held fixed, the halo
case's score fell from 13.647 to 0.0901. Its worst doubled-resolution residual
was 0.104 adopted errors, and the largest numerical refinement shift was
0.00138 errors. This intentionally capped, unconverged smoke test validates
the optimization path; it is not the final mass-recovery result.

These checks are saved under `results/df/mass_recovery_validation_20260922/`,
`results/df/mass_recovery_validation_no_dm_20260922/` and
`results/df/mass_recovery_smoke_20260922/`. The batch controller launch record
and logs are in `results/diagnostics/df_mass_recovery/`.

## Commands

```bash
python bin/run_df_mass_recovery_batch.py \
  --out results/df/new_mass_recovery_batch \
  --workers 2 --noise-realizations 8 --max-evaluations 60

python bin/report_df_mass_recovery.py \
  --batch results/df/mass_recovery_20260922T151408Z
```

The reporter writes PNGs to `plots/`, plotted data and provenance to
`results/plot_data/`, and the generated Markdown status page beside this note.
`--watch` refreshes that status once per minute while the batch is running.
