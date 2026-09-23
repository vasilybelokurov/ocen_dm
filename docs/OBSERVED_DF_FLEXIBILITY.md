# Can the compact regularized DF fit the Omega Cen data?

Status: 23 September 2026. Tested result; batches
`results/df/observed_df_20260923/` (default bounds) and
`results/df/observed_df_wide_20260923/` (wide bounds, multi-start).

## Summary

- **No.** The best fit to all 89 kinematic bins has chi2_kin/N = 4.80, against
  an acceptance threshold of 1.3. The outer Gaia proper motions miss by up to
  6 adopted errors.
- **Optimizer and search bounds are excluded as causes.** With wide bounds, all
  10 valid starts converged: five free-halo starts to one minimum (refined objective
  623.3) and five no-halo starts to another (698.1). These include eight
  Latin-hypercube starts spread across the parameter box.
- **The limitation is the anisotropy structure.** The anisotropy action scale
  J_a reaches its lower bound in every fit, whether that bound is 50 or 5.
  The fit wants one anisotropy amplitude at all actions and still cannot
  match the outer tangential and radial PM dispersions together.
- **A cored halo does not rescue the fit.** It lowers the total objective by
  75 but only trades misfit between the Gaia radial and tangential PMs. Both
  models fail, so this difference is not evidence for DM.
- **Needed, in order:** fitter robustness fixes; a fit with free distance;
  the implemented second anisotropy transition; then, if needed, separate
  mass and light DFs and a rotating DF. Model changes await approval.

## Question and scope

The question is whether the current stellar DF family can reproduce the whole
observed dataset, irrespective of whether it splits mass correctly between
stars and DM. This is a flexibility test of the model, not a mass measurement.
A model that cannot fit the data cannot support an inference of the DM
content. The mass split is examined separately in
[COMPACT_DF_RECOVERY.md](COMPACT_DF_RECOVERY.md), where noiseless mocks show a
strong star/halo degeneracy.

## Model

One spherical, self-consistent equilibrium per trial
([REGULARIZED_DF_IMPLEMENTATION.md](REGULARIZED_DF_IMPLEMENTATION.md)):

- **Stars:** one regularized exponential action DF, which serves as both the
  mass and the light DF (mass follows light). The five free parameters are
  stellar mass M_star, action scale J0, profile exponent alpha, anisotropy
  amplitude b_out, and anisotropy action scale J_a. The anisotropy bias is
  zero (isotropic) as J goes to 0 and approaches b_out for J much greater
  than J_a. b_out is a DF parameter, not the velocity anisotropy beta.
- **Dark remnants:** a Plummer sphere with free mass M_rem and scale a_rem.
- **Halo (free-halo branch only):** a cored profile (gamma = 0) with free
  density rho20 at 20 pc, free scale r_s, and a fixed Gaussian taper at
  1000 pc. The no-halo branch fixes rho20 = 0.
- **Fixed inputs:** distance 5.43 kpc, no rotation in the DF. The published
  rotation is subtracted at moment level: sigma_model = sqrt(<v^2> - vbar^2),
  and trials where this is infeasible are rejected.

## Data

| Dataset | Kind | Bins | Projected radii |
|---|---|---:|---|
| HST PM, radial (`hst_pm_radial_ours`) | PM dispersion | 21 | 0.1-9.5 pc |
| HST PM, tangential (`hst_pm_tangential_ours`) | PM dispersion | 21 | 0.1-9.5 pc |
| MUSE (`muse_los_dispersion`) | LOS dispersion | 29 | 0.15-8.7 pc |
| Gaia EDR3 PM, radial (`gaia_edr3_ours_radial`) | PM dispersion | 9 | 9.4-63 pc |
| Gaia EDR3 PM, tangential (`gaia_edr3_ours_tangential`) | PM dispersion | 9 | 9.4-63 pc |
| Surface brightness | mag | 82 | 0.05-68 pc |

Radii are converted at 0.02633 pc per arcsec (5.43 kpc). Kinematic errors are
the published asymmetric errors, used in a split-normal chi2. The photometry
splices oMEGACat star counts inside 25 arcsec with the Trager, King &
Djorgovski (1995) V-band compilation outside. Its errors (0.1-0.58 mag) are
**adopted weights**, not calibrated measurements. Ten pairs of full-weight
points at equal radius scatter by 0.13 mag per point. The photometric zero
point is profiled.

## Test design

Driver: `bin/run_compact_df_recovery.py prepare-real`. Data snapshot:
`results/df/pilot_no_dm_20260922/`.

| Setting | Default-bound batch | Wide-bound batch |
|---|---|---|
| Starts per branch | 2 hand-picked | 2 hand-picked + 4 Latin hypercube (seed 20260923, central 60%) |
| J0 bounds [pc km/s] | 100-800 | 30-800 |
| J_a bounds [pc km/s] | 50-3000 | 5-3000 |
| r_s bounds [pc] | 5-150 | 5-500 |
| Other bounds | M_star 1e6-6e6, alpha 0.6-3, b_out -1 to 1, M_rem 1e2-1e6, a_rem 0.3-20 pc, rho20 0-10 | same |
| Jacobian | serial, base-seeded | 5 parallel probe workers, base-seeded |

Both batches use bounded trust-region least squares (scipy `least_squares`,
TRF) in coordinates scaled to [0, 1]. Finite-difference steps are 0.002. The
equilibrium iteration tolerance is 5e-5, chosen after a Jacobian noise probe.
A fit also stops when the objective falls by less than 0.01 over two accepted
steps. The cap is 600 evaluations. Each selected solution is replayed cold
and rebuilt at doubled resolution. Its fast projected moments are checked
against direct projection over the data radii.

**Acceptance criteria**, fixed before the production results. The photometric
criterion was amended before any production result, because no model can
meet an unweighted RMS below 0.05 mag given the 0.13 mag duplicate scatter.
The criteria are:

1. chi2_kin/N < 1.3 over the 89 kinematic bins;
2. chi2/n < 2 in every kinematic dataset;
3. chi2_phot/N < 2 with the adopted errors (a perfect model is expected near
   1.7);
4. the optimizer terminates below the evaluation cap;
5. numerical validation passes (cold and refinement shifts < 0.1 errors;
   projection agreement < 0.5%).

## Results

### Convergence

| Batch / branch | Start | Evals | Stop | Refined objective | Parameters at a bound |
|---|---|---:|---|---:|---|
| default, free halo | 0, 1 | 140, 244 | xtol | 653.6, 653.6 | J0 low, J_a low, r_s high |
| default, no halo | 0, 1 | 104, 146 | ftol/xtol | 717.6, 717.6 | J_a low |
| wide, free halo | 0, 1, 3, 4, 5 | 199-349 | ftol/xtol/absolute | 623.3 (all five) | J_a low |
| wide, free halo | 2 | 118 | xtol (stuck) | 72870 | none (see Defects) |
| wide, no halo | 0, 1, 3, 4, 5 | 96-237 | ftol/xtol | 698.1 (all five) | J_a low |
| wide, no halo | 2 | 1 | failed | - | start infeasible (see Defects) |

Every valid start reaches its branch minimum: five of five in each wide-bound
branch and two of two in each default-bound branch. All selected solutions
pass numerical validation. Maximum cold shift is 0.0046 errors, maximum
refinement shift 0.012 errors, and maximum projection discrepancy 0.11%.

### Best fits (wide bounds)

| Quantity | Free halo | No halo |
|---|---:|---:|
| chi2_kin (N = 89) | 427.5 | 456.1 |
| chi2_kin/N | **4.80** | **5.12** |
| HST PM radial, chi2/n (n = 21) | 2.1 | 2.1 |
| HST PM tangential, chi2/n (n = 21) | 4.5 | 6.2 |
| MUSE LOS, chi2/n (n = 29) | 2.6 | 2.7 |
| Gaia PM radial, chi2/n (n = 9) | 10.1 | 2.6 |
| Gaia PM tangential, chi2/n (n = 9) | 13.6 | 20.2 |
| chi2_phot/N (N = 82) | 2.39 | 2.95 |
| M_star [Msun] | 3.21e6 | 3.33e6 |
| J0 [pc km/s] | 90 | 109 |
| alpha | 1.02 | 1.12 |
| b_out | 0.21 | 0.18 |
| J_a [pc km/s] | 5.0 (bound) | 5.0 (bound) |
| M_rem [Msun], a_rem [pc] | 1.6e5, 1.21 | 1.3e5, 0.82 |
| rho20 [Msun/pc^3], r_s [pc] | 1.50, 476 | 0 (fixed) |
| M(<24 pc) total [Msun] | 3.07e6 | 3.07e6 |

Both minima fail criteria 1-3. Widening the bounds improved the objective
only from 653.6 to 623.3 (free halo) and from 717.6 to 698.1 (no halo).

### Residual pattern

Figure: `plots/observed_df_wide_20260923_best_fits.png`. Data and provenance:
`results/plot_data/observed_df_wide_20260923_best_fits.json`. The script is
`bin/plot_observed_df_fits.py`.

- **PM tangential:** the model is too high at 3-7 pc (residuals to -6 errors)
  and too low beyond 12 pc (to +6 errors).
- **PM radial:** the model is too high beyond about 20 pc (to -5 errors with a
  halo, -3 without).
- **LOS:** the model is 3-5 errors too high at 5-8 pc, at the outer edge of
  the MUSE field.
- **Inside about 2 pc** all kinematic profiles are fitted within about
  2 errors.
- **Photometry:** fitted with chi2/N 2.4-3.0; the largest residuals are
  low-weight Trager points.

The same pattern appears in the default-bound fits
(`plots/observed_df_20260923_fits.png`).

## Conclusions

1. The current compact regularized DF, with mass following light, one
   anisotropy transition, no rotation and a fixed distance, cannot reproduce
   the combined HST, MUSE, Gaia and photometric profiles. The best attainable
   chi2_kin/N is 4.8.
2. The optimizer is not the cause: every valid start converges to the same
   branch minimum. The search bounds are not the cause either: widening them
   changes chi2_kin by only 29 (free halo) and 20 (no halo).
3. The tension is between the inner and outer velocity-ellipsoid shape. J_a
   runs to its lower bound, i.e. one anisotropy amplitude everywhere. The
   model still cannot make the tangential dispersion low at 3-7 pc and high
   beyond 12 pc while keeping the radial and LOS dispersions low there.
4. The halo's 75-unit objective gain comes from moving misfit between the Gaia
   components. Neither model is acceptable, so this gain is not a DM
   detection and must not be quoted as one.

### Caveats

- The kinematic chi2 treats bins as independent. Datasets that share stars
  (HST radial and tangential; the two Gaia components) are correlated, so
  absolute chi2/N values are approximate. The robust statements are the
  failure pattern and the reproducibility of the minima.
- The rotation subtraction is a moment-level approximation. The outer MUSE
  and Gaia bins, where rotation matters most, are where the model fails.
- Photometric errors are adopted weights. The photometric verdict is
  indicative only.
- Gaia EDR3 dispersion systematics (error underestimation, spatially
  correlated PM errors) have not been re-examined here.

## Defects found and their status

| Defect | Effect | Status |
|---|---|---|
| A rejected equilibrium inside a Jacobian enters as a 1e6 residual | Corrupts the Jacobian; the trust region collapses to a spurious xtol stop (wide free halo start2 at objective 72870) | To fix: retry the probe in the opposite direction; flag stops near rejections as unconverged |
| An infeasible starting point crashes the job ("base potential unavailable") | Start lost with a misleading error (wide no halo start2: streaming exceeded the DF second moment) | To fix: record as `start_infeasible` |
| Frozen batches still print the old unweighted photometric RMS gate | Report field only | Fixed in the working driver (`data_fit_gates`) |

Neither defect affects the converged minima above.

## What is needed

Revised on 23 September after an external review, checked in `JOURNAL.md`.
Model changes (steps 2-6) await approval; steps 0-1 do not change the model.

| Step | What | Question answered |
|---|---|---|
| 0 | Fitter robustness: retry a rejected probe in the opposite direction; flag stops near rejections; record infeasible starts | Removes the two defects above |
| 1a | Spline-beta(r) Jeans inversion (a new anisotropy class in `jeans.py`), with and without Gaia | Which beta(r) do the data demand? |
| 1b | Current DF fitted to HST+MUSE only and to Gaia only | Do the datasets need incompatible DFs? |
| 1c | Gaia systematics: HST/Gaia overlap bins, systematic error floor | Is the outer tension partly systematic? |
| 1d | Distance free | Sanity check only: D rescales both PM components equally and cannot change sigma_T/sigma_R |
| 2 | Two-transition anisotropy (`b_outer`, `J_outer`; b > 0 is radial), **no halo**, multi-start; plot the implied beta(r) | Can a radial-then-tangential ellipsoid fit the data? |
| 3 | Step 2 plus cored halo | Is rho20 still preferred once the DF is flexible? |
| 4 | Rotating DF (odd in L_z), with mean velocities fitted jointly | Removes the annular, spherical streaming approximation, which is largest where the fit fails (Gaia vbar^2/sigma^2 up to 0.38) |
| 5 | Separate tracer and mass DFs | Relaxes mass follows light (mass segregation, multiple populations) |
| 6 | Axisymmetric f(J_R, J_z, L_z) | Needed if spherical models cannot absorb flattening and rotation |

The halo stays frozen at zero while DF flexibility is being tested. Only a
model that passes the observed-data criteria should define mock truths (v3)
for the profile scan and the star/DM identifiability test.

## Reproduction

```sh
PY=/Users/vasilybelokurov/Work/venvs/.venv/bin/python
O=results/df/CHOOSE_FRESH_NAME
$PY bin/run_compact_df_recovery.py prepare-real --out $O --stop-delta 0.01 \
    --jacobian-seed base --iteration-tolerance 5e-5 --probe-workers 5 --n-starts 4 \
    --bound stellar.J0=30:800 --bound stellar.J_a=5:3000 --bound matter.r_s=5:500
$PY bin/run_compact_df_recovery.py run --out $O --workers 2 --max-calls 600 --detach
$PY bin/plot_observed_df_fits.py $O --jobs observed_free_halo_start0 observed_no_halo_start0 --suffix _best_fits
```

Batches freeze their source, configurations, tests and inputs under
`<batch>/code` and `<batch>/observed`. The wide batch ran from commit
`5486309`. With the other workloads on the machine paused, 12 jobs took
about 2 h.
