# Why the compact DF fails, and the proposed anisotropy extension

Status: 23 September 2026. The diagnostics below are complete. The model
adjustment in the final section was approved after external review and launched. This follows
[OBSERVED_DF_FLEXIBILITY.md](OBSERVED_DF_FLEXIBILITY.md). That test showed
the current DF cannot fit the combined Omega Cen data: the best fit has
chi2_kin/N = 4.80 with every valid start converged.

## Summary

- **The current DF fits each part of the data on its own.** HST+MUSE alone
  reach chi2_kin/N = 1.43, and Gaia alone 0.88. Only the combination fails.
- **The two parts demand opposite anisotropy.** The inner data (HST+MUSE)
  want radial bias switching on at J ~ 30 pc km/s. The outer data (Gaia)
  want strong tangential bias switching on at J ~ 1700 pc km/s. The current
  DF has one transition, so it cannot supply both.
- **An independent Jeans fit on the same data needs the same shape.** It
  wants radial beta ~ +0.15 near 6 pc turning tangential beyond ~20 pc,
  and reaches chi2 = 190 against the DF's 456.
- **The combined verdict is sensitive to the Gaia errors, but not only
  through them.** A 0.02 mas/yr floor lowers chi2_kin/N from 5.12 to 2.09.
  Gaia tangential still misses (chi2/n = 8.2), and J_a stays at its bound.
- **Distance is not the cause.** Freeing it gives D = 5.30 kpc and
  chi2_kin/N = 4.78.
- **Proposed next model:** the already implemented second anisotropy
  transition, with the halo frozen at zero. Starts come from the subset
  fits, and the ordering is reparametrized so that J_outer > J_a always.

## What was tested

All fits use the current compact regularized DF with no halo (rho20 = 0). The
data snapshot is the same as before (89 kinematic bins, 82 photometric
points). Photometry is always included. Each batch uses two hand-picked
starts and wide search bounds (J0 30-800, J_a 5-3000 pc km/s). Fits use
parallel base-seeded Jacobians, iteration tolerance 5e-5, an absolute stop
of 0.01, and a 600-evaluation cap. Runner: `bin/run_df_step1_diagnostics.sh`.
Batches: `results/df/dfdiag_*_20260923/`, commit `8a7e377`.

| Batch | Kinematic data | Change | Purpose |
|---|---|---|---|
| `dfdiag_hst_muse` | HST PM R, HST PM T, MUSE LOS (71 bins) | - | Can the DF fit the inner data alone? |
| `dfdiag_gaia` | Gaia PM R, Gaia PM T (18 bins) | - | Can the DF fit the outer data alone? |
| `dfdiag_gaiafloor01` | all 89 | Gaia errors ⊕ 0.01 mas/yr | Error-floor sensitivity |
| `dfdiag_gaiafloor02` | all 89 | Gaia errors ⊕ 0.02 mas/yr | Error-floor sensitivity |
| `dfdiag_freedist` | all 89 | distance free, 4.8-6.0 kpc | Distance check |

The floors are added in quadrature. They test sensitivity to spatially
correlated Gaia EDR3 systematics
([Lindegren et al. 2021](https://arxiv.org/abs/2012.03380);
[Vasiliev & Baumgardt 2021](https://arxiv.org/abs/2102.09568)); they are
not calibrated error models. At 5.43 kpc, 0.01 mas/yr is 0.26 km/s.

## Results

Both starts of every batch converged to the same solution (equal chi2 to
0.1), and all solutions passed numerical validation. The reference row is the
earlier all-data, no-halo wide-bound fit.

| Fit | chi2_kin / N | chi2_kin/N | HST R | HST T | MUSE | Gaia R | Gaia T | chi2_phot/N |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| All data (reference) | 456.1 / 89 | 5.12 | 2.1 | 6.2 | 2.7 | 2.6 | 20.2 | 2.95 |
| HST+MUSE only | 101.6 / 71 | **1.43** | 1.0 | 1.1 | 2.0 | - | - | 2.78 |
| Gaia only | 15.9 / 18 | **0.88** | - | - | - | 0.8 | 1.0 | 2.53 |
| All, Gaia floor 0.01 | 288.1 / 89 | 3.24 | 1.3 | 2.4 | 2.3 | 0.8 | 15.3 | 2.88 |
| All, Gaia floor 0.02 | 185.8 / 89 | 2.09 | 1.1 | 1.2 | 2.1 | 0.3 | 8.2 | 2.82 |
| All, distance free | 425.1 / 89 | 4.78 | 1.8 | 6.6 | 1.6 | 2.6 | 19.6 | 2.94 |
| *Jeans turnover, raw errors (for comparison)* | *190.0 / 89* | *2.13* | *2.9* | *1.7* | *1.8* | *1.9* | *2.6* | *(separate MGE)* |

The dataset columns give chi2/n (n = 21, 21, 29, 9 and 9 bins).

Fitted parameters. b_out > 0 is radial bias; the DF factor is
exp(-b sin(pi c/2)) in orbit circularity c.

| Fit | M_star [Msun] | J0 | alpha | b_out | J_a [pc km/s] | M_rem [Msun] | a_rem [pc] | D [kpc] | At a bound |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| All data (reference) | 3.33e6 | 109 | 1.12 | +0.18 | 5 | 1.3e5 | 0.82 | 5.43 | J_a |
| HST+MUSE only | 3.12e6 | 120 | 1.24 | **+0.39** | **32** | 1.7e5 | 1.31 | 5.43 | none |
| Gaia only | 2.75e6 | 129 | 1.14 | **-1.00** | **1665** | 6.6e5 | 2.61 | 5.43 | b_out (lower) |
| All, Gaia floor 0.01 | 3.25e6 | 112 | 1.16 | +0.27 | 5 | 1.4e5 | 0.98 | 5.43 | J_a |
| All, Gaia floor 0.02 | 3.19e6 | 116 | 1.20 | +0.33 | 5 | 1.5e5 | 1.14 | 5.43 | J_a |
| All, distance free | 3.10e6 | 105 | 1.13 | +0.18 | 5 | 1.2e5 | 0.80 | 5.30 | J_a |

### Intrinsic anisotropy

Figure: `plots/df_subset_beta_20260923.png` (data in
`results/plot_data/df_subset_beta_20260923.json`, script
`bin/plot_beta_profiles.py`). It shows beta = 1 - sigma_t^2/sigma_r^2 of the
refined models together with the maximum-likelihood rung-2 Jeans turnover
fits on the same data snapshot.

- **HST+MUSE only:** beta = 0 in the centre (the DF bias vanishes as J -> 0),
  rising to +0.24 at 10 pc and +0.5 at 80 pc.
- **Gaia only:** beta ~ 0 inside 10 pc, turning tangential beyond (-0.08 at
  20 pc, -0.16 at 30 pc, -0.70 at 80 pc).
- **All data:** a compromise, rising monotonically to +0.3 at 80 pc. It
  matches neither.
- **Jeans turnover:** -0.5 at the centre (at its prior limit; the ML value is
  -0.503, so the data prefer less tangential), a radial peak of +0.15 near
  6 pc, zero near 20 pc, and -0.6 to -0.7 at 80 pc. It follows the HST+MUSE
  DF inside ~5 pc and the Gaia DF outside ~25 pc.

## Analysis

1. **The failure is a conflict between datasets within one model.** Each
   subset is fitted well by the same DF family with different anisotropy
   parameters. The combined fit pins J_a at its lower bound: it applies one
   anisotropy amplitude at all actions and satisfies neither subset.
2. **The conflict is radial, not a handover artefact.** The Gaia fit's
   tangential turn and the Jeans zero crossing lie at ~20 pc, well inside
   the Gaia range (9.4-63 pc). They are not at the 9.5 pc HST/Gaia
   boundary. In the one overlapping bin (9.0-9.4 pc) the two surveys agree
   within errors (sigma_T 0.465 +- 0.027 vs 0.474 +- 0.049 mas/yr).
3. **The Gaia errors set how strongly the conflict is penalized.** Raw Gaia
   dispersion errors are 0.005-0.007 mas/yr at 17-35 pc (1-2%), so the
   18 Gaia bins dominate the combined fit. A 0.02 mas/yr floor halves
   chi2_kin/N and lets the inner fit relax (HST T 6.2 -> 1.2). Even so,
   Gaia tangential misses by ~3 of its inflated errors per bin, and J_a
   remains at its bound. The structural conflict survives any floor tested.
   Its size, and hence any final goodness-of-fit verdict, depends on an
   unresolved Gaia systematic error model.
4. **Distance cannot resolve it.** D rescales both PM components equally
   and cannot change sigma_T/sigma_R. Freeing it lowers chi2_kin by 31 with
   D = 5.30 kpc. That is 2.6 sigma below the project's adopted prior of
   5.43 +- 0.05 kpc, taken from
   [Baumgardt & Vasiliev 2021](https://arxiv.org/abs/2105.09526) (per-cluster
   value not re-checked here), while Gaia tangential stays at 19.6.
5. **Mass is not yet the question, but there is a hint.** The Gaia-only fit
   puts 6.6e5 Msun in an extended remnant component (a_rem = 2.6 pc),
   against 1.7e5 Msun for HST+MUSE. Once the anisotropy conflict is
   removed, the outer data may still require more mass at large radii. That
   question belongs to the later halo step, not to this one.

### Visual versus statistical fit

Data-model figures such as `plots/observed_df_smoke_20260923_fits.png` can
look acceptable in the top panels while failing badly. PM errors of
0.005-0.01 mas/yr are 1-2% of the dispersion, so a 5 sigma miss (~0.03
mas/yr) is invisible at the plotted scale. Fits must be judged from the
residual panels and chi2 per dataset. The smoke figure shows a 12-evaluation
unconverged run (chi2_kin/N 8.5-11.6), not a converged fit.

## Conclusion

The current DF is too rigid in one specific, identified way. It has a single
anisotropy transition, so its velocity ellipsoid can only go from isotropic
to one asymptotic state. The Omega Cen data require isotropic or mildly
radial bias at small actions, radial bias at intermediate actions
(J ~ 30-1000 pc km/s, r ~ 2-15 pc), and tangential bias at large actions
(J >~ 1700 pc km/s, r >~ 20 pc). Subset fits, an independent Jeans analysis
and the error-floor tests all point to this. Distance, the optimizer, the
search bounds and the HST/Gaia boundary do not explain it. The Gaia error
model sets how large the misfit looks, not whether it exists.

## Model adjustment: second anisotropy transition (approved)

### Model

Enable the second anisotropy transition already implemented in
`ExponentialParameters` (`src/ocen_dm/kinematics/regularized_df.py:46`). It
is covered by `test_two_transitions_have_correct_limits` and
`test_two_transition_velocity_anisotropy`:

    b(J) = b_out s(J; J_a) + (b_outer - b_out) s(J; J_outer),
    s(J; J_t) = J^2 / (J^2 + J_t^2),        J_outer > J_a,

so that b -> 0 for J << J_a, b ~ b_out between the two scales, and
b -> b_outer for J >> J_outer. The target is b_out > 0 (radial) and
b_outer < 0 (tangential). Nothing else changes: mass follows light, no
rotation, fixed distance, and **no halo** (rho20 = 0). Freezing the halo
keeps it from absorbing DF inadequacy.

### Action coverage of the data

`plots/action_coverage_observed_df_wide_20260923.png` (script
`bin/plot_action_coverage.py`) shows 200,000 tracers sampled from the
all-data no-halo DF. The bias acts on q = J_r + L. The 5/50/95% quantiles of
q per projected-radius range, in pc km/s:

| Projected radius | q 5% | q 50% | q 95% |
|---|---:|---:|---:|
| 0.1-1 pc | 22 | 101 | 320 |
| 1-3 pc | 38 | 119 | 311 |
| 3-9.5 pc | 81 | 187 | 390 |
| 9.4-20 pc | 174 | 303 | 517 |
| 20-63 pc | 313 | 476 | 734 |

Two consequences for the parametrization:

- The Gaia-only J_a = 1665 lies above the 95% point of even the outermost
  bin. That fit uses only the low-q tail b ~ b_out (q/J_a)^2, a smooth ramp
  in which only b_out/J_a^2 is constrained. A second transition must lie
  inside q ~ 20-800 to be identified.
- Neighbouring radial ranges overlap strongly in q, and one q does not map to
  one radius. Any transition in q is soft in radius. Its interpretation must
  come from p(q | R) and beta(r), not from a J-to-r conversion.

### Adopted setup (revised after external review, 23 September)

Code: coordinate `stellar.log_J_outer_ratio` = Delta = ln(J_outer/J_a),
decoded after all other coordinates (`df_fit.config_at`), so that
J_outer > J_a for every trial (unit-tested). Driver option:
`prepare-real --two-transition`.

| Coordinate | Bounds | Start 0 | Start 1 |
|---|---|---:|---:|
| M_star [Msun] | 1e6-6e6 (log) | 3.1e6 | 2.8e6 |
| J0 [pc km/s] | 30-800 (log) | 120 | 130 |
| alpha | 0.6-3 | 1.24 | 1.14 |
| b_out | -2 to +2 | +0.39 | +0.30 |
| J_a [pc km/s] | 5-500 (log) | 32 | 50 |
| b_outer | -4 to +2 | -1.0 | -1.5 |
| Delta = ln(J_outer/J_a) | ln 2 - ln 100 | ln 15 (J_outer 480) | ln 8 (J_outer 400) |
| M_rem [Msun] | 1e2-1e6 (log) | 1.7e5 | 4e5 |
| a_rem [pc] | 0.3-20 (log) | 1.3 | 2.0 |

The J_a upper bound drops from 3000 to 500 because larger values lie outside
the sampled q. Four Latin-hypercube starts (seed 20260923) are added. There
is no halo, raw errors are used first, and the 0.02 mas/yr Gaia floor is
repeated only if the raw fit is structurally successful. A 12-evaluation
smoke run from start 0 validated (cold shift 9e-5, refinement 7e-4 errors,
projection <= 0.18%) with no rejected evaluations.

### Acceptance criteria, split (revised)

The formal Gaia errors (1-2%) are known to be incomplete, so structure and
statistics are judged separately:

- **Structural pass:** independent starts converge (Delta chi2 < 1 for at
  least half of the valid starts); no coordinate at a bound; no coherent
  residual pattern; the implied beta(r) is radial at 2-15 pc and tangential
  beyond ~20 pc; no dataset with chi2/n > 3; M_rem and a_rem checked for
  collapse toward the HST+MUSE values (1.7e5 Msun, 1.3 pc).
- **Statistical pass:** chi2_kin/N < 1.3, every dataset chi2/n < 2, and
  chi2_phot/N < 2. This is required only under a justified Gaia systematics
  model. With raw errors it is reported, not required.

After the run, for every selected fit: beta(r), sigma_r(r) and sigma_t(r);
b(q) with p(q | R) and both transitions marked; and the DF in the (J_r, L)
plane, checked for artificial ridges. The Jeans comparison is a shape
target only. A positive equilibrium DF need not reproduce an arbitrary
Jeans beta(r).

### Decision after the run

- **Passes:** add the cored halo (step 3) and ask whether rho20 > 0 is still
  preferred. The passing no-halo and halo fits then define realistic mock
  truths (v3) for the fixed-rho20 profile scan and the star/DM
  identifiability test.
- **Fails only with raw Gaia errors:** the result depends on the Gaia
  systematic error model. Build a calibrated floor from the Gaia
  astrometric systematics before going further.
- **Fails in both runs:** a rotating DF with mean velocities fitted jointly
  (step 4), then separate tracer and mass DFs (step 5). An axisymmetric DF
  is the last resort (step 6).

## Result of the two-transition, no-halo fit (23 September)

Batch `results/df/dftwo_20260923/` (commit `f6ce414`), 6 starts, raw errors.
Two Latin starts were infeasible at their initial point and recorded as such.
The other four reached the same minimum (refined objective 378.3-378.4,
Delta chi2_kin <= 0.6). Starts 0 and 1 terminated by the absolute rule after
502 and 492 evaluations; starts 3 and 5 reached the cap at the same point.
No evaluation was rejected, and no coordinate is at a bound. Validation
passed (start 0: cold shift 0.0010, refinement 0.0003 errors, projection
0.03%).

| Quantity | One transition | **Two transitions** |
|---|---:|---:|
| chi2_kin (N = 89) | 456.1 | **159.2** |
| chi2_kin/N | 5.12 | **1.79** |
| HST R / HST T / MUSE, chi2/n | 2.1 / 6.2 / 2.7 | **1.0 / 1.8 / 2.0** |
| Gaia R / Gaia T, chi2/n | 2.6 / 20.2 | **0.9 / 3.7** |
| chi2_phot/N | 2.95 | 2.67 |

Best-fit parameters: M_star = 3.25e6 Msun, J0 = 65 pc km/s, alpha = 0.85,
b_out = +0.46, J_a = 27 pc km/s, b_outer = -1.26, J_outer = 812 pc km/s,
M_rem = 1.72e5 Msun, a_rem = 1.27 pc.

Figures: `plots/dftwo_20260923_best_fits.png` (data and residuals),
`plots/dftwo_beta_20260923.png` (beta vs the one-transition DF and Jeans),
and `plots/df_structure_dftwo_20260923_observed_no_halo_start0.png`
(intrinsic moments, b(q) over sampled actions, f(J_r, L); script
`bin/plot_df_structure.py`).

- **Anisotropy:** beta(r) = 0.02, 0.14, 0.16, 0.06, -0.13, -0.55 and -1.27
  at 1, 6, 10, 20, 30, 45 and 63 pc. Inside ~30 pc this reproduces the Jeans
  turnover shape (radial peak near 6-10 pc, zero near 20 pc). Beyond ~40 pc
  it falls faster (-2.1 at 80 pc, outside the data). There sigma_t exceeds
  sigma_r.
- **Bias vs coverage:** b(q) peaks at +0.41 at q = 103 pc km/s, where the
  HST/MUSE tracers sit (q 5/50/95% = 56/170/379). It crosses zero at
  q = 488, inside the Gaia tracers (183/347/650). Both transitions act on
  sampled actions. The negative asymptote b_outer is reached only beyond the
  Gaia q range, so the steep outer beta is extrapolation.
- **DF shape:** f(J_r, L) is smooth and monotonic, with no ridge at either
  transition scale.
- **Remnants:** M_rem and a_rem return to the HST+MUSE-only values
  (1.7e5 Msun, 1.3 pc). The 6.6e5 Msun extended remnant of the Gaia-only fit
  was compensating for the missing tangential outskirts.
- **Remaining residuals:** tangential PM -2.5 errors at 4-7 pc and +3.5 at
  15-20 pc (first Gaia bins). LOS is 2-3.6 errors high at 5-8 pc (MUSE outer
  edge, where the rotation subtraction is largest). The three innermost HST
  PM bins (0.1-0.2 pc) lie above the model; there is no central point mass.

**Against the criteria:** the structural criteria pass except one. Starts
converge, nothing is at a bound, the beta shape is correct, the DF is smooth,
and the remnants collapse. The exception is Gaia tangential at chi2/n = 3.7,
against a limit of 3, with a mildly coherent residual at 4-20 pc. The
statistical criterion is not met (chi2_kin/N 1.79 > 1.3); with raw Gaia
errors it is reported, not required.

### Decision on Gaia error floors (23 September)

The 0.01/0.02 mas/yr floors above were ad hoc, not a calibrated error model.
They are withdrawn as a tool for judging models and will not be used again.
All fits use the published (raw) errors. The incompleteness of the Gaia
error model stays an open caveat. The proper remedy is a separate data task:
re-derive the Gaia dispersion profiles with a literature-based per-star error
calibration, or calibrate against HST in the overlap. Its method needs
approval before implementation. Step 3 (two transitions plus free halo)
therefore runs with raw errors.

## Step 3: two transitions plus free cored halo (23 September)

Batch `results/df/dftwo_halo_20260923/`, raw errors. Three of six starts
were infeasible at their initial point. Two failed numerically (action
normalization error 5.1e-5 above the 2e-5 tolerance; a density-closure
failure) and one physically (streaming exceeded the DF second moment). All
three valid starts reached one minimum: refined objective 346.9 (no halo:
378.3), chi2_kin 143.3 (159.2). By dataset (chi2/n): HST R 1.02, HST T
1.39, MUSE 1.83, Gaia R 2.04 (was 0.93), Gaia T 2.34 (was 3.7); chi2_phot/N
is 2.48.

Not a clean DM result:
- **Two coordinates are at bounds.** r_s = 500 pc (upper), so the halo is
  effectively uniform (rho20 = 0.92 Msun/pc^3) inside the data.
  J_outer/J_a = 2 (lower), so the transitions are pushed together.
- **The anisotropy re-arranges** (b_out 0.46 -> 1.76, J_a 27 -> 101,
  b_outer -1.26 -> -0.25): a mass-anisotropy trade.
- **The halo is < 1% of the enclosed mass at 20 pc**, 8% at 42 pc and 32% at
  80 pc. It again trades misfit between the Gaia R and T components.

The gain (Delta objective 31, Delta chi2_kin 16, for two extra parameters)
cannot be read as a DM preference. A DM constraint needs a rho20 profile
likelihood with all other parameters refitted, then mock validation.

## Central point mass: is there evidence? (23 September)

`plots/central_bh_evidence_20260923.png` (script
`bin/plot_central_bh_evidence.py`, data in `results/plot_data/`) compares
the 43 kinematic bins inside 2 pc: HST PM R and T, and MUSE LOS.

| Model | central chi2 / 43 | HST R | HST T | MUSE |
|---|---:|---:|---:|---:|
| DF two transitions, no BH, no halo | **48.8** | 11.7 | 12.5 | 24.6 |
| Jeans turnover + BH (M_bh 4.4e4), no DM | 59.0 | 21.5 | 16.7 | 20.8 |
| Jeans turnover + BH (M_bh 4.4e4), cored DM | 59.3 | 23.7 | 14.1 | 21.4 |

The binned data give no evidence that the DF needs a point mass. Only the
HST bins at 0.12-0.16 pc lie 1-2 errors above the no-BH DF, and the
innermost bin (0.10 pc) lies below it. The comparison is not controlled,
since the Jeans fits also differ in anisotropy (beta_0 = -0.5 imposed) and
mass model. A matched DF with and without a BH would be the clean test. The
binned profiles start at 0.1 pc, inside r_infl ~ 0.4 pc for 4e4 Msun, and
contain none of the individual fast central stars behind the IMBH claim of
Haeberle et al. 2024 ([arXiv:2405.06015](https://arxiv.org/abs/2405.06015)).
The BH is therefore dropped from the objective-1 list.

## Best-fit model performance (objective 1: fit the data)

Figures: `plots/best_fit_performance_20260923_data.png` (all data, both
two-transition fits and the rung-2 Jeans turnover fit, residuals
(data-model)/err) and `plots/best_fit_performance_20260923_diagnostics.png`
(residual histogram, chi2/n, runs test, beta, enclosed mass). Script:
`bin/plot_best_fit_performance.py`; numbers in
`results/plot_data/best_fit_performance_20260923.json`. Published errors are
used and bins are treated as independent, so the p-values are indicative
only.

| Dataset (n) | No halo chi2/n (p) | + halo chi2/n (p) | Runs z, no halo / halo | max abs residual |
|---|---|---|---|---|
| HST PM R (21) | 1.01 (0.45) | 1.02 (0.43) | -0.9 / 0.0 | 2.0 / 2.0 |
| HST PM T (21) | 1.81 (0.013) | 1.39 (0.11) | -1.1 / -1.1 | 2.6 / 2.1 |
| MUSE LOS (29) | 2.01 (0.001) | 1.84 (0.004) | +0.5 / +0.5 | 3.6 / 3.2 |
| Gaia PM R (9) | 0.93 (0.50) | 2.04 (0.03) | -2.3 / -1.2 | 1.9 / 2.2 |
| Gaia PM T (9) | 3.70 (1e-4) | 2.34 (0.01) | -0.8 / -1.0 | 3.6 / 2.7 |
| Photometry (82) | 2.67 (2e-14) | 2.48 (3e-12) | -1.2 / -2.1 | 5.7 / 5.8 |

The runs test counts sign runs along radius; z < -2 means too few runs,
i.e. coherent residuals.

- **The kinematic fits are now better than the best Jeans fit.** chi2_kin is
  159.2 (no halo) and 143.3 (halo), against 190.0 for the rung-2 Jeans
  turnover (with BH, no DM) on the same bins. All profiles are tracked
  across 0.1-63 pc. The kinematic residual distribution is somewhat wider
  than N(0,1), with tails to about 3.6.
- **Remaining structured misfits:** MUSE LOS at its outer edge, 5-8 pc
  (model 2-3.6 errors high; largest streaming correction); PM tangential at
  4-7 pc (model 2-2.6 high) and 15-20 pc, the first Gaia bins (model up to
  3.5 low). With a halo, the outer Gaia PM R becomes about 2 errors too
  high. Only Gaia R (no halo) and the photometry (halo) pass z = -2 in the
  runs test.
- **Photometry:** chi2/N 2.5-2.7, dominated by 10-40 pc points with
  residuals of opposite sign at the same radii (-3 to -6 and +2 to +4). This
  matches the mutually inconsistent Trager-compilation points: duplicate
  full-weight pairs differ by 0.13 mag, against 0.1 mag adopted errors. This
  points to heterogeneity in the photometric data rather than the model. It
  is a hypothesis; it has not been tested.
- **The two models agree inside the data:** beta(r) and the total M(<r)
  coincide to ~25 pc. The halo adds mass only beyond ~30 pc (8% at 42 pc).
  Inside the data, the data do not distinguish the two.

Verdict on objective 1: a structural near-pass; not statistically acceptable
under the published errors. The remaining misfit is localised (MUSE outer
edge and the 4-20 pc tangential PMs), not a global failure.
