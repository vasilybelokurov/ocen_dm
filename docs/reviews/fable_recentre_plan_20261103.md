# Fable review: finishing the centre marginalisation (plan of 2026-11-03)

Read-only. Code cited as file:line in oCen_dm. Numbers below come from
`results/plot_data/centre_marginal_20261103.json` (first run) and two short Python checks
(< 1 min each); no fits were run.

## Verdict (summary)

- **The plan should change.** Steps 1-3 (re-centre all 50 bins, Monte-Carlo centres, warm-started
  `dispersion_2d`) spend ~3-5 h of CPU converging an integral that cannot be converged and is not the
  right object. The per-centre log-likelihood of the *same* model varies with s.d. 5.2 (range 25)
  across the sigma_c = 2" prior; 99.5 % of that variation is the Poisson deviance and chi^2 of
  re-binned data (resampling noise), not centre physics. The prior-weighted sum is then carried by
  1-3 centres (effective sample size 1.3-3.5 out of 36-100 nodes). That, not Gauss-Hermite
  per se, is T4.
- **Recommended replacement: marginalise on the model side.** Keep every binned product exactly as
  frozen (bins and stars about AvdM10, errors as measured) and let the centre enter the *model
  prediction*: mu_k(theta, c) = a sum over the bin's coverage cells of Sigma_theta(|x - c|);
  sigma_k^2(theta, c) = mean over the bin's stars of sigma_theta,proj^2(|x_i - c|, rotated
  components). The integrand is smooth in c, so the existing n = 6 Gauss-Hermite grid converges
  (T8's unbinned, smooth counts term already does: 9.4 vs 9.6 at n = 6 vs 10), the cost per centre
  is ~10-100 ms of interpolation instead of 265 s, and the result is a bona fide
  marginal likelihood p(d | theta) = sum_k w_k p(d | theta, c_k) of fixed data. The marginalised
  objective then costs nothing extra in a refit.
- Q4 premise is wrong: the DF fits do **not** drop `streaming2`. `df_fit.py:291-292` strips it
  from the `ProfileLikelihood` copy and `df_fit.py:306-310` re-applies it after bin averaging.
  The frozen snapshot carries it (MUSE inner bins 53.5, 11.2, 34.6 (km/s)^2). Under the current
  data-side MUSE differencing this creates a real inconsistency (issue 4).
- The decision rule's "scatter < half the value" is not usable with data-side re-centring: the
  scatter is dominated by resampling noise (BH 4e4 - no BH per-centre difference ranges
  -14.6 .. +5.2 in ln L on the n = 10 grid).
- BH 8,200 conclusion (stable, -2 Delta ln L ~ -5.7) survives any of this; the per-centre
  difference to no BH has s.d. 0.26-0.63. The BH 4e4 and the central-deficit statements do not
  yet rest on a computable quantity.

## Ranked issues

### 1. The data-side marginal is a mixture over different statistics and is dominated by resampling noise (critical; verified)

What the code computes (`centre_marg.py:232-236`, `bin/centre_marginal.py:176-182`):
ln L_m(theta) = logsumexp_k [ln w_k + ln p(d(c_k) | theta, c_k)], where d(c_k) is a different
summary of the same stars for each k. For each fixed k this is a valid binned likelihood; the
weighted sum over k is not the likelihood of any fixed dataset. It is an average over K
(correlated) resamples of the inner bins.

Measured on the first run (`centre_marginal_20261103.json`, sigma_c = 2", n = 10, no-BH model,
s.d. over the 100 nodes):

| part of ln L | s.d. across centres |
|---|---|
| total | 5.19 |
| counts: -D/2 (model-dependent) | 4.23 |
| counts: saturated + ln Gamma terms (model-independent) | 0.26 |
| kinematics: -chi^2/2 | 2.51 |
| kinematics: split-normal normalisation | 0.32 |
| `lnl:ocen_counts_hst_f625w19` | 3.95 (range 19.6) |
| `lnl:ocen_counts_gaia_g17` | 0.81 (range 3.5) |
| `lnl:hst_pm_radial_ours` / `tangential` | 2.05 / 1.73 |
| `lnl:muse_los_dispersion` | 0.72 |

The Gaia-count row is the cleanest demonstration: a 2" shift has no physical effect at
300-2400", yet re-binning the Gaia counts per centre (`centre_marg.py:133-139`) injects a +-1.7
swing in ln L. All of it is noise; it is identical for the three models to first order, so it
cancels in a per-centre comparison but not in the marginal, where it reweights the centres.

Effective sample size of the weights w_k exp(ln L_k) (per model, per grid):

| grid | no BH | 8,200 | 4e4 | top weight |
|---|---|---|---|---|
| s1_n6 (36) | 2.1 | 2.0 | 1.5 | 0.67-0.82 |
| s2_n6 (36) | 2.9 | 3.2 | 1.3 | 0.53-0.87 |
| s3_n6 (36) | 1.5 | 1.5 | 3.5 | 0.34-0.82 |
| s2_n10 (100) | 1.5 | 1.5 | 3.3 | 0.47-0.82 |

Each model's marginal is the likelihood at its own luckiest resample (for BH 4e4 on s2_n6 the
dominant node is (1.23, 1.23)", for no BH it is (-1.23, 3.78)"). The n = 6 and n = 10 grids pick
different lucky nodes, hence T4 = +3.9 for the model whose per-centre differences are largest.

Why Monte Carlo cannot fix it (answer to Q2): ln L_k across centres is, to a good approximation,
mu + noise with s = 5.2. Estimating E[exp(ln L)] of a lognormal with that s needs K >> exp(s^2)
~ 10^11 draws; with K = 400 the estimate sits at the single largest draw and keeps rising with K.
Two independent sets of 200 would disagree by several units, after ~5 h of CPU (see issue 3).
"Non-smoothness" is the mechanism that makes Gauss-Hermite nodes behave like random draws; the
ESS collapse is the reason no node set converges.

Statistical reading for Q1: the construction does not double count the measurement noise, but it
is not a marginal likelihood either. In the idealised limit where the K resamples of a bin were
independent draws n_k ~ N(0, eps^2), E_k[L(theta | mu* + n_k)] is the likelihood *convolved* with
the resampling distribution: inner-bin errors effectively inflated by up to sqrt(2) (under-stating
the information), and the finite-K estimate of that convolution is a tail statistic dominated by
the best few draws (over-stating it, randomly, per model). Both effects are present here; the
observed inner `hst_pmr` scatter across centres (0.051, 0.047, 0.027 mas/yr in the three innermost
bins vs AvdM10 errors ~0.03) is that resampling noise.

### 2. Model-side marginalisation removes issue 1 at near-zero cost (recommendation; verified feasible, not implemented)

Fix the data at AvdM10 (every bin, every error, every `r_nodes` as frozen) and move the centre
into the prediction:

- Counts: mu_k(theta, c) = a sum_{cells j in bin k about AvdM10} Sigma_theta(|x_j - c|) dA
  (field b and amplitude a profiled per theta as now, `counts.py:188-196`). The coverage point set
  already exists (`centre_marg.py:73-96`, `px, py, dA`); inside 30" use the same point set rather
  than the analytic annulus. Gaia counts: identical, full annuli. The prediction is smooth in c.
- HST PM dispersions: for bin k the frozen estimator is an (approximately) unweighted average over
  its stars, which is what `r_nodes` already encodes (`hst_profile.py:84-96`,
  `likelihood.py:437-447, 458-462`). Generalise it to two dimensions:
  sigma_R'^2 = (1/N_k) sum_i [sigma_R^2(R_i') cos^2 dphi_i + sigma_T^2(R_i') sin^2 dphi_i],
  with R_i' = |x_i - c| and dphi_i the angle between the AvdM10 and c radial directions at star i
  (sigma_T'^2 with cos and sin swapped). Star positions come from `load_hst_sample`
  (`hst_profile.py:129-153`, r and phi about AvdM10). Subsample to ~5,000 stars per bin for
  R > 15"; inside 15" use all. Everything is interpolation on one 1-D table per theta (as
  `model_lnsigma`, `bin/centre_marginal.py:97-104`).
- MUSE: same per-star average of sigma_los^2(|x_i - c|) over the `selection_hq_pm_and_los` stars
  of each published bin (`centre_marg.py:82, 94-96`); no estimator differencing, no T9.
- `streaming2`: keep the frozen values (issue 4 explains why that is now consistent).
- Marginal: the existing `marginal_loglike` over the existing Gauss-Hermite grid. T4 becomes a
  real convergence test (expect < 0.1 between n = 6 and n = 10, as for the unbinned T8 term).
- Per centre cost: one interpolation pass over ~10^5 positions, ~10-100 ms. The marginalised
  objective can therefore go straight into refits with no precomputation.

Approximation it keeps: the Gaussian mixture estimator applied to a bin whose stars have a
spread of true radii and a rotation of the (R, T) axes by dphi (up to ~35 deg in the 3.24-4.13"
bin at |c| = 2"); its ML variance is close to the mean variance when the spread is modest.
Test it once with a mock: draw velocities from a saved model about c = (2, 0)", bin about AvdM10
with `hst_profile`, compare with the model-side prediction; accept at < 0.2 of the bin error in
every bin.

Checks that carry over: T1 (trivially, data unchanged), T2/T3 (c = 0 must reproduce the frozen
`loglike_profiled` to 1e-9: the per-star average at c = 0 must equal the `r_nodes` average,
which it does only if the same stars and equal weights are used; otherwise expect a ~1e-3
difference and say so), T5/T6 become moot, T7 is replaced by "prediction shift per bin at
+-2 sigma_c" which is the physical sensitivity (expect < 0.1 sigma beyond 60").

### 3. If the data-side route is kept anyway: the cost and warm-start estimates are wrong (verified in code)

`dispersion_2d` (`outer_profile.py:485-649`): the 4-start Nelder-Mead (lines 94-102, maxiter 300)
is not the dominant cost. The two profile-likelihood interval scans (lines 105-157) run 2 sides x
up to ~8 bracketing steps x a Nelder-Mead of up to 60 iterations per component, i.e. ~1,000-2,000
`loglike` evaluations, comparable to or more than the multi-start. Dropping three starts saves
at most ~1/3, not 4x. Wall time for 300 centres at 265 s on 6 cores is ~3.7 h before the warm
start, ~2.5 h after; the brief's "~1 h" is off by 2-4x. A larger saving would be to keep the
AvdM10 errors rescaled by sqrt(N_0/N(c)) and skip the interval scans (the per-centre error change
is a few per cent; T1 shows errors reproduce to 5e-10), validated at 8 centres. But issue 1
makes the whole route pointless.

Warm start fairness (Q3): fair in principle (same objective, keep the better of warm/cold), and
the proposed checks are adequate. The risk is not the start but the sigma_max = 1.5 clip and the
field-solution runaway documented at lines 91-93, which matter only beyond 300".

### 4. MUSE rotation term is inconsistent under data-side re-centring; fine under model-side (verified by a quick check)

`centre_marg.py:150-175, 211-215`: the per-centre MUSE value is published + [ours(c) - ours(0)],
where "ours" is a dispersion about a rotation a sin phi + b cos phi re-fitted *about c*, while
`streaming2` stays at the published v_rot^2/2 about AvdM10. In the three inner bins (~50 stars
each) the re-fitted rotation amplitude swings strongly with a 2" shift:

| centre (dRA*, dDec)" | bin 1 v_rot, v^2/2 | bin 2 | bin 3 |
|---|---|---|---|
| (0, 0) | 10.9, 59 | 4.3, 9 | 7.8, 30 |
| (2, 0) | 6.4, 21 | 6.3, 20 | 7.5, 28 |
| (-2, 0) | 9.2, 42 | 8.7, 38 | 7.9, 31 |
| (0, 2) | 7.8, 31 | 6.5, 21 | 1.0, 0.5 |
| (0, -2) | 11.1, 62 | 7.5, 28 | 5.5, 15 |
| (-2, -2) | 5.4, 15 | 13.1, 85 | 6.2, 19 |

(km/s and (km/s)^2; frozen `streaming2` = 53.5, 11.2, 34.6.) The model prediction
sqrt(sigma_model^2 - streaming2) with sigma_model ~ 21 km/s changes by ~1 km/s (0.5 published
sigma) between streaming2 = 10 and 53. So under data-side re-centring the streaming term would
have to be re-fitted per centre, and it is itself noise-dominated with 50 stars. Under
model-side marginalisation the data (and their rotation correction) are fixed and
centre-independent, and the aperture shift changes <vbar^2> only at second order
((2/7.9)^2 ~ 6 % of 53 ~ 3 (km/s)^2, negligible). HST: locally corrected PMs
(`hst_profile.py:142-145` global median, `likelihood.py:186-193`) and the per-bin fitted mean in
`dispersion_2d` (lines 47-74) make the HST dispersions independent of the centre choice;
`pm_rotation_curve(r_median)` for the tangential term changes negligibly. No per-centre
recomputation needed in either case for HST.

### 5. The central-deficit statement mixes two different quantities (verified)

`bin/centre_marginal.py:163, 183-184`: P(N <= n_obs) is averaged over centres with *both* n_obs
and mu re-binned; n_inner2 ranges 14-48 over the prior. "P = 0.027 at AvdM10, 0.30-0.37
marginalised" therefore compares a fixed-data probability with an average over resampled data,
and the latter is not a tail probability of anything observed. Model-side: n_obs = 18 fixed,
mu(theta, c) marginalised over the prior, P = sum_k w_k Poisson-cdf(18 | mu_k). Report that.
(Whether a count deficit "exists" at a position is a question about the point pattern; the
unbinned T8 likelihood is the right tool and already shows the 4e4 count penalty is robust.)

### 6. Smaller points

- Decision rule: with model-side marginalisation the per-centre scatter is pure centre sensitivity
  and the "< half the value" clause is meaningful; keep it. Add: report the posterior ESS of the
  centre weights per model (should be >> 1) and the posterior mean centre per model (a model that
  needs a 3" offset to fit is telling you something).
- Model-independent normalisation terms (ln Gamma(n+1), split-normal ln norm) vary by only 0.36
  across centres: not the problem, but they vanish entirely with fixed data.
- Correlated bins: with the centre as a model parameter, neighbouring bins are correlated through c
  in the usual nuisance-parameter way; nothing extra to do.
- `bin/centre_marginal.py:131-134` rounds centres to 12 decimals to deduplicate; fine.
- The BH position = cluster centre assumption stands (earlier reviews); unchanged by this plan.

## Answers to the seven questions

1. **Soundness.** As implemented, no: the prior-weighted sum of likelihoods of centre-dependent
   summaries is a mixture over statistics, not p(d | theta). It does not double count the
   measurement noise; it convolves the inner-bin likelihood with the resampling distribution
   (information loss) and its finite-K estimate is a tail statistic of that convolution
   (random gain, per model). Measured: ln L of one model varies with s.d. 5.2 over the prior,
   model-independent parts 0.36, ESS 1.3-3.5. The right object is p(d | theta) with fixed d and
   the centre inside the model (issue 2), or the unbinned p(x | theta, c) (T8 for counts).
2. **MC/QMC.** Not viable for the data-side integrand (lognormal with s ~ 5 needs ~10^11 draws;
   200-400 draws give the maximum draw, which rises with K). The ESS of w_k L_k is the right
   diagnostic and it already says 1-3. T4 is explained by the ESS collapse; non-smoothness only
   makes the nodes behave as random draws. With model-side predictions the integrand is smooth
   and Gauss-Hermite n = 6 vs 10 is a genuine convergence test (T8 shows 9.4 vs 9.6).
3. **Warm start.** Fair, checks adequate, but the saving is <= 1/3 because the interval scans
   dominate `dispersion_2d`; "~1 h" is really 2.5-4 h on 6 cores. Cheaper: keep AvdM10 errors
   rescaled by sqrt(N_0/N). Cheapest and better: no per-centre rebuilds at all (issue 2).
4. **Rotation.** The premise "DF fits drop streaming2" is false (`df_fit.py:291-292, 306-310`;
   snapshot carries it). Data-side MUSE differencing is inconsistent with the fixed streaming2
   (inner-bin v_rot^2/2 swings 0.5-85 (km/s)^2 with 2" shifts, ~0.5 sigma effect on the
   prediction). Under model-side marginalisation keep all streaming terms fixed; HST needs
   nothing either way.
5. **Gaia.** Kinematic bins: skip (second order, (2/300)^2). Counts: with data-side re-binning
   they add +-1.7 of pure noise to ln L (measured), so do not re-bin them; with model-side they are
   free and effectively constant. One check at 2" is enough.
6. **Unbinned.** Deferring the per-star kinematic likelihood is right (model change, 1-2 days).
   The middle path is model-side marginalisation of the existing binned products: the data stay
   binned (same estimators, same errors, same frozen fits), only the centre enters the model. For
   counts the unbinned fixed-window term (T8) is already the clean version; keep it as the
   cross-check and, if wanted later, as the R < 60" count likelihood.
7. **What would invalidate the result as planned.** (a) Issue 1: no finite K gives a converged
   number; whatever is reported for BH 4e4 is the luckiest resample. (b) Issue 4 for MUSE.
   (c) Issue 5 for the deficit probability. (d) Decision-rule scatter conflating noise with
   sensitivity. All four disappear with fixed data and a model-side centre.

## Proposed plan (replacement for steps 1-4)

1. Implement `predict_offset(theta, c)` for counts (coverage cells), HST PM (per-star 2-D
   average, rotated components) and MUSE (per-star average); data products untouched. Tests:
   c = 0 reproduces `loglike_profiled` of the frozen problem to 1e-9 (or state the ~1e-3 node
   difference); mock at c = (2, 0)" from a saved model: estimator vs prediction < 0.2 bin errors;
   Gaia prediction change at 2" < 1e-3 of the error.
2. Marginal over the existing grids (n = 6, 10; sigma_c = 1, 2, 3"). Tests: T4 < 0.1; ESS per
   model > 10; per-model posterior mean centre and its spread reported.
3. Report -2 Delta ln L_m, per-centre scatter (now pure sensitivity), and the fixed-n_obs deficit
   probability. Apply the decision rule.
4. Refit no BH and 4e4 with the marginalised objective (cheap now); mocks at 0-3" offsets with and
   without a BH as the earlier reviews asked. Launch only after asking.

Compute: minutes, 1 core.

## Verified vs assumed

Verified (read or ran): `centre_marg.py` (all), `bin/centre_marginal.py` (all),
`tests/test_centre_marg.py`, `hst_profile.py:84-153, 181-229`, `outer_profile.py:485-649`,
`likelihood.py:172-203, 275-306, 437-515`, `df_fit.py:276-317`, `counts.py:178-196, 308-310`,
`docs/CENTRE_MARGINALISATION.md`, both earlier reviews; the ESS / decomposition / correlation
numbers from the first-run JSON; the MUSE rotation-amplitude table from the catalogue with the
same `_ml_sigma` form as `centre_marg.py:150-163`.

Assumed (not run): the ML mixture estimator's variance equals the mean variance over a bin to
better than ~0.2 errors under a 2" offset (test proposed); `model.projected_moments` is cheap
enough to tabulate once per theta on ~600 radii (it is used that way in
`bin/centre_marginal.py:97-104`); the 265 s per-centre figure and the "4 NM starts" cost split
(argued from the code structure, not profiled); the oMEGACat MUSE bins are about AvdM10 (as the
project's notes state).
