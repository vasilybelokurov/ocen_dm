# Rung 0: analysis of the three isotropic fits

*2026-09-21. Runs `rung0_K1`, `rung0_K2_cored`, `rung0_K2_nfw` in `results/fits/`; comparison
table `results/fits/comparison.md`; figures `plots/fit_posterior_profiles.png`,
`plots/fit_rung0_*_posterior_profiles.png`. Plan: `docs/MODELLING_PLAN.md`; data:
`docs/data_analysis.pdf`; method: `docs/mass_modelling.pdf`.*

**Baseline stage complete (2026-09-21).** These runs establish the intended reference
for the next modelling steps. All three saved maximum-likelihood values reproduce
exactly with the current code and data. Their residuals guide tests of anisotropy,
rotation and other extensions. No new posterior was sampled during the audit or
subsequent reporting repairs; see [CODE_ANALYSIS_AUDIT.md](CODE_ANALYSIS_AUDIT.md).

## 1. What was run

Rung 0 of the ladder: the simplest defensible model. Isotropic (β ≡ 0), no instrument scale
factors, composite tracer, our own Jeans solver, 89 binned points:

| dataset | points | radial range |
|---|---|---|
| `hst_pm_radial_ours` | 21 | 3.7–342″ (0.1–8.8 pc) |
| `hst_pm_tangential_ours` | 21 | 3.7–342″ |
| `muse_los_dispersion` | 29 | 5.7–290″ |
| `gaia_edr3_ours_radial` | 9 | 357–2149″ (9.2–55 pc) |
| `gaia_edr3_ours_tangential` | 9 | 357–2149″ |

Switches recorded in `run.yaml`: `gaia_errors: raw`, `gaia_rotation: published`,
`isotropic: true`, `no_scales: true`. Three families: K1 (stars + remnants + point mass, 5
parameters), K2-cored and K2-NFW (K1 plus a two-parameter halo, `M_dm_100` and `r_s`).

## 2. Results

| | K1 (no DM) | K2 cored | K2 NFW |
|---|---|---|---|
| ln Z | −190.02 ± 0.48 | **−188.94 ± 0.36** | −190.09 ± 0.35 |
| max ln L | −169.7 | −164.4 | −165.5 |
| χ² at max L (89 points) | 772 | 761 | 763 |
| M★ [M☉] | 2.88 [2.67, 3.03] ×10⁶ | 2.30 [1.68, 2.77] ×10⁶ | 2.34 [1.57, 2.85] ×10⁶ |
| M_rem [M☉] | 3.1 [2.0, 4.8] ×10⁵ | 7.7 [3.9, 13.3] ×10⁵ | 7.1 [3.2, 13.4] ×10⁵ |
| a_rem [pc] | 4.0 [3.2, 4.7] | 5.2 [4.3, 5.7] | 5.1 [4.1, 5.7] |
| M_BH [M☉] | 4.34 [3.90, 4.71] ×10⁴ | 4.65 [4.29, 4.97] ×10⁴ | 4.53 [4.11, 4.88] ×10⁴ |
| M_DM(<100 pc) [M☉] | — | 1.32 [0.19, 2.47] ×10⁶ | 5.3 [0.23, 11.9] ×10⁵ |
| r_s [pc] | — | 261 [46, 663] | 123 [7, 544] |
| D [kpc] | 5.311 [5.291, 5.332] | 5.312 [5.291, 5.333] | 5.313 [5.292, 5.334] |
| 95 % upper limit M_DM(<100 pc) | — | 3.13 ×10⁶ | 1.64 ×10⁶ |
| likelihood calls / wall time | 3.2 ×10⁶ / 4.0 h | 6.6 ×10⁶ / 10.0 h | 14.3 ×10⁶ / 20.6 h |

**Δ ln Z**: cored − K1 = +1.08 ± 0.60; NFW − K1 = −0.08 ± 0.59; cored − NFW = +1.16 ± 0.50.
The cored-to-K1 Bayes factor is about 2.9. **Neither halo family is persuasively preferred
within this baseline.** Because all three models fit poorly, these evidence differences
cannot establish the physical presence or absence of a halo.

**Input-hash discrepancy.** The report flags the runs because the Gaia product's hash
changed between launches. Re-evaluating the three saved best samples now reproduces their
stored log likelihoods exactly and all per-dataset χ² values to their recorded precision.
This supports the journal's interpretation of a metadata-only change. It is not a proof
that every original likelihood input was identical: the run directories do not contain
immutable copies of the fitted arrays. The report warning is retained. Launch-time hashes
now reduce this risk, but do not freeze inputs or cover the tracer source files.

## 3. Goodness of fit: the excess χ² is concentrated beyond 50″

χ² at the K2-cored maximum-likelihood point, split by projected radius:

| dataset | R < 50″ | 50–150″ | > 150″ |
|---|---|---|---|
| HST radial | 12.9 (11 pts, 1.2/pt) | 49.5 (5) | 90.1 (5) |
| HST tangential | 13.6 (11 pts, 1.2/pt) | 165.8 (5) | 304.4 (5, 61/pt) |
| MUSE LOS | 9.5 (13 pts, 0.7/pt) | 27.5 (10) | 38.9 (6) |

Inside 50″ the residuals are consistent with the adopted uncertainties. Most of the excess
χ² comes from larger radii; HST tangential contributes 484 of the total 761.

The misfit is small in physical terms. Fractional residuals (data − model)/data:

| | 108″ | 138″ | 176″ | 224″ | 271″ | 311″ |
|---|---|---|---|---|---|---|
| HST radial | +0.3 % | +1.4 % | +1.2 % | +1.1 % | +0.2 % | −1.2 % |
| HST tangential | −2.0 % | −3.7 % | −5.1 % | −7.5 % | −9.6 % | −10.0 % |
| statistical error | 0.26–0.30 % | 0.24–0.35 % | 0.20–0.49 % | 0.16–0.79 % | 0.21–1.19 % | 0.51–1.61 % |

χ² = 761 corresponds to a model wrong by 1 % in one component and 5–10 % in another. It looks
catastrophic only because the HST errors are 0.16–0.35 % in those bins. The fit errors also include propagated rotation uncertainty where applicable; the table
shows statistical errors alone. The large χ² establishes model/data tension. Missing
anisotropy and calibration or selection effects are possible causes, and mass estimates
can remain biased until those causes are tested.

## 4. Component ratios motivate anisotropy tests

The geometric mean of the two PM components is matched to 2–3 % at 150–300″, while their
ratio has a stronger radial residual. Matching that mean is not an independent measurement
of enclosed mass: projection, tracer density, anisotropy and streaming also enter. The isotropic model does not predict
exactly 1 for this ratio, because the tangential model has the published rotation subtracted;
the residual excess (data ratio over model ratio) is the part that needs new physics:

| R | 2.2 pc | 3.6 pc | 5.8 pc | 8.0 pc | 9.2 pc | 13.7 pc | 21.4 pc | 26.9 pc | 34.2 pc | 43.6 pc | 55.3 pc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| data σ_T/σ_R | 0.979 | 0.937 | 0.881 | 0.847 | 0.931 | 0.859 | 0.940 | 0.988 | 1.081 | 1.103 | 1.176 |
| model σ_T/σ_R | 0.995 | 0.985 | 0.958 | 0.920 | 0.901 | 0.850 | 0.907 | 0.939 | 0.967 | 0.989 | 0.998 |
| excess | −1.6 % | −4.8 % | −8.0 % | −8.0 % | +3.3 % | +1.1 % | +3.6 % | +5.3 % | **+11.8 %** | **+11.6 %** | **+17.9 %** |

(The first four columns are HST, the rest Gaia EDR3; 9.2 pc onwards is a different instrument,
so the step at 8 → 9 pc should not be over-interpreted.)

This pattern is consistent with radial anisotropy at intermediate radii and tangential
anisotropy farther out, as expected if tides preferentially remove radial orbits. The
projected ratios do not directly recover intrinsic β(r), and the radii are projected.
Within Gaia alone, the measured ratio rises from 0.86 near 11 pc to 1.18 near 55 pc,
crossing unity around 27 pc. That observation motivates a turnover model; it does not
uniquely determine one.

MUSE is also below the isotropic prediction at large projected radius (9.4% at 290″).
Its response to radial anisotropy makes this a useful independent-instrument check.
Different tracer populations and the approximate rotation treatment remain alternatives.

### The alternative reading

The tangential component is also where rotation enters, so an error in the rotation curve
produces the same signature. Quantitatively: reconciling the HST tangential data at 224″
requires an extra 5.8 km/s in quadrature on top of the 4.45 km/s the published curve already
supplies there, i.e. a ~64 % larger rotation amplitude at 5.8 pc; and reconciling the outer
Gaia points requires *less* rotation than the model subtracts at 34–55 pc. A single
normalisation error cannot do both — it would need the shape of the rotation curve to be wrong
in opposite directions at 6 and 40 pc. This makes a common amplitude error an incomplete explanation. Test the Gaia alternative
with `--gaia-rotation ours`; there is no equivalent HST switch. HST needs an external
absolute rotation estimate because local corrections removed its streaming. Reports now
replay saved Gaia options and check the stored likelihood. Anisotropy and rotation
errors can coexist.

## 5. Error-floor sensitivity at the existing best sample

χ² at the same best model with a systematic floor added in quadrature (not refitted):

| floor | 0 | 0.5 % | 1 % | 2 % |
|---|---|---|---|---|
| HST radial | 152 | 38 | 19 | 11 |
| HST tangential | 484 | 272 | 149 | 63 |
| MUSE | 76 | 72 | 62 | 42 |
| Gaia radial | 25 | 24 | 19 | 12 |
| Gaia tangential | 24 | 23 | 21 | 16 |
| **total (89)** | **761** | **428** | **271** | **144** |

A 0.5 % floor fixes the HST radial component almost completely (152 → 38): that misfit really
is error-bar-driven. The tangential component is not: it still contributes 63 of the 144 at a
2 % floor. The structured residual remains. This calculation does not refit the model and therefore
does not rule out every alternative error model. `--gaia-errors eta` changes the Gaia
measurements; it does not add an HST systematic floor or directly address this HST pattern.

## 6. What the "evidence for dark matter" is actually made of

Δχ² between K1 and K2-cored, by dataset:

| | HST radial | HST tangential | MUSE | Gaia radial | Gaia tangential | total |
|---|---|---|---|---|---|---|
| Δχ² (K2-cored − K1) | −8 | +7 | −1 | **+6** | **−14** | −10 |

The halo buys Δχ² ≈ 10 for two parameters, and the gain comes almost entirely from the outer
Gaia points: the tangential ones improve by 14 and the radial ones worsen by 6. Those are the
same 34–55 pc points that demand β < 0. An isotropic model can only raise both components
together, and raising them helps the tangential data (which sit above the model) more than it
hurts the radial data (which sit below). **The improvement is consistent with the halo compensating for missing outer tangential
anisotropy.** That mechanism needs a fit with the additional orbital freedom to establish
it. The evidence difference of +1.1 is insufficient for a DM claim.

## 7. Other parameters

**Distance.** All three runs give D = 5.311–5.313 ± 0.021 kpc against a prior of
5.43 ± 0.05 — the data pull it down by 2.4 prior σ and tighten it by a factor 2.4. This is a
kinematic distance from matching proper motions in mas/yr to line-of-sight velocities in km/s,
and it is driven mainly by the inner data where the fit is good. It is nevertheless a result
that should be re-checked at rung 1, because anisotropy changes the PM-to-LOS ratio and so can
bias D.

**Central mass.** The total M★ + M_rem + M_BH is 3.23 ×10⁶ (K1) and 3.12 ×10⁶ (K2-cored), i.e.
stable, but the split moves a lot: adding a halo moves 6 ×10⁵ M☉ from the stars to the
remnants. M_BH is 4.3–4.7 ×10⁴ M☉ with an 8 % error in all three runs; the inner fit is good
(χ²/point = 1.2), but this does not make the point-mass estimate robust. The central mass–anisotropy
degeneracy is untested. Moreover, an exactly cored tracer in a central point-mass potential
requires β(0) ≤ −1/2 for a non-negative DF; the isotropic moment baseline violates that
necessary condition in the central limit.

**Halo parameters.** r_s is unconstrained in both K2 runs (46–663 pc cored, 7–544 pc NFW),
which is expected: the data end at 55 pc and the halo is a small perturbation there.

## 8. Independent-engine check

Evaluating each best sample with JamPy instead of our Jeans solver:

| run | our solver | JamPy |
|---|---|---|
| K1 | ln L = −169.66 | −170.35 |
| K2-cored | −164.44 | −165.08 |

The total log likelihoods differ by about 0.7. Individual HST χ² terms shift by roughly
18 and largely cancel; the other terms shift by less than 1. The engines share the intended
spherical moment equations, but JamPy uses an MGE approximation to non-Gaussian mass
components and independent numerical integration. This audit has not isolated the cause
of the residual difference. Total-likelihood agreement alone is not a per-component
accuracy certificate.

## 9. Sampling cost

K1 took 4 h and 3.2 ×10⁶ likelihood calls; K2-cored 10 h / 6.6 ×10⁶; K2-NFW 20.6 h /
14.3 ×10⁶. The NFW run spent most of its time crawling along the M_DM–r_s degeneracy ridge
(one accepted draw per ~4000 calls for several hours, with 82 % of the evidence still in the
live points after 12 h). For rungs 1 and 2, which add one and three parameters, the NFW family
should test UltraNest's step sampler (`--step-sampler`) and check convergence. A prior
change requires scientific justification and changes the evidence; it is not a numerical
substitute for improving sampling.

## 10. Conclusions and next experiments

1. All three isotropic models have acceptable central residuals and a poor outer fit.
2. The strongest residual is in the relative PM components. It motivates radial variation
   in anisotropy, but does not independently certify the mass profile.
3. A model with intermediate radial and outer tangential anisotropy is a testable
   interpretation. Rotation, tracer selection and correlated errors remain relevant.
4. The small cored-halo preference may compensate for missing orbital freedom. It is
   not a DM detection; retain the table's limits only as conditional baseline outputs.

Follow [MODELLING_PLAN.md](MODELLING_PLAN.md): constant-β controls first, then a specified
turnover extension. The existing three-parameter β(r) is monotonic and is not the proposed
turnover model. Neither rung 1 nor the extension has yet superseded these runs with a
completed current-data posterior.

Gaia rotation/error variants, an independently defined HST rotation test, prior sensitivity
and injection recovery follow. The report replay and mock-generator repairs are complete;
use the source-run recipes in the modelling plan. Instrument scales remain a separate diagnostic, outside the
compared sequence. New fits require an explicit request.
