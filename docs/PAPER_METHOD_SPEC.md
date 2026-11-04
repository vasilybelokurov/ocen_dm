# Paper method specification (draft v0, 2026-11-04)

Status: **draft for discussion**. Nothing here is adopted until approved. Items marked
**[DECIDE]** need a decision; the recommendation follows each.

## 1. Objectives and the numbers the paper reports

| # | Objective | Reported quantity |
|---|---|---|
| O1 | Dark matter | Posterior (or calibrated upper limit) on M_DM(<r) at r = 10, 30, 63 pc (63 pc = outermost kinematic bin); ρ_DM(20 pc) as secondary |
| O2 | Central black hole | Posterior on M_BH including M_BH = 0; upper limit if no detection |
| O3 | Dark remnants | Posterior on the extra dark central mass (M_rem, a_rem) and the enclosed dark mass M_dark(<r) = M_BH + M_rem(<r) at r = 0.1, 0.3, 1, 3 pc |

The joint M_BH–M_rem posterior is a deliverable in its own right: it shows how far an
extended remnant cluster can replace a point mass.

## 2. Primary model

**[DECIDE D1] Primary dynamical model.** Recommendation: the self-consistent regularized
action DF (`src/ocen_dm/kinematics/regularized_df.py`) as primary. Jeans (J-counts, and JAM
where relevant) as a cross-check on the same frozen data.

- Reason: the DF is positive by construction and self-consistent, so it predicts the counts
  and the kinematics from one object. Jeans can fit kinematics with unphysical anisotropy and
  density pairs.
- Known disagreement that the paper must explain: Jeans gives M_BH ≈ 3.6–4.3×10⁴ M☉
  (`docs/kb/models/dm_fits.md`); the DF prefers 8,200 over no BH by Δ(−2 ln L) ≈ −5.5 and
  finds 4×10⁴ neutral (`docs/CENTRE_MARGINALISATION.md`). The mock suite (§6) must show
  which method recovers an injected BH, and why they differ.

Free parameters (DF):

| Block | Parameters | Current treatment |
|---|---|---|
| Stars | M_star; DF shape J0, alpha, b_out, J_a, b_outer, J_outer | free, bounds as in `results/df/dfbh8200_restart_20261015/batch.json` |
| Remnants | M_rem, a_rem (Plummer, extra to the tracer DF) | free |
| BH | M_BH ≥ 0 | **fixed** at 0 / 8,200 / 4×10⁴ so far → becomes free |
| DM | ρ20, r_s (γ = 0 cored; γ = 1 NFW as a variant), taper 1 kpc | free in one DF run only (driven to 0) |
| Distance | D, prior N(5.43, 0.05) kpc | free |
| Centre | model centre c, prior N(0, σ_c = 2″) per axis about AvdM10; data frozen | marginalised (model side, `src/ocen_dm/centre_model.py`) |
| Rotation | VB21 rotation curve subtracted (HST tangential); v_rot²/2 (MUSE) | fixed, not fitted |

## 3. Data (frozen once)

One versioned snapshot with SHA-256 hashes recorded, as now (`data_snapshot`,
`counts_sha256` in each run):

- HST PMs (our oMEGACat-based estimator, radial and tangential, OmegaCAT systematics
  treatment as now);
- MUSE line-of-sight dispersions;
- Gaia EDR3 PMs (radial and tangential);
- HST and Gaia star counts: completeness-corrected, amplitudes tied, area-weighted bins.

**[DECIDE D2] Count range and outer Gaia counts.** Cristina's wavy residuals at 300–2,400″
are unresolved. Recommendation: keep the counts in, add a pre-declared variant without
Gaia counts beyond 300″, and report both if they differ.

## 4. Likelihood

As implemented (`DFJointProblem.evaluate`): split-normal kinematic bins, Poisson counts,
distance prior, and the centre logsumexp over Gauss–Hermite nodes (n = 10, σ_c = 2″).

**[DECIDE D3] Model adequacy.** The best fits have χ²_kin/89 ≈ 1.7 and fail the data-fit gate
(`batch.json` gates: 1.3 per point). Likelihood differences are then not calibrated by
Wilks' theorem. Options:

1. calibrate all Δ ln L thresholds and intervals on mocks (§6) that reproduce the
   misfit — no change to the errors;
2. add flexibility where the residuals sit (the DF outer envelope; a second stellar
   population), then re-check the gate;
3. a pre-declared per-dataset error scale, fitted as a parameter.

Recommendation: 1 plus 2. Option 3 changes the error model and needs your explicit
approval; I do not propose it by default.

**[DECIDE D4] Inner-bin estimator bias.** The centre mock found that our dispersion
estimator is biased by −0.2 to −0.5 σ in the two innermost HST bins for a steep profile.
Recommendation: measure the bias as a function of the true profile on mocks, and report the
BH result with and without the two innermost PM bins (already computed for fixed BH:
`dm2_drop2`). Do not apply a correction unless the bias is model-independent.

## 5. Inference

The cost of one DF evaluation (a model build plus the centre grid) is ≈ 1.5–3 min on one core
(83 s per build with a BH, `JOURNAL.md` 2026-10-15). Full sampling of about 15 parameters
needs ~10⁵ evaluations ≈ 3,000–5,000 core-hours. That is not practical.

**[DECIDE D5] Inference scheme.** Recommendation: profile likelihoods in the target
quantities, calibrated on mocks.

- **Grid in M_BH:** 0, 2×10³, 5×10³, 10⁴, 2×10⁴, 4×10⁴, 8×10⁴. All other parameters are
  optimised at each point (as for the existing fixed-BH fits: 25–60 min each, warm-started
  from the neighbouring grid point plus one independent start).
- **Grid in M_DM(<63 pc):** 0, 10⁴, 10⁵, 3×10⁵, 10⁶, 3×10⁶ M☉, same scheme.
- **M_rem–M_BH:** a 2-D grid.
- **Intervals:** from Δ(−2 ln L_profile) thresholds calibrated by the mock ensemble (the
  thresholds the mocks need for 68 / 95 % coverage), not from χ²₁ = 1 / 3.84.

Alternative: an emulator of ln L plus MCMC. It gives proper posteriors but adds a layer
that itself needs validation. Keep it as an option if the profiles turn out to be
multi-modal.

Jeans cross-check: nested sampling as now (400 live points), with the same data and the
same centre treatment.

## 6. Validation: tests and pass criteria

All mocks use the real star positions, real per-star errors, real bins and the real
estimator, so the selection and the noise are those of the data.

| ID | Test | Pass criterion |
|---|---|---|
| V0 | Unit tests: likelihood pieces, centre geometry, count quadrature (exist: 12 pass) | all pass |
| V1 | Noiseless self-mock (truth = a DF of our own family) | recovered parameters within 0.1 σ of the truth; Δ ln L at the truth < 0.1 |
| V2 | Noisy self-mocks: 4 truths × 10 noise realisations — (M_BH, M_DM) = (0, 0), (10⁴, 0), (4×10⁴, 0), (0, 10⁶ in 63 pc); remnants as the best fit | median bias < 0.3 σ; calibrated 68 / 95 % intervals cover within ±10 % (in the stated number of realisations) |
| V3 | Independent-family mocks (`src/ocen_dm/kinematics/df_mock.py` energy DFs + point mass; mass-segregated multi-mass model) | the BH and DM conclusions do not flip; biases reported |
| V4 | Data-like misfit: V2 truths with the real best-fit residual pattern added | the calibrated thresholds from V2 still cover |
| V5 | Jeans on the V2 mocks | quantifies the Jeans-vs-DF BH difference on known truths |
| V6 | Real-data robustness, pre-declared variants: σ_c = 1, 3″; drop the 2 innermost PM bins; no outer Gaia counts; NFW instead of cored | report every variant; none used to choose the headline |

Before any real-data headline: V0–V2 pass, and V3–V5 are reported. Figures with the
data, the model and the residuals accompany every fit.

## 7. Compute estimate (to be confirmed before launch)

- V2 profiles: 40 mocks × 7 BH points ≈ 280 fits × 0.5–1 h ≈ 140–280 core-hours, i.e. 14–28 h
  on 10 cores. A DM grid of the same size adds the same again.
- Real-data profiles: about 20 fits ≈ 10–20 core-hours.
- Each launch is announced with its threads and wall time first.

## 8. Build list (my working time)

1. A free-BH and free-DM-grid driver with the centre-marginalised objective, built on the
   existing batch controller (~2 h).
2. A mock generator for V1/V2 from a DF truth through the real estimator (adapted from
   `bin/centre_model_mock.py`) (~2 h).
3. Profile assembly, threshold calibration and the standard figures (~2 h).
4. V3 independent-family mocks with a point mass (~2–3 h; the BH extension of `df_mock.py`
   is new).
5. Single method write-up (`docs/METHOD.md`) replacing the scattered notes; older docs stay
   as history.

## 9. Open decisions (summary)

| ID | Question | Recommendation |
|---|---|---|
| D1 | Primary model | DF primary; Jeans as cross-check |
| D2 | Outer Gaia counts | keep, plus a variant without them beyond 300″ |
| D3 | χ²/pt ≈ 1.7 | mock-calibrated thresholds plus targeted DF flexibility; no error scaling without approval |
| D4 | Inner-bin estimator bias | measure on mocks; report with and without the 2 innermost bins |
| D5 | Inference | mock-calibrated profile likelihoods; emulator plus MCMC as a fallback |
| D6 | Remnant prior | flat in M_rem, a_rem (as now), with the IMF budget (`docs/kb/` remnant page) used only for interpretation; or adopt the proposed prior |
