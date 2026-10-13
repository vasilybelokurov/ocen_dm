# Results

Current findings by strand. Each links the report it rests on; the
[journal](../journal/index.md) has the full record, and later entries can
supersede earlier ones.

## Cluster mass model

- The two-transition DF with Poisson star counts and a distance prior fits
  all kinematic data without dark matter: χ²_kin = 156 for 89 bins. It does
  not pass the project's own data-fit gate (χ²_kin per bin 1.76 against ≤ 1.3;
  HST count deviance 4.8 per bin). See the [base model](models/base_no_dm.md).
- The data identify the total mass profile, not how it splits between stars,
  remnants and dark matter. On realistic mocks \(\rho_{20}\) scatters by
  ±0.7 M☉ pc⁻³.
- With the stellar mass per counted star fixed from the oMEGACat mass
  function, \(\rho_{20}\) ≈ 2–4 at the photometric 8–9 M☉ per F625W < 19 star
  and 0 at the dynamically preferred 10.7. The photometric uncertainty
  (≳ 1 M☉ per star) is too large to decide between them.
- 2026-10-13: a bias in the count likelihood (star-weighted rather than
  area-weighted bin averages) was fixed. The refitted base model moves by
  ≤ 4% in every parameter and by −1.8% in mass per counted star (10.67 →
  10.48 M☉). The halo and mass-per-star scans have not been rerun.

Reports: [DF observed-fit write-up (PDF)](../df_observed_fit.pdf);
[data analysis (PDF)](../data_analysis.pdf).

## Dark-matter survival and tidal tails

- Restricted and live N-body runs of a no-DM and a DM model on Omega Cen's
  orbit (static McMillan17, 1.96 Gyr) give tails whose kinematics agree to
  10%. The DM model loses 1.3–1.6× fewer stars. Neither difference is
  testable with the published tail data (Kuzma+2021, Kuzma+2026, Ibata+2019).
- Fitted halos are compact (r_s at the 5 pc bound) and traded against stellar
  mass, so across the kinematically allowed range the Jacobi radius at
  pericentre changes by ≤ 6.5%. Dark matter that could affect the tails lies
  beyond ~80–200 pc and is not constrained by the kinematics.

Reports: [dynamical experiments](../DYNAMICAL_EXPERIMENTS.md),
[first batch (PDF)](../dynamical_batch.pdf),
[batch analysis](../DYNAMICAL_BATCH_ANALYSIS.md),
[progenitor orbits](../PROGENITOR_ORBITS.md),
[Section-11 tests](../SECTION11_TESTS.md),
[lifetime evolution](../LIFETIME_EVOLUTION.md) (paused),
[DM capture constraints (PDF)](../dm_capture_constraints.pdf),
[N-body code assessment](../NBODY_CODE_ASSESSMENT.md),
[N-body ICs](../NBODY_IC_PROPOSAL.md).

## Stream 54 (Fimbulthul) and the proper-motion knee

- Beyond φ1 ≈ 17° the stream bends in Gaia DR3 proper motion. No published
  model reproduces this, nor do sprays in static McMillan17 or DB98 hosts.
- A rotating Hunter+2024 bar bends the arm part-way and steepens the pmdec
  gradient. Bar angle 16° beats 28° by ~285 in ln L on held-out seeds (paired
  seed sd 101); the distal PM residual remains.
- The LMC (rigid MW + LMC, reflex frame, benchmarked against the Tango
  N-body runs) shifts distal pmra by 0.2–0.6 mas/yr, but a 0.03–0.05 mas/yr
  change in Omega Cen's PM absorbs it.
- Spiral forcing, the conditional likelihood's density dependence, the
  Chen+2024 escape law and both launch arms do not remove the mismatch.
- Backtracked in the barred host, the 13 knee RV stars connect to the cluster
  much worse than model controls (1.9–2.0 against ~1.1). This is a mismatch in
  the recent (≤ 700 Myr) fixed-potential models tested, not yet an exclusion.

Reports: [PM refinement](../STREAM_PM_REFINEMENT_2026-10-13.md) (current),
[knee investigation](../STREAM_KNEE_INVESTIGATION.md),
[diagnostic plan](../STREAM_DIAGNOSTIC_PLAN.md),
[fit campaign plan](../STREAM_FIT_CAMPAIGN_PLAN.md),
[experiments plan](../STREAM_EXPERIMENTS_PLAN.md),
[stream-fitting methods review](../STREAM_FIT_METHODS_REVIEW.md),
[stream distances literature](../literature/omega_cen_stream_distances_2026-10-10.md).
