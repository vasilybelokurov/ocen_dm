# Results

Current findings by strand. Each links the report it rests on.

## Base model fits

The reference result of the project: the [base no-DM model](models/base_no_dm.md)
(two-transition action DF for the stars, a Plummer remnant component, no dark
matter, no central black hole) fitted jointly to 89 kinematic bins and the HST
and Gaia star counts.

### Current fit (2026-10-13)

Batch `results/df/dftwo_counts_area_20261013`, refitted after the count-likelihood
correction. All four valid starts reach the same minimum (objectives within
0.01); numerical validation passes; no parameter sits at a bound.

| Quantity | 2026-10-13 (current) | 2026-09-23 (superseded) |
|---|---|---|
| \(\chi^2_{\rm kin}\) (89 bins) | 156.4 (1.76 per bin) | 158.1 (1.78) |
| HST count deviance (20 bins) | 96.9 | 98.9 |
| Gaia count deviance (12 bins) | 25.2 | 26.3 |
| \(M_\star\) | 3.00×10⁶ M☉ | 3.00×10⁶ M☉ |
| \(M_{\rm rem}\), \(a_{\rm rem}\) | 2.73×10⁵ M☉, 2.18 pc | 2.62×10⁵ M☉, 2.15 pc |
| Distance | 5.356 kpc | 5.356 kpc |
| Mass per F625W < 19 star | 10.48 M☉ | 10.67 M☉ |
| Data-fit gate | fails | fails |

![Data versus the base model, both fits](../../plots/base_area_20261013_data.png)

*Star counts, MUSE line-of-sight and HST/Gaia proper-motion dispersions
against the current (blue) and superseded (orange) fits; dashed: the earlier
Jeans model. Bottom row: residuals (signed Poisson deviance residuals for the
counts). Both fits are evaluated with the corrected count likelihood.*

**Enclosed mass of the current fit**

| r (pc) | Stars (10⁶ M☉) | Remnants (10⁶ M☉) | Total (10⁶ M☉) | Remnant fraction |
|---|---|---|---|---|
| 1 | 0.012 | 0.020 | 0.032 | 62% |
| 3 | 0.201 | 0.145 | 0.346 | 42% |
| 10 | 1.436 | 0.255 | 1.691 | 15% |
| 30 | 2.716 | 0.271 | 2.987 | 9% |
| 100 | 2.987 | 0.273 | 3.260 | 8% |

The total inside 10 pc differs from the superseded fit by −0.7%.

### What the base fits establish

- **No dark matter is needed** to fit the kinematics and counts. A free cored
  halo goes to \(\rho_{20}=0\) (2026-09-23, old count likelihood).
- **The total mass profile is constrained; its split is not.** Mocks drawn
  from the base model recover the total mass to ±3% at 20 pc and ±5% over
  2–40 pc, but stars trade against remnants plus halo through the count
  amplitude, and the recovered \(\rho_{20}\) scatters by ±0.7 M☉ pc⁻³.
- **The answer depends on the stellar mass per counted star.** Fixed at the
  photometric 8–9 M☉ per F625W < 19 star, the data prefer \(\rho_{20}\) ≈ 2–4;
  at the dynamical value (now 10.5) they prefer none. The photometric
  estimate (8.3–9.8 in the zones that dominate the kinematics) is not precise
  enough (~1 M☉ per star; ~0.3 is needed) to decide.
- **The fit is not statistically acceptable.** χ²_kin per bin is 1.76 against
  the gate of 1.3, and the HST count residuals are structured: data below the
  model near 1 pc and in the outermost HST bin (~6 pc), above it at 3–4 pc.
  Either the error bars are underestimated or the model lacks freedom
  (mass segregation, rotation, a different remnant shape); this is not yet
  separated.

### Not yet rerun with the corrected count likelihood

The free-halo fit, the \(\rho_{20}\) profile, the mass-per-star scan and the
mocks (all 2026-09-23). See [Models](models/index.md) and the
[DF observed-fit write-up (PDF)](../df_observed_fit.pdf);
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
