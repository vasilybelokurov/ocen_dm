# Results

Current findings by strand. Each links the report it rests on.

## Base model fits

The reference result of the project: the [base no-DM model](models/base_no_dm.md)
(two-transition action DF for the stars, a Plummer remnant component, no dark
matter, no central black hole) fitted jointly to 89 kinematic bins and the HST
and Gaia star counts.

Since 2026-10-16 the reference DF fit is the
[tied, completeness-corrected refit](models/base_no_dm.md#refit-with-completeness-corrected-counts-and-tied-amplitudes-2026-10-16)
(`results/df/dftie_restart_20261016`; objective 295.77, χ²_kin 159.4; best fit
M_★ 3.08×10⁶ M☉, M_rem 2.00×10⁵ M☉, a_rem 1.88 pc). Its predictions differ from
the 2026-10-13 fit below by a few per cent, but M_rem falls by 27 %: the split
between stars and remnants is weakly constrained. The black-hole fits
([DF + central black hole](models/df_bh.md)) and the figures below use the
2026-10-13 fit; the [dark-matter fits](models/dm_fits.md) use the tied refit.

### Current fit (2026-10-13)

Batch `results/df/dftwo_counts_area_20261013`, refitted after the count-likelihood
correction. All four valid starts reach the same minimum (objectives within
0.01); numerical validation passes; no parameter sits at a bound.

| Quantity | 2026-10-13 (current) | 2026-09-23 (superseded) |
|---|---|---|
| \(\chi^2_{\rm kin}\) (89 bins) | 156.4 (1.76 per bin) | 158.1 (1.78) |
| HST count deviance (20 bins) | 96.9 | 98.9 |
| Gaia count deviance (12 bins) | 25.2 | 26.3 |
| \(M_\star\) (best fit) | 3.00×10⁶ M☉ | 3.00×10⁶ M☉ |
| \(M_{\rm rem}\), \(a_{\rm rem}\) | 2.73×10⁵ M☉, 2.18 pc | 2.62×10⁵ M☉, 2.15 pc |
| Distance | 5.356 kpc | 5.356 kpc |
| Mass per F625W < 19 star | 10.48 M☉ | 10.67 M☉ |
| Data-fit gate | fails | fails |

![Data versus the base model, both fits](../../plots/base_area_20261013_data.png)

*Star counts, MUSE line-of-sight and HST/Gaia proper-motion dispersions
against the current (blue) and superseded (orange) fits; dashed: the earlier
Jeans model. Bottom row: residuals (signed Poisson deviance residuals for the
counts). Both fits are evaluated with the corrected count likelihood.*

The two count profiles do not join, by design: they are different selections in their own units. HST counts are oMEGACat stars with F625W < 19 at 2–250″; Gaia counts are stars with G < 17 at 300–2400″. The brighter Gaia limit counts fewer stars per unit stellar mass (fitted amplitude ratio 13.7 with free amplitudes and raw counts; measured 11.42 ± 0.23 for complete samples, i.e. 11.95 ± 0.3 expected for the raw counts; [count-ratio check](data/cluster.md#ratio-of-the-two-count-products-in-the-overlap-2026-10-16)), so each profile has its own fitted count amplitude (a<sub>HST</sub>, a<sub>Gaia</sub>), and the Gaia model also includes a uniform field of about 1.6 stars arcmin⁻². The model curves are the same stellar surface density scaled by each amplitude. The gap between 250″ and 300″ is a gap in coverage.

### Comparison with the Jeans model (dashed curves)

Full description, figures and posterior: [Jeans comparison model](models/jeans_comparison.md).

The dashed curves are the best rung-2 spherical Jeans fit without dark matter
(`results/fits/rung2_K1_turnover`, J-Trager), the best no-DM model of the earlier
[Jeans ladder](models/index.md#earlier-spherical-jeans-ladder-historical).
Both models are fitted to the same 89-bin kinematic snapshot, so their
kinematic χ² compare directly; the Jeans fit uses a fixed light profile (an MGE
of the Trager+1995 profile joined to HST star counts) instead of the Poisson
star counts.

| | Jeans (rung 2, no DM) | DF base model |
|---|---|---|
| Stellar tracer and mass | MGE light profile, mass-to-light free | self-consistent action DF |
| Free stellar parameters | \(M_\star\) only (density shape fixed by the light) | \(M_\star\) and six DF parameters (shape fitted to the counts) |
| Anisotropy | turnover β(r): β₀ = −0.50 (centre) → +0.20 → −0.98 (outside) | isotropic core → +0.15 near 8 pc → −0.53 at 100 pc |
| Central black hole | 4.41×10⁴ M☉ (ML); 4.25 (3.92–4.58)×10⁴ (median, 16–84%) | none |
| \(M_\star\) | 2.37×10⁶ M☉ (ML) | 3.00×10⁶ M☉ (best fit) |
| Remnants (Plummer) | 7.7×10⁵ M☉, a = 5.7 pc (ML) | 2.7×10⁵ M☉, a = 2.2 pc (best fit) |
| Total mass, M<sub>★</sub> + M<sub>rem</sub> + M<sub>BH</sub> | 3.18×10⁶ M☉ (ML); 3.24×10⁶ (median) | 3.27×10⁶ M☉ (best fit) |
| \(M_\star/L_V\), \(M_{\rm tot}/L_V\) | 2.0, 2.7 (ML); 2.3, 2.8 (median) | 2.6, 2.8 |
| Distance | 5.34 kpc | 5.36 kpc |

| χ² (per bin) | Jeans | DF |
|---|---|---|
| HST PM radial (21) | 61.8 (2.94) | 42.6 (2.03) |
| HST PM tangential (21) | 36.6 (1.74) | 39.8 (1.89) |
| MUSE LOS (29) | 51.7 (1.78) | 38.1 (1.31) |
| Gaia PM radial (9) | 16.9 (1.88) | 14.4 (1.60) |
| Gaia PM tangential (9) | 23.0 (2.56) | 21.5 (2.39) |
| **Total (89)** | **190.0** | **156.4** |

- **The DF fits better than J-Trager overall** (Δχ² = 34) and in every dataset
  except the HST tangential proper motions. This holds for the light-profile
  tracer only: the Jeans twin with the star-count tracer
  ([J-counts](models/jeans_comparison.md#twin-model-with-the-star-count-tracer-j-counts-2026-10-16))
  reaches χ² = 152.2 (ML), below the DF's 156.4, with M_BH = 3.67
  (3.30–4.02)×10⁴ M☉ (posterior median) and M_rem = 8.4×10⁵ M☉.
- **The Jeans model follows the central rise** of the HST proper-motion
  dispersions inside 0.3 pc (~11″), which the DF, with an isotropic core and
  no black hole, stays below by 1.5–2.4σ in several bins: χ² = 7.3 (Jeans)
  against 21.5 (DF) for those 10 bins. This is where the black hole matters.
- **It pays for this at intermediate radii:** at 2–20 pc (36 bins of HST,
  MUSE and Gaia) χ² = 122.2 (Jeans) against 77.8 (DF), mostly in the HST
  radial proper motions (40.2 against 20.2) and the MUSE line-of-sight
  dispersions (30.9 against 18.7).
- **The black hole is tied to an imposed anisotropy.** A cored stellar density
  around a point mass needs β₀ ≤ −½ for any non-negative DF; the prior enforces
  it and the fit sits at the bound (β₀ = −0.503). Even so, under the standard
  separable construction a non-negative DF fails to exist at 0.02–0.18 pc,
  where the black hole dominates (89.7% of posterior draws fail; no failure
  outside ~0.2 pc; a non-separable DF is not excluded;
  [DF consistency audit](../DF_CONSISTENCY_AUDIT.md)). The DF base model is
  non-negative by construction.
- **Both put about the same total mass** in the cluster (3.2–3.3×10⁶ M☉) but
  divide it differently: the Jeans fit has 20% less stellar mass and nearly
  three times more remnant mass, spread over a 2.5× larger scale. This is the
  same stars-versus-dark-mass degeneracy seen within the DF fits.
- **Why \(M_\star\) differs by 11% at the same total mass.** In neither model
  is the stellar mass fixed by the light or the counts; it is the
  normalisation the kinematics leave after the remnants (and, for Jeans, the
  black hole and the imposed central anisotropy) take their share. The
  remnant mass is close to the total minus \(M_\star\), so with remnants at
  8–17% of the total, a 1% change in \(M_\star\) moves \(M_{\rm rem}\) by ~10%.
- **A central black hole in the DF** ([DF + central black hole](models/df_bh.md),
  2026-10-15, fixed masses, 2026-10-13 count likelihood): 8,200 M☉ lowers the
  objective by 5.5; 4×10⁴ M☉ lowers χ²_kin by 12 but loses 15 in the HST
  counts, because its cusp over-predicts the central counts. The differences
  are not calibrated with mocks. Whether the inner rise needs a black hole,
  radially biased central orbits or mass segregation is not settled.

### Mass components

![Mass components of the base model](../../plots/base_model_20261013_components.png)

*Current fit, by component. (a) Density; (b) enclosed mass with the
half-mass radii dotted; (c) circular velocity; (d) remnant share, local
(\(\rho_{\rm rem}/\rho_{\rm tot}\)) and cumulative (\(M_{\rm rem}/M_{\rm tot}\));
(e) dynamical-to-stellar mass ratio. Grey: outside the radial range of the
kinematic data (0.1–63 pc).*

| r (pc) | Stars (10⁶ M☉) | Remnants (10⁶ M☉) | Total (10⁶ M☉) | Remnant share of M(<r) |
|---|---|---|---|---|
| 1 | 0.012 | 0.020 | 0.032 | 62% |
| 3 | 0.201 | 0.145 | 0.346 | 42% |
| 10 | 1.436 | 0.255 | 1.691 | 15% |
| 30 | 2.716 | 0.271 | 2.987 | 9% |
| 100 | 2.987 | 0.273 | 3.260 | 8% |

- **Two components with different scales.** The stars have a 3D half-mass
  radius of 10.4 pc (projected 7.7 pc; the Baumgardt database lists
  r_h = 7.56 pc at 5.43 kpc, i.e. 7.46 pc at the fitted distance); the remnant
  Plummer has 2.8 pc.
- **The remnants dominate the core.** They are 66% of the central density
  (6.0×10³ against 3.2×10³ M☉ pc⁻³ for the stars at 0.3 pc). The stars take
  over the local density near 1.5 pc (remnants 59% at 1 pc, 26% at 3 pc) and
  the enclosed mass near 2 pc.
- **Mass-to-light-following mass.** \(M_{\rm tot}(<r)/M_\star(<r)\) falls from
  2.9 in the centre to 1.09 overall. This is the quantity a mass-segregated
  stellar population could also produce.
- **Circular velocity** peaks at 27.0 km/s at 9.0 pc; the stars alone peak
  at ~25 km/s near 10 pc, the remnants at ~14 km/s near 3 pc. The central
  escape speed is 68 km/s.
- The total inside 10 pc differs from the superseded fit by −0.7%.
- **The remnant component matches published fits.** Its best-fit mass and scale (2.73×10⁵ M☉, 2.18 pc; the tied refit gives 2.00×10⁵ M☉, 1.88 pc) agree with the extended central mass of Bañares-Hernández+2025 (2.65×10⁵ M☉ in a 1.88 pc Plummer); black-hole-only estimates are 1.35–1.82×10⁵ M☉. See [Remnant mass and spatial scale](literature_comparison.md#remnant-mass-and-spatial-scale).

### Internal kinematics

![Internal kinematics of the base model](../../plots/base_model_20261013_kinematics.png)

*(a) Intrinsic radial and tangential (one component) dispersions, with the
local escape speed/√12 for scale; (b) anisotropy β(r); (c) stellar surface
density with the projected and 3D half-mass radii; (d) projected
line-of-sight and proper-motion dispersions from the DF, before the
published rotation is subtracted in the likelihood.*

- **Isotropic core, radial bias at intermediate radii, tangential outskirts.**
  β = 0.00 at 0.3 pc, 0.09 at 3 pc, peaks at +0.15 near 8 pc, crosses zero
  near 22 pc, and reaches −0.08 at 30 pc and −0.53 at 100 pc.
- **Central dispersion** σ_r ≈ σ_t ≈ 21.4 km/s; σ_r = 14.7 and σ_t = 13.5 km/s
  at 10 pc.
- **Projected dispersions** in the core: 19.8 km/s along the line of sight
  and 19.3 km/s in each proper-motion component. The radial and tangential
  proper-motion profiles differ by ≲ 0.5 km/s inside 30 pc and separate in the
  tangentially biased outskirts.

### Distribution function

![Structure of the base-model DF](../../plots/df_structure_dftwo_counts_area_20261013_observed_no_halo_start3.png)

*(a) Intrinsic moments and β(r); (b) the anisotropy bias \(b(q)\),
\(q=J_r+L\), over the action distributions of the stars in the HST/MUSE
(R ≤ 9.5 pc) and Gaia (9.4–63 pc) ranges, with the two transitions
\(J_a\) and \(J_{\rm outer}\) dotted; (c) \(f(J_r, L)\).*

- The radial bias (b ≈ +0.45 at its peak) sits on the actions of the
  HST/MUSE stars; the tangential tail (b → −0.33) sits on the Gaia stars.
  The data constrain the DF over \(q\) ≈ 10–1000 pc km/s; outside that range
  its shape is extrapolated.
- \(f\) is smooth and monotonic in both actions, with no features at the
  transitions.

### What the base fits establish

- **No dark matter is needed** to fit the kinematics and counts. A free cored
  halo goes to \(\rho_{20}=0\), with the old count likelihood (2026-09-23) and
  again with the tied, completeness-corrected counts (2026-10-16); the Jeans
  fits show no preference for a halo either ([dark-matter fits](models/dm_fits.md)).
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

### Comparison with published models

Total mass, M/L, central dark mass and anisotropy against the literature (van de Ven+2006 to Bañares-Hernández+2025): [Comparison with published models](literature_comparison.md). The remnant masses and scales are compared in [Remnant mass and spatial scale](literature_comparison.md#remnant-mass-and-spatial-scale).

### Not yet rerun with the corrected count likelihood

The \(\rho_{20}\) profile, the mass-per-star scan and the mocks (all
2026-09-23). The free-halo fit has been rerun with the tied counts
([dark-matter fits](models/dm_fits.md)). See [Models](models/index.md) and the
[DF observed-fit write-up (PDF)](../df_observed_fit.pdf);
[data analysis (PDF)](../data_analysis.pdf).

## Dark-matter survival and tidal tails

- Restricted and live N-body runs of a no-DM and a DM model on Omega Cen's
  orbit (static McMillan17, 1.96 Gyr) give tails whose kinematics agree to
  10%. The DM model loses 1.3–1.6× fewer stars. Neither difference is
  testable with the published tail data (Kuzma+2021,
  [2021MNRAS.507.1127K](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1127K/abstract); Kuzma+2026,
  [2026MNRAS.550g1147K](https://ui.adsabs.harvard.edu/abs/2026MNRAS.550g1147K/abstract); Ibata+2019,
  [2019NatAs...3..667I](https://ui.adsabs.harvard.edu/abs/2019NatAs...3..667I/abstract)).
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

## Stream 54 (Fimbulthul) and the Milky Way bar

Since 2026-10-25 the stream work has its own project and site, [oCen_stream](https://vasilybelokurov.github.io/ocen_stream/) (code: https://github.com/vasilybelokurov/ocen_stream). Its status, reports and reviews moved there.
