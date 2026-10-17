# Models

## Current: action-based DF models

| Model | What it adds to the base model | Purpose | Status |
|---|---|---|---|
| [Base no-DM model](base_no_dm.md) | — | Reference DF fit until 2026-10-16; parent of the black-hole fits | Refitted 2026-10-13 with the corrected count likelihood; free count amplitudes, no completeness. Objective 280.7, χ²_kin 156.4 |
| [Tied, completeness-corrected DF refit](base_no_dm.md#refit-with-completeness-corrected-counts-and-tied-amplitudes-2026-10-16) | per-bin count completeness; HST/Gaia amplitude ratio tied to the measured 11.42 ± 0.23 | **Reference DF fit** (used for the dark-matter fits) | Fitted 2026-10-16. Objective 295.77, χ²_kin 159.4; best fit (ML) M_rem 2.00×10⁵ M☉ |
| [Jeans comparison model, J-Trager](jeans_comparison.md) | spherical Jeans, fixed MGE light tracer (Trager+1995), turnover anisotropy, free black hole | Earlier best model; tests a central black hole | χ²_kin 190.0 (ML); M_BH 4.25 (3.92–4.58)×10⁴ M☉ (posterior median, 16–84 %) |
| [Jeans twin, J-counts](jeans_comparison.md#twin-model-with-the-star-count-tracer-j-counts-2026-10-16) | as J-Trager, tracer an MGE of the completeness-corrected star counts (tied ratio) | Tests the tracer dependence | χ²_kin 152.2 (ML), below the DF's 156.4; M_BH 3.67 (3.30–4.02)×10⁴ M☉ (posterior median) |
| [DF + central black hole](df_bh.md) | exact point mass M_BH = 8,200 or 4×10⁴ M☉ (fixed) in the DF potential | Is a central black hole compatible with the DF fit? | Objective 275.2 / 279.7 (base 280.7); 4×10⁴ over-predicts central counts. Free amplitudes, no completeness: not refitted with the tied counts. [Cusp-depletion term](bh_cusp_depletion.md) **proposed, not adopted** |
| Free cored halo | gNFW halo with γ = 0, free \(\rho_{20}\) and \(r_s\) | Do the data want extended dark matter? | Free fit goes to \(\rho_{20}=0\): 2026-09-23 (old counts) and again with the tied counts (2026-10-16, `results/df/dftie_halo_20261016`; [dark-matter fits](dm_fits.md)) |
| \(\rho_{20}\) profile | \(\rho_{20}\) fixed at 0.25, 0.5, 1, 2, 4 M☉ pc⁻³, all else refitted | How much DM is allowed at fixed fit quality | Δ objective +0.5 … +10.6; a sensitivity curve, not a limit (mocks scatter by ±0.7) |
| Mass-per-star scan | HST count amplitude pinned at 8–12 M☉ per F625W < 19 star, free halo | How the star/dark split depends on the assumed stellar mass | \(\rho_{20}\) ≈ 2–4 at the photometric 8–9, 0 at the dynamical 10.7 |

Main results, figures and mock tests: [DF observed-fit write-up](../../df_observed_fit.pdf)
([LaTeX](../../df_observed_fit.tex)). The \(\rho_{20}\) profile and the
mass-per-star scan above use the count products from before the 2026-10-13 fix
and have not been rerun; the free cored halo has been rerun with the tied counts.

Which DF fit is the reference: the tied, completeness-corrected refit
(2026-10-16), on which the dark-matter fits are built; it changes the
predictions of the 2026-10-13 base fit by a few per cent but moves M_rem by
−27 %. The black-hole fits and most figures on the base-model page still use
the 2026-10-13 fit.

### Parameters of each model

The [base model](base_no_dm.md#all-parameters) has 10 free parameters (7
stellar, 2 remnant, the distance) and 3 profiled count nuisances. The
dark-matter variants add a halo; the [Jeans model](jeans_comparison.md#parameters)
has its own 10 parameters.

| Model | Parameters beyond the base model | Range | Fixed |
|---|---|---|---|
| Free cored halo | \(\rho_{20}\) (DM density at 20 pc), \(r_s\) (scale) | 0–10 M☉ pc⁻³; 5–150 pc, log | \(\gamma=0\) (core), taper \(r_t=1\) kpc |
| \(\rho_{20}\) profile | \(r_s\) | 5–150 pc, log | \(\rho_{20}\) = 0.25, 0.5, 1, 2, 4 |
| Mass-per-star scan | \(\rho_{20}\), \(r_s\) | as the free halo | HST count amplitude pinned: 8–12 M☉ per counted star |

Halo density: \(\rho(r)=\rho_{20}\,(r/20\,{\rm pc})^{-\gamma}\,[(1+r/r_s)/(1+20/r_s)]^{\gamma-3}\,e^{-(r^2-20^2)/r_t^2}\).

**What the data constrain.** Twelve mocks drawn from the base model recover
the total mass to ±3% at 20 pc and ±5% over 2–40 pc, but the recovered
\(\rho_{20}\) scatters by 0.5–0.9 M☉ pc⁻³. The leading degeneracy is stellar
mass against dark mass (remnants plus halo) through the count amplitude.
Pinning the mass per counted star removes it but needs that number to ~3%;
the photometric mass function gives 8.3–9.8 M☉ per star against the
dynamical 10.7.

**Implementation and history.**

- [Regularized DF implementation](../../REGULARIZED_DF_IMPLEMENTATION.md): the
  DF family, numerics and tests.
- [Anisotropy diagnostics](../../DF_ANISOTROPY_DIAGNOSTICS.md): why the model
  needs two anisotropy transitions (inner data want radial, outer data
  tangential anisotropy).
- [Observed-data flexibility](../../OBSERVED_DF_FLEXIBILITY.md): the
  one-transition DF fails (χ²/N ≈ 4.8).
- [Compact DF proposal (PDF)](../../compact_df_model.pdf),
  [compact DF recovery](../../COMPACT_DF_RECOVERY.md),
  [frequency-mapper diagnostic](../../COMPACT_DF_FREQUENCY_DIAGNOSTIC.md).
- Earlier mixture DF (superseded): [AGAMA DF models](../../AGAMA_DF_MODELS.md),
  [capacity challenge](../../DF_CAPACITY_CHALLENGE.md),
  [mass recovery](../../DF_MASS_RECOVERY.md).
- [DF model literature review](../../DF_MODEL_LITERATURE_REVIEW.md).
- Next steps under discussion: [Codex review, 2026-10-13](../../reviews/codex_mass_model_next_2026-10-13.md)
  (multi-mass stellar populations, remnant treatment).

**Dark matter (2026-10-16).** All three reference models were refitted with a dark-matter halo: J-Trager and J-counts each with an NFW and a cored halo, and the tied DF with a cored halo. None prefers one, as the [dark-matter fits](dm_fits.md) show.

## Earlier: spherical Jeans ladder (historical)

The first modelling sequence used the project's spherical anisotropic Jeans
solver, with JamPy as a cross-check, and nested sampling. It is superseded by
the DF models, which guarantee a positive distribution function.

| Rung | Anisotropy | χ² (89 bins): K1 no halo / K2 cored / K2 NFW |
|---|---|---|
| 0 | isotropic | 772 / 761 / 763 |
| 1 | constant β | 342.5 / 315.3 / 316.0 |
| 2 | varying β | see [model comparison](../../RUNG_MODEL_COMPARISON.md) |

- [Modelling plan](../../MODELLING_PLAN.md) and [mass-modelling procedure (PDF)](../../mass_modelling.pdf).
- [Rung 0](../../RUNG0_ANALYSIS.md), [rung 1](../../RUNG1_COMPARISON.md),
  [all rungs](../../RUNG_MODEL_COMPARISON.md).
- Dark-matter density at 20 pc from rung 2 (superseded: the DF mocks show
  \(\rho_{20}\) is not identified): [posterior](../../DM_DENSITY_POSTERIOR.md),
  [data constraints](../../DM_DENSITY_DATA_CONSTRAINTS.md),
  [profile limit](../../DM_DENSITY_PROFILE_LIMIT.md),
  [checks](../../DM_DENSITY_POSTERIOR_CHECKS.md).
- [DF consistency audit](../../DF_CONSISTENCY_AUDIT.md) of the Jeans fits;
  [literature baselines](../../LITERATURE_BASELINES.md);
  [code and analysis audit](../../CODE_ANALYSIS_AUDIT.md).
