# Models

## Current: action-based DF models

| Model | What it adds to the base model | Purpose | Status |
|---|---|---|---|
| [Base no-DM model](base_no_dm.md) | — | Reference fit | Refitted 2026-10-13 with the corrected count likelihood |
| Free cored halo | gNFW halo with γ = 0, free \(\rho_{20}\) and \(r_s\) | Do the data want extended dark matter? | Free fit goes to \(\rho_{20}=0\) (2026-09-23, old counts) |
| \(\rho_{20}\) profile | \(\rho_{20}\) fixed at 0.25, 0.5, 1, 2, 4 M☉ pc⁻³, all else refitted | How much DM is allowed at fixed fit quality | Δ objective +0.5 … +10.6; a sensitivity curve, not a limit (mocks scatter by ±0.7) |
| Mass-per-star scan | HST count amplitude pinned at 8–12 M☉ per F625W < 19 star, free halo | How the star/dark split depends on the assumed stellar mass | \(\rho_{20}\) ≈ 2–4 at the photometric 8–9, 0 at the dynamical 10.7 |

Main results, figures and mock tests: [DF observed-fit write-up](../../df_observed_fit.pdf)
([LaTeX](../../df_observed_fit.tex)). The halo and scan numbers above use the
count products from before the 2026-10-13 fix and have not been rerun.

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
