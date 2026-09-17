# Literature baselines for the Milestone 3 dynamical fits

*Compiled 2026-09-17 from ADS. Every entry below has a resolved bibcode; the numbers are
taken from the papers' own abstracts (marked **A**) or from a catalogue product on disk
(marked **C**). Nothing here is from memory. Where a paper's number lives only in its body
text, it is not quoted.*

The purpose is a list of quantities our K1/K2 posteriors must be checked against, and
the published disagreements they must be able to express.

## Distance

| Value | Method | Source |
|---|---|---|
| **5494 ± 61 pc** | kinematic (σ_PM vs σ_LOS, r < 100″) | oMEGACat VI, Häberle+ 2025, [2025ApJ...983...95H](https://ui.adsabs.harvard.edu/abs/2025ApJ...983...95H) — read from the paper source |
| **5.43 ± 0.05 kpc** | combined Gaia EDR3 + HST + literature | Baumgardt & Vasiliev 2021, [2021MNRAS.505.5957B](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.5957B) — **C**, via `gc_catalog_updated.fits` |
| 4.59 ± 0.08 kpc | discrete axisymmetric Jeans models | Watkins+ 2013, [2013MNRAS.436.2598W](https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2598W) — **A** |

The 2013 value is 15% lower than the two recent ones; a fit with D free should recover
the recent ones. D is a parameter in K2.

## Luminous mass, M/L, anisotropy

| Quantity | Value | Source |
|---|---|---|
| Total mass | 3.94 × 10⁶ M☉ | Baumgardt database compilation — **C** (`baumgardt_ocen_parameters.mass`; version unrecorded, UNVERIFIED against upstream) |
| M/L_V | **2.71 ± 0.05** | Watkins+ 2013 — **A** |
| Anisotropy β | **0.10 ± 0.02** (mildly radial), inclination 50 ± 1° | Watkins+ 2013 — **A** |
| Core isotropy (22 GCs) | σ_t/σ_r = 0.992 ± 0.005 in cores; 0.8–1.0 near r_h | Watkins+ 2015, [2015ApJ...803...29W](https://ui.adsabs.harvard.edu/abs/2015ApJ...803...29W) — **A** |
| Our own data | σ_R/σ_T rises from ≈1 (core) to 1.18 at 300″ | oMEGACat VI products, `plots/omegacat_vi_profiles.png` |

## Energy equipartition (relevant to the HST–Gaia σ discrepancy)

| Quantity | Value | Source |
|---|---|---|
| η (σ ∝ m^−η), 1G stars | −0.007 ± 0.026 | Bellini+ 2018, [2018ApJ...853...86B](https://ui.adsabs.harvard.edu/abs/2018ApJ...853...86B) — **A** |
| η, 2G stars | 0.074 ± 0.029 | same |
| Theoretical ceiling | η_max ≈ 0.15 ± 0.03 (centre), η_∞ ≈ 0.08 ± 0.02 | Trenti & van der Marel 2013, [2013MNRAS.435.3272T](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.3272T) — **A** |

Implication: with η ≲ 0.1, a factor ~3 in tracer mass changes σ by ≲ 12%. The 20–25%
HST–EDR3 σ_PM difference we measured at 1–3′ is therefore **too large to be equipartition
alone**; crowding systematics in the Gaia profile inside ~3′ must carry part of it. This
sharpens the Milestone 3 decision to use the Gaia profiles only outside the overlap.

## The central dark mass: IMBH versus stellar-mass black holes

The long-standing disagreement our K1 model must be able to express.

| Claim | Value | Source |
|---|---|---|
| IMBH (isotropic, spherical, Gemini+HST) | 4.0 (+0.75/−1.0) × 10⁴ M☉; σ rises 18.6 → 23 km/s inside 14″ | Noyola+ 2008, [2008ApJ...676.1008N](https://ui.adsabs.harvard.edu/abs/2008ApJ...676.1008N) — **A** |
| IMBH (VLT) | (4.7 ± 1.0) × 10⁴ M☉; central σ 22.8 ± 1.2 km/s | Noyola+ 2010, [2010ApJ...719L..60N](https://ui.adsabs.harvard.edu/abs/2010ApJ...719L..60N) — **A** |
| Upper limit (HST PMs, anisotropic models) | M_BH ≲ 1.2 × 10⁴ (1σ), ≲ 1.8 × 10⁴ (3σ); core (γ=0) models fit with M_BH = 0; cusp models need (8.7 ± 2.9) × 10³ or a dark cluster ≲ 0.16 pc | van der Marel & Anderson 2010, [2010ApJ...710.1063V](https://ui.adsabs.harvard.edu/abs/2010ApJ...710.1063V) — **A** |
| IMBH from σ profile (N-body grid) | ~4 × 10⁴ M☉ | Baumgardt 2017, [2017MNRAS.464.2174B](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.2174B) — **A** |
| No IMBH; stellar-mass BHs instead | 4.6% of the cluster mass in a central BH cluster fits all data; a 4.5 × 10⁴ IMBH would produce ~20 fast stars not seen (at the time) | Baumgardt+ 2019, [2019MNRAS.488.5340B](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.5340B) — **A** |
| Stellar-mass BHs | ~5% of mass in BHs reproduces the central kinematics | Zocchi+ 2019, [2019MNRAS.482.4713Z](https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4713Z) — **A** |
| Multimass models | BH mass fraction > 5% for ω Cen (0–1% typical elsewhere) | Dickson+ 2024, [2024MNRAS.529..331D](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529..331D) — **A** |
| **IMBH from fast stars** | **M_BH ≳ 8200 M☉** (firm lower limit from velocities alone) | Häberle+ 2024, [2024Natur.631..285H](https://ui.adsabs.harvard.edu/abs/2024Natur.631..285H) — **A** |

Reading: the σ-profile-based IMBH masses (4–5 × 10⁴) are degenerate with a ~5% remnant
population; the fast-star lower limit (≳ 8200 M☉) is the one constraint that is not. For
K1 this means: **the remnant component and the IMBH are expected to be degenerate in the
second moments alone**, which is exactly the specification's point (section 1) — and the
K1 posterior should show it rather than pick one.

## The direct precursor: a non-luminous extended component

**Evans, Strigari & Zivick 2022**, *Dark and luminous mass components of Omega Centauri from
stellar kinematics*, [2022MNRAS.511.4251E](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4251E),
[arXiv:2109.10998](https://arxiv.org/abs/2109.10998) — **A**:

- Gaia EDR3 + HST proper motions + line-of-sight velocities; steady-state axisymmetric model;
  Gaussian and NFW profiles for a dark component distinct from the stars.
- **Dark mass ≲ 10⁶ M☉ inside the half-light radius.**
- For the HST + RV data, models with a non-luminous component are *strongly preferred* over
  a constant-M/L stellar-only model.
- A compact remnant core "may account for up to ~5 × 10⁵ M☉" but "likely cannot explain the
  higher end of the range" — leaving particle dark matter open.
- The inferred dark component is much more centrally concentrated than dSph haloes;
  J-factors ~10²²–10²⁴ GeV² cm⁻⁵.

This is the closest published analysis to Milestone 4 and the natural first comparison:
our K2 posterior on `M_DM(<r_h)` should be confronted with their ≲ 10⁶ M☉, and our
remnant/DM degeneracy with their remnant ceiling of ~5 × 10⁵ M☉. Their use of the same
Gaia + HST combination also means they faced the σ discrepancy we measured; how they
handled it is worth reading in the body text before we fix our own treatment.

## What is *not* in the literature (as far as this search found)

No published attempt to use ω Cen's **tidal tails** to constrain an extended dark
component — the specification's central idea. The tails literature (Ibata+ 2019;
Kuzma+ 2021, 2025, 2026; Youakim+ 2023) is morphological and spectroscopic, not a mass
constraint on the progenitor's outer potential.

## How these get used

- **K0 check:** our ingested profiles against Watkins+ 2015 (isotropic core) and the
  oMEGACat VI numbers — already visible in `plots/omegacat_vi_profiles.png`.
- **K1 targets:** M/L_V ≈ 2.7, β ≈ 0.1, remnant fraction ~5%, M_BH between ≳ 8.2 × 10³
  (fast stars) and ≲ 1.2–1.8 × 10⁴ (PM limits); D ≈ 5.43–5.49 kpc.
- **K2 comparison:** Evans+ 2022's ≲ 10⁶ M☉ dark mass within r_h.
- A K1 posterior that cannot accommodate the published disagreements is wrong, not decisive.
