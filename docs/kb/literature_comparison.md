# Comparison with published mass models

How our two no-DM models (the [DF base model](models/base_no_dm.md) and the
[Jeans comparison model](models/jeans_comparison.md)) compare with published
models of Omega Centauri's mass distribution. Published numbers are taken
from the papers' abstracts (checked in ADS, 2026-10-14) unless marked
otherwise; all papers are in the [SciX library](literature.md).

!!! abstract "Summary"
    - **Total mass and M/L agree within the spread of published values.** Our
      total M/L_V of 2.8 is close to Watkins+2013 (2.71); our total mass,
      3.26×10⁶ M☉, is 13–22% below the recent estimates of 3.7–4.2×10⁶ M☉
      (D'Souza & Rix 2013 scaled to our distance; the Baumgardt database).
    - **The central dark mass is where the literature splits, and our two
      models fall on the two sides.** The Jeans model's 4.3×10⁴ M☉ black hole
      reproduces the classic IMBH claims (Noyola+2008/2010, Baumgardt 2017,
      Jalali+2012). The DF model's compact remnant component (2.7×10⁵ M☉,
      2.2 pc scale, 8% of the mass, no black hole) matches the "extended
      central mass" picture: 2–3×10⁵ M☉ favoured by kinematics plus pulsar
      timing (Bañares-Hernández+2025: 2.65×10⁵ M☉ in a Plummer of 1.88 pc,
      the same form, mass and nearly the same scale as ours) and ~5% of the
      mass in stellar-mass black holes (Zocchi+2019, Baumgardt+2019,
      Dickson+2024).
    - **Anisotropy:** the DF profile (isotropic core, mildly radial at
      intermediate radii, tangential outside) agrees with oMEGACat VI and
      van de Ven+2006; the Jeans model's tangential core does not.
    - **The fast-moving stars require a point mass of at least 8,200 M☉**
      (Häberle+2024, robust bound). Their model-dependent estimates go
      higher (≥ 21,100 M☉ with accelerations; 39,000–47,000 M☉ preferred by
      an N-body comparison, which brackets our Jeans black hole). Whether
      the DF model with such a point mass still fits is not yet tested.

## Total mass, mass-to-light ratio and distance

Proper-motion-based masses scale as \(M\propto D^3\) and M/L as \(\propto D\);
published values are given at their own distance, and scaled to ours
(5.356 kpc) only where the paper states the scaling.

| Study | Method | D [kpc] | Total mass [10⁶ M☉] | M/L_V |
|---|---|---|---|---|
| **This work, DF base model** | spherical action DF, HST+MUSE+Gaia kinematics, star counts | 5.356 (prior 5.43 ± 0.05) | **3.26** | **2.8** (stellar 2.6) |
| **This work, Jeans model** | spherical Jeans, black hole, same kinematics | 5.34–5.36 | **3.18–3.24** | **2.7–2.8** (stellar 2.0–2.3) |
| [van de Ven+2006](https://ui.adsabs.harvard.edu/abs/2006A%26A...445..513V) | axisymmetric Schwarzschild, ground-based PMs + LOS | 4.8 ± 0.3 (fitted) | 2.5 ± 0.3 | 2.5 ± 0.1 |
| [Watkins+2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2598W) | discrete axisymmetric Jeans | 4.59 ± 0.08 (fitted) | — | 2.71 ± 0.05 |
| [D'Souza & Rix 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.429.1887D) | discrete PM likelihood | scaled to 5.356 | 4.20 (spherical isotropic), 3.74 (flattened) | — |
| Baumgardt GC database | N-body grid (catalogue value; version unrecorded, UNVERIFIED against the live database) | 5.43 | 3.94 | — |
| [Merritt+1997](https://ui.adsabs.harvard.edu/abs/1997AJ....114.1074M) | non-parametric, 469 LOS velocities | — | — | mass "slightly more extended than" the light |

D'Souza & Rix give \(4.55\times10^6\,(D/5.5\,{\rm kpc})^3\) M☉ (spherical) and
\(4.05\times10^6\,(D/5.5\,{\rm kpc})^3\) M☉ (flattened); at 5.356 kpc the factor
is 0.923.

- **M/L is close to published values:** our total M/L_V of 2.8 is slightly
  above Watkins (2.71) and van de Ven (2.5). Scaled to our distance with
  M/L ∝ D their values would be 3.2 and 2.8; their distances were fitted
  jointly, so this scaling is only indicative.
- **Our total mass is 13–22% below the recent estimates** at our distance
  (3.26 against 3.74–4.20×10⁶ M☉ from D'Souza & Rix and 3.94×10⁶ M☉ from
  the Baumgardt database). van de Ven+2006 found 2.5×10⁶ M☉ at 4.8 kpc. Possible reasons, not yet tested: our model
  is spherical (D'Souza & Rix find flattening lowers the mass by 10%, in the
  same direction, so it does not explain it); our kinematics reach ~60 pc and
  the mass beyond is extrapolated; the Gaia data may have systematics in the
  outer bins; and the Baumgardt value comes from mass-segregated N-body
  models.
- **Distance:** ours (5.356 kpc) is pulled by the prior; the oMEGACat VI
  kinematic distance is 5.494 ± 0.061 kpc
  ([Häberle+2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...983...95H)),
  2.2σ larger. At 5.494 kpc our proper-motion masses would be ~8% higher.

## The central dark mass

| Study | Method | Central dark mass |
|---|---|---|
| **This work, Jeans model** | spherical Jeans, β₀ ≤ −½ imposed | **IMBH 4.3 (3.6–4.9) ×10⁴ M☉** (95%) |
| **This work, DF base model** | self-consistent DF, isotropic core | **no IMBH; remnants 2.7×10⁵ M☉ in a 2.2 pc Plummer (8% of the mass)** |
| [Noyola+2008](https://ui.adsabs.harvard.edu/abs/2008ApJ...676.1008N) | isotropic spherical, Gemini + HST | IMBH 4.0 (+0.75/−1.0) ×10⁴ M☉ |
| [Noyola+2010](https://ui.adsabs.harvard.edu/abs/2010ApJ...719L..60N) | VLT kinematics | IMBH (4.7 ± 1.0) ×10⁴ M☉ |
| [Jalali+2012](https://ui.adsabs.harvard.edu/abs/2012A%26A...538A..19J) | N-body, isotropic multimass King ICs | IMBH 5×10⁴ M☉ (2% of the stellar mass) |
| [Baumgardt 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.2174B) | N-body grid | IMBH ~4×10⁴ M☉ |
| [van der Marel & Anderson 2010](https://ui.adsabs.harvard.edu/abs/2010ApJ...710.1063V) | anisotropic spherical Jeans, HST PMs | core models fit with M_BH = 0; M_BH ≲ 1.2×10⁴ (1σ), ≲ 1.8×10⁴ M☉ (3σ) |
| [Zocchi+2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.468.4429Z) | radially anisotropic limepy | central kinematics fitted without an IMBH |
| [Zocchi+2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4713Z) | multimass limepy | ~5% of the mass in stellar-mass black holes |
| [Baumgardt+2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.5340B) | N-body grid + fast-star test | no IMBH; 4.6% of the mass in a central cluster of stellar-mass black holes |
| [Evans+2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4251E) | axisymmetric Jeans, dark + luminous components | non-luminous component ≲ 10⁶ M☉ within r_h preferred over constant M/L; a compact remnant core can supply up to ~5×10⁵ M☉ |
| [Dickson+2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529..331D) | multimass models | black-hole mass fraction above 5% |
| [Häberle+2024](https://ui.adsabs.harvard.edu/abs/2024Natur.631..285H) | seven fast-moving stars in the central 3″ (0.08 pc), multi-epoch HST; projected velocities 66–113 km/s against a no-IMBH escape velocity of 62 km/s | **IMBH ≥ 8,200 M☉** (robust bound: velocities and escape velocity only). Model-dependent: ≥ 21,100 M☉ (99%) with the measured accelerations (all consistent with zero within 3σ); an N-body comparison of the inner 10″ prefers 39,000–47,000 M☉ and disfavours no IMBH, a stellar-mass black-hole cluster and > 50,000 M☉, with the authors' caveat of low number statistics (full text, main text p. 3) |
| [Bañares-Hernández+2025](https://ui.adsabs.harvard.edu/abs/2025A%26A...693A.104B) | Jeans (GravSphere) + millisecond-pulsar accelerations | extended central mass ~2–3×10⁵ M☉ favoured over an IMBH; IMBH < 6×10³ M☉ (3σ) |

- **Our two models reproduce the two sides of the debate with the same
  data.** The black-hole mass of our Jeans model is the classic value
  (4–5×10⁴ M☉). That value comes from models with isotropic or imposed
  central anisotropy; in ours the black hole is tied to β₀ at its −½ bound
  and gains Δχ² ≈ 6 from a few bins inside 0.3 pc.
- **The DF model's dark mass is the "extended central mass".** Its remnant
  component, 2.7×10⁵ M☉ with a 2.2 pc scale (half-mass radius 2.8 pc),
  matches the 2–3×10⁵ M☉ that Bañares-Hernández+2025 favour from kinematics
  and pulsar timing, the ~5% stellar-mass black-hole fraction of
  Zocchi+2019, Baumgardt+2019 and Dickson+2024, and the ≲ 5×10⁵ M☉ remnant
  core that Evans+2022 allow. Our component is a generic Plummer, so it can
  also contain white dwarfs and neutron stars.
- **The fast-moving stars set a robust lower bound on a point mass.**
  Häberle+2024 find ≥ 8,200 M☉ from the velocities and escape velocity
  alone. This already exceeds the 3σ upper limit of Bañares-Hernández+2025
  (6×10³ M☉): the published results are in tension. Their model-dependent
  estimates are higher: ≥ 21,100 M☉ (99%) when the accelerations (all
  within 3σ of zero) are modelled, and 39,000–47,000 M☉ preferred by a
  comparison with N-body models of the inner 10″, which also disfavours a
  stellar-mass black-hole cluster without an IMBH; the authors caution that
  low number statistics limit that comparison. The 21,100 M☉ limit would
  also exceed the 3σ limit of van der Marel & Anderson 2010 (1.8×10⁴ M☉).
- **Our two models map onto that tension.** The Jeans black hole
  (4.3×10⁴ M☉) lies inside Häberle's preferred range; the DF model's
  extended remnants match Bañares-Hernández. The DF model has no point mass;
  refitting it with a black hole of 2–5×10⁴ M☉ is the natural next test.

### Remnant mass and spatial scale

Numbers from the full texts (arXiv versions, read 2026-10-14), with their
location in each paper. Values marked *computed* are ours, derived from the
paper's numbers as shown.

| Study | Component (definition) | Mass [M☉] | Fraction | Form and scale | Total mass, distance | Where |
|---|---|---|---|---|---|---|
| **This work, DF base model** | remnant Plummer: all mass not following the counted light | **2.73×10⁵** | 8.4% | Plummer, **a = 2.18 pc** (3D; half-mass 2.85 pc) | 3.26×10⁶, 5.356 kpc | [base model](models/base_no_dm.md) |
| **This work, Jeans model** | remnant Plummer + black hole | **5.4×10⁵ + 4.3×10⁴** | 16.5% + 1.3% | Plummer, **a = 5.3 pc** (3D) | 3.24×10⁶, 5.36 kpc | [Jeans model](models/jeans_comparison.md) |
| [Bañares-Hernández+2025](https://ui.adsabs.harvard.edu/abs/2025A%26A...693A.104B) | "extended central mass" on top of a light-traced M★ (which includes remnants that follow the light); interpreted as stellar-mass black holes | **2.65 ± 0.19×10⁵**; 2.20×10⁵ without the pulsar accelerations; 2.1–3.1×10⁵ across centres and profile shapes | 9.0% (*computed*: 2.65/(26.8+2.65)) | Plummer, **r_cen = 1.88 ± 0.10 pc** (3D; 1.59 pc without pulsars); Gaussian alternative σ = 1.37 pc | M★ = 2.68 ± 0.05×10⁶; 5.2 kpc fixed | Fig. 2, Table 2, Secs. 6, 8–9, App. B–D |
| [Baumgardt+2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.5340B) | stellar-mass black holes only | **1.65×10⁵** | 4.6 ± 0.5% | mass-segregated N-body; radius not given ("centrally concentrated") | 5.24 kpc adopted | Secs. 4–5 |
| [Zocchi+2019](https://ui.adsabs.harvard.edu/abs/2019MNRAS.482.4713Z) | black holes as the dark component of two-component limepy models (white dwarfs and neutron stars in the visible component) | 1.35×10⁵ (*computed*: F₂ × M) | 4.5 (+4.0/−2.5)% anisotropic; 5.2% isotropic | not fully segregated, extending to about the half-mass radius (r_h = 8.57 pc); ~1.4×10³ M☉ within 0.2 pc | 3.01×10⁶, 5.14 kpc (fitted) | Table 2, Secs. 3–4 |
| [Dickson+2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529..331D) | (i) black holes; (ii) all dark remnants (white dwarfs, neutron stars, black holes) | (i) **1.82×10⁵** (12,400 black holes, mean 14.8 M☉) | (i) 5.9 ± 0.2% (statistical only); (ii) 50.3% | multimass limepy; radius not given ("concentrated within the central regions") | in their Paper I | Table 1, Sec. 5.3 |
| [Evans+2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4251E) | generic non-luminous component (remnants or dark matter) | 10⁴–10⁶ within the ~6 pc half-light radius (posterior contours only) | — | NFW with **r_max ≲ 1 pc** (r_s = r_max/2.16), or Gaussian | 5.24 kpc fixed | Secs. 6–7, Figs. 5, 7 |

**Like-for-like with the closest model (Bañares-Hernández+2025, same Plummer form):**

| | DF base model | Bañares-Hernández+2025 |
|---|---|---|
| Mass | 2.73×10⁵ M☉ | 2.65 ± 0.19×10⁵ M☉ |
| Plummer scale | 2.18 pc | 1.88 ± 0.10 pc at 5.2 kpc; 1.94 pc scaled to 5.356 kpc (*computed*) |
| Mass fraction | 8.4% | 9.0% |
| M(<1 pc), M(<3 pc) (*computed*) | 2.0×10⁴, 1.45×10⁵ M☉ | 2.7×10⁴, 1.61×10⁵ M☉ |
| Data | HST, MUSE, Gaia kinematics; star counts | stellar kinematics, surface density, millisecond-pulsar positions and accelerations |
| Method | self-consistent action DF | Jeans (GravSphere) |

- **The DF remnant component reproduces the Bañares-Hernández central mass:**
  the same mass within 3%, a scale 12% larger once both are at the same
  distance, and the same mass fraction. They reach it with a different method
  and with pulsar timing; from stellar kinematics alone they find a smaller
  and more compact component (2.20×10⁵ M☉, 1.59 pc).
- **The black-hole-only estimates are lower,** 1.35–1.82×10⁵ M☉ (4.5–5.9%).
  That is expected: our component also holds whatever neutron stars and
  segregated white dwarfs do not follow the light, whereas those models put
  white dwarfs and neutron stars in the visible or light-traced mass
  (in Dickson+2024 all dark remnants make up 50% of the mass).
- **Scales:** besides Bañares-Hernández, only Evans+2022 constrain one
  (r_max ≲ 1 pc for an NFW), and Zocchi+2019 give a qualitative extent (to
  about the half-mass radius, 8.6 pc). The concentration differs between
  models: within 0.2 pc Zocchi's black-hole component holds ~1.4×10³ M☉,
  against ~300 M☉ for the Bañares-Hernández Plummer and ~200 M☉ for ours
  (*computed* for the Plummer models).
- **The Jeans model's remnant component** (5.4×10⁵ M☉ at 5.3 pc) is twice as
  massive and 2.5 times more extended than any of these; its central mass is
  carried by the black hole instead.
- **Our remnant mass is the least constrained of our numbers.** It trades
  against the stellar mass (±23% on mocks), and pinning the mass per counted
  star between 8 and 12 M☉ moves it from 1.6 to 4.4×10⁵ M☉.

How these remnant masses relate to the stellar mass and the IMF, and which are physically allowed: [Remnant mass and the stellar population](remnant_mass.md).

Caveat: the Bañares-Hernández values are from arXiv:2408.00939v6, which may
differ in detail from the journal version; in it, the Sec. 8.1 statement on
the centre dependence does not match the Fig. C.1 numbers.

## Anisotropy

| Study | Anisotropy |
|---|---|
| **This work, DF base model** | isotropic core; radial, peak β ≈ +0.15 near 8 pc; tangential beyond ~22 pc (−0.53 at 100 pc) |
| **This work, Jeans model** | tangential core (β₀ ≈ −0.52, imposed ≤ −½); radial +0.15 at 3–15 pc; tangential outside |
| [Häberle+2025](https://ui.adsabs.harvard.edu/abs/2025ApJ...983...95H) (oMEGACat VI) | isotropic in the core, increasing radial anisotropy outwards (HST field) |
| [Watkins+2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2598W) | mildly radial, β = 0.10 ± 0.02 (single value) |
| [van de Ven+2006](https://ui.adsabs.harvard.edu/abs/2006A%26A...445..513V) | close to isotropic inside ~10′, increasingly tangential outside |
| [Zocchi+2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.468.4429Z) | radially anisotropic models fit the central kinematics |

- **The DF profile matches the published picture:** isotropic core,
  radially biased within the HST field, tangential in the outskirts. Our
  zero crossing (~22 pc, 14′ at our distance) is near van de Ven's ~10′
  transition.
- **The Jeans model's tangential core disagrees** with the isotropic core
  measured directly by oMEGACat VI. It is a consequence of the black hole
  (a cored tracer around a point mass needs β₀ ≤ −½), which is a further
  argument against taking its black-hole mass at face value.

## Dark matter

[Evans+2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4251E) find a
non-luminous component of ≲ 10⁶ M☉ within the half-light radius strongly
preferred over a constant-M/L stellar model and leave open that part of it is
dark matter. Our models agree that a centrally concentrated non-luminous mass
is needed (2.7×10⁵ M☉ in the DF model, ~6×10⁵ M☉ of remnants plus black hole
in the Jeans model), but find that it can be supplied by remnants
without dark matter. Our dark-matter variants show that the data constrain the
total mass profile, not how it splits between stars, remnants and dark matter
([Results](results.md#base-model-fits)).

## Gaps

- No comparison yet with published axisymmetric or rotating models on our
  data; ours are spherical.
- The published total masses use different data and extents; a like-for-like
  test (our pipeline on their data, or theirs on ours) has not been done.
- Earlier literature numbers assembled for the Jeans fits:
  [Literature baselines](../LITERATURE_BASELINES.md).
