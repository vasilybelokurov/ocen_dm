# Omega Centauri Dark-Matter/Tidal-Tail Inference Project

> **For Codex CLI / Claude Code:** Treat this file as the project specification and implementation roadmap. Work incrementally, keep all scientific calculations in importable/testable modules, and do not replace missing data or metadata with invented values. Every external dataset must have a provenance record, URL/DOI, checksum when available, and a script that reproduces the local processed product from the raw download.

## Goal

Build a reproducible inference pipeline that asks:

> **Can the observed tidal tails of Omega Centauri break the mass-decomposition degeneracy left by its internal stellar kinematics, and thereby constrain a possible surviving dark-matter halo?**

The central experiment is intentionally hierarchical:

1. Fit the **bound stellar kinematics** of Omega Cen with a family of mass models containing luminous stars, centrally concentrated dark remnants, an optional IMBH, and a possible extended DM component.
2. Use the kinematic posterior to identify **iso-kinematic solutions**: models that fit the internal data comparably well but contain different amounts and radial distributions of DM.
3. Evolve those models on the same Galactic orbit and predict their tidal-tail morphology and phase-space structure.
4. Compare the predictions with the observed Omega Cen tails, Fimbulthul, and tail spectroscopy.
5. Reweight or jointly refit the kinematic posterior using the tail likelihood.
6. Only after establishing a robust posterior on the extended DM distribution, derive quantities such as `rho_DM(r)`, `M_DM(<r)`, escape-speed profiles, decay `D` factors, and—more cautiously—annihilation `J` factors.

The scientific focus is **not** to force a best-fit NFW halo. The primary target is the radial distribution of extra mass, especially on scales where tidal stripping is sensitive to the gravitational potential.

---

# 1. Scientific hypothesis

Internal stellar kinematics mainly constrain the total potential in the central/bound cluster. Several physically distinct decompositions may produce nearly identical observables:

\[
\Phi_{\rm tot}
= \Phi_\star
+ \Phi_{\rm rem}
+ \Phi_{\rm IMBH}
+ \Phi_{\rm DM}.
\]

A centrally concentrated remnant population or an IMBH can mimic extra central mass, while an extended surviving DM halo can contribute additional mass at tens to hundreds of parsecs.

The tails respond to the **total potential near the escape/Jacobi region**, not to the microscopic identity of the mass. Therefore the experiment has leverage when two kinematically acceptable solutions satisfy approximately

\[
M_A(<5\,{\rm pc}) \simeq M_B(<5\,{\rm pc})
\]

but differ at larger radii,

\[
M_A(<50{-}200\,{\rm pc}) \neq M_B(<50{-}200\,{\rm pc}).
\]

The tails are not expected to measure an unresolved central cusp slope directly. They should be treated primarily as a constraint on the **radial extent and enclosed mass of any extra component**.

## Primary derived parameters

Report posterior constraints on quantities that are close to observables:

- `M_DM(<5 pc)`
- `M_DM(<20 pc)`
- `M_DM(<60 pc)`
- `M_DM(<100 pc)`
- `M_DM(<r_J)` where `r_J` is the present-day Jacobi scale for each model
- `rho_DM(20 pc)` as a useful derived quantity for comparison with the deep WD field
- `v_esc(r)`
- total bound mass as a function of radius

Avoid presenting `rho_s`, `r_s`, or a formal inner slope `gamma` as the main result unless injection/recovery tests show they are identifiable.

---

# 2. Key observational datasets

## 2.1 Tidal-tail / outer-cluster data

### A. Kuzma & Ishigaki (2025): Pristine + Gaia Omega Cen periphery

**Paper:** *A Pristine View of Galactic Globular Clusters and their Peripheries: Omega Centauri*  
ArXiv: https://arxiv.org/abs/2502.01135  
MNRAS: https://academic.oup.com/mnras/article/537/3/2752/8002861  
Dataset DOI: https://doi.org/10.5281/zenodo.14791430

Zenodo contains:

- `wCen_table.fits`
- size: 26.5 MB
- published MD5: `3dcd58e69767901fbbdf69385a14c768`

Use this as the primary wide-field tail catalogue. It combines Gaia astrometry with photometric/chemical information and probabilistic membership information out to several degrees.

**Important likelihood rule:** do not simply use the published membership probability as a likelihood weight and then reuse the same proper-motion/photometric information in the scientific likelihood. That risks double-counting. Use published membership probabilities for diagnostics and initial sample exploration; the final likelihood should either reconstruct the relevant mixture model or clearly condition on a defined selection.

### B. Kuzma et al. (2026): spectroscopy of the periphery and tails

**Paper:** *Disruption of a Giant: Spectroscopic Identification of Members in the Periphery and Tidal Tails of Omega Centauri*  
ArXiv: https://arxiv.org/abs/2605.23474  
MNRAS DOI: https://doi.org/10.1093/mnras/stag1147  
Publisher page: https://academic.oup.com/mnras/article/550/2/stag1147/8709286

The paper provides the full star list and measurements as online supplementary material. The underlying ESO observations are available under:

- programme ID: `108.22MM.001`
- PI: Kuzma

Use the published table first. Raw reduction of FLAMES spectra is not required for the first project version.

**Important likelihood rule:** the spectroscopy is targeted. Use velocities and metallicities **conditional on the observed target positions**. Do not use the number of spectroscopic detections as an unbiased measurement of the tail surface density unless the target-selection function is explicitly reconstructed.

### C. Fimbulthul stream

**Paper:** Ibata et al. (2019), *Identification of the long stellar stream of the prototypical massive globular cluster Omega Centauri*  
DOI: https://doi.org/10.1038/s41550-019-0751-x  
Publisher: https://www.nature.com/articles/s41550-019-0751-x

The Fimbulthul candidate stars are supplied in Supplementary Table 1.

Use Fimbulthul mainly as a **long-baseline stream-track constraint**. Because the original sample comes from STREAMFINDER and has a complex selection function, do not initially interpret its absolute number density as an unbiased tail density.

### D. Kuzma et al. (2021): Gaia EDR3 tail morphology

**Paper:** *Detecting globular cluster tidal extensions with Bayesian inference - I. Analysis of Omega Centauri with Gaia EDR3*  
MNRAS: https://academic.oup.com/mnras/article/507/1/1127/6347368

Useful for reproducing earlier tail morphology, radial star-count profiles, ellipticity, and membership methodology.

### E. Youakim, Lind & Kushniruk (2023): long-baseline debris and stream simulations

**Paper:** *Tidal debris from Omega Centauri discovered with unsupervised machine learning*  
ArXiv: https://arxiv.org/abs/2307.03035  
DOI: https://doi.org/10.1093/mnras/stad1952

This is an important validation target because it already used particle-spray stream models for Omega Cen-like orbits. Reproduce at least one of their no-DM stellar-progenitor calculations as an integration test before introducing an extended DM component.

---

## 2.2 Bound/internal kinematics

### A. oMEGACat VI: kinematic profiles and maps

**Paper:** *oMEGACat. VI. Analysis of the overall kinematics of Omega Centauri in 3D: velocity dispersion, kinematic distance, anisotropy, and energy equipartition*  
ArXiv: https://arxiv.org/abs/2503.04903  
DOI: https://doi.org/10.3847/1538-4357/adbe67  
Data: https://zenodo.org/records/14978551

The public data products include:

- proper-motion dispersion profiles
- radial and tangential PM components
- LOS rotation/dispersion profiles
- kinematic maps
- mass-binned PM profiles
- quality-selection definitions

**Use these binned products for v1 of the kinematic inference.** They make the initial model and likelihood tractable and avoid immediately fitting millions of individual stars.

### B. oMEGACat II: individual HST astrometry

**Paper:** https://arxiv.org/abs/2404.03722  
Data DOI: https://doi.org/10.5281/zenodo.11104046

The astrometric catalogue contains approximately 1.4 million proper-motion measurements. Use only after the binned-profile pipeline is validated.

### C. oMEGACat I / MUSE spectroscopy

Survey homepage and data access: https://omegacatalog.github.io/

Use individual LOS velocities only in a later likelihood version after the profile-based model is established.

---

# 3. Modelling philosophy

## 3.1 Stage 1: fit the equilibrium/bound cluster

Start from an axisymmetric or spherical-equivalent model of the bound cluster. The mass decomposition is

\[
\Phi_{\rm tot}
= \Phi_\star
+ \Phi_{\rm rem}
+ \Phi_{\rm IMBH}
+ \Phi_{\rm DM}.
\]

### Luminous stars

Represent the projected light/star-count distribution using an MGE or another smooth basis that can be deprojected consistently. Prefer deriving the luminous model from public Omega Cen photometric/star-count data or a reproducible literature profile rather than hard-coding a single analytic Plummer model.

Free nuisance parameters should include a stellar mass normalization or effective `M/L` and, if required, a modest radial variation.

### Dark remnants

Use a compact extended component rather than forcing remnants into the luminous profile. Initial simple option:

- Plummer or multi-Gaussian remnant density
- free total remnant mass `M_rem`
- free scale radius `a_rem`

This component is phenomenological in v1. Later it may be replaced by a multi-mass collisional cluster prediction.

### IMBH

Point mass with parameter `M_bh >= 0`.

### DM

Use at least two families in separate inference runs:

1. truncated generalized NFW / NFW-like cusp
2. cored profile such as Burkert or cored-gNFW

A suitable generic form is

\[
\rho_{\rm DM}(r)
= \rho_s
\left(\frac{r}{r_s}\right)^{-\gamma}
\left(1+\frac{r}{r_s}\right)^{\gamma-3}
T(r;r_t),
\]

where `T` is a smooth outer truncation.

Do not let the first inference have more halo-shape freedom than the data can support. Start with fixed `gamma=0` and `gamma=1` families. Promote `gamma` to a free parameter only after injection tests show it is recoverable.

### Anisotropy

Use a smooth low-dimensional radial anisotropy parameterization, for example

\[
\beta(r)
= \beta_0
+ (\beta_\infty-\beta_0)
\frac{r^2}{r^2+r_\beta^2}.
\]

Do not use many independent anisotropy bins in the first fit; that can absorb the mass-profile signal.

### Rotation

Include rotation consistently with the chosen Jeans formalism and the oMEGACat rotation constraints. Keep the initial rotation model low dimensional.

---

# 4. Recommended software stack

## Core Python

Recommended baseline: Python 3.12 in a dedicated environment.

Required:

- `numpy`
- `scipy`
- `astropy`
- `pandas` or `polars`
- `pyarrow`
- `h5py` and/or `zarr`
- `matplotlib`
- `corner`
- `pooch`
- `pyyaml`
- `pytest`
- `ruff`

## Kinematic modelling

### Primary: JamPy

PyPI: https://pypi.org/project/jampy/

JamPy 9.x includes the 2026 spectral axisymmetric Jeans solver and supports projected proper-motion and LOS moments with general anisotropy.

**Licensing note:** JamPy is not an ordinary permissively licensed dependency; verify its current research/non-commercial terms before packaging or redistribution. Do not vendor it into this repository.

### Cross-check / alternative

AGAMA: https://github.com/GalacticDynamics-Oxford/AGAMA

AGAMA supports potentials, orbit integration, DFs, moments, self-consistent multi-component models, and construction of equilibrium realizations. It will also be used later for equilibrium N-body initial conditions.

## Stream modelling

### Primary fast engine: galpy particle spray

Documentation: https://docs.galpy.org/en/latest/tutorials/streams/streamspraydf.html

Use the current Chen/Fardal spray implementations. Current galpy supports:

- `chen24spraydf`
- `fardal15spraydf`
- a progenitor potential via `progpot`
- time-dependent/callable progenitor mass
- leading/trailing or both tails
- smooth stream-track fitting with phase-space covariance

This makes it suitable for quickly propagating many posterior mass models.

### Independent cross-check: Gala

`MockStreamGenerator` documentation: https://gala-astro.readthedocs.io/en/latest/api/gala.dynamics.mockstream.MockStreamGenerator.html

Gala supports an explicit progenitor potential and is useful for reproducing earlier Omega Cen stream calculations.

## Inference

Primary:

- `ultranest`: https://johannesbuchner.github.io/UltraNest/

Cross-check:

- `dynesty`

Nested sampling is useful because the kinematic mass decomposition may be multimodal and because model-family evidences can be informative, although evidence ratios are not the main scientific product.

## Workflow

- `snakemake` for reproducible data -> inference -> simulations -> figures
- `pre-commit` optional
- `uv` or `conda-lock` for environment reproducibility

---

# 5. Repository layout

```text
ocen-dm-tails/
├── README.md
├── PROJECT.md
├── pyproject.toml
├── environment.yml
├── configs/
│   ├── data.yaml
│   ├── kinematics_baseline.yaml
│   ├── dm_nfw.yaml
│   ├── dm_cored.yaml
│   ├── mw_axisym.yaml
│   ├── mw_barred.yaml
│   └── experiments/
├── data/
│   ├── raw/          # never edited manually; gitignored
│   ├── interim/      # reproducible intermediate products
│   └── processed/    # analysis-ready products
├── provenance/
│   ├── datasets.yaml
│   └── checksums.txt
├── src/ocen_dm/
│   ├── data/
│   │   ├── download.py
│   │   ├── kuzma2025.py
│   │   ├── kuzma2026.py
│   │   ├── fimbulthul.py
│   │   └── omegacat.py
│   ├── coordinates.py
│   ├── light_model.py
│   ├── mass_models/
│   │   ├── stellar.py
│   │   ├── remnants.py
│   │   ├── imbh.py
│   │   ├── dark_matter.py
│   │   └── composite.py
│   ├── kinematics/
│   │   ├── jam_model.py
│   │   ├── likelihood.py
│   │   └── derived.py
│   ├── inference/
│   │   ├── priors.py
│   │   ├── ultranest_runner.py
│   │   └── posterior.py
│   ├── stream/
│   │   ├── orbit.py
│   │   ├── galpy_spray.py
│   │   ├── gala_spray.py
│   │   ├── observables.py
│   │   └── likelihood.py
│   ├── selection/
│   │   ├── kuzma2025.py
│   │   ├── spectroscopy.py
│   │   └── fimbulthul.py
│   ├── experiments/
│   │   ├── select_isokinematic.py
│   │   ├── injections.py
│   │   ├── reweight.py
│   │   └── posterior_predictive.py
│   ├── nbody/
│   │   ├── agama_ic.py
│   │   └── export.py
│   ├── plotting/
│   └── cli.py
├── tests/
├── notebooks/        # inspection only; no unique science logic
├── workflows/
│   └── Snakefile
└── results/
```

## CLI target

Expose a small command surface, for example:

```bash
ocen fetch-data
ocen preprocess
ocen fit-kinematics --config configs/kinematics_baseline.yaml
ocen select-isokinematic --posterior results/kinematics/posterior.h5
ocen simulate-tails --ensemble results/isokinematic/models.parquet
ocen fit-tails
ocen reweight-posterior
ocen run-injections
ocen compare-mw-potentials
ocen posterior-predictive
ocen export-nbody
```

---

# 6. Experiment sequence

## Experiment K0: reproduce published kinematic observables

Before inference, load the oMEGACat VI products and reproduce the published profile figures numerically:

- PM radial dispersion
- PM tangential dispersion
- LOS dispersion
- rotation
- anisotropy diagnostic

Deliverable:

`results/validation/omegacat_vi_profiles.*`

Failure here blocks subsequent inference.

## Experiment K1: no-DM baseline fit

Fit a model containing:

- luminous stars
- compact remnant component
- optional IMBH
- anisotropy
- rotation
- no particle DM

This establishes whether the implementation can reproduce the equilibrium kinematics without using DM as a catch-all component.

Outputs:

- posterior samples
- posterior predictive kinematic profiles
- residual maps/profiles
- derived enclosed-mass posterior

## Experiment K2: DM-enabled kinematic fits

Run at least two model families:

1. cored DM
2. cuspy DM

Keep the remainder of the model identical.

Store derived quantities for every posterior sample:

```text
M_star(<r)
M_rem(<r)
M_DM(<r)
M_total(<r)
f_DM(<r)
r_J
v_esc(r)
```

at a fixed common radial grid.

The key product is not a best fit. It is the posterior degeneracy between the compact and extended components.

## Experiment K3: construct the iso-kinematic ensemble

Do **not** invent arbitrary halo models at this stage.

Sample directly from `p(theta | D_kin)` and stratify the accepted models by a physically relevant extended-DM variable such as

```text
M_DM(<100 pc)
M_DM(<r_J)
```

Suggested v1 selection:

- divide the posterior into quantile bins in `M_DM(<r_J)`
- draw several posterior samples from each bin
- preserve posterior weights
- optionally require posterior-predictive residual metrics to pass a common quality threshold

The result should contain models spanning low-DM through high-DM solutions while all remaining valid fits to the same kinematic data.

Deliverable:

`results/isokinematic/models.parquet`

Each row must contain all mass-model, anisotropy, orbit, and nuisance parameters needed to regenerate the model exactly.

---

# 7. Tail forward modelling

## Experiment T0: no-DM literature validation

Before using the isokinematic ensemble:

1. integrate the present-day Omega Cen orbit in the baseline Milky Way potential;
2. generate a no-DM particle-spray stream;
3. reproduce the approximate observed orientation/track of the local tails;
4. reproduce at least one configuration comparable to Youakim et al. (2023) or Ibata et al. (2019).

This validates coordinates, orbit conventions, time direction, units, and stream-release machinery.

## Experiment T1: controlled compact-vs-extended mass test

Construct pairs of models with similar central enclosed mass but different outer mass profiles.

Hold fixed:

- present-day cluster phase space
- Milky Way potential
- stellar mass model where possible
- random numbers used in particle release
- integration settings

Vary only the mass decomposition.

Compare:

- surface-density morphology
- tail width
- leading/trailing asymmetry
- PM track
- LOS velocity track
- stripping-time distribution
- bound mass evolution
- stellar escape rate

Use common random numbers across nearby models to suppress Monte Carlo noise in model comparisons.

## Experiment T2: propagate the kinematic posterior

For every selected isokinematic model:

1. construct the corresponding progenitor potential;
2. generate both leading and trailing tails;
3. transform to heliocentric observables;
4. apply the relevant observational footprint/selection;
5. store a compact summary plus the mock particles needed for likelihood evaluation.

This experiment maps

\[
\theta_i \sim p(\theta|D_{\rm kin})
\quad\longrightarrow\quad
D_{\rm tail}^{\rm mock}(\theta_i).
\]

---

# 8. Tail likelihood

The first scientifically useful tail likelihood should be modular:

\[
\ln {\cal L}_{\rm tail}
=
\ln {\cal L}_{\rm sky}
+
\ln {\cal L}_{\rm PM}
+
\ln {\cal L}_{v_{\rm los}}
+
\ln {\cal L}_{\rm Fimbulthul}.
\]

## 8.1 Sky morphology

Preferred long-term solution: a Poisson point-process likelihood containing cluster/tail and field components plus the angular selection function.

For a fast v1, compare smoothed tail-density fields only in regions where the published sample selection is understood. Clearly label this as an approximation.

## 8.2 Proper motions

Compare the Gaia PM data to the local mock-stream phase-space distribution convolved with the full Gaia astrometric covariance where available.

Do not assume independent `pmra` and `pmdec` errors if covariance information exists.

## 8.3 Kuzma et al. 2026 radial velocities

Condition on target positions and compare observed `v_los` with the model distribution at those positions. Include measurement errors and an outlier/field component if required.

## 8.4 Fimbulthul

Initially use the stream as a track constraint in sky position / PM / LOS velocity where available. Avoid an absolute surface-density likelihood until STREAMFINDER selection is modeled.

---

# 9. Posterior updating

The conceptual update is

\[
p(\theta|D_{\rm kin},D_{\rm tail})
\propto
{\cal L}_{\rm tail}(\theta)
\,p(\theta|D_{\rm kin}).
\]

The first implementation should perform **importance reweighting** of the kinematic posterior:

\[
w_i \propto {\cal L}_{\rm tail}(\theta_i).
\]

Compute effective sample size. If the tail likelihood collapses the posterior so strongly that the effective sample size is poor, switch to a joint sampler or sequential Monte Carlo initialized from the kinematic posterior.

Primary comparison plots:

- `p[M_DM(<r) | kinematics]`
- `p[M_DM(<r) | kinematics + tails]`

and similarly for

- `f_DM(<r)`
- `v_esc(r)`
- `M_rem`
- `M_bh`

The central scientific diagnostic is whether the tails narrow the posterior on **extended** DM while leaving compact-remnant/IMBH degeneracies comparatively unchanged.

---

# 10. Injection/recovery experiments

These are mandatory before interpreting the real-data posterior.

Generate synthetic observations using the real angular footprint and error model for at least:

1. no DM;
2. low extended-DM mass;
3. intermediate extended-DM mass;
4. high extended-DM mass;
5. compact extra mass with little extended DM;
6. cored and cuspy halos with matched `M_DM(<r_J)`.

For each injection:

- generate internal kinematic observables;
- fit the kinematics without revealing the truth to the recovery stage;
- select the resulting iso-kinematic ensemble;
- generate tail predictions;
- reweight with the mock tail data;
- evaluate recovery of `M_DM(<r)`.

Questions the injection suite must answer:

- At what radius is `M_DM(<r)` best constrained?
- Can tails distinguish compact from extended extra mass?
- Can they distinguish cored from cuspy halos at fixed outer enclosed mass?
- How much bias is introduced by particle spray?
- How strongly do results depend on the Milky Way potential?
- Does the inference spuriously detect DM when the truth is no DM?

Do not quote a DM constraint from real data until the no-DM false-positive test passes.

---

# 11. Milky Way potential systematics

Use a hierarchy:

## MW-A: baseline static axisymmetric potential

Use for development and the first posterior.

## MW-B: alternative static potential

Repeat the full inference to measure sensitivity to reasonable mass-model changes.

## MW-C: barred/time-dependent model

Use AGAMA or another validated integrator for a rotating/time-dependent bar. This is a later systematic, not part of the minimal viable experiment.

The nearby/current tails should be emphasized before attempting to infer the entire ancient accretion history, because the latter introduces much larger uncertainties in Galactic evolution and dynamical friction.

---

# 12. From particle spray to live N-body

Particle spray is the **inference engine**, not the final validation.

After the tail-reweighted posterior is obtained, select representative models:

- no/very-low DM
- median allowed extended DM
- upper allowed extended DM
- cored and cuspy variants
- compact-remnant-dominated solution

Use AGAMA to construct equilibrium multi-component DFs / particle realizations. Evolve each model in isolation first.

Required equilibrium check:

- density profile remains stable
- velocity-dispersion profile remains stable
- no artificial rapid expansion/contraction

Then evolve representative systems through the same Galactic orbit with a live N-body code or appropriate collisionless solver.

Compare particle-spray and live results for:

- stripping rate
- tail width
- stream track
- PM/LOS gradients
- bound-mass evolution

Use the comparison to calibrate or bracket particle-spray modelling error.

Do not begin expensive N-body parameter inference before the spray-level injection tests establish useful sensitivity.

---

# 13. Testing requirements

Every implementation task must add tests. At minimum:

1. **Units:** all internal mass/radius/velocity conversions are tested.
2. **Coordinates:** sky -> Galactocentric -> sky round trip agrees within numerical tolerance.
3. **Mass profiles:** analytic/numerical enclosed masses agree for Plummer, point mass, NFW/gNFW, and cored models.
4. **Potential gradients:** numerical finite-difference forces agree with analytic/package forces.
5. **Kinematic limiting cases:** no-DM/remnant-only models reproduce expected qualitative profiles.
6. **Orbit cross-check:** galpy and Gala/AGAMA agree for a simple static potential to a documented tolerance.
7. **Stream reproducibility:** fixed seed produces bitwise or numerically stable mock summaries where feasible.
8. **No-DM injection:** does not produce a false statistically significant extended-DM signal.
9. **DM injection:** recovers the injected outer enclosed mass within calibrated coverage.
10. **Timestep convergence:** stream observables are stable to timestep refinement.
11. **Particle-number convergence:** likelihood summaries are stable to increased spray-particle count.
12. **Selection-function test:** the same selection function applied to real and mock catalogues produces identical bookkeeping/columns.
13. **Posterior reweighting:** a constant tail likelihood leaves the kinematic posterior unchanged.
14. **Importance weights:** effective sample size and weight normalization are unit tested.

---

# 14. Implementation milestones for the coding agent

## Milestone 1: repository + provenance + data ingestion

Deliver a working repository where

```bash
ocen fetch-data
ocen preprocess
```

produces analysis-ready tables for:

- Kuzma & Ishigaki 2025
- Kuzma et al. 2026 supplementary spectroscopy
- oMEGACat VI
- Fimbulthul supplementary table

No science fit yet.

Acceptance criteria:

- raw files are not committed;
- checksums/provenance are saved;
- processed products have explicit schemas;
- loaders have tests.

## Milestone 2: mass-profile library

Implement and test:

- luminous MGE / baseline light model
- compact remnant component
- point-mass IMBH
- truncated cored DM
- truncated cuspy DM
- composite enclosed mass / potential / force / escape speed

Acceptance criteria:

- all profiles produce finite, physically consistent values over the configured radial range;
- enclosed-mass curves are independently numerically verified.

## Milestone 3: kinematic likelihood

Implement the oMEGACat VI profile likelihood using JamPy or the agreed solver.

Acceptance criteria:

- no-DM baseline fit completes;
- posterior-predictive profiles reproduce the observed data within the expected residual distribution;
- a synthetic injected model is recoverable.

## Milestone 4: DM kinematic posterior

Run cored and cuspy model families and save derived radial mass profiles for every sample.

Acceptance criteria:

- posterior includes compact-remnant/IMBH/DM degeneracies;
- results are stored in a restartable machine-readable format;
- no result depends on notebook-only code.

## Milestone 5: iso-kinematic model selection

Create a stratified ensemble spanning the posterior in `M_DM(<r_J)` and/or `M_DM(<100 pc)`.

Acceptance criteria:

- each selected model regenerates its kinematic prediction exactly;
- the ensemble covers low-to-high allowed extended DM without manually modifying parameters.

## Milestone 6: stream baseline

Implement Omega Cen orbit integration and no-DM particle spray.

Acceptance criteria:

- leading/trailing tails appear in the correct broad orientation;
- a literature no-DM stream test is reproduced at the expected qualitative/quantitative level;
- galpy/Gala comparison is documented.

## Milestone 7: propagate iso-kinematic models into tails

Run the same stream pipeline for the full selected ensemble.

Acceptance criteria:

- simulations use common random numbers where appropriate;
- outputs include sky, PM, LOS and stripping-time information;
- all simulations carry model IDs back to the kinematic posterior.

## Milestone 8: tail likelihood and posterior reweighting

Implement modular tail likelihoods and importance reweighting.

Acceptance criteria:

- constant-likelihood test preserves the original posterior;
- effective sample size is reported;
- before/after radial DM credible bands are produced.

## Milestone 9: injection/recovery suite

Run the mandatory synthetic experiments.

Acceptance criteria:

- false-positive rate and coverage are measured;
- the radius at which DM enclosed mass is actually identifiable is reported;
- particle-spray bias is quantified at least approximately.

## Milestone 10: Milky Way systematics

Repeat the key result under at least one alternative Galactic potential.

Acceptance criteria:

- the DM constraint is presented with potential-model sensitivity explicitly separated from statistical uncertainty.

## Milestone 11: live N-body validation

Only after Milestones 1-10 are scientifically stable.

Acceptance criteria:

- representative equilibrium models remain stable in isolation;
- live simulations quantify the error made by the spray approximation.

---

# 15. Expected science outputs

The first paper-quality result should emphasize plots like:

1. **Kinematic degeneracy:** posterior samples showing equally good internal fits with different `M_DM(<r)`.
2. **Mass profiles:** `M_star`, `M_rem`, `M_DM`, and `M_total` as functions of radius for representative iso-kinematic models.
3. **Tail response:** mock tail maps for low-, medium-, and high-extended-DM solutions.
4. **Discriminator plot:** `ln L_tail` versus `M_DM(<r_J)` or `M_DM(<100 pc)` for models already acceptable to the internal kinematics.
5. **Posterior update:** `p[M_DM(<r) | kinematics]` versus `p[M_DM(<r) | kinematics + tails]`.
6. **Compact-vs-extended mass:** joint posterior in `M_rem`, `M_bh`, and `M_DM(<100 pc)`.
7. **Model-systematics panel:** baseline versus alternate Milky Way potential.
8. **Posterior predictive tails:** sky positions, PM trends, and LOS velocity trends with real data overlaid.

A null result is scientifically valid: if tail likelihoods are insensitive across the entire kinematically allowed DM range, report an upper bound on the information content of present tail data rather than forcing a DM constraint.

---

# 16. Non-goals for the first version

Do not initially:

- infer a DM annihilation cross section;
- interpret a central kinematic dark component as particle DM by definition;
- fit the entire ancient Omega Cen progenitor/debris history simultaneously;
- use the absolute Fimbulthul star count without a selection model;
- run a giant live-N-body grid before particle-spray identifiability is demonstrated;
- overfit anisotropy or DM inner slope with highly flexible functions;
- use published membership probabilities twice in the same likelihood;
- hide unit conversions or coordinate-frame conventions in notebooks.

---

# 17. Agent working rules

1. **Reproducibility first.** Every external file must be fetched by code or have a documented manual retrieval step plus checksum.
2. **No invented metadata.** If a catalogue column meaning, covariance convention, or selection rule is unclear, inspect the paper/data documentation and record the uncertainty.
3. **Tests before expensive runs.** Unit/injection tests must pass before launching posterior or stream grids.
4. **Keep models regenerable.** Every posterior/stream product must store the exact config and package/version metadata used to create it.
5. **No scientific logic only in notebooks.** Notebooks may call package functions for inspection and plotting only.
6. **Separate observational likelihoods.** Internal kinematic likelihood and tail likelihood must be callable independently.
7. **Preserve the sequential Bayesian structure.** The tail experiment is explicitly a test of models already allowed by internal kinematics.
8. **Report derived radial mass quantities.** Do not make raw NFW parameters the headline output.
9. **Treat the Milky Way potential as a systematic.** Do not silently fix it and present the result as purely statistical.
10. **Escalate fidelity only when justified.** Particle spray -> live collisionless N-body -> collisional modelling only if the previous level demonstrates sensitivity.

---

# 18. Immediate first task for Codex CLI / Claude Code

Start with **Milestone 1 only**.

The first implementation PR/commit series should:

1. scaffold the Python package and CLI;
2. create `provenance/datasets.yaml` with the four primary public datasets and their source URLs;
3. implement deterministic download/parsing for the Zenodo datasets;
4. document retrieval of publisher supplementary tables that cannot be fetched automatically;
5. implement typed/validated loaders returning standardized tables;
6. add tests for file checksums, required columns, units, and unique IDs;
7. create one validation command that prints a concise dataset inventory and produces no scientific inference.

Do **not** implement the Jeans model or stream model in the same first task. Get data provenance and schemas correct first.

Once Milestone 1 passes, proceed to the mass-profile library, then the kinematic fit, then the iso-kinematic ensemble, and only then tidal-tail simulation.

---

# 19. Compact statement of the experiment

The project tests the following chain:

\[
D_{\rm kin}
\rightarrow
p(\theta|D_{\rm kin})
\rightarrow
\{\text{iso-kinematic models spanning extended DM}\}
\rightarrow
\text{tidal-tail forward models}
\rightarrow
{\cal L}_{\rm tail}(\theta)
\rightarrow
p(\theta|D_{\rm kin},D_{\rm tail}).
\]

The central scientific question is:

> **Do the observed Omega Cen tidal tails constrain the radial extent and mass of a possible surviving DM component beyond what can be learned from the bound stellar kinematics alone?**

That question should remain the organizing principle for implementation choices, model complexity, validation tests, and final outputs.
