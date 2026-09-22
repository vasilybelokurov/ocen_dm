# oCen_dm — dark matter in Omega Centauri

Can Omega Centauri's tidal debris distinguish an extended dark-matter halo from
the stars, remnants and central mass that fit its internal kinematics?

The project has two connected approaches: infer the present mass distribution
from stellar kinematics, and estimate how much dark matter could survive the
progenitor dwarf's disruption. The eventual joint kinematic-and-tail inference
is still to be built.

## Current state — 21 September 2026

Data ingestion, mass components, a spherical Jeans likelihood and nested sampling
are implemented. The current comparison uses **89 binned measurements**: HST
radial and tangential proper motions (21 each), MUSE line-of-sight dispersions
(29), and our Gaia EDR3 radial and tangential profiles (9 each). Pristine,
periphery spectroscopy and Fimbulthul support the outer-cluster work; they do
not enter this likelihood.

Three isotropic, no-instrument-scale fits (rung 0) are complete:

| Family | Halo added to stars, remnants and a point mass | ln Z | χ² / number of points |
|---|---|---:|---:|
| K1 | none | −190.02 ± 0.48 | 772 / 89 |
| K2-cored | truncated gNFW, γ = 0 | −188.94 ± 0.36 | 761 / 89 |
| K2-NFW | truncated gNFW, γ = 1 | −190.09 ± 0.35 | 763 / 89 |

Rung 0 establishes the intended reference for the modelling sequence. There is
no persuasive halo preference within this baseline, and the structured residuals
outside the central 50 arcsec identify where additional freedom should be tested.
The component ratios motivate radial anisotropy at intermediate radii and tangential
anisotropy farther out, with rotation as a competing explanation. The baseline
posteriors and these diagnostics guide the next controlled modelling step.

Orbital histories, bar migration, friction experiments and approximate stripping
calculations are available. The Kepler energy-truncation formula has passed an
independent check. The [dynamical experiment pipeline](docs/DYNAMICAL_EXPERIMENTS.md)
now prepares joint spherical equilibria and runs live or frozen-potential remnant
and progenitor experiments on prescribed orbits. The first remnant single-passage
pilots, numerical controls and progenitor isolation runs have completed. The
[batch analysis](docs/DYNAMICAL_BATCH_ANALYSIS.md) finds mild inner remnant depletion
after one passage and stable progenitor profiles on well-sampled scales; the
progenitor core needs better central particle sampling. Longer tidal experiments, converged
survival predictions, particle-spray predictions and a tail likelihood remain
outstanding. The [lifetime batch](docs/LIFETIME_BATCH_REPORT.md) is running full-passage
controls, with responsive stellar-nucleus experiments and multi-Gyr frozen-potential
comparisons queued behind numerical checks. All three [rung-1 fits](docs/RUNG1_COMPARISON.md)
have completed: χ² = 342.51 (K1), 316.01 (NFW), 315.29 (core), for 89 measurements.
The constant-anisotropy extension improves substantially on rung 0; the core–NFW
evidence difference remains small.

## Reading guide

| Document | Purpose |
|---|---|
| [Code and analysis audit](docs/CODE_ANALYSIS_AUDIT.md) | Verified implementation, numerical checks, and unresolved issues |
| [Rung-0 analysis](docs/RUNG0_ANALYSIS.md) | Latest fit results and their interpretation |
| [Rung-1 comparison](docs/RUNG1_COMPARISON.md) | All three completed constant-anisotropy fits, rung-0 comparison and plots |
| [Data and models across all rungs](docs/RUNG_MODEL_COMPARISON.md) | Best-fit no-DM and DM profiles, including both rung-2 variants, with residuals |
| [DM density probability at 20 pc](docs/DM_DENSITY_POSTERIOR.md) | Evidence-weighted rung-2 density distributions, explicit model priors and a probability at zero for no DM |
| [Which data constrain DM density?](docs/DM_DENSITY_DATA_CONSTRAINTS.md) | Dataset and radial-bin sensitivity, Gaia residuals and the tradeoff with stellar/remnant mass |
| [Testing densities above the DM limit](docs/DM_DENSITY_PROFILE_LIMIT.md) | Fixed-density refits distinguish the Bayesian upper limit from fit deterioration and stellar-mass bounds |
| [High-density posterior checks](docs/DM_DENSITY_POSTERIOR_CHECKS.md) | Independent 800-live-point core repeat and rho20 = 2 conditional run; prior derivation, launch records and validation |
| [Density and DF consistency audit](docs/DF_CONSISTENCY_AUDIT.md) | All runs and reruns: positive spatial densities, central DF condition, and finite-radius tests with explicit separability assumptions |
| [Positive AGAMA DF models](docs/AGAMA_DF_MODELS.md) / [PDF](docs/agama_df_models.pdf) / [LaTeX](docs/agama_df_models.tex) | Self-consistent stellar action DFs, joint photometric/kinematic pilot fitting, numerical controls and scope |
| [Independent DF flexibility challenge](docs/DF_CAPACITY_CHALLENGE.md) | Three physical mocks, one/two/three positive DFs in known gravity, multiple starts, wider-bound controls and PNG comparisons |
| [Free-potential DF mass recovery](docs/DF_MASS_RECOVERY.md) / [live status](docs/DF_MASS_RECOVERY_STATUS.md) | Independent stellar mass/light weights, self-consistent gravity, noiseless recovery and gated noisy mock experiments |
| [Lifetime batch](docs/LIFETIME_BATCH_REPORT.md) | Active queue, live stellar response, numerical gates, long orbital forcing and physical heating estimate |
| [Modelling plan](docs/MODELLING_PLAN.md) | Adopted fit sequence and exact CLI choices |
| [Data preparation](docs/data_analysis.pdf) / [source](docs/data_analysis.tex) | Selections, estimators and remaining systematics |
| [Mass modelling](docs/mass_modelling.pdf) / [source](docs/mass_modelling.tex) | Model, priors, likelihood and limitations |
| [Progenitor orbits](docs/PROGENITOR_ORBITS.md) | Orbital-history classes and assumptions |
| [DM survival and WD note](docs/dm_capture_constraints.pdf) | Conditional density estimates and white-dwarf applications |
| [Section-11 tests](docs/SECTION11_TESTS.md) | Numerical work needed to validate survival predictions |
| [Dynamical experiments](docs/DYNAMICAL_EXPERIMENTS.md) | Initial conditions, live/frozen drivers, controls and diagnostics |
| [Dynamical batch analysis](docs/DYNAMICAL_BATCH_ANALYSIS.md) / [PDF](docs/dynamical_batch.pdf) / [LaTeX](docs/dynamical_batch.tex) | First-passage response, isolation controls, resolution and next experiments |
| [Original specification](OCEN_DM_TAILS_PROJECT.md) | Long-term scope; proposed interfaces are not implemented commands |
| [Journal](JOURNAL.md) | Chronological record; later entries can supersede earlier results |

Keep compiled PDFs in the same directory as their LaTeX sources.
Save all analysis plots as PNGs in `plots/`; do not export plots as PDFs or place
them in `docs/` or `results/`.
Keep `plots/` image-only. Plotted tables and JSON provenance belong in
`results/plot_data/`; migration and housekeeping records belong in
`results/maintenance/`.

## Use the existing checkout

The working environment is `~/Work/venvs/.venv`. From the repository root:

```bash
source ~/Work/venvs/.venv/bin/activate
export PYTHONPATH=src
python -m ocen_dm.cli --help
python -m ocen_dm.cli inventory
python -m ocen_dm.cli fit --help
```

The installed console entry point is **`ocen`**, declared in `pyproject.toml`.
`python -m ocen_dm.cli` works without installing the package. `ocen-dm` is not a
declared command.

Data-building commands write products:

```bash
python -m ocen_dm.cli fetch-data
python -m ocen_dm.cli preprocess
python -m ocen_dm.cli plot-data
```

Journal supplements need [manual retrieval](docs/MANUAL_DOWNLOADS.md).
`preprocess` standardises catalogues and published profiles; our measured HST
and Gaia products have separate builders in `kinematics/hst_profile.py` and
`kinematics/outer_gaia.py`. Their loaders build missing products automatically,
so fitting or plotting can write data on an incomplete checkout.

**Fits are launched only when requested, one rung at a time.** A bare `fit`
command still selects the older varying-anisotropy model with instrument scales.
Use explicit switches and unique output labels as shown in the
[modelling plan](docs/MODELLING_PLAN.md). Fitting refuses an existing run directory,
including an interrupted run; choose a fresh label.

## Code and products

```text
configs/                 dataset locations and column mappings
provenance/              source registry, download manifest and checksums
data/raw/                source downloads (not committed)
data/processed/          standardised and measured products (not committed)
src/ocen_dm/data/        retrieval, schema validation and catalogue loaders
src/ocen_dm/selection/   covariance joins, field templates and crossmatches
src/ocen_dm/kinematics/  estimators, dynamical engines, fits and reports
src/ocen_dm/mass_models/ stellar, remnant, point-mass and halo components
src/ocen_dm/tails/       orbital histories and exploratory survival models
src/ocen_dm/dynamics/    controlled remnant/progenitor initial conditions and evolution
results/                local fits, orbit tables and diagnostics (not committed)
results/plot_data/      plotted tables and figure provenance
results/maintenance/    migration and housekeeping records
plots/                  PNG figures only
docs/                   current notes, plans and dated reviews
tests/                  synthetic tests plus tests requiring local inputs
```

`workflows/` and `notebooks/` are placeholders. There is no complete Snakemake
workflow, environment lockfile, axisymmetric inference, or joint tail fit.

Internal dynamics use pc, M☉ and km/s; Galactic-orbit modules use their own
explicit conversions. The main engine is the local spherical Jeans solver.
JamPy supplies a numerical cross-check. The legacy AGAMA inversion backend is
a restricted diagnostic; the separate [positive stellar DF branch](docs/AGAMA_DF_MODELS.md)
derives density and anisotropy from action DFs and fits photometry jointly.
JamPy is an external non-commercial dependency and is not
redistributed here. UltraNest drives production sampling; a dynesty driver
is not implemented.

## Validation and provenance

Column roles resolve against documented candidates and fail on ambiguity.
Downloads record source information and hashes; processed tables carry their
likelihood-use restrictions. Fit outputs record posterior samples, priors,
settings and summary statistics. New runs save the actual likelihood arrays and
resolved model, including the tracer MGE, before sampling. Reports verify the
snapshot and reproduce the saved best-sample likelihood before plotting. Older
runs use recorded options and current products, with an explicit warning.
The [audit](docs/CODE_ANALYSIS_AUDIT.md) distinguishes its original findings from
the completed reporting, mock-handling and run-preservation repairs.

```bash
python -m pytest tests -q
```

The full suite includes local-data and optional-engine tests and can take tens
of minutes. Some tests skip without those inputs. See the dated audit for the
checks actually run; an old test count or saved log is not a current pass claim.

Repository: <https://github.com/vasilybelokurov/ocen_dm>
