# Code and analysis audit — 21 September 2026

## Scope and status

**Follow-up:** the reporting, mock-handling and run-preservation repairs described at
the end of this document are now implemented. The audit findings and verification
below record the state before that repair pass. Rung 0 completed its intended purpose:
establish a baseline and identify the residuals that guide the next research step.
The planned model extensions are separate from the software repairs.

This audit starts from commit `6e41a86` and compares the current CLI, data loaders,
light/mass models, Jeans likelihood, fit/report machinery and DM-survival calculations
with saved products and the analysis notes. The WD literature compilation and bibliography
were not re-audited in full; several bibliography entries remain incomplete. It updates
documentation and code comments;
it does not change numerical algorithms, rebuild data, resample posteriors or integrate
new orbital histories. Dated reviews and earlier journal entries remain historical records.

The three 89-point isotropic fits provide the intended baseline for staged modelling.
Several explanations and operational claims needed correction.
The most consequential ones are below, with confirmed software gaps kept separate from
scientific hypotheses.

## What the current pipeline actually does

| Stage | Implemented behaviour | Important boundary |
|---|---|---|
| Data | Registry, retrieval, checksums, schema validation; published and measured profiles | `preprocess` does not build every later analysis product |
| Tracer | HST counts inside 25 arcsec joined to Trager light, fitted by spherical MGE | Same shape is used for tracer and stellar mass; tracer construction is outside fit input hashes |
| Default observations | HST R/T: 21 + 21; MUSE: 29; Gaia EDR3 R/T: 9 + 9 | Pristine, DR2 and tail spectroscopy are outside the current 89-point fit |
| Family K1 | Stars + remnant Plummer + point mass + distance; selectable anisotropy | Baseline isotropy is not a positive-DF guarantee |
| Family K2 | K1 plus gNFW γ = 0 or 1; log-uniform `M_dm_100`, `r_s`; fixed `r_t=1000 pc` | Cored means γ = 0 gNFW, not Burkert; taper is not an orbit-derived tidal radius |
| Likelihood | Star-radius quantiles for our PMs; annular averaging for MUSE; streaming correction; split normal | Diagonal likelihood; coherent calibration and rotation errors are not marginalised |
| Engines | Local spherical Jeans; JamPy spherical comparison; AGAMA QuasiSpherical diagnostic | No axisymmetric fit or automatic positive-DF filter |
| Inference | UltraNest, saved posterior, summary, profile subset and YAML metadata | No dynesty driver; no complete immutable run snapshot |
| Tails | Orbits, friction/bar experiments, approximate tracks and analytic truncation | No particle-spray likelihood, tail posterior reweighting or production disruption run |

The installed entry point is `ocen`, not `ocen-dm`. A bare `fit` still selects the
legacy varying-β family with three scale parameters, including an unused Gaia DR2 scale.
The adopted baseline requires `--isotropic --no-scales`. The Python model constructor
also retains `tracer="trager"`, whereas the CLI selects `composite`.

## Numerical checks performed

### Saved fits reproduce at their best samples

For each run, reconstruct its family from the recorded ladder switches, read the saved
maximum-likelihood vector and load the current default profiles. No optimiser or sampler
is needed.

| Run | Stored and recomputed ln L | Recomputed χ² | ln Z from saved run |
|---|---:|---:|---:|
| `rung0_K1` | −169.65805562516505 | 771.5338604 | −190.0162626 ± 0.4767496 |
| `rung0_K2_cored` | −164.4418542662639 | 761.1014577 | −188.9384154 ± 0.3591258 |
| `rung0_K2_nfw` | −165.45250672809564 | 763.1227626 | −190.0937139 ± 0.3523136 |

The log likelihoods agree exactly in this environment; all per-dataset χ² values match
the summaries to their stored two decimal places. From the unrounded evidences,
Δln Z is +1.08 ± 0.60 (cored minus K1), −0.08 ± 0.59 (NFW minus K1), and
+1.16 ± 0.50 (cored minus NFW), combining sampler errors in quadrature.

The Gaia input hash differs between the archived runs. This reproduction supports the
journal's metadata-only explanation, but evaluating one point per run cannot prove
identity of all original likelihood inputs. The arrays were not archived separately;
the generated report's warning remains intact.

The K2-cored radius-split residuals and the error-floor experiment also reproduce:
χ² = 761.10, 428.05, 270.71 and 144.42 for added fractional floors of 0, 0.5%, 1%
and 2%, at the same parameter vector. This is a sensitivity calculation, not a refit.

A minimal reproduction, after activating the existing environment:

```bash
PYTHONPATH=src python - <<'PY'
import numpy as np
from ocen_dm.kinematics.fit import FitProblem
from ocen_dm.kinematics.likelihood import KinematicData
from ocen_dm.kinematics.report import load_run, _family_for
for label in ('rung0_K1', 'rung0_K2_cored', 'rung0_K2_nfw'):
    r = load_run(label)
    s = r['summary']
    f = _family_for(s, options=r['run']['dataset_options'])
    # Valid here because these runs all use raw errors and published rotation.
    d = KinematicData.load(list(s['datasets']))
    p = FitProblem(f, d)
    x = np.array([s['parameters'][n]['ml'] for n in f.names])
    print(label, s['lnL_max'], p.loglike_vector(x), p.chi2(x))
PY
```

The recipe requires the original local inputs. Some loaders build missing products, so
check their presence before treating it as a read-only operation on another checkout.

### DM-survival arithmetic

Using the existing `shock_heating` function and its stated passage counts:

| Local radius | Class 1: 114 × ΔE/|E| per passage | Class 3: 66 × ΔE/|E| per passage |
|---|---:|---:|
| 20 pc | 0.00122 | 0.11228 |
| 35 pc | 0.06698 | 1.52175 |
| 50 pc | 0.67188 | 5.45810 |
| 70 pc | 4.22219 | 16.15328 |

The former claim of order-unity class-3 heating **at 20 pc in 66 passages** was an
arithmetic error. The estimate requires about 588 identical passages there. It suggests
outer erosion but does not establish a shock-limited final 20-pc density. These local,
fixed-structure sums cease to be reliable physical predictions once the structure changes.

The energy-cut factor multiplies the **initial** density. The plotted density already
contains `exp(-r/r_J)`, so replacing that taper requires `F_gamma * exp(r/r_J)` times
the plotted value, not an additional factor `F_gamma` alone.

For the source grid `M200 = 10^9, 10^10, 10^11 Msun`, `c = 5, 10, 15`, `z = 2`:

| At r = 20 pc | r_J = 35 pc | r_J = 70 pc |
|---|---:|---:|
| Initial NFW density range [Msun/pc³] | 1.019–22.390 | same |
| Exponentially tapered range | 0.575–12.644 | 0.766–16.826 |
| Kepler γ = 1 retention factor | 0.139326 | 0.353387 |
| Energy-cut estimate from the initial density | 0.142–3.120 | 0.360–7.912 |
| Multiplier replacing the existing exponential | 0.246718 | 0.470257 |

Independent quadrature reproduces the retention factors to better than 10⁻¹⁴. Applying
a pure-power-law Kepler factor to the local NFW density is still an approximation;
none of these ranges is an observational constraint or a self-consistent evolved halo.
The `m_dm_in_rj` table field integrates the **untruncated** NFW profile, so it is not the
mass integral of the plotted density. The plotted stellar comparison borrows the K1 mass
but uses a fixed distance and a different MGE width grid.

### Focused test suite

**104 passed in 9.37 s**, using the existing Python environment:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src MPLCONFIGDIR=/private/tmp/ocen-doc-audit/mpl \
python -m pytest -q -p no:cacheprovider \
  tests/test_likelihood_normalisation.py tests/test_gaia_options.py \
  tests/test_modelling_ladder.py tests/test_truncated_equilibrium.py \
  tests/test_master_plot.py::test_point_count_matches_the_likelihood \
  tests/test_mass_models.py tests/test_jeans.py
```

The old full-suite log records 392 passes and one stale 98-versus-89-point assertion.
The current point-count test expects 89 and passed here. This audit did not rerun the
full suite or the long orbital grids. Passing tests do not validate the software paths
identified below; the narrow reproductions demonstrate gaps in their coverage.

### Documentation and preservation checks

All four companion PDFs were rebuilt and their rendered pages inspected. The final
LaTeX logs contain no overfull/underfull boxes, undefined references or errors.
The documented fit recipes parse with the current CLI. Local links in the main
Markdown guides resolve, and `git diff --check` passes.

Removing comments and docstrings leaves identical Python syntax trees in all eleven
modified source/test files. SHA-256 checks confirm that all 39 snapshotted processed-data
and rung-0 result files are unchanged. These checks establish the scope of this edit;
they do not supply the missing historical run snapshots discussed below.

## Software gaps identified at the audit

These findings describe code at this audit's starting commit. The initial documentation
pass did not change their behaviour; see the repair status below for subsequent changes.

| Issue | Evidence and consequence | Required repair |
|---|---|---|
| Report replays default Gaia options | `report.data_for` ignores recorded `gaia_errors` and `gaia_rotation`. A synthetic summary with eta/ours still reloaded raw/published arrays exactly; the correct eta values differ by up to 0.003911 mas/yr and the streaming variance by 0.026520 (mas/yr)². | Reconstruct data using stored options; compare arrays and likelihood against the run. Default rung-0 replay is unaffected. |
| Mock fingerprint uses the wrong schema | `comparison_table` reads `mock_from`, `family`, `parameters`; the writer stores `source_file`, `generating_family`, `truth`. Two otherwise identical mock records with different generating families/truth and the same seed received a numerical Δln Z without a mismatch warning. | Fingerprint the actual generating model, full truth and realisation, preferably the fitted arrays. |
| CLI mock generator omits ladder switches | `cmd_fit` builds the generator with only family and tracer. A rung-0 K1 vector has five entries; the default generator expects eleven. | Reconstruct the generating run rather than infer it from a bare `.npy` vector. |
| Run snapshot is incomplete | `run_nested` hashes all kinematic ECSV files, including unused ones; omits the raw HST/Trager tracer inputs; reads git commit at completion; omits `step_sampler`, dirty-tree state, some custom family settings and the legacy Gaia cut. | Snapshot actual inputs, resolved family/settings and code state before computation. Capture enough information to reproduce custom and preset runs. |
| Family replay is partial | `_family_for` recognises ladder switches, but does not restore every preset, fixed parameter, distance prior, custom halo taper or saved engine. | Serialize and replay the resolved model configuration. |
| Mock noise differs from fitted likelihood | `mock_data` uses symmetric Gaussian errors from the mean of each error pair and clips low values. | Define what the recovery test should validate; use an appropriate generative model for asymmetric data. |
| Legacy defaults and automatic writes | Bare `fit` differs from the adopted sequence; `resume="overwrite"`; missing measured products are built by their loaders. | Keep explicit recipes now; make any future default/overwrite policy change deliberately. |

Some legacy metadata text is also stale, notably the combined Gaia loader's midpoint
note. Current component loaders and product columns define the adopted likelihood.
Fixing metadata during a future product rebuild should be recorded rather than silently
changing archived run fingerprints.

## Scientific interpretation that needed tightening

- **Anisotropy:** projected PM ratios motivate a radial trend, but do not directly measure
  intrinsic β(r). Matching the geometric mean of the two components does not establish
  enclosed mass to the same fractional precision. A turnover is a proposed test, not a
  result recovered by the present monotonic parametrisation.
- **Central mass:** a good central moment fit does not establish a robust IMBH measurement.
  For a cored tracer in a central point-mass potential, the necessary condition
  `beta(0) <= -1/2` also applies to constant β. Rungs 0/1 deliberately allow diagnostic
  models outside that condition. See [An & Evans](https://arxiv.org/abs/astro-ph/0511686).
- **Engine agreement:** the JamPy comparison shares the intended equations but approximates
  non-Gaussian masses with MGEs. Opposing HST χ² shifts can cancel in total ln L. Their
  cause has not been isolated; agreement in the total is not an exact identity test.
- **Tidal boundary:** the circular spherical reference is
  `r_J³ = G M_sat / (Omega² - Phi'')`, hence a denominator coefficient
  `3 - d ln M_host/d ln R`, not `2 - d ln M_host/d ln R`.
  [galpy documents this circular-orbit definition](https://docs.galpy.org/en/v1.5.0/reference/potentialrtide.html).
  Eccentric barred histories need a stated instantaneous approximation, not a claim of
  a conserved Jacobi integral. The historical alternative radius of 59 pc is unverified.
- **Energy versus apocentre:** `E < Phi(r_J)` is stricter than an apocentre cut at nonzero
  angular momentum. The analytic retention formula cannot be described as removing
  exactly and only those orbits that reach the boundary.
- **Mass accounting:** the original N2 plan added the cluster MGE to a Plummer nucleus
  already representing that cluster. These are alternative baryonic models. Seven pc is
  the assumed half-mass radius, not a radius containing the entire 3.55-million-Msun mass.
- **Self-gravity:** the high-DM grid is not everywhere subdominant to the nucleus. Test
  enclosed-mass ratios before adopting massless particles or a small-adjustment estimate.
- **Orbital experiments:** constant-mass friction scans are not universal upper bounds on
  arbitrary mass-loss histories. Class-3 contour overlap is strongly pattern-speed-dependent;
  the 10⁸-Msun, Ω = 24 trial has overlap 0.34, contrary to the old figure caption's zero.
  The first two selected class-2 histories end above the present nucleus mass.
- **Tidal tracks:** extrapolations outside calibration and initial aperture budgets are
  not validated remnant masses. Reducing the problem to the inner 100 pc remains a proposed
  approximation that must retain or account for centre-crossing outer orbits.

## Next work

1. Reporting, provenance and rung-specific mock generation: repaired as described below.
2. When requested, run the constant-β controls with fixed observations and no scales;
   specify and validate the turnover extension before fitting it.
3. For survival, establish the tidal-radius convention and a non-double-counted baryonic
   model, then solve the truncated equilibrium and test tidal heating (N0/N2/N3).
4. Keep numerical convergence criteria separate from expected scientific outcomes.
   A converged result that contradicts the analytic estimate is a finding, not a failed test.

## Repair pass — 21 September 2026

The user approved a bounded repair of reporting, mock handling and run preservation.
They also clarified the purpose of the research sequence: rung 0 establishes a useful
baseline. Its residuals motivate the next step; provisional scientific scope is not a
software defect. The README, modelling plan and current analysis now lead with that purpose.

Implemented in `kinematics/run_io.py`, `kinematics/report.py`, `kinematics/fit.py` and the CLI:

- New runs save the exact likelihood arrays and resolved model before sampler construction.
  The model includes MGE coefficients, parameter order and priors, fixed values, distance,
  halo taper, instruments and backend. Launch metadata includes source hashes, git commit
  and dirty state, available package versions, Gaia options and sampler settings.
- Reports verify snapshot checksums and reproduce the saved maximum log likelihood before
  plotting. Legacy runs load recorded Gaia options, warn about using current products and
  undergo the same likelihood check. Missing new snapshots never trigger a silent reload.
- New evidence comparisons identify observations by their saved arrays. Legacy comparisons
  use the actual mock schema, including generator, full truth and seed, and withhold a
  numerical difference when provenance is incomplete. Model-stage switches are separated
  from observational choices.
- Injection from a run directory or its parameter file restores and checks the source model
  independently of the fitted model. A bare vector now needs an explicit `--mock-family`;
  the tracer, backend and ladder flags then define its generator. Generated arrays and the
  full generating model are saved.
- Both CLI and driver refuse an occupied run directory. Launch snapshots survive an
  interrupted sampler. Reusing or resuming such a directory is deliberately not implemented.

The historical mock-noise prescription (symmetric Gaussian using averaged errors, with
a positive floor) remains unchanged and is explicitly named in new provenance. Numerical
solvers, priors, modelling stages and production results were not changed. Missing-product
builders retain their previous behaviour. Snapshots preserve resolved fit inputs and code
identity, not a complete software installation or raw-catalogue reduction archive. Old
custom configurations cannot be recovered when the necessary metadata was never recorded;
inconsistent replay raises an error instead of substituting defaults.

Validation: **78 focused tests passed**, including **36 new regression cases** and short
CLI sampler/report tests. Two assertions were then strengthened (distinct new snapshot
comparisons and preset replay), and all 40 replay/CLI tests passed again. Coverage includes
all three source rungs for K1/K2 mocks, fixed-distance presets, custom halo settings,
backend configuration, changed observations, broken snapshots and interrupted runs.
The full suite was not rerun; a broader exploratory plotting run was stopped after a
fixed intermediate report error and is not counted as a passing suite.

All three archived rung-0 maximum log likelihoods still reproduce exactly; their original
hash-mismatch warning is retained. All 39 snapshotted processed-data and fit files are
unchanged. The two updated PDF notes were rebuilt and inspected. No production fit was launched.
