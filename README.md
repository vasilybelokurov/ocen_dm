# oCen_dm — dark matter in Omega Centauri

Can Omega Centauri's tidal debris distinguish an extended dark-matter halo from
the stars, remnants and central mass that fit its internal kinematics?

**Knowledge base: <https://vasilybelokurov.github.io/ocen_dm/>.** It describes
every data product (source, processing, use), the base model and its variants,
the current results and the literature. Its sources are in
`docs/`; preview locally with `mkdocs serve`.

**Literature:** public SciX library
[omega Cen dark matter (oCen_dm)](https://scixplorer.org/public-libraries/RSXq3-NtTTCApvfPUbdwrg).

## Status (2026-10-17)

- **Cluster mass model.** A two-transition action DF with Poisson star counts
  fits the 89 kinematic bins without dark matter (χ²_kin = 156) but fails its
  own data-fit gate. The data fix the total mass profile, not the split
  between stars and dark mass. [Base model](docs/kb/models/base_no_dm.md).
  - Fixed central black holes of 8,200 and 4×10⁴ M☉ were fitted; neither
    passes the gate and the mass is not measured. [DF + black hole](docs/kb/models/df_bh.md).
  - A refit with completeness-corrected counts and a tied HST/Gaia amplitude
    ratio changes the predictions by a few per cent (χ²_kin = 159).
  - Dark matter: the DF with a free cored halo and tied counts drives the halo
    to zero; Jeans fits with NFW or cored halos show no preference for dark
    matter with either tracer. [DM fits](docs/kb/models/dm_fits.md).
  - The ρ₂₀ profile, mass-per-star and mock scans still use the old
    (schema-1) count products.
- **DM survival and tails.** N-body tails of no-DM and DM models agree to 10%;
  the differences are not testable with published tail data.
- **Stream 54.** A rotating bar is needed: barred models beat the best
  axisymmetric model by about 3,300–4,100 in ln L, but none reproduces the
  steep proper motions of the turn. The LMC shifts the distal proper motions
  by a few tenths of a mas/yr but is absorbed by the cluster's PM. Fits to the
  binned tracks disagree with the star-by-star likelihood. That likelihood did
  not converge with particle number; it now uses fixed χ bins with a jackknife
  bias correction, and a 9-D emulator training set is running.
  [Stream status](docs/kb/stream_status.md),
  [likelihood versions](docs/kb/models/stream_likelihood_versions.md).

Details: [results](docs/kb/results.md).

## Quick start

```bash
source <venv>/bin/activate
export PYTHONPATH=src
python -m ocen_dm.cli --help        # installed entry point: `ocen`
python -m pytest tests -q           # full suite takes tens of minutes
mkdocs serve                        # knowledge base at http://127.0.0.1:8000
```

Fits are launched only on request, each with a fresh output label; every batch
freezes its source, configuration and inputs. See
[getting started](docs/kb/getting_started.md).

## Layout

```text
configs/                 dataset locations, column mappings, potentials
provenance/              source registry, download manifest and checksums
data/raw/, data/processed/  downloads and products (not committed)
src/ocen_dm/data/        retrieval, schema validation and catalogue loaders
src/ocen_dm/selection/   covariance joins, field templates and crossmatches
src/ocen_dm/kinematics/  estimators, Jeans and DF engines, likelihoods, fits
src/ocen_dm/mass_models/ stellar, remnant, point-mass and halo components
src/ocen_dm/tails/       orbital histories and survival models
src/ocen_dm/dynamics/    remnant/progenitor initial conditions and evolution
src/ocen_dm/streams/     restricted N-body, particle spray, barred/LMC hosts, stream likelihoods
bin/, bin/streams/       drivers, experiments and plotting scripts
results/                 fits and diagnostics (not committed)
plots/                   PNG figures
docs/                    knowledge base, reports, LaTeX write-ups, reviews
tests/                   synthetic tests plus tests that need local inputs
```

Internal dynamics use pc, M☉ and km/s. JamPy is an external non-commercial
dependency and is not redistributed here.

Repository: <https://github.com/vasilybelokurov/ocen_dm>
