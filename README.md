# oCen_dm — dark matter in Omega Centauri

Can Omega Centauri's tidal debris distinguish an extended dark-matter halo from
the stars, remnants and central mass that fit its internal kinematics?

**Knowledge base: <https://vasilybelokurov.github.io/ocen_dm/>.** It describes
every data product (source, processing, use), the base model and its variants,
the current results, the literature and the journal. Its sources are in
`docs/`; preview locally with `mkdocs serve`.

**Literature:** public SciX library
[omega Cen dark matter (oCen_dm)](https://scixplorer.org/public-libraries/RSXq3-NtTTCApvfPUbdwrg).

## Status (2026-10-13)

- **Cluster mass model.** A two-transition action DF with Poisson star counts
  fits the 89 kinematic bins without dark matter (χ²_kin = 156) but fails its
  own data-fit gate. The data fix the total mass profile, not the split
  between stars and dark mass. A count-likelihood bias was fixed on 2026-10-13
  and the base model refitted; the DM scans are not yet rerun. [Base model](docs/kb/models/base_no_dm.md).
- **DM survival and tails.** N-body tails of no-DM and DM models agree to 10%;
  the differences are not testable with published tail data.
- **Stream 54.** A rotating bar and the LMC bend the stream part-way; the
  proper-motion knee is not reproduced. [Current note](docs/STREAM_PM_REFINEMENT_2026-10-13.md).

Details: [results](docs/kb/results.md) and the chronological [journal](JOURNAL.md).

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
