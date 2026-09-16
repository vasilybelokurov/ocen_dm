# oCen_dm — dark matter in Omega Centauri from its tidal tails

Can the observed tidal tails of Omega Centauri break the mass-decomposition
degeneracy left by its internal stellar kinematics, and so constrain a surviving
dark-matter halo?

The full specification, experiment sequence and milestone list live in
[`OCEN_DM_TAILS_PROJECT.md`](OCEN_DM_TAILS_PROJECT.md). Running notes are in
[`JOURNAL.md`](JOURNAL.md).

Repository: https://github.com/vasilybelokurov/ocen_dm

**Current state: Milestones 1–2 done** — all four primary datasets ingested and
validated; the mass-profile library (`src/ocen_dm/mass_models/`) verified against
quadrature and against AGAMA. No scientific inference is implemented yet; the
kinematic fit (Milestone 3) is next. Overview figures of every dataset are in
[`plots/`](plots/) and are regenerated with `ocen plot-data`.

## Quick start

```bash
source ~/Work/venvs/.venv/bin/activate
cd "$HOME/Work/Code/oCen_dm"

export PYTHONPATH=src            # the package is used in place, not installed
python -m ocen_dm.cli inventory              # dataset inventory (no inference)
python -m ocen_dm.cli fetch-data             # download + checksum + manifest
python -m ocen_dm.cli preprocess             # standardize -> data/processed
python -m ocen_dm.cli inspect-omegacat --columns
python -m ocen_dm.cli plot-data              # PNG figures -> plots/
python -m pytest tests -q                    # 175 tests, no network needed
```

Two datasets are journal supplementary material and must be fetched by hand:
see [`docs/MANUAL_DOWNLOADS.md`](docs/MANUAL_DOWNLOADS.md).

## Layout

```text
configs/       data.yaml (which file holds which product), column_maps.yaml
provenance/    datasets.yaml (registry), manifest.json + checksums.txt (written)
data/raw/      downloads, never edited by hand, never committed
data/processed/ validated ECSV products
src/ocen_dm/   package: paths, provenance, data loaders, mass_models, plotting, CLI
plots/         PNG overview figures of the ingested data
tests/         pytest suite on synthetic fixtures only
docs/          manual retrieval instructions
```

## Design rules enforced in code

- **No invented metadata.** A checksum without a stated source is rejected by
  the registry loader. Column roles resolve only against documented candidate
  names; an unresolved or ambiguous role raises and prints the table's real
  columns instead of guessing (`src/ocen_dm/data/schema.py`).
- **Provenance for every file.** `fetch-data` verifies the published MD5,
  records SHA-256 and byte counts in `provenance/manifest.json`, and stores the
  interpreter and package versions used.
- **Nothing is faked on failure.** A failed download leaves no partial file; a
  manual dataset is reported as `manual`, not substituted.
- **Likelihood rules travel with the data.** Each processed table carries
  `meta['ocen_likelihood_rule']` recording the spec's usage restriction
  (membership probabilities are diagnostic only; spectroscopy is conditional on
  target positions; Fimbulthul is a track constraint, not a density).

## Environment

Python 3.13.5 in `~/Work/venvs/.venv`. Present: numpy 1.26.4, scipy 1.16.1,
astropy 7.1.0, pandas 2.3.2, pyarrow 21.0.0, h5py 3.14.0, matplotlib 3.10.3,
corner 2.2.3, pyyaml 6.0.2, emcee 3.1.6, agama 1.0.152, galpy 1.9.2, gala 1.9.1,
ultranest 4.5.0, dynesty 3.1.0, jampy 9.0.2.

**JamPy licence.** JamPy (Cappellari) is non-commercial and **may not be
redistributed**. It is a dependency installed into the environment, never copied
into this repository.

## Mass models

Internal units are pc, Msun and km/s (`G = 4.300917270036e-3`, asserted against
`astropy.constants`). Components with closed forms use them; the truncated halos
are tabulated once per parameter set on a 2000-point log grid and splined,
which costs ~0.5 ms per composite model (adaptive quadrature took 256 ms). The
quadrature paths remain as the independent reference that `verify()` checks
against.
