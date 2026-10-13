# Getting started

## Environment

The working environment is a Python virtual environment with AGAMA, JamPy
(non-commercial licence, imported, never vendored), UltraNest and the usual
scientific stack. From the repository root:

```bash
source <venv>/bin/activate
export PYTHONPATH=src
python -m ocen_dm.cli --help        # or the installed entry point `ocen`
python -m ocen_dm.cli inventory     # what data are present
```

`ocen` is declared in `pyproject.toml`; `ocen-dm` is not a command.

## Data

```bash
python -m ocen_dm.cli fetch-data    # registered downloads, checksums verified
python -m ocen_dm.cli preprocess    # standardised catalogues and published profiles
python -m ocen_dm.cli plot-data
```

- Two journal supplements need [manual retrieval](../MANUAL_DOWNLOADS.md).
- Our own HST and Gaia PM profiles and the count products have separate
  builders (`kinematics/hst_profile.py`, `kinematics/outer_gaia.py`,
  `bin/build_count_profiles.py`). Loaders build missing products
  automatically, so a fresh checkout can write data on first use.
- External catalogues live in a shared store (`$OCEN_CATALOGUES`) and are
  referenced, not copied; some model inputs come from a clone of
  [oCen_bar](https://github.com/adllmr/oCen_bar). See
  [Data products](data/index.md).

## Fitting

- Fits are launched only when requested, with a fresh output label; a run
  directory is never reused. Every batch freezes its source, configuration
  and inputs (with sha256), so it can be replayed exactly.
- The base model is fitted with `bin/run_compact_df_recovery.py`; the exact
  command is on the [base model](models/base_no_dm.md) page.
- The older Jeans fits use `python -m ocen_dm.cli fit` with explicit switches
  (see the [modelling plan](../MODELLING_PLAN.md)); a bare `fit` selects an
  outdated default.

## Tests

```bash
python -m pytest tests -q
```

The full suite includes local-data and optional-engine tests and takes tens
of minutes; some tests skip without those inputs.

## Conventions

- Units: internal dynamics in pc, M☉ and km/s; Galactic-orbit modules carry
  their own explicit conversions (AGAMA time unit = 977.8 Myr).
- Plots are PNG files in `plots/`; plotted tables and JSON provenance go to
  `results/plot_data/`. Compiled PDFs sit next to their LaTeX in `docs/`.
- `data/` and `results/` are not committed.
- Results are recorded in `JOURNAL.md` as they happen; corrections are new
  entries, not edits of old ones.
- External reviews go to `docs/reviews/`.

## Building this site

```bash
source <venv>/bin/activate
mkdocs serve            # local preview at http://127.0.0.1:8000
mkdocs build --strict   # the link check used before publishing
mkdocs gh-deploy        # publish to GitHub Pages
```

The build hook (`docs/hooks/kb.py`) serves `plots/`, splits the journal by
month, generates the [registered sources](data/sources.md) and
[literature](literature.md) pages from committed files, and points code links
to GitHub.
