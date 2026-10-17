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

## External locations

Inputs that are not in the repository are found through environment
variables; nothing in the code assumes a machine-specific path.

| Variable | Contents | Needed by |
|---|---|---|
| `OCEN_CATALOGUES` | shared catalogue store: Vasiliev & Baumgardt 2021 members and cluster zip, Baumgardt GC tables, STREAMFINDER, HST Large Programme, Tango archive | provenance checks, Gaia products, all stream scripts |
| `OCEN_ISOCHRONES` | MIST grids prepared with minimint (`minimint_hst/` with the HST filters) | mass function, HST CMD plots |
| `OCEN_BAR_DIR` | clone of [oCen_bar](https://github.com/adllmr/oCen_bar) (Hunter+2024 and Portail+2017 potentials, GSE contours) | barred hosts, bar migration, stream fits |
| `NEMO_DIR` | NEMO/falcON build | `bin/nemo_env.sh` (gyrfalcON runs) |
| `PYTHON` | interpreter for the shell drivers (default `python`) | `bin/*.sh` |

Put them in a file named `.env` in the repository root. It is gitignored, so the locations
never enter the repository, and `import ocen_dm` reads it, so scripts, tests and fits find the
inputs without any shell setup (`bin/nemo_env.sh` reads `NEMO_DIR` from it too):

```bash
# .env -- one NAME=value per line; ~ and $VAR are expanded
OCEN_CATALOGUES=/path/to/catalogues
OCEN_ISOCHRONES=/path/to/isochrones
OCEN_BAR_DIR=/path/to/oCen_bar
NEMO_DIR=/path/to/nemo
```

A variable exported in the shell always wins over the file. An unset variable shows up
literally in the error message, for example
`FileNotFoundError: ... '$OCEN_CATALOGUES/streamfinder_ibata2024_dr3.fits'`.

## Data

```bash
python -m ocen_dm.cli fetch-data    # registered downloads, checksums verified
python -m ocen_dm.cli preprocess    # standardised catalogues and published profiles
python -m ocen_dm.cli plot-data
```

- One journal supplement (Kuzma et al. 2026) needs [manual retrieval](../MANUAL_DOWNLOADS.md);
  Fimbulthul is now fetched from VizieR.
- `data/raw/omegacat_vi_kinematics` may be a local link to an external copy of
  the oMEGACat VI data rather than a directory in the checkout.
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

The full suite includes local-data and optional-engine tests and takes about
13-25 min; some tests skip without those inputs. `pytest-xdist` (in the `dev`
extra of `pyproject.toml`) runs it in parallel:

```bash
python -m pytest tests -q -n 6
```

## Conventions

- Units: internal dynamics in pc, M☉ and km/s; Galactic-orbit modules carry
  their own explicit conversions (AGAMA time unit = 977.8 Myr).
- Plots are PNG files in `plots/`; plotted tables and JSON provenance go to
  `results/plot_data/`. Compiled PDFs sit next to their LaTeX in `docs/`.
- `data/` and `results/` are not committed.
- Results are recorded in a project journal that is kept locally (not in
  git); the reports in `docs/` carry the results that are shared.
- External reviews go to `docs/reviews/`.

## Building this site

```bash
source <venv>/bin/activate
mkdocs serve            # local preview at http://127.0.0.1:8000
mkdocs build --strict   # the link check used before publishing
mkdocs gh-deploy        # publish to GitHub Pages
```

The build hook (`docs/hooks/kb.py`) serves `plots/`, generates the [registered sources](data/sources.md) and
[literature](literature.md) pages from committed files, and points code links
to GitHub.
