# oCen_dm Journal

Running log of what was done, what was verified, and what remains open.
Newest entries at the bottom of each day. Dates are absolute.

---

## 2026-09-16 — Milestone 1: repository, provenance, data ingestion

### Environment audit

- Interpreter: `~/Work/venvs/.venv`, Python 3.13.5 (macOS 15.7.9, arm64).
- Available: numpy 1.26.4, scipy 1.16.1, astropy 7.1.0, pandas 2.3.2,
  pyarrow 21.0.0, h5py 3.14.0, matplotlib 3.10.3, corner 2.2.3, pyyaml 6.0.2,
  emcee 3.1.6, astroquery 0.4.10, requests 2.32.5, **agama 1.0.152**,
  **galpy 1.9.2**, **gala 1.9.1**.
- Missing, needed later: `jampy` (Milestone 3), `ultranest`, `dynesty`
  (Milestone 4), `pooch`, `polars`, `snakemake`. Nothing was installed; the
  Milestone 1 code deliberately uses only the standard library plus astropy and
  pyyaml, so no install is required yet.
- `/Users/vasilybelokurov/Work/Code` is a symlink to the Dropbox `Code`
  directory, so the two declared working directories are the same tree.
- Neighbouring project `../satellite_experiment/` uses AGAMA for host
  potentials/orbits and NEMO `gyrfalcON` for live N-body, driven by numbered
  `scripts/` stages over reusable `src/` modules. That is the pattern to follow
  for Milestone 11.

### Reference verification (spec section 9 policy)

Checked live on 2026-09-16:

- arXiv abstract pages HTTP 200, titles match the spec: 2502.01135
  (Kuzma & Ishigaki 2025), 2605.23474 (Kuzma et al. 2026), 2503.04903
  (oMEGACat VI), 2404.03722 (oMEGACat II), 2307.03035 (Youakim et al. 2023).
- Crossref resolves 10.1038/s41550-019-0751-x (Ibata et al. 2019, Nature
  Astronomy), 10.1093/mnras/stag1147, 10.3847/1538-4357/adbe67,
  10.1093/mnras/stad1952 — all titles match.
- DataCite resolves the three Zenodo DOIs 10.5281/zenodo.14791430 (2025),
  14978551 (2025), 11104046 (2024) — titles match.
- Kuzma et al. (2021) MNRAS URL is recorded **UNVERIFIED**: OUP blocks
  automated requests. It is not needed for any data product.

### Blocker: Zenodo outage

`zenodo.org` resolves and answers ICMP but returned **HTTP 504** for both
`/api/records/<id>` and the file endpoints throughout the session (checked
~17:00–17:15 UTC). Consequently:

- no Zenodo file was downloaded;
- the oMEGACat VI file list and all real column names are **unknown**;
- the published MD5 for `wCen_table.fits` (`3dcd58e6…c768`) is recorded from
  the project specification and marked as such, not confirmed against Zenodo.

Nothing was invented to work around this. `ocen fetch-data` is re-runnable and
will discover, download and verify the Zenodo files when the service returns.

### Implemented

- Repository scaffold per spec section 5 (`configs/`, `provenance/`, `data/`,
  `src/ocen_dm/`, `tests/`, `docs/`, `results/`, `notebooks/`, `workflows/`),
  `pyproject.toml`, `.gitignore` (raw/interim/processed/results untracked).
- `provenance/datasets.yaml`: 7 datasets (4 primary + 3 optional/reference),
  each with role, paper links, verification status and verification note.
- `src/ocen_dm/provenance.py`: registry parsing that **refuses a checksum
  without a `checksum_source`**; MD5/SHA-256 digests; `verify_file`;
  `manifest.json` + `checksums.txt` writers recording interpreter and package
  versions.
- `src/ocen_dm/data/download.py`: atomic downloads (`.part` then rename),
  exponential-backoff retries, cached-copy reuse, `--force`, Zenodo file
  discovery, manual-dataset reporting. A failure never leaves a partial file
  and never substitutes a placeholder.
- `src/ocen_dm/data/schema.py`: `ColumnRole`/`TableSchema`, role resolution
  against documented candidate names (Gaia archive data model), explicit
  `configs/column_maps.yaml` overrides, and validation of units, dtypes,
  non-finite counts and unique IDs.
- Loaders `kuzma2025.py`, `kuzma2026.py`, `fimbulthul.py`, `omegacat.py`. Each
  attaches the spec's likelihood restriction to the product metadata. The
  oMEGACat loader refuses to load a profile until `configs/data.yaml` names the
  file, because the real file/column names are not yet known.
- `src/ocen_dm/cli.py`: `fetch-data`, `preprocess`, `inventory`,
  `inspect-omegacat`.
- `docs/MANUAL_DOWNLOADS.md` for the two publisher-supplementary datasets.

### Tests (all run, all passing)

`python -m pytest tests -q` → **55 passed** in 0.35 s, no network access.

Coverage of the Milestone 1 acceptance criteria:

- raw data gitignored — `test_raw_data_is_not_committed`;
- checksums/provenance saved — manifest/checksum writer tests, corrupted- and
  missing-file detection;
- no invented metadata — registry rejects an unsourced checksum; unresolved and
  ambiguous column roles raise with the real column list;
- processed products have explicit schemas — unit, dtype, duplicate-ID and
  non-finite validation, ECSV round-trip preserving units and metadata;
- loaders have tests — synthetic FITS/CSV/ECSV fixtures for all four datasets;
- downloader — atomicity, retry, give-up-cleanly, cache reuse, `--force`,
  checksum-mismatch detection, Zenodo 504 surfaced rather than swallowed;
- CLI inventory runs with no data and produces no science output.

Two defects were found by the tests and fixed: candidate-name case variants
(`source_id` / `SOURCE_ID`) were counted as an ambiguous double match, and the
column-role deduplication needed to key on the resolved source column.

### Live command run

`ocen fetch-data --timeout 30 --retries 2` behaved as designed under the
outage: `wCen_table.fits` reported `failed` with the timeout message, the two
publisher datasets reported `manual`, oMEGACat VI discovery reported the 504,
the manifest was written with an empty file list, and no partial or placeholder
file was left on disk.

### Open items

1. Re-run `ocen fetch-data` when Zenodo is back; then
   `ocen inspect-omegacat --columns` and fill `configs/data.yaml` and
   `configs/column_maps.yaml` with the **real** file and column names.
2. Retrieve the Kuzma et al. (2026) and Ibata et al. (2019) supplementary
   tables by hand (`docs/MANUAL_DOWNLOADS.md`) and record the download date and
   exact file names here.
3. Confirm the `wCen_table.fits` MD5 against Zenodo's own published checksum
   once the API answers.
4. Then Milestone 2 (mass-profile library): stellar MGE, compact remnants,
   IMBH point mass, cored and cuspy truncated DM, composite enclosed
   mass/potential/force/escape speed, with independent numerical verification
   of every enclosed-mass curve.

---

## 2026-09-16 (later) — Codex second opinion + independent verification

Ran the `ask-codex` consultation (Codex CLI, read-only sandbox, `codex exec`
directly because this project is not a git repository so the skill's diff-based
script does not apply). Sent: the question, and whatever project files Codex
opened (spec, `src/`, `configs/`, `provenance/`, `tests/`). `JOURNAL.md` was
excluded. Review saved outside the repo in the session scratchpad.

**Codex's verdict: do not sign off Milestone 1.** Twelve findings, of which I
reproduced every one I tested. Codex ran its own read-only probes and found the
same three criticals I found independently.

### Confirmed defects (my own probes, not taken on trust)

| # | Defect | Evidence |
|---|---|---|
| 1 | A schema-INVALID table is still written to `data/processed/` and the CLI prints success; the error is buried in `preprocess_report.json` | RA in km → product written with `ra` in km, exit status hides it |
| 2 | Zenodo discovery discards the publisher's MD5 and size; `FileRecord` gets `published_md5=None` | discovered file with `md5:aaa…` → `published_md5=None` |
| 2b | `verify_file` returns `ok=True` when there is nothing to check — including a 0-byte file | `bytes: 0, ok: True` |
| 3 | Misspelled `fetch-data --dataset`, failed discovery and missing required manual files all exit 0 | `--dataset omegacat_typo` → exit 0 |
| 4 | `size_ok=False` still reports `status='ok'`; a manual file failing its MD5 reports `cached`, `.ok=True` | reproduced both |
| 5 | Units are **asserted, not converted**: a radius really in arcmin, unitless in CSV, silently becomes arcsec with unchanged values | `[1,2] arcmin → [1,2] arcsec` |
| 6 | 64-bit identifiers read as float lose precision and still pass the uniqueness check | `6000000000000000001 → 6000000000000000000`, `unique_id_ok=True` |
| 9 | Truncated download published: `Content-Length` is read but never enforced, and the `.part` file is renamed **before** checksum verification | 14 bytes delivered against a claimed 100000 → published |
| 10 | Masked values are not counted as non-finite (`np.asarray` drops the mask) | masked RA → `n_nonfinite = 0` |
| 11 | A misspelled override key is silently ignored; two roles may map to the same column | `{'raa': …}` ignored; `dec ← ra` accepted |
| 12 | A per-dataset fetch overwrites the whole manifest, erasing other datasets' provenance | fetch B → manifest holds only B |

Findings 7 (processed products carry no raw-file hash, so lineage is broken)
and 8 (the profile schema has no bin boundaries and only one symmetric error)
are confirmed by inspection; 8 is now proven by the published tables, below.

The unifying theme is the failure mode I said I cared about most: **several
paths report success while accepting unverified data or writing an invalid
scientific product.** `resolve_columns` fails loudly, but almost everything
downstream of it does not.

### Independent literature/data work done while Codex ran

- **Fimbulthul is in VizieR** as `J/other/NatAs/3.667/tables1` — Codex could not
  confirm any VizieR deposit; TAP on `METAcat` found it. Downloaded 309 rows
  (sha256 `9d06aca4…c705f`, 71578 B) and recorded in the registry as a
  CDS-curated *representation*, not the publisher's bytes. This removes one of
  the two manual downloads. Real columns: `RAJ2000, DEJ2000, pmRA, e_pmRA,
  pmDE, e_pmDE, plx, e_plx, G0mag, __BP-RP_0, SimbadName, recno` — note there
  is **no Gaia `source_id` column**; 306/309 carry a `Gaia DR3 …` SimbadName and
  3 carry CRTS/UCAC4 names. My `fimbulthul` schema, which requires `source_id`,
  is therefore wrong and would fail on the real data.
- **Kuzma & Ishigaki 2025 is NOT in VizieR**: `J/MNRAS/537/2752` does not exist
  (7 catalogues from MNRAS vol. 537, none of them this one). Nor is oMEGACat VI.
  Kuzma et al. 2026 (MNRAS 550) is not there either.
- **oMEGACat VI file names verified from the paper itself** (arXiv HTML
  2503.04903v2, Appendix A), independently of Zenodo:
  `proper_motion_dispersion_{log,lin,equaln,voronoi}_bins.fits`,
  `los_profile.fits`, `los_dispersion_and_rotation_voronoi_bins.fits`,
  `energy_equipartition_{profiles,massbins}.fits`, `catalog_and_selections.fits`.
  Codex quoted MD5s for some of these from an indexed copy of the Zenodo page;
  those are **not recorded** anywhere in this repository, because I cannot check
  them while Zenodo is down.
- **The published profile tables are obtainable without Zenodo**: the arXiv
  source package `https://arxiv.org/src/2503.04903v2` (5.58 MB) contains
  `revision1_dispersion_profile_table_pm.tex` (Table 2) and
  `dispersion_profile_table_los.tex` (Table 3), with columns
  `r_lower, r_median, r_upper, N_stars, sigma_PM,c/rad/tan` in both mas/yr and
  km/s plus an anisotropy ratio, and `v_LOS, sigma_LOS, theta_0` — all with
  **asymmetric** errors. This is a usable v1 kinematics dataset with verifiable
  provenance, and it proves finding 8: a schema with one symmetric error column
  and no bin boundaries cannot represent these data.
- Refuted Codex's guessed IOP machine-readable-table URL (`…apjadbe67t2_mrt.txt`
  → HTTP 404). Codex had already flagged it as unconfirmed.
- Zenodo still 504 at 18:40 UTC. Internet Archive was itself offline
  ("Temporarily Offline"), so no Wayback corroboration was possible. DataCite
  metadata carries no file list or checksum, so it cannot confirm the published
  MD5; it does reveal the concept DOI `10.5281/zenodo.14791429`.

### Not yet fixed

No code was changed in response to the review. The fix list is in the chat
summary and awaits a decision on scope/ordering.

---

## 2026-09-16 (later still) — all 12 findings fixed, plus 2 found by self-audit

All twelve confirmed defects were fixed, each with a regression test that fails
against the pre-fix code. Test suite: **104 passed**, no network needed.

### What changed

| # | Fix | Where |
|---|---|---|
| 1 | `schema_report` no longer swallows `SchemaError`; the new `_product.build_product` validates **before** writing, so an invalid table is never published and a previously valid product is never replaced by a bad one | `data/_product.py`, `data/base.py` |
| 2 | Zenodo discovery parses `'<alg>:<hex>'` and carries the publisher checksum, size and required-ness into the `FileRecord`; `verify_file` reports a verification level and rejects 0-byte files | `data/download.py`, `provenance.py` |
| 3 | Unknown `--dataset` key → exit 2; failed discovery → a `failed` result; missing/unverified required file or skipped product → exit 1 | `data/download.py`, `cli.py` |
| 4 | `DownloadResult.ok` now requires `verification['ok']`; manual files failing their checksum are `failed`; ambiguous glob matches are `failed`; the real filename is recorded, not the wildcard | `data/download.py` |
| 5 | Units are **converted**, never asserted. A unitless source column requires a declared input unit (`{column: X, unit: Y}`); validation requires the product to be exactly in the schema unit | `data/base.py`, `data/schema.py` |
| 6 | `ColumnRole.kind` / `is_identifier`; float identifiers rejected; generic `ID` removed from `source_id` candidates | `data/schema.py`, loaders |
| 7 | Products record raw path/sha256/bytes/HDU/column map/source units/package version; `check_against_manifest` raises if the raw file changed since it was verified | `data/base.py`, `data/_product.py` |
| 8 | Profile schemas carry `r_lower/r_median/r_upper` and asymmetric `*_err_lo/_err_hi` (following published Tables 2 and 3); unclaimed columns kept as `src_*` or reported in `meta` | `data/omegacat.py`, `data/base.py` |
| 9 | `tempfile.mkstemp`, `Content-Length` enforced, verification on the **temp** file, then `os.replace` — a corrupt refetch cannot destroy a good copy | `data/download.py` |
| 10 | Masked entries and NaNs fail any role without `allow_missing`; empty tables rejected; `valid_range` enforced | `data/schema.py` |
| 11 | Unknown override keys and two-roles-one-column rejected; the column-map path resolves dynamically from `OCEN_DM_ROOT` | `data/schema.py`, `data/base.py` |
| 12 | `write_manifest` merges by `(dataset, name)`, records failures, writes atomically; processed ECSV written atomically | `provenance.py`, `data/base.py` |

### Two further defects found by my own audit of the fixes

- **`Column.to()` drops the mask.** Converting a masked column silently turned a
  missing measurement into a real-looking number (masked `2000 m/s` → `2.0 km/s`,
  unmasked). Conversion now rebuilds the `MaskedColumn` with its mask.
- **Missing values were miscounted.** One masked row plus one NaN row was
  reported as one missing value, and the `valid_range` check used the wrong
  "usable" mask. Both now use the union of masked and non-finite.

Both have regression tests (`test_selfaudit_*`).

### Consequences worth remembering

- Plain CSV/ASCII carries no units, so the Kuzma 2026 supplementary table will
  require declared input units in `configs/column_maps.yaml` before it loads.
  This is intended (fail closed) and is documented in `docs/MANUAL_DOWNLOADS.md`.
- The `fimbulthul` product now builds from the real VizieR table end to end:
  309 rows, lineage sha256 `9d06aca4…`, written to
  `data/processed/tails/fimbulthul_members.ecsv`.

---

## 2026-09-16 (round 2) — Codex re-review of the fixed tree

Codex re-reviewed the fixed code: **5 fixes genuine, 7 partial, none wrong**,
plus 10 remaining/new defects and a criticism of two of my regression tests.
Everything below was reproduced independently before being fixed.

### Confirmed and fixed

| Sev | Defect | My evidence |
|-----|--------|-------------|
| P1 | `check_against_manifest` returned `"match"` for a file whose recorded verification **failed** — a checksum-rejected file could fail `fetch-data` and then pass `preprocess` | reproduced: entry with `ok=False, checksum_ok=False` → `"match"` |
| P1 | A failed refetch **erased the verified anchor**: the manifest entry became `sha256=None` while the good file stayed on disk, so later preprocessing saw `"unregistered"` | reproduced: `sha256 ee4d5e2c… → None` |
| P1 | Dimensionless roles bypassed unit handling: `0.5 percent` passed a `[0,1]` probability role as `0.5` instead of `0.005` | reproduced |
| P2 | `load()` / `load_profile()` bypassed validation entirely — `dec = 120°` and a float Gaia id came straight back | reproduced |
| P2 | Empty discovery, and `--no-discover` on a dataset with no declared files, both exited **0** | reproduced |
| P2 | **Regression I introduced:** `load_registry` dropped `published_algorithm`, silently turning a declared SHA-256 into MD5 | reproduced |
| P2 | `is_identifier` rejected float dtype but not an *object* array of floats; blank identifiers passed | reproduced |
| P2 | Profile validation accepted `r_lower > r_upper`, overlapping bins and negative dispersions | reproduced |
| P2 | Retry cleanup could `os.close()` a descriptor `os.fdopen()` already owned; `KeyboardInterrupt` bypassed temp-file cleanup | fixed by explicit ownership tracking + `except BaseException` |

Fixes: manifest entries now carry `last_verified`, and consumption rejects a
recorded checksum failure; `ColumnRole.dimensionless` converts percent → pure
number; `load()` validates by default (`validate=False` only for inspection);
a required dataset with nothing to fetch is an explicit `failed` result;
`TableSchema.ordered_triples` enforces `r_lower <= r_median <= r_upper` and
non-overlap; dispersions and uncertainties carry non-negative bounds.

### Two of my regression tests passed for the wrong reason

- `test_f1_invalid_table_is_not_written` used RA in km, which fails in
  `standardize()` **before** the code I fixed. Now uses duplicate identifiers,
  which pass resolution and conversion and fail only in `validate_table`.
- `test_f9_temporary_files_are_unique` collected names but never asserted
  uniqueness; `test_f12_processed_write_is_atomic` never exercised a failed
  write. Both now assert the real property (distinct temp names; previous
  product byte-identical after an injected mid-write `OSError`).

Codex was right about both. Test suite now **118 passed**.

### Accepted, not fixed

- **Concurrent manifest merging is not transactional** (two `fetch-data`
  processes can lose each other's entries). Single-user pipeline; documented
  rather than locked.
- Codex's remaining milestone gaps — enforced coordinate frame/epoch, Gaia
  release, velocity frame, tracer selection, and the fact that the combined and
  radial/tangential PM profiles **share stars and are not independent
  likelihood factors** — are recorded as Milestone 3 requirements, not Milestone
  1 work. The last point matters: treating all produced profiles as independent
  would overstate the kinematic information and bias the DM posterior.

### State

`ocen fetch-data --dataset fimbulthul_ibata2019` → cached, size-checked;
`ocen preprocess --dataset fimbulthul` → 309 rows, lineage `manifest: verified`.
Zenodo was still returning 504 at the end of the session, so the two Zenodo
datasets remain unfetched.

---

## 2026-09-16 (evening) — Zenodo recovered: Milestone 1 data is in

Zenodo returned HTTP 200 again at ~20:35 UTC after ~3.5 h of 504s. `ocen fetch-data`
was re-run unchanged and the blocked datasets came down.

### Acquired and verified

| Dataset | Result |
|---|---|
| Kuzma & Ishigaki 2025 | `wCen_table.fits`, 25.2 MB, **157481 rows**, MD5 `3dcd58e6…c768` — **matches the published checksum**, which was UNVERIFIED this morning |
| oMEGACat VI | 21 files, every one **checksum-verified** against the Zenodo API digests |
| Fimbulthul | 309 rows from VizieR (unchanged) |
| Kuzma et al. 2026 | still missing — see below |

The oMEGACat discovery path worked exactly as rebuilt this afternoon: it parsed each
`md5:<digest>` from the API, carried it into the `FileRecord`, and enforced it on
download. Before the fix it would have discarded all 21 checksums.

### Real oMEGACat VI structure (was unknowable while Zenodo was down)

`proper_motion_dispersion_log_bins.fits` (40 bins) and `los_profile.fits` (29 bins) carry
exactly the structure inferred this morning from the published Tables 2–3:
`r_lower, r_median, r_upper, nstars`, then values with **asymmetric** errors
`e_lower_*` / `e_upper_*`. The redesign under finding 8 was right; only the names differ.

Two things checked rather than assumed:

- `vlos` in `los_profile.fits` ranges 0.76–10.34 km/s, so it is the **rotation amplitude
  per bin**, not the ~232 km/s systemic velocity. Mapped to the `v_rot` role.
- The three PM profiles come from **one file** and share bins and stars; the same holds for
  the two LOS profiles. Recorded in `configs/data.yaml`: fitting more than one binning of
  the same stars would double-count the data. This was Codex's Milestone 3 warning and it
  is now a fact about the files, not a hypothetical.

### Kuzma & Ishigaki 2025 columns

157481 rows, and **no column carries a unit**, so the strict unit rule fired and refused to
load until each input unit was declared. Declared in `configs/column_maps.yaml` with the
justification recorded there: RA 194–209°, Dec −52 to −42° around ωCen; PMs ±60 with errors
0.01–0.1 (Gaia mas/yr); `G_0` 6.5–16 mag; `FeH_CaHKsyn` −4 to 0 dex; `P_mem` ∈ [0,1].
2737 stars have `P_mem > 0.5`. `source_id` is full-precision int64.

Schema correction: the catalogue's magnitudes are **dereddened** (`G_0`, `BP_0`, `RP_0`,
`CaHK_0`), so the roles are now `g0_mag`, `bp0_mag`, `rp0_mag`, `cahk0_mag`. Mapping `G_0`
onto a role called `phot_g_mean_mag` would have labelled a corrected quantity as a raw one.

### Kuzma et al. 2026 — researched, still not obtainable programmatically

Used the new `ask-codex --research` path, then verified its claims:

- **ADS** (real bibcode `2026MNRAS.550g1147K`, found with the local token) lists exactly one
  data product: `ESO:1`. No CDS link.
- **VizieR** has no catalogue for this paper (MNRAS 550 holds only two, both unrelated).
- **arXiv source** downloads fine (2.5 MB; Codex wrongly reported this as failing) but
  contains only `wCen_paper.tex` with **five sample rows** and the note "This table is fully
  available as online material".
- **OUP** returns HTTP 403 to automated requests; the supplement needs a browser session.

So the only route is the publisher's supplementary file, via an authenticated browser.

The paper's Table 3 column structure is now known from the LaTeX: `Star ID, DR3 Source ID,
R.A., Dec., G, V_h, σ_Vh, [Fe/H], σ_[Fe/H], P_mem, Member`, with 616 observed targets, 593
after quality cuts and 157 adopted members.

**Data-quality warning for when the table arrives:** the DR3 source IDs printed in the
LaTeX are float-damaged — e.g. `6086552125451670000`, a 19-digit identifier ending in four
zeros. If the supplementary file carries the same values, the IDs cannot be used for
cross-matching and must be recovered by position. The loader's `is_identifier` rule rejects
float identifiers, so this will fail loudly rather than silently mis-join.

### State

`ocen preprocess` builds 7 validated products: 5 oMEGACat profiles, the Kuzma 2025
periphery catalogue and Fimbulthul. Tests: **118 passed** (one test-isolation bug fixed —
it read the now-populated real column map instead of a clean one).

Milestone 1 is complete for three of the four primary datasets.

---

## 2026-09-16 (evening, cont.) — Kuzma et al. 2026 delivered: Milestone 1 complete

The user supplied the extended version of Table 3 directly, which was the one route
research had identified as the only one available (publisher supplementary material;
`academic.oup.com` returns HTTP 403 to automated requests). Saved verbatim as
`data/raw/kuzma2026_spectroscopy/kuzma2026_members.txt`.

### Integrity checks on the delivered file

- **592 data rows**, every one with exactly 11 fields.
- **157 rows flagged `Member = True`** — matches the 157 adopted members stated in the
  paper. An independent check that the table is complete, not a fragment.
- The paper reportedly quotes 593 stars surviving quality cuts against the 592 rows here.
  That 593 came from an unverified secondary summary, so the one-row difference is noted,
  not resolved.
- **The DR3 source IDs are full precision** — `6086552125451679744` for `OT_T_325`, and
  **none of the 592 ends in `0000`**. The float damage seen this afternoon
  (`6086552125451670000`) is confined to the typeset LaTeX in the arXiv source; the
  supplementary data is clean. The warning was right about the paper and wrong as a
  prediction about this file, which is the better outcome: the IDs are usable for
  cross-matching as delivered.

### Ingestion

The file is whitespace-separated ASCII whose column-name line is the **last** comment line,
which astropy's default reader does not find (it returned `col1…col11`). Rather than
hard-code that in the loader, `configs/data.yaml` now carries `format:
ascii.commented_header` and `reader_kwargs: {header_start: -1}`, and the settings are
recorded in the product's lineage so the exact read is reproducible.

Units come from the file's own header block (degrees, mag, km/s, dex) and are declared in
`configs/column_maps.yaml`; the strict unit rule refused to load it otherwise. Three roles
were added to the schema for columns the file actually carries: `star_id`, `g_mag` and
`membership_prob`.

### Scientific sanity check on the 157 members

- mean heliocentric velocity **234.8 km/s**, dispersion **7.3 km/s** — consistent with the
  systemic velocity of ω Cen (~232 km/s), so the flagged members are genuinely associated;
- median **[Fe/H] = −1.21**.

Neither number was assumed; both were computed from the delivered table.

### State — Milestone 1 complete

All four primary datasets are ingested, validated and written:

| Product | Rows |
|---|---|
| `kuzma2025_periphery` | 157481 |
| `kuzma2026_spectroscopy` | 592 (157 members) |
| `fimbulthul_members` | 309 |
| `omegacat_vi_pm_radial` / `pm_tangential` / `pm_combined` | 40 bins each |
| `omegacat_vi_los_dispersion` / `los_rotation` | 29 bins each |

`ocen fetch-data` exits 0, `ocen preprocess` builds 8 products, tests: **118 passed**.

Also fixed: `ocen inventory` reported a present wildcard file as MISSING because it did not
resolve the glob the way `fetch-data` does.

### Next

Milestone 2, the mass-profile library: stellar MGE from the light profile, compact remnant
component, IMBH point mass, cored and cuspy truncated DM, and a composite enclosed
mass/potential/force/escape-speed interface, with every enclosed-mass curve verified
numerically against an independent integration.

---

## 2026-09-16 (night) — Milestone 2: mass-profile library

`src/ocen_dm/mass_models/` (1044 lines) implementing
`Phi_tot = Phi_star + Phi_rem + Phi_IMBH + Phi_DM`.

| Module | Contents |
|---|---|
| `base.py` | `MassComponent` interface, `G`, numerical fallbacks, `verify()` |
| `stellar.py` | `Plummer`, `MGE` (spherical multi-Gaussian, deprojects consistently) |
| `remnants.py` | `RemnantPlummer` — free `M_rem`, `a_rem`, its own name in the posterior |
| `imbh.py` | `PointMass`, `M_bh >= 0` |
| `dark_matter.py` | `TruncatedGNFW` (`.nfw()` γ=1, `.cored()` γ=0) and `Burkert` |
| `composite.py` | `CompositeMassModel` + `M_DM(<r)`, `f_DM(<r)`, `v_esc(r)`, `r_J` |

**Units.** pc, Msun, km/s internally; `G = 4.300917270036e-3` asserted against
`astropy.constants` (CODATA 2018) by a test, so it cannot drift.

**Truncation.** `T(r) = (1 + (r/r_t)^2)^-2`: smooth everywhere, keeps the total mass
finite, falls as r^-4 outside r_t. It makes the enclosed mass non-analytic, so truncated
profiles integrate numerically while the untruncated limits keep their closed forms — which
is what tests the quadrature.

### Verification

Every component's analytic enclosed mass agrees with quadrature of its own density to
**4e-13 (Plummer), 3.8e-9 (MGE), 1.9e-11 (NFW), 2.2e-7 (Burkert)**, and every analytic force
agrees with the derivative of the potential to **2e-9 – 4e-5**.

**Independent cross-check against AGAMA** (no shared code): Plummer enclosed mass to 1.2e-15,
potential and force to 1.5e-10 — which is exactly the difference between AGAMA's value of G
and CODATA 2018, i.e. the floor of the comparison, not an error in either. NFW density,
enclosed mass and potential to 2.4e-6, 2.0e-7 and 8.9e-9.

### Two defects found and fixed, both in my own checks rather than the physics

- **The force test was arithmetically broken.** A central difference of the potential with a
  10^-5 step failed for every cored profile. The cause was catastrophic cancellation: at
  r = 10^-3 pc the Plummer potential is −860 (km/s)² and varies by 7e-10 across the step,
  only ~3600× double-precision noise. For components whose potential is itself computed by
  quadrature the effect was worse, because quadrature error is not smooth in r and
  differencing amplifies it — hence force "errors" of 2–5%. Replaced by a five-point stencil
  with a *larger* step (`dphi_dr`): truncation falls as step⁴ while roundoff falls as 1/step.
  Errors dropped by three to four orders of magnitude and all eight components now pass.
- **My first AGAMA comparison was wrong in the units,** reporting a 99.9% potential
  discrepancy. `agama.setUnits(length=1e-3)` makes the length unit 1 pc, so radii must be
  passed in pc; I was also scaling them by 1e-3. The enclosed mass agreed anyway — because
  scaling r and the scale radius identically leaves `M(<r)` invariant — which masked the
  error until the potential exposed it. A reminder that agreement in a scale-invariant
  quantity is not evidence of correct units.

### Sanity check against the Milestone 1 data

A deliberately simple, **unfitted** fiducial (MGE 2.2e6 + 1.3e6 Msun at 3 and 12 pc, plus
1.5e5 Msun of remnants) against the oMEGACat VI LOS profile now on disk (29 bins,
0.15–7.6 pc, σ 21.4 → 12.5 km/s at the kinematic distance 5.43 kpc):

- total mass 3.65e6 Msun, against the ~3.5e6 Msun usually quoted for ω Cen;
- median `v_circ/σ_los` = **0.91** across the measured bins — order unity, as a bound system
  requires;
- `M(<5 pc) = 1.40e6`, `M(<20 pc) = 3.09e6` Msun; `v_esc(5 pc) = 66.2` km/s;
- Jacobi radius 174 pc for R = 6 kpc and `M_G(<R) = 5e10` Msun; the escape speed to that
  radius is 64.8 km/s rather than the 66.2 km/s to infinity.

The innermost bin is the informative failure: the fiducial gives `v_circ = 2.0` km/s at
0.15 pc against a measured σ of 20.9 km/s, because a Gaussian core has `M ∝ r³`. The real
cluster keeps σ ≈ 20 km/s inwards, so an MGE fitted to the actual light profile will need
narrow components. That is Milestone 3's work, and the parameters above were invented for a
scale check, not fitted to anything.

### Tests

`tests/test_mass_models.py`: **48 tests**, covering the specification's section 13 items 1–5
(units, mass profiles, potential gradients, limiting cases). Whole suite: **166 passed**.

### Next

Milestone 3, the kinematic likelihood: an MGE fitted to a real ω Cen light/star-count
profile, a Jeans solver (JamPy is not yet installed in the venv — AGAMA can serve as the
cross-check), and the oMEGACat VI profile likelihood, respecting the recorded fact that the
combined/radial/tangential PM profiles share stars and are not independent.
