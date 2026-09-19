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

---

## 2026-09-16 (late) — public repo, fast path, solvers, figures

### Repository

`git init`; public repo **https://github.com/vasilybelokurov/ocen_dm**. Pre-flight before
publishing: secret scan of every committable file (clean), `.gitignore` covers `data/`,
`results/`, `.DS_Store`, review packets and temp files. Raw and processed data (528 + 55 MB)
are not committed; `ocen fetch-data` regenerates them. Two commits pushed so far.

### Fast path (was the top-ranked problem)

Measured before: a composite with a truncated NFW cost **256 ms** per evaluation (adaptive
quadrature for both M(<r) and Φ), i.e. 7 h per 10⁵ likelihood calls. A nested sampler changes
the halo parameters on every call, so caching across calls cannot help; the tables themselves
had to be cheap to build. `ProfileTables`: vectorised Simpson on a fixed 2000-point log grid
over 10⁻⁵–10⁶ pc, inner power-law and outer-tail extensions from the local log slopes, cubic
splines in ln r. Build + evaluate: **0.54 ms** (470×). Agreement with quadrature 1.5–7.5e-8.
The quadrature paths survive as `quad_enclosed_mass` / `quad_potential`, and `verify()` now
compares against them explicitly rather than against a copy of itself. Nine new tests,
including a speed budget (<5 ms per fresh-parameter evaluation). Side effect: the mass-model
suite runs in 2.3 s instead of 9.7 s.

### Solvers

Installed into the venv: **ultranest 4.5.0, dynesty 3.1.0, jampy 9.0.2**. JamPy's licence
(read from the installed `LICENSE.txt`): non-commercial use, modification for internal use,
**no redistribution**. It stays a dependency; nothing from it enters the public repository.

### Figures (`plots/`, regenerated by `ocen plot-data`)

Seven PNGs from `src/ocen_dm/plotting/data_overview.py`, built with the validated reference
palette (three categorical slots — blue Kuzma 2025, orange Kuzma 2026, aqua Fimbulthul — pass
every colour-vision check for scatter; neutral grey for field; one-hue blue ramp for density).

| File | Shows |
|---|---|
| `sky_overview.png` | wide field (Fimbulthul ~20° north) + periphery zoom with the FLAMES fields |
| `kuzma2025_proper_motions.png` | field density (hexbin) and the member clump at (−3.2, −6.7) mas/yr |
| `kuzma2025_cmd_metallicity.png` | dereddened CMD; members' [Fe/H] peaks at ≈ −1.75 against a field near 0 |
| `kuzma2026_spectroscopy.png` | v–[Fe/H] plane and velocity histogram; members at ⟨v⟩ = 234.8, σ = 7.3 km/s |
| `fimbulthul.png` | positions with coherent PM vectors (4 mas/yr key), and the CMD |
| `omegacat_vi_profiles.png` | σ_PM (3 series, asymmetric errors), σ_R/σ_T, σ_LOS, v_rot |
| `omegacat_vi_rotation_axis.png` | θ₀(r), wrapped to a window on the outer median 105° |

Rendered and inspected before acceptance; four were wrong on first render and fixed: the sky
map was squashed by the 20° offset to Fimbulthul (split into two panels), PM arrows were too
long to read, θ₀ split across ±180°, and two legends sat on data. Physics read off the
figures matches the literature: the member PM clump, the [Fe/H] ≈ −1.75 peak, σ_R/σ_T rising
to 1.18 outward, the FLAMES fields strung along the NE–SW tail axis.

### Codex

Usage limit reached; resets **19 Sep 2026, 09:10**. The Milestone 2 critical review is
pending until then.

---

## 2026-09-16 (late, cont.) — the velocity figures that were missing

The user asked where the v_los and 3D-velocity plots were; there were none. `v_los` had
appeared only in the Kuzma 2026 v–[Fe/H] plane, and the 1.48M-star oMEGACat
`catalog_and_selections.fits` — 24,928 stars with high-quality **proper motions and v_los** —
had been downloaded but never plotted. Two figures added.

### `kuzma2026_vlos_along_tails.png`

v_los against signed angular distance along the tails, sign from the target-name field code
(I/M/O × L/T; leading drawn positive), members highlighted, per-field member means labelled.
Inner fields: **IT (trailing) 236.8 ± 0.8 km/s (n=75)** versus **IL (leading) 233.0 ± 0.8
km/s (n=77)** — a 3.8 ± 1.1 km/s offset between the two sides at ~0.7° from the centre. The
middle and outer fields hold only 0–2 members each, so the tails' velocity structure is
measured essentially at the inner fields alone. Recorded as data, not interpreted.

### `omegacat_vi_3d_velocities.png`

Four panels from the high-quality PM+LOS sample: the LOS velocity field (a clean rotation
dipole, ±5 km/s); rotation curves in the plane of the sky (from PM) and along the line of
sight (paper); the three velocity components of the same stars; and σ_PM converted at 5.43
kpc against σ_LOS.

Facts read off it:

- **σ_R = 16.8, σ_T = 15.4, σ_LOS = 16.9 km/s** for the same 24,928 stars — the three
  components agree once the PMs are converted at the paper's kinematic distance, and the
  σ_PM(r) and σ_LOS(r) profiles track each other over two decades in radius. That agreement is
  the kinematic distance; the plot shows it directly.
- σ_T < σ_R in the plane of the sky: mild radial anisotropy, matching σ_R/σ_T rising outward in
  the profile figure.

### The θ₀ convention, settled empirically

The tabulated θ₀ (paper: "position angle of the rotation axis", 104.3 ± 1.4° for r > 30″) did
not match the map, whose zero-velocity line runs nearly N–S. Rather than guess, I fitted
`v = v0 + A sin(φ − φ0)` to the 24,928 stars with φ the PA from North through East, r > 30″:
**v0 = 232.8 km/s, A = 5.7 km/s, φ0 = 12.9°**. The paper's value equals φ0 + 90° (102.9° vs
104.3°, within the errors): the tabulated angle marks the direction of *maximum* rotation
velocity in the North-through-East convention, i.e. it is the rotation-axis angle in a frame
whose origin differs by 90°. Recorded in the role description in `omegacat.py` and drawn on
the map. The paper's PMs are stated to be measured relative to the bulk motion; the
plane-of-sky ⟨v_tan⟩(r) panel shows empirically whether differential rotation survives that.

Two figures were wrong on first render and fixed: field labels colliding with a title, and a
per-cell PM rotation *map* that was noise-dominated (~2 km/s per-cell error against a
few-km/s signal) — replaced by a binned rotation *curve* with errors.

---

## 2026-09-16 (night) — the canonical catalogue store; Vasiliev & Baumgardt 2021 ingested

### Correction of approach

I had started copying a 119 MB Zenodo zip from another project into `data/raw/`. The user
stopped that twice: external catalogues live in **`~/data/catalogues/`**, the canonical
store, and are referenced there, not duplicated. Saved as a memory. The pipeline gained a
retrieval kind **`local`**: a file at its canonical path, verified in place (sha256, size,
published checksum where one exists), recorded in the manifest, never copied — and
`fetch-data` no longer creates an empty `data/raw/<dataset>/` for such datasets (a defect a
test caught).

### What was found there

- **`gc_members_gaia_vasiliev.fits`** — the user's FITS compilation of Vasiliev & Baumgardt
  (2021, Zenodo 10.5281/zenodo.4891252) for 170 clusters. Its ω Cen block was checked against
  the zip's `catalogues/NGC_5139_oCen.txt`: **identical `source_id` set (228,055)**, PMs equal
  to 1.1e-5 mas/yr (float32 storage), membership probabilities to 3e-8. The zip itself
  (`apogee_halo_rotation/data/external/`) matches the published MD5 `3a07093b…`.
- **`gc_catalog_updated.fits`** — Baumgardt & Hilker database compilation (mtime 2025-05-16;
  database version not recorded → flagged UNVERIFIED against upstream). ω Cen: μ = (−3.236,
  −6.731) ± 0.011 mas/yr, D = 5.43 ± 0.05 kpc, RV = 232.78 ± 0.21 km/s, M = 3.94e6 Msun,
  r_h = 7.56 pc.

Paper reference verified live: arXiv 2102.09568 → "Gaia EDR3 view on Galactic globular
clusters", DOI 10.1093/mnras/stab1475, MNRAS 505, 5978.

### Products (3 new, 12 total)

| Product | Rows | Note |
|---|---|---|
| `vasiliev2021_ocen_members` | 228,055 | to G = 21, 0.67° ≈ 63 pc; 156,530 at P > 0.9; per-star PM covariance |
| `vasiliev2021_ocen_pm_profiles` | 101 radii | authors' σ_PM and v_rot percentiles; **file has 12 columns, readme documents 11** — the 12th carried as `undocumented_col12`, not interpreted |
| `baumgardt_ocen_parameters` | 1 | systemic phase space; UNVERIFIED flag propagated into the likelihood rule |

`build_product` gained `row_filter` (rows in file / rows selected recorded in the lineage).
Units for both FITS files, which carry none, are declared from the zip readme and the
database's conventions in `configs/column_maps.yaml`.

### Checks against what we already had

- Systemic PM: catalogue (−3.236, −6.731) vs our Kuzma 2025 member centroid (−3.259, −6.726):
  offsets **23 and 5 µas/yr**. RV 232.78 vs our 234.8 (Kuzma 2026 members) and 232.6
  (oMEGACat median). D 5.43 in both catalogue and oMEGACat — the "distance tension" I flagged
  in the plan does not exist in these files.
- Reach: HST PM profile 0.05–8.2 pc; **EDR3 profile 0.6–63 pc** — into the spec's 50–200 pc
  window but not through it. Beyond 40′ only the Kuzma 2025 giants (2737) remain, plus the
  Baumgardt 2019 DR2 profile still to fetch.
- **EDR3 σ_PM sits systematically below HST σ_PM where both exist**, converging outward:

  | r | HST (oMEGACat) | EDR3 (VB21) |
  |---|---|---|
  | 60″ | 18.9 km/s | 14.1 (+2.0/−1.3) |
  | 120″ | 17.2 | 13.8 (+1.6/−1.1) |
  | 200″ | 15.5 | 13.5 (+1.1/−0.8) |
  | 300″ | 13.4 | 12.9 (+0.6/−0.5) |

  Recorded as a fact; two candidate explanations, **neither established**: (i) the EDR3
  profile inside ~1′ is affected by crowding-driven systematics in a way its error scaling
  does not fully capture; (ii) **energy equipartition** — inside the core Gaia sees only
  bright (massive) giants, which move more slowly than the fainter stars HST measures.
  oMEGACat VI measured equipartition directly, and its `energy_equipartition_*.fits` products
  are on disk. Either way the two datasets **cannot be combined naively**: the Milestone 3
  likelihood must either model σ(m) or restrict both to a common magnitude range. This is
  the most consequential thing learned today about the modelling.

Tests: **188 passed**.

### WSDB timing, settled (WP4 feasibility)

The 162 s `local_join` earlier today was a transient stall, not the cost of the query: on
re-measurement the same 157,481-id join to `gaia_dr3.gaia_source` took **3.5 s** (1,000 ids
0.4 s; 10,000 ids 1.1 s), and an upload-free `unnest(%s::bigint[]) LEFT JOIN` form 2.2 s.
A 1° cone count on Gaia DR3 (401,654 stars) takes 1.2 s. One run also produced a `COPY`
timeout at 10k rows that vanished on retry. Conclusion for WP4 and for the `wsdb` skill,
which now carries a measured-performance table and the covariance-join recipe: one call,
seconds, no chunking; minutes or a `COPY` error mean "retry", not "split the input".

---

## 2026-09-17 — WP1/WP2/WP5: light profile, MGE, outer kinematics

### Fetched through the pipeline (VizieR)

- **Trager, King & Djorgovski 1995** `J/AJ/109/218/tables`: 73 V-band points for `ngc5139`,
  10.5–2588″ (0.28–68 pc), μ_V 16.15–28.17, with the authors' weights and data-set flags.
  Units come from the VOTable itself (`logr` in log(arcsec), μ in mag/arcsec²).
- **Baumgardt+ 2019** `J/MNRAS/482/5138/table4`: the Gaia DR2 σ_PM profile of NGC 5139 —
  9 bins, 179–1747″ (4.7–46 pc), asymmetric errors. The outer kinematics the plan flagged.
- **Noyola & Gebhardt 2006** `J/AJ/132/447`: checked and **not applicable** — 0 rows for
  NGC 5139 in either table. Recorded in the registry so the check isn't repeated.

Two ingestion rules had to be refined, both because the strictness caught real structure:
CDS uses `E_sigma`/`e_sigma` (upper/lower error) — a case collision the resolver refused —
so matching is now exact-first with a case-insensitive fallback that refuses ambiguity; and
Trager's `logr` carries a `log(arcsec)` unit, so `build_product` gained a recorded
`transform` step (`r_arcsec = 10**logr`, noted in the lineage as a re-expression, not new
information).

### `src/ocen_dm/light_model.py` — the MGE

Spherical MGE fitted **in projection** by weighted non-negative least squares with the
widths on a log grid; deprojection is analytic and keeps the widths, mapping directly onto
`mass_models.MGE`. Only the profile's *shape* is used (M/L is a free parameter, so the
zero-point and extinction cancel). Rows are scaled by `sqrt(weight)/I`, i.e. the fit is in
magnitudes with the authors' weights as inverse variances.

Result on the real profile (default 64-node grid): **11 Gaussians, weighted rms 0.18 mag**
against **0.21 mag for Trager's own Chebyshev fit**; the worst residuals are Trager's own
worst points (ours −1.62 vs theirs −1.67 at 2148″, weight 0.10; the SB-2342 points at weight
0.03). **Projected half-light radius 280″ = 7.37 pc** — against Harris's 300″ and the
Baumgardt catalogue's 7.56 pc. That agreement also suggests the catalogue's `RH` is the
projected half-light radius rather than the 3D half-mass (9.70 pc for this MGE) — recorded
as suggestive, not settled. `plots/light_profile_mge.png`.

Grid density measured rather than assumed: a noiseless single Gaussian is recovered to
0.36 / 0.10 / 0.027 / 0.007 mag rms with 16 / 32 / 64 / 128 nodes (fixed-grid NNLS cannot
place a node exactly on the true width), while the real data plateau at 0.184 mag from 32
nodes on — the floor is the data scatter, not the grid. Default set to 64.

One test I wrote was wrong, not the code: a +3 mag outlier at weight 0.01 legitimately bends
a noiseless fit by ~0.1 mag (χ² = 0.09 is not zero); the test now asserts what matters —
the unit-weight points stay fitted and the half-light radius is unchanged when the outlier's
weight sits at the floor.

### State

14 products, **201 tests**, 10 figures. Milestone 3 now has its luminous model, its
systemic phase space, and internal kinematics spanning 0.05–46 pc from three instruments
(HST, MUSE, Gaia) — with the recorded caveat that the EDR3 and HST σ_PM disagree by 20–25%
where they overlap and must not be stacked.

Left from the plan: WP4 (Gaia covariance product; recipe written, trivial), WP6 (Kuzma 2025
footprint), and the star-count cross-check of the light profile from our own catalogues.

### WP4 and the star-count cross-check (same day)

**WP4 done.** `selection/gaia_covariance.py`: one `sqlutilpy.local_join` fetches the full Gaia
DR3 astrometric covariance for all 157,481 Kuzma 2025 stars in **3.7 s**. The product
`kuzma2025_periphery_gaia_covariance` carries `pmra_pmdec_corr` (median |corr| 0.11),
the parallax terms and RUWE, and the build refuses to write if any star is unmatched, if the
catalogue's PMs differ from Gaia's beyond their rounding (the release check), or if a
correlation is out of range. Five tests with a mocked join.

**Light profile validated by star counts.** `plots/light_profile_star_counts.png` compares
the Trager MGE with surface-density profiles from the three catalogues we hold, each anchored
only where it is complete:

| Sample | Anchored | rms vs MGE |
|---|---|---|
| HST oMEGACat, F625W < 20, hq (222,792 stars) | 0.5–4.5′ | **0.09 mag** |
| Gaia EDR3 members P > 0.9, G < 19 (76,118) | 8–30′ | **0.19 mag** |
| Gaia+Pristine members P > 0.5, G₀ < 16 (2,737) | 15–42′ | **0.09 mag** |

So the MGE shape holds from 0.5′ to 42′ (0.8–66 pc) against three independent tracers. Gaia
is crowding-incomplete inside ~6–10′ (residuals rise to +4 mag inward) and HST inside ~0.3′;
beyond the Trager data (43′) the MGE is extrapolation and the member counts sit above it.

The first version of this figure anchored EDR3 over 2–20′, inside its incomplete core, which
shifted the whole sample 1.4 mag high and made the MGE look wrong at 5–11′. The slopes had
agreed all along (6.1 vs 6.3 mag over 7–36′); the anchoring was the error. Fixed and the
anchoring ranges are now drawn on the figure.

**Consequence for the model:** the luminous profile beyond ~43′ (> 68 pc) is not constrained
by the light data and would need the member star counts as a tracer — with their own
selection effects — if the halo's outer extent is ever to be tied to the stars there. Inside
that radius it is solid.

State: 15 products, **206 tests**, 11 figures. Remaining from the plan: WP6 (Kuzma 2025
selection footprint). Milestone 3 can start.

### Correction: the distance I attributed to oMEGACat VI

I had been calling 5.43 kpc "the oMEGACat kinematic distance" (plots, journal 2026-09-16, plan)
and concluded there was no distance tension. Wrong: reading the oMEGACat VI source, its
kinematic distance is **5494 ± 61 pc** (section 4.4). The 5.43 ± 0.05 kpc is the Baumgardt &
Vasiliev (2021) catalogue value (`baumgardt_ocen_parameters.distance`). The two differ by
0.8σ — mild, and a 1.2% effect on pc scales and on the PM→km/s conversion, so no result above
changes materially; but the record was wrong and is now corrected in the plotting constants,
the plan and here. D stays a parameter in Milestone 3 with both values as priors/checks.

---

## 2026-09-17 — Milestone 3 begins: the spherical Jeans solver, verified three ways

Decision (approved): our own spherical Jeans solver is primary; JamPy and AGAMA are the
independent checks; JamPy stays available for the axisymmetric extension.

### `src/ocen_dm/kinematics/`

- `anisotropy.py` — the spec's β(r) = β₀ + (β∞−β₀) r²/(r²+r_β²), with its integrating
  factor g(r) = r^{2β₀}(r²+r_β²)^{β∞−β₀} in closed form. **This is exactly JamPy's
  `logistic` anisotropy with α = 2**, which is what lets the JamPy comparison be exact.
- `jeans.py` — `SphericalJeans(mass, tracer, anisotropy)`: νσ_r² by cumulative Simpson on a
  600-point log grid (outside in), splined in the log; the three projected second moments
  (LOS, PM radial, PM tangential) by the substitution r = R cosh u, which removes the
  singularity at r = R and leaves smooth integrals that 96 Gauss–Legendre nodes handle,
  vectorised over R. `dispersions_observed(R_arcsec, D)` returns σ_LOS in km/s and the PM
  dispersions in mas/yr. **5 ms per model** (build + three projections at 70 radii).

### Verification

| Check | Result |
|---|---|
| Isotropic Plummer, σ_r(r) = √(GM/6√(r²+a²)) | 1.3e-5 |
| Isotropic Plummer, σ_p(R) = √(3πGM/64√(R²+a²)) | 1.1e-5 |
| Isotropy identity (pmr = pmt = los) | exact |
| **JamPy** `jam_sph_proj`, 3-Gaussian MGE + 10⁴ M☉ BH, β 0.05→0.5, 25 radii, all three projections | **≤ 1.0e-3** (JamPy's own interpolation level) |
| **AGAMA** QuasiSpherical DF, Plummer, β = 0 / −0.5 | **2.8e-5 / 3.6e-5** |
| AGAMA, cuspy NFW tracer, β = 0 / +0.3 | 3e-2 (limited by AGAMA's Multipole+DF numerics; same at β = 0) |

One defect found and fixed: with a Gaussian tracer the density underflows to zero far
inside `r_max`; a spline through ln(0) developed a kink that the projection nodes sampled
as a 10¹⁰ km/s value at one radius. The table now stops where ν has fallen 20 decades below
its peak (`r_cut`, 131 pc for the test MGE).

### A physical constraint on the K1 prior, from the AGAMA check

For a cored Plummer tracer with constant β₀ = +0.1 / +0.3, AGAMA's realised density exceeds
the input by ×1.27 / ×2.58 at 0.3 pc: **no positive distribution function exists** for
radial anisotropy in a core (An & Evans 2006: the tracer cusp must satisfy γ ≥ 2β). Our
solver still returns a Jeans solution for such models — the equation has one — but the
model is unphysical. The ω Cen MGE is cored (μ_V flat inside ~100″), so **β₀ ≤ 0** in K1;
radial anisotropy may only develop outward (β∞ > 0), consistent with Watkins+ 2015's
isotropic cores and Watkins+ 2013's global β = 0.10. Encoded as a test so it is not lost.

Tests: `tests/test_jeans.py`, 17 tests (two external, skipped where JamPy/AGAMA absent).

## 2026-09-17 — Milestone 3, part 2: the likelihood, first fits, and what the data said back

### `kinematics/likelihood.py`, `kinematics/fit.py`, `plotting/fits.py`

The likelihood connects the solver to the ingested profiles under three rules enforced
in code: (1) a published dispersion is a statistic of the stars between `r_lower` and
`r_upper`, so the model enters as the tracer-weighted bin mean of σ²
(`∫Σσ²R dR / ∫ΣR dR`, 6-node Gauss–Legendre per bin; profiles published without edges are
point-evaluated and say so); (2) asymmetric errors are a split normal; (3) the oMEGACat
combined-PM profile shares its stars with the radial/tangential ones and is refused
alongside either (`KinematicData` raises). Datasets: `hst_pm_radial`, `hst_pm_tangential`
(40 bins each), `muse_los_dispersion` (29), `gaia_dr2_pm` (Baumgardt+ 2019, 9 points, no
edges), `gaia_edr3_pm` (Vasiliev & Baumgardt 2021, fine grid from overlapping samples —
**thinned to 8 log-spaced points beyond 300″** because adjacent points are not independent;
using all 88 would give one dataset 8× the weight of HST). Both Gaia profiles are
*assumed* to be 1-D PM dispersions — flagged in the product notes, to verify against the
papers before publication.

Speed: sharing the projection nodes between Σ and the three second moments
(`SphericalJeans.projected_moments`) and batching every profile's radii into one solver
call took the K1 likelihood from 33 to **4.1 ms** (results identical to 1e-12).

`fit.py`: unit-cube priors (`uniform`, `loguniform`, `normal`, `truncnormal`), families
**K1** (stars = Trager MGE × M★, remnant Plummer (M_rem, a_rem), point mass M• with a
log-uniform prior down to 100 M☉ ≡ zero, β₀ ∈ [−1, 0] (An & Evans), β∞, r_β, nuisances)
and **K2** = K1 + truncated gNFW (γ = 0 or 1) parametrised by **M_DM(<100 pc)** and r_s so
the prior is flat in the quantity the data constrain (ρ_s derived; exact by construction,
tested). `FitProblem` wraps a family and data for ultranest/dynesty; `maximum_likelihood`
is a multi-start Nelder–Mead for smoke tests; `run_nested` writes posterior, summary,
per-sample `radial_profile()` on a 0.1–500 pc grid and a `run.yaml` with git commit and
input sha256s. Tests: `test_likelihood.py` (13), `test_fit.py` (8) — including a synthetic
injection recovered by the optimiser (total dynamical mass to 5 %, Gaia scale to 0.05).
Suite: **243 passing**.

### First real fit, and three defects it exposed

The first K1 maximum-likelihood fit (fixed D = 5.43) gave χ² = 879 for 126 points and a
model nobody would believe: M• = 3.4e4 M☉ making a cusp that overshoots the innermost HST
bin by 4.6σ, a dip at 8″ absent from the data, MUSE 2–3σ below the model everywhere
outside 30″. Each had a cause:

1. **Cold spike in the light model.** Unconstrained NNLS put 0.24 % of the light in two
   Gaussians of σ = 5–6″ — *inside* the innermost Trager datum (10.5″). A compact tracer
   component sitting in a harmonic core is dynamically cold: without a BH the model
   dispersion at 1–8″ collapsed to 0.41–0.53 mas/yr (data 0.84), and the optimiser bought
   it back with 3e4 M☉ of point mass. Fix: `MGE_SIGMA_RANGE_ARCSEC = (10.5, 3000)` — the
   smallest Gaussian is the innermost datum; rms unchanged (0.186 mag), R_h unchanged
   (280″). lnL improved by 26 from this alone. Limitation noted: inside 10″ the tracer is
   a flat core by construction; the HST star counts in the oMEGACat catalogue reach ~1″
   and should replace Trager there (follow-up).
2. **Second moments vs dispersions.** oMEGACat's `sigma_los` is fitted jointly with a
   rotation curve (`vlos`, `theta_0` per annulus), i.e. it is a dispersion about the
   rotating mean; the spherical Jeans model predicts the full second moment. Added
   `BinnedProfile.streaming2 = ⟨v̄²⟩ = v_rot²/2` (annulus average of v_rot sin φ), read
   from our `los_rotation` product on the same bins; the model compares
   `√(σ²_model − ⟨v̄²⟩)`. Size: 12 % of σ² at 6–16″ (v_rot ≈ 10 km/s), 10–12 % at
   190–260″, <5 % between. **PMs: the catalogue PMs are "locally corrected", and the mean
   tangential PM per annulus is 0 ± 0.005 mas/yr at every radius** — the local correction
   removed the rotation, so the PM streaming term cannot be measured from these data.
   Expected size ~2 % in σ_T at 100–300″ (v_rot,PM ≈ 3–4 km/s); an external rotation
   curve is needed. Open item.
3. **MUSE and HST do not measure the same stars.** MUSE targets giants (median F625W 17.4,
   90th percentile 17.9); the HST PM sample is 3 mag fainter (median 20.7). Restricting the
   PM sample to the MUSE magnitude range gives σ_PM lower by 4.4 %, 2.7 %, 2.6 % at
   20–60″, 60–120″, 120–250″ (raw error-subtracted std; N = 2839/7207/15118 bright stars):
   energy equipartition. This is the LOS/PM offset the first fit tried to absorb with
   D = 5.17 kpc (2.6σ below both literature distances). Added `s_MUSE ∈ [0.85, 1.15]` as a
   nuisance relative to HST; **distance is now a parameter with the literature prior
   N(5.43, 0.05)** in every family (the spec fixed D in K1, but D is degenerate with the
   LOS/PM ratio that the two samples disagree on, so fixing it would bake the equipartition
   offset into the mass).

Second K1 ML fit (11 parameters): D = 5.41, s_MUSE = 0.955, s_DR2 = 0.88, s_EDR3 = 1.05,
M★ = 3.0e6, M_rem = 1.5e5, M• = 4.2e4, β₀ = −0.14, β∞ = 0.43, r_β = 3.9 pc; χ² still 879
with Nelder–Mead — the 11-D search is unreliable, so a differential-evolution run for K1,
K2-cored and K2-NFW is going in the background before anything is concluded.

### The residual pattern, recorded before the sampler runs

With the corrected model the inner profiles are reproduced to the 2σ level except a
shallow dip at 8–15″ (model 0.73, data 0.78–0.80 mas/yr). **Outside ~150″ all four
instruments say the same thing: the no-DM model falls off too fast.** σ_T and σ_LOS lie
3–4σ above the model at 200–300″ while σ_R lies 3σ below it (the fit is pushing β∞ up
to lower σ_T and σ_LOS, and overshoots σ_R); Gaia DR2 sits 3–4σ above at 700–2000″ and
Gaia EDR3 (1 % errors) 8–12σ above at 1300–2500″. Whether that is mass beyond the light
(the K2 question), a tracer profile that is too steep beyond 500″ (Trager's outer points
vs the star counts — the light profile was validated only to 66 pc = 2500″), or a
mass-segregation/M/L gradient, is exactly what the evidence comparison and the
injection tests must decide. No conclusion is drawn here.

Figures: `plots/fit_K1_ml_profiles.png` (current model), `plots/fit_K1_ml_smin10_profiles.png`
(intermediate). Results: `results/fits/*.npy` (ML vectors), `results/fits/K1_de.log`.

### Outer tracer checked against star counts; runs launched

Could the outer residuals be a tracer profile that is too steep? The Vasiliev & Baumgardt
EDR3 members (P > 0.9, G < 19, 70k stars) anchored to the MGE over 8–30′ track it to
**±0.1 mag from 12′ to 36′ (700–2200″)** — the radii where both Gaia dispersion profiles
lie above the K1 model. Inside 10′ Gaia counts fall below the MGE (crowding
incompleteness, known); the member catalogue ends at ~40′. So the light model is not what
makes the no-DM model fall too fast out there. (`fit_mge_projected(..., sigma_range_arcsec
= (10.5, 3000))`; check script inline, numbers above.)

`ocen fit --family {K1,K2-cored,K2-nfw} [--mock-from x.npy] [--n-live 400 --dlogz 0.5]`
writes `results/fits/<label>/{posterior.ecsv, summary.json, profiles.npz, run.yaml,
ml_x.npy}` (`results/` is git-ignored; numbers go in this journal). End-to-end smoke test
on mock data with a 1200-call budget: `tests/test_cli_fit.py` (5.7 s). Launched in the
background: **K1** (11 parameters) and **K2-cored** (13), 400 live points, dlogz 0.5, seed
42; K2-NFW and the injection runs (`--mock-from results/fits/K1_noDM/ml_x.npy` fitted with
K2: the no-DM false-positive test of spec §13 item 8) follow once K1 finishes.

### The rotation term, properly: ω Cen's PM rotation is not small

Looking for the definition of the Gaia dispersion profiles in the Vasiliev & Baumgardt
(2021) readme turned up what the PM streaming term needed: their `profiles/` hold "the
rotational PM component" alongside the dispersion (percentiles vs radius; our processed
product already carried it as `vrot_pm`). The mean tangential PM is **0 at the centre,
0.12 mas/yr at 144″, 0.21 at 288″, peaks at 0.25 mas/yr near 430–580″ and is down to
0.02 mas/yr by 2000″** — i.e. ⟨μ̄_T⟩²/σ² = 5 % at 144″, 18 % at 288″, 29–30 % at
430–580″, 25 % at 720″, 10 % at 1150″. Every published dispersion we fit is measured
about the rotating mean (HST: locally corrected PMs; Vasiliev: joint rotation+dispersion
fit; MUSE: rotation curve per annulus), so the σ² + v̄² correction is a first-order effect
in the very region where the K1 model was "falling off too fast".

Implemented (`likelihood.py`): `pm_rotation_curve(R)` interpolates the Vasiliev curve;
`hst_pm_tangential` gets `streaming2 = ⟨μ̄_T⟩²` at the bin median radius, `gaia_edr3_pm`
(1-D combined) gets `⟨μ̄_T⟩²/2`, `hst_pm_radial` none (rotation is tangential), MUSE keeps
v_rot²/2. **Baumgardt+ 2019 (Gaia DR2) is left uncorrected** because whether its
dispersions were fitted about a rotating mean cannot be settled from the local files; if
they were, the model is biased high there by up to v_rot²/2 ≈ 10 % of σ² at 700″ — flagged
in the product note. Using ⟨μ̄_T⟩² for ⟨μ̄_T²⟩ is a lower bound when the rotation axis is
inclined. On the 1-D question: the fitted Gaia scales are 0.94 (DR2) and 1.04 (EDR3), not
≈1.4, so both profiles are 1-D dispersions as assumed.

Effect on the K1 maximum likelihood (Nelder–Mead, 10 starts): **lnL −88 → +82.5**, total
χ² 879 → 536 for 126 points: MUSE 137 → 51/29, Gaia DR2 82 → 12/9, HST tangential
179 → 103/40; HST radial 148/40 and Gaia EDR3 222/8 remain poor. Parameters moved to
M★ = 3.35e6, M_rem = 4.5e4 at a_rem = 0.8 pc, M• = 4.0e4, β₀ = −0.14, β∞ = 0.27,
r_β = 3.4 pc, s_MUSE = 0.98, D = 5.44 kpc. The stale nested runs (no PM term) were killed
and their outputs deleted; K1, K2-cored and K2-NFW relaunched with the corrected
likelihood (400 live points, dlogz 0.5, seed 42). Figure: `plots/fit_K1_ml_profiles.png`.

### The tracer is stars, not light: composite profile from HST counts

The persistent 10–20″ dip (model 0.74, data 0.78–0.80 mas/yr in both PM components) sent
me back to the light model. Star counts from the oMEGACat catalogue (`selection_hq_f625w`,
pixel-based centre, 1–120″) anchored to the Trager MGE over 30–100″ lie **0.2–0.3 mag
below it inside 20″ for every cut F625W < 18, 19, 20** (N = 47–560 per bin at < 19). The
agreement across cuts rules out incompleteness; the V-band light is boosted by a few bright
giants in the core, while the Jeans equation needs the number density of the stars whose
velocities are measured. A flatter tracer core raises the predicted central dispersion —
the sign of the dip.

`light_model.star_count_profile(mag_cut=19)` (Poisson weights relative to 100 stars),
`composite_profile(outer, inner, r_switch=25″, anchor=(30, 100)″)` (zero-point by weighted
mean offset; tested on a synthetic Plummer that splices exactly), `load_tracer_profile
('trager'|'composite')`. The composite MGE (12 Gaussians, smallest 14.5″ although 1.7″ was
allowed — the counts *are* flat inside 15″; rms 0.178 mag; R_h 279.8″ vs 280.4″) has
I(1″)/I(30″) = 1.30 against Trager's 1.61.

K1 maximum likelihood with the composite tracer: **lnL 82.5 → 150.7**, χ² 536 → 399/126;
HST radial 148 → 75/40, tangential 103 → 50/40, MUSE 51 → 44/29, Gaia DR2 13/9; Gaia EDR3
still 217/8. The dip is gone (figure `plots/fit_K1_composite_ml_profiles.png`). The
parameters moved where the physics says they should: **M• → 9e2 (prior floor: no IMBH
needed by the dispersion profile)**, M_rem = 8.9e4 at a_rem = 0.3 pc (the prior's lower
edge — a compact dark-remnant core is the fit's preferred way to supply the central mass;
this M•/M_rem degeneracy is the well-known one), M★ = 3.46e6, β₀ = 0, β∞ = 0.28,
r_β = 5 pc, s_MUSE = 0.97, s_DR2 = 0.94, s_EDR3 = 1.05, **D = 5.51 kpc** (now 1.6σ above the
5.43 prior and close to oMEGACat's kinematic 5.494 ± 0.061). Remaining structure: HST σ_R
falls below the model by 3–4σ in the outermost three bins (200–300″), and the EDR3 profile
disagrees in shape at the 1 % level of its percentile errors (below the model at 300–700″,
above beyond 1300″).

`ocen fit --tracer {composite,trager}` (composite is the default). Six nested runs are now
in progress: K1/K2-cored/K2-NFW with the Trager tracer (kept as the systematics variant)
and with the composite tracer.

**Correction (same day, before any run finished):** `composite_profile` interpolated the
Trager profile through *all* its points when setting the splice zero-point, including the
weight-0.03 outliers that sit 0.5–0.9 mag off the curve at 42″, 67″, 93″; the offset was
biased by ≈0.25 mag (inner counts placed too bright). Caught by the regression test
(`test_star_counts_below_trager_light_in_the_core` measured a 0.00 mag deficit where the
diagnostic had 0.2–0.3). Fix: only outer points with weight ≥ 0.5 enter the anchor; a
synthetic test with planted low-weight outliers now guards it. The three composite nested
runs were killed, their outputs deleted, and relaunched. With the corrected splice the
composite MGE is flatter still (I(1″)/I(30″) = 1.13; smallest Gaussian 7.1″) and the K1
Nelder–Mead maximum is lnL = 109.5 (HST radial 110/40, tangential 94/40, MUSE 46/29,
DR2 13/9, EDR3 218/8; M• = 4e3, M_rem = 7.5e4 at 0.3 pc, D = 5.54). The mis-spliced
profile, 0.25 mag brighter in the core, had fit better (150.7) — the ML numbers from an
11-D Nelder–Mead differ by tens between runs, so the tracer question (which stars, how
segregated) is deferred to the nested-sampling maxima and to a per-sample tracer profile
(the HST PM stars are 2–3 mag fainter than the F625W < 19 count sample). Six runs in
progress; nothing concluded yet.

**Which stars does the composite trace?** Count profiles by magnitude, each anchored to the
composite MGE over 30–100″ (residuals in mag inside 30″): F625W 17–19 within ±0.1;
19–21 within ±0.14; **the actual HST PM sample (hq astrometry + membership, 18.5–21)
within ±0.16 and mostly ±0.05**; only the faintest bin 21–23 shows a +0.2–0.35 deficit
inside 30″ (crowding incompleteness and/or mass segregation — not separable here, and not
used). So the composite tracer describes the stars whose proper motions we fit; the
segregation worry is retired for the PM data. The MUSE giants (F625W < 18) also match to
±0.1 outside 3″.

Sampler progress at 35 min: Trager-tracer runs at 130–150k likelihood calls (Lmax still
below the Nelder–Mead maxima — the 11–13-D posteriors are broad); composite runs just
started. Expect 1–2 h per run with six sharing the machine.

### Sampler efficiency (same day)

After 70 min the MLFriends runs had used 300–500k likelihood calls for ~3000 iterations
(≈150 calls per iteration and rising) and their best likelihoods were still 250–450 below
the optimiser's maxima: the region-based rejection sampler collapses on the needle-like
M_rem–a_rem–M• degeneracy in 11–13 dimensions. Benchmark on the real K1-composite problem,
200 live points, fixed 40k-call budget: MLFriends reached Lmax = −487 (logZ −499) in 1472
iterations; a **slice step sampler** (`ultranest.stepsampler.SliceSampler`, 2·ndim steps,
mixture random directions) reached **Lmax = −152 (logZ −163) in 1305 iterations at the same
cost**. `ocen fit --step-sampler` added (recorded in `run.yaml`); the four runs (K1-Trager,
K1/K2-cored/K2-NFW composite) were restarted with it. The two Trager-tracer K2 runs were
dropped — the Trager tracer is kept only as the K1 systematics variant.

## 2026-09-18 — First posterior: K1 with the Trager tracer (systematics variant)

Finished after 1.55M likelihood calls / 124 min (slice sampler, 400 live points, dlogz 0.5).
**ln Z = 86.76 ± 0.29; max ln L = 126.8; χ² at the best sample 446 / 126 points.** Medians
[16–84 %]: M★ = 2.52e6 [2.33, 2.72], M_rem = 6.5e5 [4.9, 8.3] at a_rem = 5.15 pc [4.78,
5.42], **M• = 4.16e4 [3.90, 4.42]**, β₀ = −0.033 [−0.056, −0.014], β∞ = 0.255 [0.245,
0.266], r_β = 4.0 pc, s_MUSE = 0.982, s_DR2 = 0.953, s_EDR3 = 1.072, D = 5.436 [5.389,
5.483] kpc.

Reading: with the light-weighted (Trager) tracer the no-DM model needs a 4e4 M☉ point mass
*and* a 6.5e5 M☉ extended remnant sphere (26 % of the stellar mass, a_rem ≈ 5 pc ≈ 0.7 R_h,
i.e. not segregated at all), and still leaves χ²/N = 3.5. The ±6 % width on M• is the width
of a misspecified model's posterior, not a measurement: the value is set by the depth of
the artificial 10–20″ dip (the giant-boosted core), exactly the failure mode diagnosed
yesterday. This run is retained as the "tracer = V-band light" systematics variant, not
as a result. Figures: `plots/fit_K1_noDM_posterior_profiles.png`,
`plots/fit_posterior_profiles.png` (to be overwritten by the multi-run comparison).
Composite-tracer runs at the time of writing: K1 remainder 70 % (Lmax 204), K2-cored
remainder 99.9 % (Lmax 224.6), K2-NFW remainder 99.9 % (Lmax 210).

### K1 with the composite tracer (the reference no-DM model)

2.0M calls / 161 min. **ln Z = 160.17 ± 0.54; max ln L = 206.5; χ² 287 / 126** (HST R
64/40, T 77/40, MUSE 50/29, DR2 13/9, EDR3 84/8). Medians [16–84 %]: M★ = 2.61e6 [2.55,
2.97], M_rem = 6.4e5 [4.2, 7.1] at a_rem = 5.1 pc [4.6, 5.3], **M• = 4.65e4 [4.27, 4.98]**,
**β₀ = −0.72 [−0.74, −0.71]**, β∞ = 0.200 [0.195, 0.205], r_β = 1.09 pc [1.03, 1.16],
s_MUSE = 0.973, s_DR2 = 0.947, s_EDR3 = 1.056, **D = 5.489 [5.476, 5.530] kpc**.

Two things to take seriously and one to distrust. (i) The evidence prefers the composite
tracer over the light-weighted one by **Δ ln Z = +73** with the same K1 physics — the
tracer choice is settled by the data, not by argument. (ii) D = 5.49 kpc reproduces
oMEGACat VI's kinematic distance (5.494 ± 0.061) independently of their modelling, with the
equipartition nuisance absorbing the MUSE/HST offset (s_MUSE = 0.973 ± 0.008). (iii) The
sampler's maximum (206.5) is far above every Nelder–Mead maximum I quoted yesterday
(≤ 150.7): those optimiser numbers were not maxima, and the tracer comparison made from
them (mis-spliced vs corrected) is void. The posterior itself sits in a corner the
optimiser never found: **a 4.6e4 M☉ point mass with a strongly tangential core
(β₀ = −0.72) turning radial by 1 pc**, plus an unsegregated 6.4e5 M☉ remnant sphere.
Tangential bias near a central mass is what the cusp needs to keep σ_R flat while the
point mass lifts σ_T and σ_LOS — a legitimate solution of the Jeans equation, but M• is
4× the literature upper limit (Baumgardt+ 2019 ≲ 1.2e4; Häberle+ 2024 ≥ 8.2e3 from fast
stars), and the ±7 % width comes from a model with χ²/N = 2.3. I do not read it as a BH
measurement. It is the K1 family saying: the inner 10″ want more central mass than the
flat-cored tracer plus a diffuse remnant population can give it, and given only a point
mass and an anisotropy profile, this is the best it can do.

Launched: the **no-DM false-positive injection** (K2-cored fitted to a mock drawn from the
K1-composite best sample with the real bins and errors) and a **K1 variant without the
Gaia EDR3 profile** (its 8 points with 1 % model-percentile errors carry 84 of the 287 χ²).
K2-cored and K2-NFW composite runs at 35–38 % remainder (ln Z 184.9 and 171.8 and rising —
already above K1's 160.2; not final).

### Solver defect found from the posterior figure — composite runs invalidated and restarted

The posterior data-vs-model figure for K1-composite showed the model curve spiking to
0.44 and then 0 mas/yr just beyond the last Gaia point (2400–2600″). Probing the solution:
ν σ_r² fell from 0.17 at 65 pc to 6e-12 at 70 pc while ν itself was only 4e-6 of its
central value and r_cut was 208 pc. Cause: the outer Jeans integral was formed as
`cumulative[-1] − cumulative`; with an integrand spanning 20 decades the remainder is
below the 1e-16 floor of the total long before the tracer edge, so σ_r collapsed to zero
around 3 σ_MGE — inside the last EDR3 datum for some parameter values. Fix: accumulate the
integral from the outside inward (`cumulative_simpson` on the reversed grid), which needs
no subtraction. Regression tests: Plummer σ_r at 20–150 a (ν/ν₀ down to 1e-9) to 2e-3 of
the closed form, and a two-Gaussian MGE whose σ_r must decline smoothly to 0.9 r_cut.
All 18 Jeans tests pass.

Consequence, measured on the saved best samples: **K1-composite ln L 207.4 → 138.7** — the
sampler had found and exploited the artefact (the spurious upturn put the model at 0.239
mas/yr on the 2400″ EDR3 point whose datum is 0.228 ± 0.003; the correct model gives
0.171, a 19σ shortfall). K1-Trager: 127.26 → 127.26, unaffected at its maximum. So the
K1-composite posterior, its Δ ln Z = +73 over Trager, the comparison figures and the
injection/no-EDR3 runs launched from its best sample are **void**; they are quarantined in
`results/fits/_invalid_prefix_solver/` and the four composite runs restarted with the
fixed solver (K1, K2-cored, K2-NFW, K1-noEDR3). The Trager K1 run is retained; its
posterior will be regenerated too once the machine is free, since other samples than the
maximum could have touched the artefact.

Lesson recorded: the evidence comparison was about to be made on a solver that passed
its unit tests (closed forms at ≤ 10 a, JamPy and AGAMA at ≤ 30 pc) but had never been
exercised where ν σ_r² is 12 decades below its peak. The figure caught it; the test suite
now does.

## 2026-09-18 — Toolbox: JamPy and AGAMA-DF engines behind the same likelihood

Decision (user, this morning): the modelling must keep several dynamical methods side by
side and be tested like-for-like against published ω Cen results before any conclusion.
Implemented `kinematics/backends.py`; every engine exposes `projected_moments(R)` and plugs
into `ProfileLikelihood` unchanged, selectable with `NoDarkMatterModel(backend=...)`,
`DarkMatterModel(backend=...)` and `ocen fit --backend {jeans,jam,agama}`.

* **JamBackend** — JamPy 9.0.2 `jam_sph_proj` (imported, never vendored). Non-Gaussian
  mass components (remnant Plummer, gNFW/Burkert haloes) are projected numerically
  (`project_density`, exact on a Plummer to 1e-6) and fitted with a projected NNLS MGE
  (`mge_from_component`; enclosed mass to 0.3 % from 0.1 to 1000 pc). Agreement with our
  solver on stars + BH + remnants + cored halo, β(r) = −0.2 → 0.4: **8×10⁻⁴** on all three
  projections. Cost 0.6 s per model (3 tensors) — a cross-check engine, not a sampler engine.
* **AgamaDFBackend** — AGAMA 1.0.152 QuasiSpherical (Cuddeford–Osipkov–Merritt) DF for the
  tracer density in the total potential (Multipole from our density; agrees with our M(<r)
  and v_circ to 1e-4), projected moments from `GalaxyModel.moments` on 2-D points. Constant
  β: **1.4×10⁻³** agreement with the Jeans solver. Its anisotropy family differs (β → 1
  beyond r_a), so the K1/K2 families get parameters (β₀, r_a) in place of (β₀, β∞, r_β)
  when this engine is chosen. 0.2 s per model.
* **`density_check`** — the DF engine's realised/input tracer density. Where a positive DF
  exists it is 1 (constant β ≤ 0: within 2 %); for Osipkov–Merritt with r_a = 5 pc in the
  Gaussian tracer it departs by 3–35 % between 0.3 and 20 pc — the Jeans equation still
  returns a solution there but no positive DF reproduces the tracer (An & Evans 2006). This
  is precisely why the 2–6 % Jeans-vs-DF difference appeared in that case: physics, not a
  bug, and now a test and a diagnostic.

Two bugs caught on the way: the numerical projection carried a spurious factor R (MGE of
the Plummer was nonsense, 2 % of the mass inside 1 pc) and AGAMA's `moments` returns
⟨v²⟩, not Σ⟨v²⟩. `report.engine_crosscheck(label)` evaluates a run's best sample under the
Jeans and JamPy engines and tabulates χ² per dataset. Tests: `tests/test_backends.py` (5)
plus two in `test_fit.py`; suite green.

Next in this thread: the like-for-like literature presets (Watkins+ 2013, oMEGACat VI,
Baumgardt & Hilker 2018, Baumgardt+ 2019) as `ocen fit --preset ...` with the published
numbers printed beside ours in the comparison report; an AGAMA-DF K1 fit as the
positive-DF counterpart of the Jeans K1.

### Like-for-like literature presets

`ocen fit --preset {watkins2013, omegacat6, baumgardt2018, imbh_limit}` runs our pipeline
under a published analysis's assumptions and prints our posterior beside the published
number with a pull in σ. Family options added for this: `fixed={...}` (parameters held and
removed from the vector), `constant_beta`, `beta0_max`, `distance_kpc` with `fix_distance`.

| preset | their assumptions we impose | compared quantity |
|---|---|---|
| watkins2013 | HST PMs only, D = 4.59 kpc, constant β (prior up to +0.5), constant M/L, no remnants/BH, Trager tracer | M/L_V = 2.71 ± 0.05 (L_V from Harris V_t = 3.68, E(B−V) = 0.12, R_V = 3.1 → 8.6e5 L☉ at 4.59 kpc), β = 0.10 ± 0.02 |
| omegacat6 | HST + MUSE, composite tracer, **flat** D prior 4.5–6.5 | D = 5.494 ± 0.061 |
| baumgardt2018 | HST + MUSE + Gaia DR2, D fixed 5.24 kpc (**UNVERIFIED** — to check in their Table 1) | M_total = 3.55e6 |
| imbh_limit | HST + MUSE, remnants free | M• 95 % upper limit vs 1.2e4 (van der Marel & Anderson 2010, 1σ) |

Caveat stated in the code: our PMs are oMEGACat's, not the Watkins+ 2013 sample, so that
preset tests the modelling under their assumptions, not a re-reduction. An attribution slip
caught while checking `docs/LITERATURE_BASELINES.md`: the 1.2e4 IMBH limit is van der Marel
& Anderson 2010, not Baumgardt+ 2019 (whose result is "no IMBH, 4.6 % in stellar BHs");
preset renamed. The presets will run once the five posterior jobs release the CPUs.

**Engine cross-check on real data (K1-Trager best sample):** our solver vs JamPy, same
parameter vector — χ² per dataset 78.6/77.7/46.5/14.4/228.5 (ours) vs 80.5/76.9/46.6/14.4/
228.3 (JamPy); ln L 127.26 vs 126.50. An independent implementation reproduces the
likelihood on the real profiles to within one unit of ln L; the report now includes this
table for every finished run.

## 2026-09-18 — First fixed-solver posterior: K1 composite without Gaia EDR3

1.18M calls / 97 min. **ln Z = 210.30 ± 0.42; χ² at the best sample 137 / 118 points**
(HST R 38/40, T 47/40, MUSE 39/29, Gaia DR2 14/9). Medians [16–84 %]: M★ = 2.81e6 [2.62,
2.98], M_rem = 4.0e5 [2.8, 5.4] at a_rem = 4.3 pc [3.7, 4.8], **M• = 3.8e4 [3.4, 4.2]**,
β₀ = −0.035 [−0.065, −0.012], β∞ = 0.240 [0.230, 0.252], r_β = 3.6 pc, s_MUSE = 0.983,
s_DR2 = 0.950, **D = 5.435 [5.389, 5.481] kpc** (s_EDR3 unconstrained, as it must be).

Reading, carefully: **without the Vasiliev EDR3 profile, the no-dark-matter model is an
acceptable description of everything else** — HST PMs to 300″, MUSE LOS, and the Gaia DR2
profile to 1800″ (46 pc) — with χ²/N = 1.16 and no residual pattern of note. What it still
wants is a central mass of ~4e4 M☉ in a point (M• 3.4–4.2e4 at 68 %) on top of a 4e5 M☉
remnant sphere; the composite tracer did not remove that requirement, so yesterday's
attribution of the point mass to the light-weighted tracer was only part of the story.
The literature disagreement (van der Marel & Anderson 2010 ≲ 1.2e4; Häberle+ 2024
≥ 8.2e3) stands and is now the first item for the like-for-like presets and the engine
cross-checks. Distance: 5.435 kpc with the N(5.43, 0.05) prior — consistent with both
literature values, not independent of the prior.

Consequence for the DM question: whatever preference for K2 the full-data runs show will
rest on the EDR3 profile (8 thinned points, 1 % model-percentile errors, 300–2400″). A
**K2-cored run without EDR3** was launched to test exactly that. Figures:
`plots/fit_K1_noDM_composite_noEDR3_posterior_profiles.png`; JamPy cross-check of the best
sample in `results/fits/comparison.md`.

**K1-Trager v2 (fixed solver) = K1-Trager v1.** ln Z = 87.21 ± 0.32 vs 86.75 ± 0.16 before
the fix; every parameter median agrees to well within a third of its 68 % interval (M• 4.15e4
vs 4.16e4, D 5.435 vs 5.436, s_EDR3 1.072 vs 1.072); χ² 446.0 vs 445.6. Two conclusions: the
Trager-tracer posterior never touched the far-out artefact, and two independent slice-sampler
runs of the same 11-D problem reproduce ln Z to 0.5 — the run-to-run scatter to keep in
mind when reading Δ ln Z between families. `K1_noDM_trager_v2` is the canonical Trager
run from here; `K1_noDM` is kept as the reproducibility twin.

## 2026-09-18 — K1 composite (full data): the no-DM reference posterior

1.62M calls / 136 min. **ln Z = 134.91 ± 0.47; χ² 355 / 126**, of which **Gaia EDR3 alone
218 / 8**; the rest is the clean fit of the no-EDR3 run (HST R 38/40, T 46/40, MUSE 39/29,
DR2 14/9). Parameters are the no-EDR3 ones to within their intervals: M★ = 2.81e6 [2.55,
2.99], M_rem = 3.9e5 [2.7, 6.0] at 4.2 pc, **M• = 3.8e4 [3.4, 4.2]**, β₀ = −0.046,
β∞ = 0.233, r_β = 3.3 pc, s_MUSE = 0.984, s_DR2 = 0.950, **s_EDR3 = 1.065 ± 0.007**,
D = 5.432 [5.387, 5.479].

So the EDR3 profile does not move the no-DM solution — it cannot be accommodated by it. The
K1 family has one lever for the outskirts (β∞ and the remnant sphere) and it is already used
by the HST outer bins; the EDR3 points at 1300–2400″ demand more dispersion than any K1
model gives at those radii (data 0.23–0.28 mas/yr against model 0.17–0.24 after the
rotation term), and the nuisance scale pinned at 1.065 says the same thing over the whole
300–2400″ range. Tracer comparison at equal physics: **Δ ln Z (composite − Trager) =
+47.7 ± 0.6** — the composite tracer is decisively preferred, a real result this time
(both runs on the fixed solver, reproducibility 0.5).

Launched the **false-positive injection**: mock data drawn from this run's best sample
(M★ = 2.76e6, M_rem = 5.5e5 at 4.7 pc, M• = 4.0e4, β₀ = −0.045, β∞ = 0.23, real bins and
errors, seed 7) fitted with K2-cored. If K2 reports a non-zero M_DM(<100 pc) here, the
real-data K2 preference is not to be believed at that level. K2-cored and K2-NFW on the real
data are in the bulk of their posteriors (ln Z 171 and 163 and rising against K1's 134.9),
K2-cored-noEDR3 has started. Report: `results/fits/comparison.md`; figures
`plots/fit_posterior_profiles.png`, `plots/fit_K1_noDM_composite_posterior_profiles.png`.

### Like-for-like preset 1: Watkins et al. 2013 assumptions on the oMEGACat PMs

`ocen fit --preset watkins2013` (HST PMs only, D = 4.59 kpc, constant β with prior up to
+0.5, constant M/L, no remnants, no BH, Trager tracer; 2 parameters; 153k calls, 11 min):
**M/L_V = 2.539 [2.536, 2.542] vs 2.71 ± 0.05 published; β = 0.203 [0.201, 0.205] vs
0.10 ± 0.02.** But χ² = 1638 for 80 points: a two-parameter constant-anisotropy model
cannot describe the 2023 HST profiles (1 % errors, 40 bins per component), so the
posterior widths and the "σ pulls" the report printed (−3.4σ, +5.1σ) mean nothing — the
report now says so whenever χ²/dof > 2. What the comparison does establish: under their
assumptions our pipeline lands within 6 % of their M/L_V (2.54 vs 2.71; part of that is
D-dependence, M/L ∝ D at fixed σ_PM in mas/yr, and their tracer/geometry — axisymmetric
with inclination 50° and rotation — is not ours), and finds mildly radial anisotropy as
they did, twice as strong (a spherical β with rotation removed is not their global
axisymmetric β). No pipeline defect indicated; no precision claim possible with this
preset. Output `results/fits/preset_watkins2013/`.

## 2026-09-18 — K2-NFW (composite tracer, full data): first dark-matter posterior — read with the controls pending

2.23M calls / 228 min. **ln Z = 176.65 ± 0.63 vs K1 134.91 ± 0.47: Δ ln Z = +41.7** for
two extra parameters. χ² 260 / 126 (HST R 35/40, T 50/40, MUSE 39/29, DR2 11/9, **EDR3
124/8** — down from 218 but still 15σ-equivalent for 8 points). Posterior: **M_DM(<100 pc)
= 2.35e6 [2.09, 2.60] M☉**, r_s = 680 pc [430, 900] (the prior's upper decade; the halo is
effectively a constant-density background over the data), M★ = 1.35e6 [1.08, 2.15],
M_rem = 1.62e6 [0.92, 1.83] at 5.8 pc, M• = 4.2e4 [3.8, 4.6], β as in K1, s_EDR3 = 1.045,
D = 5.443. Zero samples below M_DM = 1e4: the posterior does not touch the no-DM corner.

What this is and is not. It is the Jeans models saying that the EDR3 dispersion profile at
300–2400″ needs ~2e6 M☉ more mass than the light-plus-remnants can supply, spread over the
whole cluster (r_s ≫ data), and that they will happily take that mass out of M★ (1.35e6 is a
V-band M/L well below 1) and put it into remnants and a halo — the three components are
degenerate where only the outer profile constrains them. It is **not yet** a detection:
(i) the K1-no-EDR3 fit was clean, so the entire preference rests on the 8 EDR3 points whose
1 % errors are model percentiles and whose shape neither family reproduces (χ² 124/8 even
here); (ii) the **false-positive injection** (K2 on a K1 mock) and (iii) the **K2 without
EDR3** are running precisely to test this and finish within the hour; (iv) the M★–M_rem–M_DM
degeneracy needs an external prior (a stellar M/L from the CMD, a remnant fraction from
N-body) before any component can be quoted on its own. Figures:
`plots/fit_K2_nfw_composite_posterior_profiles.png`, `plots/fit_posterior_profiles.png`.

## 2026-09-18 — K2-cored (composite, full data) — and what the halo it wants looks like

2.33M calls / 238 min. **ln Z = 185.92 ± 0.51**: +51.0 over K1 (134.91), +9.3 over K2-NFW
(176.65). χ² 237 / 126 (HST R 34/40, T 52/40, MUSE 38/29, DR2 11/9, **EDR3 102/8**).
Posterior: **M_DM(<100 pc) = 3.7e6 [3.3, 4.1] M☉**, r_s = 730 pc [490, 910], M★ = 2.09e6
[1.72, 2.54], M_rem = 9.8e5 [5.9, 13.0] at 5.4 pc, M• = 4.1e4 [3.7, 4.5], β and nuisances as
in K1, s_EDR3 = 1.044, D = 5.429. Dark fraction 0.3 % at 10 pc, 3 % at 28 pc, 21 % at 57 pc,
56 % at 100 pc.

The shape of this "halo" is the point. With r_s = 730 pc and γ = 0 it is a **uniform
background of ρ ≈ 3.7e6 / (4π/3 · 100³) ≈ 0.9 M☉ pc⁻³** across the entire cluster — two
orders of magnitude above any plausible dark-matter density at ω Cen's position (the local
Galactic halo is ~0.01 M☉ pc⁻³), and, continued to its own scale radius (our fixed
truncation is r_t = 1000 pc), it would hold ~10⁹ M☉. No bound remnant halo can extend beyond
the Jacobi radius (~100–200 pc for a 3–4e6 M☉ cluster at 6.5 kpc from the Galactic centre).
So the sampler is not measuring dark matter here: it is using the only component free to
add mass at large radius without changing the inner profile — a near-harmonic potential — to
lift the model dispersion in the 300–2400″ range where the EDR3 profile refuses to fall.
The cored family beats NFW because a flat core does this with less damage to the inner
profile. Both still fail EDR3 by 100+ χ² for 8 points.

Two things follow. (a) A physical prior is missing in K2: the truncation must be tied to the
Jacobi radius (`CompositeMassModel.jacobi_radius` exists; the spec asks for it) and the halo
mass must be bounded by what the cluster can hold — otherwise "K2" is a fit to systematics.
(b) The EDR3 profile's outer plateau (0.23–0.28 mas/yr at 1300–2400″ against models at
0.17–0.24) is either real (extra-tidal heating, potential escapers, an unbound envelope —
all of which spherical Jeans of a bound tracer cannot describe) or a contamination/systematic
floor in the Gaia dispersion. Either way it is not a bound halo signal. The controls now
finishing (K2 without EDR3, K2 on a K1 mock) will quantify how much of the Δ ln Z is EDR3.

## 2026-09-18 — Control 1: K2-cored without EDR3 — no dark-matter preference

1.36M calls / 142 min. **ln Z = 211.64 ± 0.40 vs K1-no-EDR3 210.30 ± 0.42: Δ ln Z = +1.3**,
i.e. nothing (run-to-run scatter is 0.5; two extra parameters). χ² 134 / 118, the same
clean fit. **M_DM(<100 pc) is unconstrained: median 1.3e5, 16–84 % [4e3, 1.6e6], 47 % of
the posterior below 1e5 M☉**, r_s anywhere in [16, 500] pc; every other parameter is the
K1 value. So the entire Δ ln Z = +51 (cored) / +42 (NFW) of the full-data K2 runs comes
from the 8 thinned Gaia EDR3 points. Without them, the HST + MUSE + Gaia DR2 kinematics
(to 46 pc) are fully described by stars + remnants + a ~4e4 M☉ central mass and put no
lower bound on dark mass; the 95 % upper limit is M_DM(<100 pc) < 3.2e6 M☉ (weak, because
the data end at 46 pc and the halo family is unbounded outward).

**Proposal recorded, NOT implemented (user asked that model changes be discussed first):**
tie the K2 truncation to the cluster's Jacobi radius instead of the fixed r_t = 1000 pc,
and cap the r_s prior accordingly. Numbers: with the Baumgardt orbit (R_peri 1.47 kpc,
R_gc now 6.37 kpc) and a flat 220 km/s rotation curve, r_J = 160 pc today (59 pc at
pericentre) for the K1-like model; the K2-cored posterior halo would push its own r_J to
420 pc. The diff is saved in the session scratchpad (`jacobi_truncation_proposal.diff`),
tests included; it changes the K2 family and therefore waits for a decision.

### Like-for-like preset 2: oMEGACat VI kinematic distance

`ocen fit --preset omegacat6` (HST + MUSE only, composite tracer, **flat** D prior
4.5–6.5 kpc, s_MUSE free; 9 parameters; 1.10M calls / 119 min): **D = 5.71 [5.26, 6.07] kpc
vs 5.494 ± 0.061 published (+0.5σ)**; χ²/dof = 1.2, a good fit. The interval is wide
because the LOS/PM dispersion ratio fixes only the product of the distance and the
MUSE/HST equipartition scale: s_MUSE came out 0.93 [0.88, 1.02], fully degenerate with D
(the K1 fits with the N(5.43, 0.05) prior had s_MUSE = 0.984 ± 0.010). Consistent with
their value; not an independent precision measurement without an external constraint on
the equipartition offset (their analysis treats the samples differently). Remaining
parameters as in K1 (M• = 4.5e4 [3.4, 5.5]). Output `results/fits/preset_omegacat6/`.
Waiting for the last run (injection) before any comparison or decision.

## 2026-09-18 — Control 2: false-positive injection, and the complete-set comparison

**Injection.** K2-cored fitted to a mock drawn from the K1-composite best sample (no dark
matter; real bins and errors; seed 7): 2.04M calls / 224 min, ln Z = 267.58 ± 0.45,
χ² 87 / 126. **Recovered M_DM(<100 pc) = 1.0e4 [2.3e3, 5.9e4]; 95 % upper limit 1.4e5 M☉;
100 % of the posterior below 1e6.** The K2 machinery does not invent a halo from noise:
the real-data K2 values (3.7e6 cored, 2.3e6 NFW) lie a factor 25 above the mock's 95 %
limit. Recovery of the other parameters exposes the expected degeneracy — total mass
inside the data recovered to 7 % (3.07e6 vs 3.31e6 truth) but split differently between
M★ (1.64e6 vs 2.76e6) and M_rem (1.43e6 vs 0.55e6); β₀ −0.15 vs −0.045 (2σ), M• 4.5e4 vs
4.0e4 (1.4σ), D 5.433 vs 5.505 (the N(5.43, 0.05) prior pulls), nuisances within 1σ.

**The complete set** (`results/fits/comparison.md`, `comparison_noEDR3.md`,
`comparison_injection.md`; figures `plots/fit_posterior_profiles_main.png` and one
data-vs-model panel per run; JamPy cross-checks of every best sample within ≤ 2.4 in ln L):

| run | data | ln Z | χ²/N | M_DM(<100 pc) |
|---|---|---|---|---|
| K1 Trager | all | 87.2 ± 0.3 | 446/126 | — |
| K1 composite | all | 134.9 ± 0.5 | 355/126 | — |
| K2 NFW composite | all | 176.7 ± 0.6 | 260/126 | 2.35e6 [2.09, 2.60] |
| **K2 cored composite** | all | **185.9 ± 0.5** | 237/126 | 3.69e6 [3.28, 4.08], r_s 730 pc |
| K1 composite | no EDR3 | 210.3 ± 0.4 | 137/118 | — |
| K2 cored composite | no EDR3 | 211.6 ± 0.4 | 134/118 | unconstrained; 95 % < 3.2e6 |
| K2 cored on K1 mock | all (mock) | 267.6 ± 0.5 | 87/126 | 1.0e4; 95 % < 1.4e5 |

What the set says, without a decision attached:
1. The composite (star-count) tracer is preferred over the light-weighted one by
   Δ ln Z = +48 at equal physics; every run with it fits HST + MUSE + Gaia DR2 with
   χ²/N ≈ 1.1–1.2 and structureless residuals.
2. The only dataset no model fits is the Gaia EDR3 profile (χ² 218/8 for K1, 102/8 for
   the best K2): its 8 thinned points with 1 % model-percentile errors alternate ±5σ
   around even the best K2 curve. All of the K2 preference (+51 / +42) comes from it; with
   EDR3 removed the K1/K2 evidences are equal (Δ = +1.3) and the DM mass is unconstrained.
3. The halo the full-data K2 runs choose (r_s ≈ 700 pc, uniform ρ ≈ 0.9 M☉ pc⁻³, M(<500 pc)
   > 1e8) is not a component bound to a 4e6 M☉ cluster whose Jacobi radius is ~160 pc
   today and ~60 pc at pericentre; in the enclosed-mass figure it is the curve that leaves
   the data region at 50 pc and never turns over.
4. The injection shows the pipeline does not produce such a halo from a no-DM truth with
   these bins and errors; the preference is therefore a genuine statement about the EDR3
   profile's shape relative to a bound-tracer Jeans model — not a false positive of the
   method, and not evidence of a bound dark halo either.
5. Every run wants M• ≈ 4e4 M☉ (3.4–4.6e4) independently of tracer, dataset and family,
   against the ≲ 1.2e4 literature limit; unresolved, and the item the imbh_limit preset
   and the AGAMA-DF engine are for.
6. Distance 5.43 ± 0.05 with the literature prior in every run; 5.71 [5.26, 6.07] with a
   flat prior (degenerate with the equipartition scale). Watkins-2013 preset: M/L_V 2.54 vs
   2.71, β 0.20 vs 0.10, χ²/dof 20 — order-of-magnitude agreement only.

Decisions on the model (K2 truncation/prior; how to treat the EDR3 profile) are the
user's and are pending; the Jacobi-truncation proposal diff remains parked.

**Figure update.** The dark-fraction panel of `plots/fit_posterior_profiles_main.png` is now
logarithmic (`plot_posterior_profiles(..., f_dm_floor=1e-5)`); the K1 runs have f_DM ≡ 0 and
are named in the panel rather than drawn at an arbitrary floor. Colour now follows the model
family and the no-EDR3 variants are dash-dotted, so the same colour never means two models.
The log axis shows what the linear one hid: **the two halo shapes differ by two orders of
magnitude in the core** — at 0.3–1 pc the NFW run carries f_DM = 7e-4–5e-3 against the cored
run's 4e-6–9e-5, and the NFW curve is already flat by 3 pc. Both reach ~0.45–0.56 at 100 pc.
So the cusp/core choice is decided inside 1 pc, where HST gives the tightest data and where
the ~4e4 M☉ point mass sits — a degeneracy to keep in view when the M• question is taken up.

**Dual radius axes on every figure (user request).** `plotting/style.py` now holds
`add_pc_axis(ax, distance_kpc, per_unit_arcsec=1)` and
`add_arcsec_axis(ax, distance_kpc, unit='arcsec'|'arcmin')`; every plot with an angular
radius axis carries a pc scale on top and every plot in pc carries an arcsec scale, each
labelled with the distance used (5.43 kpc for the data figures, the model's own fitted
distance for the fit panels — 5.40 kpc in the K2-cored figure, which is why its top axis
differs slightly from the data figures'). Added where it was missing: the MUSE dispersion
and rotation panels of `omegacat_vi_profiles`, all five panels of every data-vs-model
figure, and the three panels of the posterior-profile figure. Tests
(`tests/test_plotting_axes.py`, 4) check both conversions, the arcmin variant, the
round trip and the distance dependence. Suite 266 passing.

### Defect found while regenerating the figures: mock runs were drawn against the real data

The injection figure showed the mock-fitted model against the **real** profiles, with
χ² = 227 on Gaia EDR3 where the run itself had recorded 4. Cause: nothing in
`summary.json`/`run.yaml` said the run had been fitted to a mock realisation, so
`write_report` reloaded the real data. Fixed:

* `run_nested(..., data_provenance=...)` writes a `data` block — `{"kind": "real"}` or
  `{"kind": "mock", generating_family, seed, source_file, tracer, truth: {...}}` — into both
  `summary.json` and `run.yaml`; the CLI fills it whenever `--mock-from` is used.
* `report.data_for(summary)` rebuilds the mock realisation exactly (same generating family,
  same truth vector, same seed) and every figure title carries `[MOCK data from ..., seed n]`.
  Verified on the existing run: rebuilt χ² = 86.5 against the recorded 86.5, per dataset
  22.9 / 34.4 / 22.2 / 3.1 / 3.9. The 2026-09-18 injection run's provenance was backfilled
  from its launch command and marked `backfilled` in the file.
* Root cause of the wrong data being silently plausible: **`NoDarkMatterModel()` defaults to
  `tracer="trager"` while `ocen fit` defaults to `--tracer composite`**, so a report that
  constructs a family instead of reading the run's label gets a different light model (9 vs
  11 Gaussians). `_family_for` already derives the tracer from the run label; a test now
  pins both defaults so the mismatch cannot drift silently. Whether to make the two defaults
  agree is a model-default change and is left for discussion.

The corrected injection figure shows what the numbers always said: residuals scatter within
±2σ on all five datasets (χ² 23 / 34 / 22 / 3 / 4 = 87 for 126 points) — K2 fitted to no-DM
data reproduces it without inventing a halo. Suite 267 passing.

## 2026-09-18 — Where the constraints come from, and an independent measurement of the outer dispersion

Two questions from the user: make plots that show directly what each dataset constrains,
and check whether the outer signal rests on secure, genuinely hot distant tracers.

### `kinematics/outer_profile.py` — our own PM dispersion from the member catalogue

Per-star radial/tangential PMs about the centre with the **error covariance projected onto
the same directions** (`pmra_error`, `pmdec_error`, `pmra_pmdec_corr`), and a 1-D
maximum-likelihood dispersion per annulus with the per-star errors deconvolved and the mean
fitted simultaneously. Tests cover recovery when the errors exceed the signal, the bias from
a 10 % error underestimate (0.25 → 0.305 mas/yr) and from 3 % contamination (0.25 → >0.4).

Two things had to be right before the numbers meant anything:
1. **Subtract the systemic PM before projecting.** Projecting the absolute PM onto the
   radial/tangential directions turns the cluster's own motion into an azimuthal pattern of
   amplitude |μ_sys|, read as a spurious dispersion of |μ_sys|/√2 ≈ 5.3 mas/yr — which is
   exactly what the first run produced. Our error-weighted systemic PM from members inside
   600″: **(−3.2480, −6.7462) mas/yr**, matching the published (−3.25, −6.75).
2. **Use the authors' quality flag** (bit 2). Without it the inner annuli are inflated by
   35 % at 350″ falling to 0 by 1400″ — the Gaia crowding systematic. With it our
   measurement reproduces the published profile to **≤5 % median, ≤10 % everywhere beyond
   500″** — an independent reproduction of Vasiliev & Baumgardt's profile from their
   catalogue, with our own error budget (0.3–3 % statistical, against their ~1 %).

### Answering the question

| annulus | N (P>0.9 + quality) | σ measured | per-star error (all / G<18.5) | K1 model | K2 model |
|---|---|---|---|---|---|
| 300–369″ | 39 | 0.503 ± 0.057 | 0.25 / 0.03 | 0.574 | 0.563 |
| 689–849″ | 14224 | 0.394 ± 0.004 | 0.35 / 0.13 | 0.400 | 0.394 |
| 1045–1286″ | 14520 | 0.313 ± 0.004 | 0.38 / 0.13 | 0.317 | 0.319 |
| 1286–1583″ | 9414 | 0.281 ± 0.004 | 0.38 / 0.13 | 0.272 | 0.281 |
| 1583–1949″ | 4884 | 0.247 ± 0.005 | 0.37 / 0.13 | 0.234 | 0.251 |
| 1949–2400″ | 2337 | 0.221 ± 0.007 | 0.36 / 0.13 | 0.196 | 0.222 |

(model values are bin-averaged and carry each run's fitted Gaia scale, 1.068 / 1.050.)
χ² over these 10 bins: **K1 = 52, K2 = 16**.

* **Inside ~1200″ (32 pc) the two families are indistinguishable**: the K2/K1 model ratio is
  within 2 % from 2″ to 1000″ (K2 is in fact marginally *lower*, having moved mass from
  stars to halo), far below the data errors. No dark matter is needed or detectable there.
* **Beyond ~1300″ (35 pc) they separate**: +5 % at 1400″, +9 % at 1700″, +16 % at 2100″,
  against data errors of 1.4–3 %. The measured dispersions follow K2 and sit 2.2σ, 2.6σ,
  3.5σ above K1.
* **Secure distant tracers exist**: 9414, 4884 and 2337 members with P > 0.9 and the quality
  flag in the outer three annuli; the expected contamination (Σ(1−P)) is 1–3 %.
* **They are genuinely hot, not deconvolved artefacts**: restricting to G < 18.5, whose
  per-star errors are **0.13 mas/yr — three times smaller than the 0.25–0.28 mas/yr signal**,
  gives the same dispersion to within 10 % (and to within 4 % at 1000–2000″). Inflating all
  errors by 20 % lowers σ by ~12 %, tightening the P > 0.99 cut by ~10 %: neither erases a
  16 % model difference.

So the user's reading is confirmed, with one refinement: the transition is at ~1300″ ≈ 35 pc
≈ 5 half-light radii, not "a few hundred arcsec", and inside it the two models agree to 2 %
rather than merely "fit equally well". What remains open is whether hot tracers at 35–64 pc
mean *bound* mass: at those radii the tidal and escaping populations are exactly what a
bound-tracer Jeans model cannot represent (JOURNAL, K2-cored entry).

Figures: `plots/constraint_map.png` (all datasets in km/s, K1 and K2 curves, deviations from
K1, tracer counts per bin) and `plots/outer_tracer_audit.png` (counts, errors vs signal,
the measurement under six different cuts, robustness). CLI: `ocen plot-constraints`.
Measurement saved to `data/processed/kinematics/ocen_pm_dispersion_ours.ecsv` — **not** used
in any fit; using it in place of the published profile is decision A(d), still open.
Suite 275 passing.

### Perspective and projection effects, done properly (user's concern; Vasiliev & Belokurov 2020 §5.1)

The first measurement subtracted a constant systemic PM. For a body with one 3-D velocity
observed over a degree that is wrong in three ways (Vasiliev & Belokurov 2020, MNRAS 497,
4162, §5.1; van de Ven+ 2006 for ω Cen): the systemic velocity projects onto each star's
*own* tangent basis (perspective contraction −(v_los/D)θ plus the rotation of the
east/north directions across the field), and a star at unknown depth ζ appears slower by
μ_sys ζ, which survives as an apparent dispersion |μ_sys| σ_z(R)/D along the systemic-PM
direction. `kinematics/perspective.py`: `systemic_pm_field` (exact in angle; agrees with
astropy's cartesian-velocity transform to 1e-15 mas/yr), `depth_dispersion` (σ_z from the
tracer MGE along the line of sight). `load_members(exact=True)` now subtracts the projected
field star by star, with the centre value re-estimated once after the field is removed;
`binned_dispersion(depth_tracer=...)` adds the depth term as a per-star variance
σ_depth² cos²(φ−φ_sys) (radial) / sin² (tangential).

Sizes for ω Cen: perspective term 0.104 mas/yr at 0.66°, basis-rotation residual ≲ 0.003,
depth-induced dispersion 0.011–0.031 mas/yr (the tracer's depth at 2000″ is ~23 pc rms).
**Mock test** (`tests/test_perspective.py`): a Plummer sphere moving rigidly with ω Cen's
systemic velocity, projected with exact geometry and true depths — the exact treatment
plus depth term reads σ < 0.01 mas/yr (zero within the 0.02 errors), the exact treatment
without the depth term reads exactly the predicted |μ_sys| σ_z/D, and the naive constant
subtraction reads a spurious signal growing outward; with a 0.25 mas/yr internal
dispersion added it is recovered to 0.02.

**Real data**: naive → exact+depth changes the outer bins by ≤ 0.003 mas/yr (2122″: 0.221 →
0.218; 1722″: 0.247 → 0.245), less than half a statistical error. Why so little: the
perspective term is radial and constant around an annulus, so the free per-annulus mean of
the naive fit absorbed it (mean_R naive −0.062 at 2122″ ≈ the perspective term; after the
exact removal mean_R = −0.007 ± 0.007, i.e. **no net expansion or contraction — a
consistency check the exact treatment passes**), and the depth term is ≤ 0.03 in quadrature.
The systemic PM moves from (−3.2480, −6.7462) to (−3.2479, −6.7451). Conclusions unchanged:
χ² over the 10 bins K1 = 48, K2 = 18; outer three bins +1.9σ, +2.1σ, +3.1σ above K1. The
audit figure now carries the naive treatment as a seventh variant. Saved profile
`ocen_pm_dispersion_ours.ecsv` is the exact+depth version. Suite green.

## 2026-09-18 — Field contamination of the outer members, and rotation

### Estimating and modelling contamination (user's question)

Three handles, all in the member catalogue: Σ(1−P) from the published probabilities
(0.1–1.1 % of the P > 0.9 stars, rising outward); the field-to-member ratio per annulus
(field stars with P < 0.05 outnumber members 6:1 at 1800–2400″ — 22 178 vs 3842); and the
field PM distribution at the same radii (σ ≈ 5–7 mas/yr in two broad components).

Guarding against it: `mixture_dispersion_free` fits, per annulus and per component, a
**cluster + two-Gaussian field mixture to all quality stars with no membership probability
used at all** — cluster N(μ̄, σ² + e_i² + depth term), field widths floored at 1.5 mas/yr
(without the floor EM carves a narrow "field" out of the members' wings and biases σ low,
seen as 0.4–0.5 mas/yr "field" components at 500–1000″). Tests: heavy-contamination mock
(85 % field, errors 0.35, σ = 0.20) recovered without bias over four seeds (0.196–0.222,
mean 0.205) and identical to a full Nelder–Mead maximisation of the 8-parameter likelihood;
a pure-cluster mock returns σ to 0.01 with f < 0.02. A first version that took the field
template from P < 0.05 stars was discarded: that template has a hole at the cluster PM.

Result (quality stars, exact systemic field, depth term):

| annulus | N all | P>0.9 cut | P>0.99 cut | **mixture** | field frac. |
|---|---|---|---|---|---|
| 300–500″ | 1200 | 0.463 | 0.450 | 0.460 ± 0.008 | 0.06 |
| 500–700″ | 11366 | 0.430 | 0.418 | 0.433 ± 0.003 | 0.05 |
| 700–1000″ | 29513 | 0.376 | 0.358 | 0.375 ± 0.002 | 0.11 |
| 1000–1400″ | 30545 | 0.313 | 0.294 | 0.308 ± 0.002 | 0.28 |
| 1400–1800″ | 20101 | 0.263 | 0.242 | 0.252 ± 0.003 | 0.62 |
| 1800–2400″ | 26519 | 0.229 | 0.201 | **0.195 ± 0.004** | 0.89 |

Inside 1400″ the P > 0.9 sample is clean; at 1400–1800″ contamination inflates it by 4 %;
**at 1800–2400″ by 15 % (0.229 → 0.195)**. The bright (G < 18.5) subsample, whose cluster
peak is resolved, agrees with the mixture to 1–3 % (0.207 in the last annulus); a 3 mas/yr
field-width floor changes nothing.

### Rotation (user's question)

In the **fits**: taken into account at the second-moment level. Every published dispersion
is about the rotating mean (MUSE: per-annulus rotation curve; HST: locally corrected PMs;
Vasiliev: joint fit), so the likelihood adds ⟨v̄²⟩ back — v_rot²/2 for MUSE, μ_rot² for the
HST tangential component, μ_rot²/2 for the 1-D Gaia profiles, from the Vasiliev & Baumgardt
rotation curve — and compares total second moments, which is what a spherical Jeans model
predicts regardless of how they split into ordered and random motion. **Not** modelled:
rotation's dynamical role (the ε ≈ 0.17 flattening, axisymmetry, inclination); that is what
the JamPy axisymmetric back-end is for. In **our own measurement**: the free per-annulus
mean tangential PM is the rotation — it reproduces the published curve (−0.238, −0.222,
−0.155, −0.091, −0.033, −0.020 vs 0.249, 0.218, 0.158, 0.098, 0.051, 0.020 mas/yr) — and its
azimuthal variation is ≤ 0.05 mas/yr, adding ≤ 0.8 % to σ_T². The comparison with the
models therefore needs μ_rot²/2 removed from the model second moment; the first
constraint-map comparison omitted this and was biased at 300–700″, where μ_rot²/2 is
13–15 % of σ².

### What the corrected comparison says (6 annuli, models bin-averaged, rotation term in)

| annulus | mixture σ | K1 (s = 1) | K2-cored (s = 1) |
|---|---|---|---|
| 300–500″ | 0.460 ± 0.008 | 0.476 | 0.474 |
| 500–700″ | 0.433 ± 0.003 | 0.386 | 0.384 |
| 700–1000″ | 0.375 ± 0.002 | 0.345 | 0.346 |
| 1000–1400″ | 0.308 ± 0.002 | 0.287 | 0.295 |
| 1400–1800″ | 0.252 ± 0.003 | 0.235 | 0.251 |
| 1800–2400″ | 0.195 ± 0.004 | 0.192 | 0.219 |

Best single Gaia scale: K1 s = 1.082 → χ² = 86 / 6; K2 s = 1.063 → χ² = 245 / 6. With
contamination modelled, **the outer two annuli follow the no-DM curve and lie 5–6σ below
the cored-halo curve** that the full-data K2 fit had put through the published profile's
outer points. What remains is a *uniform* 8 % excess of the Gaia dispersions over the
HST-anchored K1 at 500–1800″ (13–47 pc), robust to the bright-star cut, the field-width
floor and the error scale, plus a 7σ deficit in the innermost, crowded annulus (300–500″,
1200 stars). Equipartition predicts the opposite sign for bright Gaia members, so this 8 %
is either a Gaia systematic in the overlap region or real extra mass at 2–7 r_h — a
compact component, not the r_s ≈ 700 pc background the unconstrained K2 chose.

Saved: `data/processed/kinematics/ocen_pm_dispersion_mixture.ecsv` (with the rotation
term per annulus). Not used in any fit. The figures now carry the mixture profile as "our
measurement (field modelled)" and the audit shows it beside the P-cut variants. Decision
A(d) — replace the published EDR3 profile by this measurement in the fits — is now
supported by the evidence and would change the K1/K2 comparison; it remains the user's.

### Figures: both versions, and the contamination model made visible (user request)

`ocen plot-constraints` now writes five figures: `constraint_map.png` (Gaia points from the
cluster+field mixture) and `constraint_map_pcut.png` (Gaia points from the P > 0.9 + quality
members, no contamination model); `outer_tracer_audit.png` (reference = P > 0.9) and
`outer_tracer_audit_mixture.png` (reference = mixture, P-cut shown as the variant); and
`contamination_model.png`, which shows the modelling itself: for 700–1000″, 1400–1800″ and
1800–2400″ the radial-PM histogram of all quality stars with the fitted cluster and field
components (wide view, log; zoom on the peak, linear), the field stars expected under the
cluster peak within ±3σ, and the field fraction and P-cut bias versus radius.

The zoom panels are the argument in one picture: under the cluster peak there are 684
field stars at 700–1000″ (2.7 % of the peak), 2651 at 1400–1800″ (35 %), and **5031 at
1800–2400″ — 174 % of the cluster stars in the same PM window**. A membership cut selects
those stars as members; only a model of the field can remove their contribution to the
dispersion. Side-by-side, the two constraint maps show the outermost Gaia point moving
from +12 % above the no-DM curve (P-cut) to −6 % below it (field modelled).

## 2026-09-18 — The background model was wrong; two dimensions fix it

The user rejected the contamination result on sight: "the model of background+cluster
reproduces the data poorly around the peak, especially in the outer regions". Correct, and
the first response was a plotting artefact hiding a real one. Three findings, in order.

**1. The drawn curve was not the fitted model.** I plotted the cluster as one Gaussian of
width sqrt(sigma^2 + <e^2>), while the model is a sum over stars of Gaussians of width
sqrt(sigma^2 + e_i^2) with e_i spanning 0.07-0.86 mas/yr. The exact per-star prediction
(each star's Gaussian integrated over the bin) fits far better than the drawn curve:
central bin at 1800-2400 arcsec, data 252, exact model 252, drawn approximation 186; chi2
over 60 bins 162 vs 460. Figures now draw the exact expectation, with residual panels.

**2. The two-Gaussian field really was a bad model, and it biased sigma low.** Even exactly
drawn it left chi2/bin = 8.2 over the full range at 1800-2400 arcsec. The Galactic field
there is not two Gaussians: measured in 0.5 degree cells it has **kurtosis 3-17** and
**unequal widths, 3.1 mas/yr in alpha* and 2.0 in delta** (MAD-based, beyond 1 deg), with
only 0.13 mas/yr rms cell-to-cell variation in its centre -- so the user's suggestion of
building the annulus from locally Gaussian pieces would not have worked either: the field
is non-Gaussian *locally*, not by spatial mixing.

**3. Projection was manufacturing the shoulders.** Projecting a field with two unequal
widths onto each star's own radial direction mixes them around the annulus. The fix is not
to project at all: fit in two dimensions.

**What is now implemented.** `selection/field_template.py` fetches an independent Gaia DR3
field from an annulus **outside** the cluster (0.75-1.6 deg, 226,474 stars, 192,054 after
quality cuts, 8.5 per arcmin^2; `ocen fetch-field-template`, WSDB, provenance and query in
the product's metadata) and builds `field_density_2d`, the smoothed empirical density of
(mu_alpha*, mu_delta) after the same exact systemic-field subtraction, using only stars
beyond 1 deg (inside that the cluster's own outskirts raise the core fraction by 20 %).
`kinematics/outer_profile.dispersion_2d` fits, per annulus and with **no membership cut**:
cluster = 2-D Gaussian with covariance R(phi_i) diag(sigma_R^2 + depth, sigma_T^2) R(phi_i)^T
plus the star's error covariance; field = the empirical density with only its normalisation
free. Maximised by multi-start simplex over (sigma_R, sigma_T) with EM for the means and
the field fraction inside, profile-likelihood intervals, bounded so the fit cannot run away
to the field solution (it did, in the outermost bin, before the bound).

**Fit quality** (projected onto mu_R for display, chi2 per bin, full range / peak):
700-1000 arcsec 1.71 / 1.39; 1400-1800 1.87 / 1.27; **1800-2400 1.73 / 1.21**, against
**8.22 / 2.71** for the two-Gaussian field. Tests: a mock annulus with 85 % field of the
right shape recovers sigma_R = 0.25 and sigma_T = 0.30 to 0.03; with the cluster at 4 % of
the stars it still recovers them and does not run away; and a 1-D two-Gaussian fit to the
same mock is biased low by more than 10 %, reproducing the artefact.

**The measurement** (quality stars, no P cut, exact perspective, depth term):

| r [arcsec] | N | cluster N | f_field | sigma_R | sigma_T | sigma_1D | P>0.9 value |
|---|---|---|---|---|---|---|---|
| 521 | 2702 | 2605 | 0.04 | 0.491 ± 0.010 | 0.441 ± 0.009 | 0.466 ± 0.007 | 0.460 |
| 773 | 15526 | 14341 | 0.08 | 0.412 ± 0.008 | 0.382 ± 0.008 | 0.397 ± 0.006 | 0.393 |
| 1156 | 19185 | 14668 | 0.24 | 0.312 ± 0.006 | 0.323 ± 0.006 | 0.317 ± 0.004 | 0.313 |
| 1423 | 17152 | 9527 | 0.45 | 0.279 ± 0.006 | 0.293 ± 0.006 | 0.286 ± 0.004 | 0.280 |
| 1760 | 16707 | 4961 | 0.70 | 0.235 ± 0.009 | 0.266 ± 0.011 | 0.251 ± 0.007 | 0.245 |
| 2174 | 19968 | 2337 | 0.88 | 0.211 ± 0.008 | 0.231 ± 0.009 | 0.221 ± 0.006 | 0.218 |

**Corrections to what I wrote earlier today.** (i) The claim that contamination inflates the
outermost dispersion by 15 % is **withdrawn**: it came from the two-Gaussian field absorbing
part of the cluster peak. Properly modelled, contamination changes the P > 0.9 values by
**1-3 %** (upward, because a membership cut also truncates the cluster's own velocity
wings). (ii) The subsequent claim that the outer points then follow the no-DM curve and lie
5-6 sigma below the cored halo is also withdrawn -- that was the same artefact.

**What the corrected measurement says.** Against the bin-averaged models with each run's
fitted Gaia scale: residuals +4.1, +4.5, +3.6, +1.9, +2.2, +4.6, +2.7, +4.1 sigma from 520
to 2170 arcsec for K1, chi2 = 107/10; for K2-cored 103/10. **With one free Gaia scale the
no-DM model fits best: s = 1.052, chi2 = 25/10, against 34/10 for the cored halo.** So the
Gaia dispersions sit ~5 % above the HST-anchored models at every radius from 13 to 57 pc,
and that offset is flat, not rising: it is an instrumental/selection scale, not a mass
gradient. No preference for dark matter survives in this dataset once the field is modelled
properly.

**A new physical result from the same fit.** The anisotropy reverses: sigma_T/sigma_R =
0.85-0.93 (radial) at 430-780 arcsec, 1.03-1.13 (tangential) beyond 1150 arcsec, with the
turnover at ~1000 arcsec = 27 pc = 3.6 r_h. Tangential anisotropy in the outskirts is the
expected signature of tidal stripping preferentially removing radial orbits, and it is
measured here on the same stars, in the same fit, without a membership cut.

Products: `ocen_field_template_dr3.ecsv` (the field, with its WSDB query and cuts),
`ocen_pm_dispersion_mixture.ecsv` (the 2-D measurement, per-component). Figures:
`contamination_model.png` (now: exact model, residual panels, the 2-D template itself, and
the anisotropy), plus both versions of the constraint map and the tracer audit. Suite 284.

**Which frame is the field template built in?** Equatorial: the residual proper motions
``(mu_alpha*, mu_delta)`` after the exact systemic-field subtraction. The field is a
Galactic population, so the physically natural frame is Galactic, and the two differ: the
angle between equatorial and Galactic north varies from -6.0 to -10.3 degrees across the
field, so ignoring the rotation displaces a 7.5 mas/yr field star by up to 0.53 mas/yr --
more than the cluster's dispersion. Tested by rebuilding the template in ``(mu_l*, mu_b)``
and rotating each member's residual and radial direction to match: **sigma changes by
0.0002-0.0005 mas/yr, 0.0-0.2 per cent**, against errors of 0.006-0.010. The reason is that
the decomposition only uses the field density *in the neighbourhood of the cluster peak*,
where the field is a smooth, slowly varying floor; rotating the bulk of the field blob
barely changes it, and the free normalisation absorbs what is left. Equatorial is kept
(fewer transformations, identical answer).

**"The background has no systemic velocity" (user, correctly).** The subtraction of
``mu_sys(position)`` is a coordinate shift applied to every star so that members scatter
about zero; for field stars it is not a physical statement, and the field's distribution is
position-independent in **absolute** proper motion, not in the shifted frame. Building the
template in the shifted frame therefore displaces it by the difference in ``mu_sys`` between
the template region and the target annulus. Measured: that difference is
**(0.0003, 0.014) mas/yr** -- negligible, because the perspective term is radial and cancels
when averaged around an annulus; only its scatter within each region survives (0.03 in
pmra*, 0.16 in pmdec for the template ring). Refitting with the template built and scored in
absolute proper motion changes sigma by **<= 0.1 per cent**. The pipeline now does it the
correct way regardless: `MemberSample` carries the systemic field per star (`sys_a`,
`sys_d`) and exposes `absolute_pm`; `field_density_2d(absolute=True)` is the default and
`dispersion_2d(field_at=...)` scores the field there. A test builds a mock whose field sits
at a large offset and shows that scoring in the wrong frame mis-assigns the field fraction
by more than 0.1.

### Per-annulus fits: figures, residuals and recorded fit quality (user request)

`ocen plot-constraints` now also writes `plots/outer_fit_annuli_radial.png` and
`..._tangential.png`: one panel per annulus over |mu| < 15 mas/yr on a logarithmic count
axis, showing the data, the fitted cluster component, the field component and their sum,
with a residual strip (data − model)/sqrt(model) beneath each and the chi2 per bin printed
for the whole range and for the cluster peak separately. The cluster peak region is shaded.
The same routine writes the fit-quality product
`data/processed/kinematics/ocen_outer_fit_quality.ecsv` (per annulus: N, f_field, N_cluster,
sigma_R, sigma_T with errors, mean_R, mean_T, and four chi2 values).

| annulus ["] | N | f_field | sigma_R | sigma_T | chi2/bin R (all / peak) | T (all / peak) |
|---|---|---|---|---|---|---|
| 300-369 | 41 | 0.049 | 0.501 ± 0.065 | 0.505 ± 0.066 | 0.05 / 0.57 | 0.06 / 0.48 |
| 455-560 | 2702 | 0.036 | 0.491 ± 0.010 | 0.441 ± 0.009 | 0.66 / 0.97 | 0.73 / 1.19 |
| 689-849 | 15526 | 0.076 | 0.412 ± 0.008 | 0.382 ± 0.008 | 1.32 / 1.04 | 1.58 / 1.62 |
| 849-1045 | 19057 | 0.127 | 0.363 ± 0.007 | 0.354 ± 0.007 | 1.66 / 1.56 | 2.00 / 1.96 |
| 1286-1583 | 17152 | 0.445 | 0.280 ± 0.006 | 0.292 ± 0.006 | 1.43 / 1.35 | 1.58 / 1.47 |
| 1949-2400 | 19968 | 0.883 | 0.212 ± 0.008 | 0.230 ± 0.009 | 1.38 / 0.97 | 1.51 / 0.96 |

The model tracks the data over three to four decades in counts in every annulus. Fit
quality is **chi2/bin = 0.05-1.7 over the full range and 0.5-2.0 at the peak**; it is worst
(1.6-2.0) at 849-1045 arcsec, where the cluster and field contribute comparably and the
statistics are largest, and best in the outermost annulus where the field dominates and is
measured directly. Nothing in the residuals is systematic at the level that moves sigma:
the peak residuals alternate in sign bin to bin rather than showing a coherent excess or
deficit.

**Also fixed while making these**: the line-of-sight depth term was being added entirely to
the radial variance in the 2-D fit, whereas it acts along the direction of the systemic
proper motion. It now enters as a rank-1 covariance along that fixed direction (and its
projection is used when drawing the model curves). The effect on sigma is at the fourth
decimal -- the term is 0.01-0.03 mas/yr against dispersions of 0.2-0.5 -- but the model is
now the one described in the docstring. Suite 288.

### Membership cut versus decomposition, component by component (user request)

`plots/outer_method_comparison.png` (also written by `ocen plot-constraints`): the radial
and tangential dispersion profiles measured both ways, their fractional difference with the
field fraction overlaid, and the anisotropy each method implies.

| r ["] | f_field | sigma_R: cut / decomp | sigma_T: cut / decomp | sigma_T/sigma_R: cut / decomp |
|---|---|---|---|---|
| 435 | 0.045 | 0.508 / 0.509 (+0.2 %) | 0.419 / 0.434 (+3.7 %) | 0.824 / 0.853 |
| 635 | 0.048 | 0.451 / 0.452 (+0.2 %) | 0.394 / 0.402 (+2.0 %) | 0.874 / 0.890 |
| 941 | 0.127 | 0.360 / 0.363 (+1.0 %) | 0.347 / 0.354 (+2.0 %) | 0.964 / 0.974 |
| 1152 | 0.235 | 0.308 / 0.313 (+1.6 %) | 0.318 / 0.322 (+1.4 %) | 1.032 / 1.030 |
| 1408 | 0.445 | 0.274 / 0.280 (+2.3 %) | 0.286 / 0.292 (+2.2 %) | 1.043 / 1.042 |
| 1722 | 0.703 | 0.231 / 0.236 (+2.3 %) | 0.257 / 0.265 (+3.1 %) | 1.115 / 1.124 |
| 2122 | 0.883 | 0.210 / 0.212 (+0.9 %) | 0.226 / 0.230 (+1.4 %) | 1.076 / 1.082 |

The two agree to **0-4 per cent in every annulus and both components**, with the
decomposition always the higher of the two (a membership cut is a proper-motion cut: it
removes the cluster's own velocity wings along with the field). The difference does **not**
grow with the field fraction -- it is +2 to +3 per cent at f = 0.05 and the same at
f = 0.88 -- so it is the wing truncation, not contamination, that separates the methods.
Contamination itself, properly modelled, moves nothing by more than its error bar. The
anisotropy profiles are identical within errors: both give sigma_T/sigma_R rising from
0.82-0.85 at 435 arcsec to 1.08-1.12 beyond 1700 arcsec, crossing unity at ~1000 arcsec
(27 pc). A test now asserts the two methods agree to 6 per cent and within 2 sigma. Suite 288.

### The "wiggle" at 10-20 pc: an artefact of my own comparison, not a feature of the cluster

The user spotted a wiggle in the residual panel of the constraint map between 10 and 20 pc:
the Gaia points dipped below the model near 350-450 arcsec and jumped ~9 per cent above it
at 500-700 arcsec. Traced and removed. The dispersions on that panel are ours, measured as
scatter **about the mean motion fitted in the same annulus**, but the streaming term I
removed from the model was the **published** Vasiliev & Baumgardt rotation curve. The two
disagree exactly where the rotation peaks:

| r ["] | our fitted mean tangential motion | published v_rot | ratio |
|---|---|---|---|
| 348 | 0.139 | 0.236 | 0.59 |
| 434 | 0.226 | 0.249 | 0.90 |
| 521 | 0.218 | 0.243 | 0.90 |
| 635 | 0.189 | 0.216 | 0.87 |
| 943 | 0.137 | 0.137 | 1.00 |
| 1760 | 0.037 | 0.033 | 1.11 |

Since the term enters as ``mu_rot^2 / 2`` and the rotation peaks at 430-580 arcsec,
over-subtracting it by 10-40 per cent there carved a dip into the model curve and produced
the apparent excess in the data. The constraint map now uses **our own measured mean
motions** for our own points (the published curve is still used for datasets whose means we
do not have), and the wiggle is gone.

With the comparison made self-consistently, in total second moments, the deviation from the
no-dark-matter model is smooth and rises monotonically outward:

| r [pc] | 9.5 | 11.7 | 13.6 | 15.8 | 18.2 | 21.1 | 24.7 | 29.2 | 35.7 | 45.1 | 57.4 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| data / K1 − 1 | −11 %* | −1 % | +5 % | +4 % | +5 % | +3 % | +2 % | +4 % | +5 % | +10 % | +17 % |

(*53 stars, ±8 per cent). So there is no structure at 10-20 pc; there is a ~4 per cent
offset from 13 to 36 pc and a genuine rise beyond 40 pc. Whether that rise is mass or a
Gaia scale error is the open question, unchanged.

Worth noting for later: our fitted mean tangential motion is **10-13 per cent below** the
published rotation amplitude between 430 and 800 arcsec, and 40 per cent below it at
350 arcsec. Both are measured from the same catalogue, so the difference is in the
definition (a sinusoid amplitude fitted per annulus versus the azimuthal mean of the
tangential component) or in the treatment of contamination. It does not affect the
dispersions, but it should be understood before any rotation result is quoted.

**Both versions kept (user request).** `plot_constraint_map(..., streaming=...)` selects how
the rotation term is removed from the model for our own points, and `ocen plot-constraints`
writes all four combinations:

| file | Gaia measurement | rotation term |
|---|---|---|
| `constraint_map_selfconsistent.png` | field contamination modelled | our own fitted mean motions |
| `constraint_map.png` | field contamination modelled | published Vasiliev & Baumgardt curve |
| `constraint_map_pcut_selfconsistent.png` | P > 0.9 members | our own fitted mean motions |
| `constraint_map_pcut.png` | P > 0.9 members | published curve |

The titles state which treatment each uses. The self-consistent pair is the one to read: in
the other two the model has 10-40 per cent too much rotation removed between 350 and 800
arcsec, which is what produced the apparent wiggle at 10-20 pc.

### Is the wiggle real? No -- but my "it is gone" was wrong too

The user looked at the new figures and said the wiggle was still there. Correct: changing
the rotation term reduced its amplitude but did not remove the feature, and the honest
answer is not that it disappeared but that **it was never significant**. Quantified, in
total second moments so that no rotation curve enters the comparison
(`plots/outer_residual_significance.png`, `ocen plot-constraints`):

| r ["] | r [pc] | stars in the annulus | after the quality flag | residual vs K1 | significance |
|---|---|---|---|---|---|
| 357 | 9.5 | 6888 | **53 (0.8 %)** | −11.4 ± 6.7 % | −1.7 sigma |
| 439 | 11.7 | 10402 | **441 (4.2 %)** | −0.9 ± 3.0 % | −0.3 sigma |
| 509 | 13.6 | 13501 | 1875 (13.9 %) | +4.7 ± 2.8 % | +1.7 |
| 592 | 15.8 | 18206 | 4665 (25.6 %) | +4.4 ± 1.4 % | +3.1 |
| 683 | 18.2 | 22390 | 8247 (36.8 %) | +5.3 ± 1.4 % | +3.7 |
| 791 | 21.1 | 24885 | 12050 (48.4 %) | +3.2 ± 1.4 % | +2.3 |
| 924 | 24.7 | 24367 | 14748 (60.5 %) | +2.1 ± 1.4 % | +1.5 |
| 1094 | 29.2 | | | +3.9 ± 1.4 % | +2.7 |
| 1336 | 35.7 | | | +4.7 ± 1.5 % | +3.2 |
| 1688 | 45.1 | | | +10.0 ± 2.4 % | +4.2 |
| 2151 | 57.4 | | | +16.8 ± 3.3 % | +5.1 |

**A constant offset of +3.9 per cent fits 460-1500 arcsec perfectly: chi2 = 3.5 for 6
degrees of freedom, p = 0.74.** There is no structure between 13 and 36 pc. The apparent
wiggle is the two innermost annuli, at −1.7 and −0.3 sigma, and they are low for a reason
that has nothing to do with the cluster: **inside 460 arcsec the Vasiliev & Baumgardt
astrometric quality flag keeps 0.8-4 per cent of the stars** (53 of 6888 at 300-380
arcsec). Gaia cannot measure a dispersion in the crowded inner region, and those points
should not be plotted as if it could -- they are a tiny, non-randomly selected subsample.
The quality fraction reaches 14 per cent at 500 arcsec and 80 per cent beyond 1500.

So the structure in the Gaia residuals is: nothing inside 460 arcsec (no usable data), a
flat +4 per cent offset from 13 to 36 pc, and a real rise to +10 and +17 per cent (4.2 and
5.1 sigma) beyond 40 pc. The flat offset is a scale difference between Gaia and the
HST-anchored model; the rise is the open question. A test now asserts the flatness
(p > 0.05), the insignificance of the inner points and the significance of the outer rise.

**Proposal, not implemented:** restrict our Gaia profile to r > 460 arcsec in everything
downstream, and mark the inner region as "no usable Gaia data" rather than plotting two
meaningless points. That is a data-selection change, so it waits for the user.

### Corrected again: the feature is a STEP between HST and Gaia, not a bump inside Gaia

The user pushed back: the eye-catching feature is the bump after the dip, and the dip lines
up with the other measurements -- "are you fudging?". They were right about the alignment,
and my flatness test was too narrow to see what matters. Redone with the fitted Gaia scale
**removed**, so both datasets are compared with the same unscaled model
(`plots/hst_gaia_step.png`):

| dataset | r ["] | deviation from the same model |
|---|---|---|
| HST | 160 | −1.0 % |
| HST | 252 | −3.2 % |
| HST | 311 (last bin, edge 346) | **−5.3 %** |
| Gaia | 357 (53 usable stars) | −5.4 ± 7.1 % |
| Gaia | 439 (441 stars) | +5.8 ± 3.2 % |
| Gaia | 509 (1875 stars) | **+11.8 ± 3.0 %** |
| Gaia | 592-1336 | +9 to +12.5 % |
| Gaia | 1688, 2151 | +17.4, +24.7 % |

So: **HST declines monotonically to −5.3 per cent of the model at its outer edge, and Gaia
sits at +10 to +12 per cent as soon as it has enough stars to measure anything.** The
Gaia point at 357 arcsec (−5.4 %) continues the HST trend exactly, as the user said -- it
is not an outlier, it is the last point where the two agree.

Is there a bump *within* the Gaia data? No: 460-730 arcsec gives +4.87 ± 0.94 per cent and
730-1500 gives +3.43 ± 0.71 (scaled units), a difference of **1.2 sigma**; each of the
"bump" points is 0.3-1.0 sigma from the plateau. My earlier p = 0.74 test was correct but
answered the wrong question -- it can only see structure inside one dataset, never a step
between two.

**The step is the real feature**, and it falls exactly in the gap where neither instrument
works: HST's outermost bin ends at 346 arcsec, and Gaia's quality-flagged sample only
passes 15 per cent of stars by ~500 arcsec. Between 346 and 500 arcsec there is no reliable
measurement from either, and the K1 fit bridges it with a single constant, s_EDR3 = 1.068 --
which is why the residuals look like a dip followed by a bump. Three readings remain open:
a Gaia crowding systematic that persists to 500 arcsec, an HST systematic at its own edge
(its last bin spans 300-346 arcsec with 18,485 stars against 60,000-79,000 in the previous
ones), or a genuine feature in the cluster at 9-13 pc. Nothing in the present data
distinguishes them, and a constant instrument scale is the wrong model for it either way.

Tests added: HST declines below −3 per cent at its edge while Gaia exceeds +5 per cent
where it becomes usable, with a step of more than 10 per cent between them.

### Why the Gaia dispersions sit above the model beyond 500 arcsec

The offset decomposes into two unrelated things (`plots/outer_offset_explained.png`).
Worked example at 683 arcsec (18 pc), all in 1-D proper-motion dispersion:

| step | value | meaning |
|---|---|---|
| model total second moment | 0.3889 | the Jeans prediction |
| minus the published rotation term | 0.3617 | what the fit actually compared with |
| times the fitted scale s_EDR3 = 1.068 | 0.3862 | |
| published EDR3 datum | 0.3925 | **fit residual +1.6 %** -- the fit matches the data it was given |
| our dispersion | 0.4175 | **+6.4 % above the published profile** |
| our total second moment | 0.4374 | the mean motions add +4.8 % |
| versus model total x scale | 0.4152 | **+5.3 %** |

**Flat part, 500-1350 arcsec (13-36 pc).** Roughly +11 per cent against the unscaled model,
of which ~6.8 per cent is the fitted instrument scale (the fit had already decided Gaia runs
hotter than the HST-anchored model) and ~5-6 per cent is our measurement sitting above the
published profile -- the decomposition keeps the cluster's velocity wings that a
membership-weighted estimate truncates. The fit reproduces the published profile to 1.6 per
cent, so this part is method and calibration, not physics, and it is flat with radius.

**Rising part, beyond 1500 arcsec (40 pc).** Here our measurement *agrees* with the published
profile (+3.1 per cent at 1688, −2.4 per cent at 2151), so it is not a measurement effect:
both lie 10-17 per cent above the model. This is a **shape** mismatch, and the likely cause
is now visible: the K1 model has beta_inf = +0.23, predicting sigma_T/sigma_R = 0.89 at every
radius outside a few pc, while the measurement turns **tangential** -- 0.98, 1.03, 1.04, 1.12,
1.08 beyond 950 arcsec. The two agree at 460-800 arcsec (0.89-0.90 both) and diverge outward.
A radially anisotropic model declines too steeply in projection, so it must fall below data
whose orbits are actually tangential.

That matters for the whole project: the K2 dark-halo preference was driven by the outer
excess, and the outer excess may be an anisotropy error rather than missing mass. The K1 fit
had freedom in beta but was given only a **1-D combined** Gaia profile, which carries no
anisotropy information; it extrapolated the radial anisotropy measured by HST at 10-20 pc
out to 60 pc.

**Proposal, not implemented:** feed the measured sigma_R and sigma_T as two separate datasets
beyond 460 arcsec instead of the 1-D combined published profile, so beta(r) is constrained by
the data at large radius, and refit K1 and K2. That is a change to the data the model is
fitted to, so it waits for the user. A test now pins the disagreement (model < 0.92,
measured > 1.0 beyond 1300 arcsec).

## 2026-09-19 — No rescaling: how the datasets actually behave

The user objected to the per-instrument scales -- "are you proposing some arbitrary
rescaling of the Gaia dispersion???" -- and asked to start again without any. The scales
were mine (introduced 2026-09-17, priors uniform 0.85-1.15 for MUSE and 0.7-1.3 for the two
Gaia sets) and I had stopped flagging them as the substantive assumption they are. Three
things make them indefensible as they stood: they are **degenerate with the signal**
(in the K2 posterior s_GaiaEDR3 and log M_DM correlate at **-0.39**; the median M_DM falls
from 3.8e6 to 3.45e6 across the fitted range of the scale), the two Gaia releases pull in
**opposite directions** (0.950 for DR2, 1.065 for EDR3 -- same telescope, same cluster), and
a radius-independent constant cannot represent either the HST/Gaia step or the outer shape
mismatch.

`ocen fit --no-scales` added; K1 and K2-cored relaunched without any instrument nuisances.
Maximum-likelihood reference (K1, 8 parameters, no scales): chi2 = 629 / 126 against 355 with
the scales free. Per dataset, with **no rescaling anywhere**
(`plots/datasets_unscaled.png`):

| dataset | 0-100" | 100-350" | 350-1000" | 1000-2500" |
|---|---|---|---|---|
| HST PM (oMEGACat) | −0.5 % | −1.9 % | | |
| MUSE line of sight | −4.7 % | −5.4 % | | |
| Gaia DR2 (Baumgardt+ 2019) | | −1.7 % | −3.5 % | −1.6 % |
| Gaia EDR3, our measurement | | −4.1 % | **+7.7 %** | **+12.8 %** |
| Gaia EDR3, published profile | | −6.2 % | +2.6 % | **+13.3 %** |

So: HST is fitted to within 2 per cent; MUSE sits 5 per cent low, the sign and size expected
from energy equipartition (its giants are 3 magnitudes brighter than the HST proper-motion
sample); **Gaia DR2 is consistent with the model everywhere, to within 4 per cent**; and
**only Gaia EDR3 shows the outer excess**.

### The two Gaia releases disagree with each other, and with a radial trend

| r ["] | 179 | 266 | 310 | 338 | 386 | 448 | 670 | 1137 | 1747 |
|---|---|---|---|---|---|---|---|---|---|
| EDR3 / DR2 | −8.5 % | −9.8 % | −9.7 % | −2.2 % | +9.4 % | +5.4 % | +2.6 % | +5.1 % | **+11.9 %** |

The weighted mean beyond 300 arcsec is 1.017 ± 0.015 -- consistent with unity -- but the
ratio runs from −10 per cent at 200-300 arcsec to +12 per cent at 1750 arcsec. Two analyses
of the same cluster with the same telescope, differing by 20 per cent across the radial
range, with a trend rather than an offset. **The outer excess that drives the entire
dark-matter preference is present in EDR3 and absent in DR2.**

That does not tell us which release is right (EDR3 has far better astrometry; DR2's profile
comes from a different membership and error treatment, and we still do not know whether its
dispersions are quoted about a rotating mean). It does tell us the excess cannot yet be
attributed to mass, and that no constant scale is the right way to absorb the difference.
Nested runs without scales (K1, K2-cored) are in progress; their evidence comparison will
say what the data prefer when nothing is free to soak up a 7 per cent offset.
