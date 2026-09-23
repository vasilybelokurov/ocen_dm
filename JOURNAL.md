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

## 2026-09-19 — HST versus Gaia on the same stars, and why EDR3 stops at ~460 arcsec

Three points from the user, all taken: line-of-sight and proper-motion measurements may
legitimately differ; two proper-motion catalogues may not; and EDR3 is the best Gaia data
available. The middle point is testable directly, so it was tested.

### The same stars, measured twice (`plots/hst_gaia_star_by_star.png`)

`selection/hst_gaia_match.py` matches the oMEGACat catalogue to the Vasiliev & Baumgardt
EDR3 members within 0.3 arcsec: **6423 matches**, of which 2211 have P > 0.9 and good HST
astrometry. HST proper motions are locally corrected, so only differences and dispersions
are comparable -- never the zero point.

| sample | N | rms(Gaia − HST) | quoted errors allow | ratio |
|---|---|---|---|---|
| G < 17 | 1800 | 0.99 | 0.13 | **7.4** |
| 17 < G < 18 | 411 | 1.15 | 0.29 | 3.9 |
| G < 18, r < 200" | 925 | 1.07 | 0.15 | 7.0 |
| G < 18, 200-300" | 1158 | 0.99 | 0.17 | 5.9 |
| G < 18, 300-460" | 128 | 0.91 | 0.18 | 5.1 |

And the dispersion each instrument reports **from the identical stars**:

| annulus | N | HST | Gaia | Gaia / HST |
|---|---|---|---|---|
| 100-200" | 761 | 0.649 | 0.766 | **+18 %** |
| 200-300" | 1158 | 0.577 | 0.714 | **+24 %** |
| 300-460" | 128 | 0.534 | 0.633 | **+19 %** |

So the two catalogues do **not** agree where they overlap: Gaia carries about 1 mas/yr of
per-star scatter that its formal errors do not describe -- five to seven times the quoted
precision -- and consequently reports a dispersion 18-24 per cent too high. HST resolves
35,300 stars per square arcminute inside 100 arcsec against Gaia's 62: this is crowding, and
it is Gaia's problem, not HST's.

### Why EDR3 does not reach the radii DR2 does

| annulus | in the EDR3 catalogue | P > 0.9 | passing the quality flag |
|---|---|---|---|
| 0-100" | 543 | 392 | **0** |
| 100-200" | 2528 | 1922 | **0** |
| 200-300" | 5171 | 3973 | **5** |
| 300-460" | 17290 | 14672 | 494 |
| 460-700" | 47369 | 42773 | 12072 |

Vasiliev & Baumgardt's astrometric quality filter removes **every** Gaia star inside
200 arcsec and all but five inside 300 -- exactly the stars the star-by-star test shows to
be wrong. Yet their published profile is tabulated from r = 0 (0.561 mas/yr at the centre,
0.556 at 24 arcsec, ...): **inside ~400 arcsec that profile is a smooth model continued
inward, not a measurement.** Our fits have been fed thinned points from it starting at
300 arcsec, with 1 per cent errors, in a range where no usable Gaia star exists. That is a
real defect in the likelihood, not a matter of taste.

Baumgardt+ 2019's DR2 profile does have points at 179-338 arcsec, and they agree with HST
(−1.7 per cent mean) where the EDR3 profile does not (−6.2 per cent). Since DR2's astrometry
is worse than EDR3's, the agreement most likely reflects their different (brighter, more
conservative) selection and larger quoted errors rather than better data.

**Proposal, not implemented:** use Gaia EDR3 only beyond 460 arcsec, from our own binned
measurement rather than the published smooth profile, and drop the Gaia DR2 profile or
restrict it likewise. With that, "PM datasets must agree" becomes testable rather than
assumed, since HST and Gaia would no longer be asked to describe the same radii.

### Correction: the match above was wrong, and so were its numbers

Everything in the two tables above is **withdrawn**. The cross-match that produced them was
broken in two independent ways, and both broke it in the direction of manufacturing
disagreement.

1. **Pixel coordinates.** The oMEGACat `x, y` columns were used as if they were sky offsets
   in arcsec. Symptom: the separation distribution peaked at 0.18 arcsec against a 0.3
   arcsec tolerance, and the match rate collapsed from 93 per cent at the centre to 1 per
   cent at 460 arcsec. Most "matches" outside the core were therefore *different stars*,
   and the difference between two unrelated cluster members is the cluster dispersion
   itself, not a measurement error.
2. **Epoch.** Matching on RA/Dec instead left a uniform 0.101 arcsec displacement between
   the catalogues -- the cluster's own systemic proper motion, 7.5 mas/yr over the ~13 yr
   between the HST and Gaia epochs. It is now measured from a generous first pass,
   (-0.0447, -0.0901) arcsec, removed, and the tolerance tightened to 0.06 arcsec.

The match now has **8064 pairs with a median separation of 0.005 arcsec**, a hundred times
inside the tolerance. `tests/test_hst_gaia_match.py` pins that first, before any physics.

### What the corrected match says

Same stars, HST high-quality astrometry, P > 0.9, G < 17, unflagged Gaia. Differences are
taken about their own median in each annulus, since the HST astrometry is relative.

| annulus | N | HST | Gaia | Gaia / HST | Gaia scatter beyond its errors | quoted |
|---|---|---|---|---|---|---|
| 0-100" | 68 | 0.654 | 0.934 | 1.43 ± 0.17 | 0.90 | 0.25 |
| 100-200" | 660 | 0.588 | 0.840 | 1.43 ± 0.06 | 0.79 | 0.20 |
| 200-300" | 1178 | 0.535 | 0.708 | 1.32 ± 0.04 | 0.59 | 0.15 |
| 300-460" | 52 | 0.466 | 0.584 | 1.25 ± 0.17 | 0.44 | 0.12 |

The qualitative conclusion survives, at about half the amplitude: the per-star difference is
0.65 mas/yr rms where the quoted errors allow 0.17, a factor 3.8 rather than the 5-7 claimed
above, and the undeclared scatter **falls monotonically outwards** exactly as crowding
should. The earlier inflation came from the bad pairs.

### The decisive test: the one annulus where both catalogues are usable

HST's high-quality astrometry ends at 380 arcsec (8145 stars have proper motions at
380-460 arcsec, **none** of them high-quality), and Gaia's quality flag passes no star
inside 300 arcsec. The entire overlap is **300-380 arcsec**. There, in the same annulus:

| sample | N | sigma (mas/yr) | ratio to HST |
|---|---|---|---|
| HST, high-quality astrometry | 15629 | 0.5238 ± 0.0021 | 1 |
| **Gaia, quality flag** | **51** | **0.482 ± 0.034** | **0.92 ± 0.07** |
| Gaia, P > 0.99 and G < 17 | 1871 | 0.5315 ± 0.0065 | 1.015 ± 0.013 |
| Gaia, G < 16 | 1124 | 0.5385 ± 0.0083 | 1.028 ± 0.016 |
| Gaia, G < 17 | 2022 | 0.5671 ± 0.0068 | 1.083 ± 0.014 |
| Gaia, no quality cut | 5714 | 0.6785 ± 0.0055 | 1.295 ± 0.012 |

**The two instruments agree once the astrometric quality flag is applied.** They disagree by
30 per cent without it. The disagreement is entirely carried by the stars the flag rejects,
which is what the flag exists for. No star-by-star test can probe this directly, because
zero quality-flagged Gaia stars have an HST counterpart -- the flag and HST's usable
footprint are almost disjoint.

The figure `plots/hst_gaia_star_by_star.png` now shows all four pieces: the per-star
difference against the quoted errors, the undeclared scatter falling outwards, the
dispersion from identical stars with the flagged point overlaid, and the quality-flag
survival fraction.

### Residual crowding inside the flagged sample

A faint star is measured worse in a crowded field, so a brightness trend in the dispersion
at fixed radius is a crowding signature. Ratio of sigma(G > 19) to sigma(G < 17), both with
P > 0.9:

| annulus | quality-flagged | no quality cut |
|---|---|---|
| 460-700" | 1.114 ± 0.023 | 1.417 ± 0.017 |
| 700-1000" | 1.135 ± 0.019 | 1.284 ± 0.017 |
| 1000-1500" | 1.077 ± 0.022 | 1.101 ± 0.020 |
| 1500-2400" | 0.985 ± 0.039 | -- |

The flag removes most of the trend but not all of it at 460-1000 arcsec, where about 10 per
cent in sigma remains between the faint and bright ends and vanishes beyond 1500 arcsec.
Mass segregation predicts a trend of the same sign, so this is an upper limit on residual
crowding rather than a detection of it; either way it is a 10 per cent effect where the
unflagged catalogue has a 30-40 per cent one.

### Consequence for the modelling

Unchanged from the proposal above, and now on firmer ground: Gaia EDR3 should enter the
likelihood only where its quality-flagged sample exists in quantity (r > 460 arcsec), from
our own binned measurement, and HST should carry everything inside. The published EDR3
profile inside ~400 arcsec remains an inward continuation over a region containing at most
53 usable stars. The earlier "HST and Gaia disagree" framing was too strong: **with the flag
applied they agree to 8 ± 7 per cent, and the apparent step between the datasets came from
comparing flagged Gaia at large radius with a profile extrapolated inwards.**

## 2026-09-19 (later) -- what exactly we were comparing, and why we disagreed

The question was whether our own outer PM dispersion and the published Gaia EDR3 profile are
apples to apples. They were not. The difference is now understood and reproduced.

### The two sides

**Theirs.** Vasiliev & Baumgardt (2021), *Gaia EDR3 view on Galactic globular clusters*,
MNRAS 505, 5978, https://arxiv.org/abs/2102.09568. The file is `profiles/NGC_5139_oCen.txt`
inside their cluster archive: 101 rows on a uniform 24 arcsec grid from 0 to 2400 arcsec,
five percentiles of the PM dispersion and five of the PM rotation. Their method, from their
Section 2:

* all EDR3 sources in the field, a 'clean' subset defined by eight quality cuts (G > 13,
  RUWE < 1.15, astrometric excess noise, IPD harmonic amplitude and multi-peak, visibility
  periods, BP/RP excess, duplicated source) -- this is the bit-2 quality flag in their
  catalogue;
* a mixture model fitted by MCMC: cluster surface density a **Plummer** profile with free
  scale radius, field **two Gaussians**, membership probabilities marginalised;
* the intrinsic PM dispersion a **cubic spline in radius with 2-5 nodes**, isotropic by
  construction ("we used isotropic PM dispersion in the mixture model"); anisotropy is a
  separate post-processing step (their Section 6);
* rotation a **fixed functional form** with one free amplitude,
  `mu_t(R) = mu_rot * 2 (R/R0) / (1 + (R/R0)^2)`, `R0` the Plummer scale radius;
* the radial PM component **fixed** to perspective expansion, no free parameter.

**Ours.** Independent maximum-likelihood fits in ten log-spaced annuli from 300 to 2400
arcsec, on the same catalogue and the same quality flag, with free mean and dispersion in
the radial and tangential directions separately, per-star error covariance deconvolved,
exact perspective removal, the line-of-sight depth term removed, and contamination handled
either by a membership cut or by an empirical two-dimensional field template.

### Not apples to apples: four structural differences

1. Their profile is a 2-5 node spline spanning 0 to 2400 arcsec. It borrows strength across
   all radii, which is why its quoted uncertainty falls to **0.002 mas/yr (0.5 per cent) at
   1050 arcsec** while ours is 1.4 per cent from 19000 stars in one annulus. That number is
   the posterior width of a stiff global model, not the precision of a local measurement.
2. Their dispersion is isotropic; ours is fitted as sigma_R and sigma_T.
3. Their rotation is one amplitude on a fixed shape; ours is a free mean per annulus. This
   one turns out not to matter: our fitted mean tangential motion tracks their rotation
   profile closely (0.226 vs 0.250 at 435 arcsec, 0.159 vs 0.177 at 773, 0.037 vs 0.036 at
   1722), so both remove the same streaming.
4. Their inner points are not measurements. Five of the eight points our likelihood takes
   sit inside 380 arcsec, where their own quality flag passes no star inside 200 and 53
   inside 380.

### The real cause of the 4-6 per cent offset: we were using raw errors

Their readme lists a `source density` column "used to determine parallax/PM uncertainty
scaling factors". Their Section 3 and Table 1 give it: `eta = (1 + Sigma/Sigma_0)^zeta` with
`zeta = 0.04` and `Sigma_0 = 10` (5-parameter) or `5` (6-parameter) for the clean subset,
reaching `eta ~ 1.1-1.15` in the densest regions. They state they apply the parallax
prescription to the proper motions as well.

**The released catalogue carries the raw, unscaled Gaia errors.** Verified directly:
`gaia_edr3.gaia_source` on WSDB joined to all 228055 members gives a median
`catalogue / raw` ratio of **1.0000** (1st to 99th percentile 0.9992 to 1.0009). The scaling
is theirs to apply at fit time, and we were not applying it.

Applying it in our annuli (`kinematics/vb2021_replication.py`, tests in
`tests/test_vb2021_replication.py`):

| r (arcsec) | theirs | ours, raw errors | ours, with eta | median eta |
|---|---|---|---|---|
| 347 | 0.4869 | 0.5032 (+3.4 %) | 0.5024 (+3.2 %) | 1.116 |
| 435 | 0.4633 | 0.4655 (+0.5 %) | 0.4568 (-1.4 %) | 1.126 |
| 521 | 0.4380 | 0.4601 (+5.0 %) | 0.4446 (+1.5 %) | 1.129 |
| 635 | 0.4054 | 0.4239 (+4.6 %) | 0.4030 (-0.6 %) | 1.129 |
| 773 | 0.3707 | 0.3932 (+6.1 %) | 0.3707 (+0.0 %) | 1.124 |
| 941 | 0.3375 | 0.3535 (+4.8 %) | 0.3342 (-1.0 %) | 1.103 |
| 1152 | 0.3061 | 0.3133 (+2.4 %) | 0.2980 (-2.6 %) | 1.080 |
| 1408 | 0.2766 | 0.2805 (+1.4 %) | 0.2701 (-2.4 %) | 1.059 |
| 1722 | 0.2501 | 0.2454 (-1.9 %) | 0.2381 (-4.8 %) | 1.044 |
| 2122 | 0.2314 | 0.2195 (-5.1 %) | 0.2132 (-7.8 %) | 1.038 |

One documented factor removes the entire discrepancy between 455 and 1000 arcsec. Their
second prescription -- use only stars with `err < 0.2 sigma(R) / (eta - 0.9)` to constrain
the dispersion, which drops about half the sample beyond 700 arcsec -- changes the answer by
under 1 per cent.

**So we can replicate their procedure, and when we do we agree with them to about 2 per cent
out to 1400 arcsec.** Beyond that we fall below them (-6 per cent at 2122 arcsec), where the
field fraction is 0.88 and their spline has few nodes; that one is not resolved and may
belong to either side.

### Consequence, and a proposal (not implemented)

**Our own outer dispersion measurement is biased high by 4-6 per cent between 460 and 1000
arcsec** because it deconvolves raw Gaia errors. A 5 per cent dispersion bias at large radius
propagates to roughly 10 per cent in enclosed mass, in exactly the region where the dark
halo is supposed to show itself. Both `our_outer_profile` and `our_mixture_profile` are
affected, and so is anything plotted from them.

Proposed, for agreement before implementation:

1. apply `eta` in `load_members` (or as an explicit option consumed by both profile
   builders), with the raw behaviour kept for comparison;
2. keep the low-error selection off by default, since it changes < 1 per cent and costs half
   the sample;
3. rebuild the outer profile products and the plots that use them.

Until then, every outer-dispersion number in this journal before today stands 4-6 per cent
high in the 460-1000 arcsec range.

### Follow-up: is the error inflation actually the right fix?

Applying their `eta` makes our numbers match theirs, but matching a curve is not evidence.
Their own non-circular validation (their Figure 6) is that a correct error model makes the
**error-deconvolved dispersion independent of magnitude**. Running that test on our annuli:

| annulus | median eta | sigma(G>19)/sigma(G<17), raw errors | with eta |
|---|---|---|---|
| 460-700" | 1.129 | 1.114 ± 0.023 | 0.954 ± 0.023 |
| 700-1000" | 1.118 | 1.135 ± 0.019 | 0.980 ± 0.018 |
| 1000-1500" | 1.078 | 1.077 ± 0.022 | 0.951 ± 0.021 |
| 1500-2400" | 1.045 | 0.985 ± 0.039 | 0.878 ± 0.038 |

Raw errors are **under**-estimated for faint stars, as they say. But their density-only
`eta` **over**-corrects: the ratio goes from 1.11 to 0.95, overshooting 1.00. Their own
Figure 5 shows `eta` depends on magnitude as well as density, and they adopted the
density-only form for simplicity; the residual shows where that costs something. So `eta` is
a better error model than the raw one, not a correct one.

Split by magnitude at 460-1000 arcsec, the whole argument becomes visible:

| G | N | median error | sigma, raw | sigma, with eta |
|---|---|---|---|---|
| 13-16 | 1386 | 0.031 | 0.3864 | 0.3860 |
| 16-17.5 | 2339 | 0.084 | 0.3903 | 0.3878 |
| 17.5-19 | 13438 | 0.188 | 0.3911 | 0.3786 |
| 19-20 | 13880 | 0.397 | 0.4127 | 0.3632 |
| 20-21 | 6819 | 0.805 | 0.4598 | 0.2703 |

For bright stars the error is a tenth of the signal and the error model is irrelevant: 0.386
either way. For G > 19 the error exceeds the dispersion and the answer is whatever the error
model says, swinging from 0.46 to 0.27. **The 4-6 per cent discrepancy never existed in the
well-measured stars; it lived entirely in the faint majority, where the measurement is a
statement about the error model rather than about the cluster.**

### The measurement to trust

`low_noise_profile()` keeps only stars with `err < 0.4 sigma(R)`, where an `eta` rescaling can
move the answer by at most about 2 per cent by construction, and reports both versions:

| r (arcsec) | N | theirs | ours, raw | ours, with eta | shift |
|---|---|---|---|---|---|
| 518 | 1101 | 0.4390 | 0.4437 (+1.1 %) | 0.4403 (+0.3 %) | 0.8 % |
| 626 | 2027 | 0.4078 | 0.4099 (+0.5 %) | 0.4060 (-0.4 %) | 0.9 % |
| 766 | 2420 | 0.3722 | 0.3761 (+1.1 %) | 0.3726 (+0.1 %) | 1.0 % |
| 935 | 2045 | 0.3384 | 0.3392 (+0.2 %) | 0.3365 (-0.6 %) | 0.8 % |
| 1152 | 1373 | 0.3060 | 0.3083 (+0.7 %) | 0.3062 (+0.1 %) | 0.7 % |

Ours and theirs agree to about 1 per cent from 520 to 1150 arcsec **whichever error model is
used**, at the cost of keeping roughly 15 per cent of the stars. Statistical errors are
0.004-0.007 mas/yr, still far below their spline's quoted 0.002-0.003 but honest.

Revised proposal, replacing yesterday's: rather than adopting `eta` (which over-corrects) or
keeping raw errors (which under-correct), build the outer profile from the low-noise subset
and quote the `eta`-on/off spread as a systematic. That keeps every outer point independent
of an error model we have just shown to be imperfect.

## 2026-09-19 (implementation) -- Gaia EDR3 rebuilt as our own measurement

Agreed with the user and implemented. Nothing about the physical model changed; this is
entirely about which numbers the likelihood is shown.

### What was built

`kinematics/outer_gaia.py` produces `data/processed/kinematics/ocen_pm_dispersion_edr3_ours.ecsv`:

* **quality-flagged EDR3 stars only**, no membership cut, contamination handled by the
  two-dimensional empirical field template;
* **only stars with `err < 0.4 sigma(R)`**. An error rescaling by `eta` can then move the
  deconvolved dispersion by at most `0.4^2 (eta^2 - 1)/2`, about 2 per cent, by construction;
* the value is the **midpoint of the raw-error and eta-inflated fits**, and half their
  separation is carried as a systematic added in quadrature to the statistical error. The
  magnitude-consistency test shows raw errors under-correct and `eta` over-corrects, so the
  truth is between them;
* **seven independent annuli from 460 to 2400 arcsec**. Nothing inside 460, where the quality
  flag leaves too few stars and HST has astrometry the two instruments agree on;
* exact perspective removal, line-of-sight depth removal, and the fitted mean radial and
  tangential motions returned as `streaming2` so the Jeans model is compared with the full
  second moment.

| r (arcsec) | N | sigma | stat | sys | sigma_T/sigma_R | field fraction |
|---|---|---|---|---|---|---|
| 533 | 1477 | 0.4350 | 0.0062 | 0.0017 | 0.897 | 0.06 |
| 662 | 2754 | 0.4016 | 0.0057 | 0.0019 | 0.890 | 0.10 |
| 830 | 3176 | 0.3614 | 0.0051 | 0.0016 | 0.953 | 0.16 |
| 1043 | 2895 | 0.3204 | 0.0045 | 0.0012 | 0.988 | 0.34 |
| 1328 | 2552 | 0.2934 | 0.0083 | 0.0008 | 1.083 | 0.59 |
| 1690 | 2615 | 0.2523 | 0.0089 | 0.0005 | 1.076 | 0.85 |
| 2148 | 3261 | 0.2314 | 0.0138 | 0.0004 | 1.171 | 0.95 |

The error-model systematic is 0.1-0.4 per cent, far below the statistical 1.4-6 per cent, so
the cut did its job. The anisotropy runs radial inside about 1000 arcsec and tangential
outside, which is what Vasiliev & Baumgardt report independently for omega Cen ("NGC 5139
transitions from being radially anisotropic in the inner part to tangentially anisotropic in
the outskirts", their Section 6).

### Wiring

* new dataset key `gaia_edr3_ours` in `kinematics/likelihood.py`, instrument label
  `GaiaEDR3` so the existing nuisance scale still applies;
* it and `gaia_edr3_pm` are forbidden together, since they are the same stars;
* `DEFAULT_DATASETS` now uses it. The published spline stays reachable with `--datasets`;
* `MemberSample.scale_errors` added so the two error models are one call apart.

The default likelihood went from 126 to 125 points: eight published EDR3 points replaced by
seven of ours, five of the eight discarded ones having been inward extrapolation.

### The plot

`plots/pm_datasets_edr3_rebuild.png` (the earlier `plots/datasets_unscaled.png` is kept).
Top panel: HST, Gaia DR2, the published EDR3 spline with its band, our old measurement and
our new one, plus the starred pair at 300-380 arcsec which is the only annulus where HST's
high-quality astrometry and Gaia's quality flag both have usable stars. Bottom panel: the
same divided by the published spline.

Two things are visible at a glance. Inside about 400 arcsec the HST points sit 10-25 per cent
**above** the published EDR3 spline, which is the extrapolation showing itself. And our new
measurement lies on the spline to about 1 per cent from 460 to 2400 arcsec, where our old one
sat 4-7 per cent high.

### Runs

K1 and K2-cored relaunched on the new default datasets as `K1_edr3ours` and
`K2_cored_edr3ours`, same sampler settings as the reference runs (400 live points,
dlogz 0.5, seed 42, composite tracer, instrument scales on). Reference values to beat:
K1 = 134.91 +- 0.47, K2-cored = 185.92 +- 0.51, so the current Delta ln Z is +51.0. Note the
reference runs took 8158 s and 14251 s, not the hour I quoted to the user.

Gaia DR2 was left untouched, although four of its nine points also sit inside 380 arcsec and
are open to the same objection. That was not part of the agreed scope.

### Correction: 300-460 arcsec is usable, and I had thrown it away

The user: "you say 'no usable data here' but then quote less than 500 stars. 500 stars is a
lot!!!" Correct, and the label was conflating two different regions. Counts:

| annulus | in the EDR3 catalogue | passing the quality flag | flagged + err < 0.4 sigma | median error |
|---|---|---|---|---|
| 100-200" | 2528 | 0 | 0 | -- |
| 200-300" | 5171 | 5 | 5 | 0.024 |
| 300-340" | 3169 | 17 | 17 | 0.026 |
| 340-380" | 3719 | 36 | 35 | 0.025 |
| 380-420" | 4659 | 116 | 80 | 0.032 |
| 420-460" | 5743 | 325 | 196 | 0.044 |
| 460-520" | 9971 | 1249 | 581 | 0.090 |

Inside 300 arcsec there really is nothing: five stars. But 300-460 arcsec holds 494 flagged
stars, 328 of them in the low-noise subset, and those are the **cleanest stars in the entire
Gaia sample** -- errors of 0.025 to 0.044 mas/yr against a dispersion near 0.47, a ratio of
5 to 9 per cent. Where the rest of the profile has to argue about error models, here the
error is negligible by a wide margin. They also fill the gap between HST's last point at 311
arcsec and the bulk of the Gaia sample, which is exactly the range where the old step between
datasets appeared.

The inner edge is now **300 arcsec**, with two narrow bins carrying that range:

| r (arcsec) | N | sigma | stat | sys | median G |
|---|---|---|---|---|---|
| 356 | 52 | 0.4789 | 0.0373 | 0.0001 | 13.9 |
| 437 | 276 | 0.4497 | 0.0191 | 0.0007 | 15.1 |

The error-model systematic on these is 0.02 and 0.15 per cent, ten times smaller than
anywhere else in the profile. The 356 arcsec point also reproduces the instrument test
independently: 0.479 +- 0.037 against HST's 0.524 +- 0.002 in the same range is a ratio of
0.91 +- 0.07, matching the 0.92 +- 0.07 found from the annulus comparison.

The profile now has nine points from 300 to 2400 arcsec and the likelihood 127. Both fits
were stopped about half an hour in and relaunched on the corrected dataset, then stopped
again and deleted: the user's "ok, implement, run" meant run the data comparison, not the
nested fits. **No K1/K2 run on the rebuilt dataset exists yet, and none should be started
without being asked for by name.** The plot's
shading and annotation were wrong in the same way and are fixed.

**Lesson for the record:** "too few stars to be worth using" was asserted from a count
without checking the precision those stars carry. The sparse bins turned out to be the most
reliable in the profile, because star count and measurement quality run in opposite
directions here -- the flag keeps only the brightest stars where crowding is worst.

## 2026-09-19 -- where does the cluster stop?

User: "what happens if we go to even larger distance from the cluster's centre ... at some
stage the cluster should stop, no?"

### The Gaia data stop before the cluster does

The Vasiliev & Baumgardt member catalogue ends at **0.670 deg = 2413 arcsec = 63 pc**. That
is their retrieval radius for this cluster, stated in their Section 2 ("a given radius from
the cluster centre ... adjusted individually"), not a physical boundary.

The Jacobi radius, computed with AGAMA in the McMillan (2017) Milky Way potential for
M = 3.55e6 Msun (Baumgardt & Hilker 2018), with the cluster at R = 6.30 kpc, z = +1.42 kpc:

| evaluated at | r_gal | r_J |
|---|---|---|
| present position | 6.46 kpc | 185 pc = 1.95 deg |
| pericentre | 1.59 kpc | 90 pc = 0.95 deg |

The orbit is eccentric (r_peri 1.59, r_apo 6.99 kpc, e = 0.63), so the pericentric value is
the one that limits the bound cluster. **The Gaia sample stops at 0.70 r_J(peri).** Nothing
in our fits has ever seen the boundary.

### What lies beyond, and why the proper motions cannot answer it alone

Two catalogues already in the project reach further. `kuzma2025_periphery` (Pristine CaHK +
Gaia, 157481 stars to 5.1 deg) and `kuzma2026_spectroscopy` (592 line-of-sight velocities to
3.15 deg, 157 flagged members).

The Pristine membership uses the proper motions, so measuring a PM dispersion from it is
circular. Demonstrated directly, at two fixed windows about the systemic motion:

| r (pc) | N (0.8 window) | sigma (0.8) | N (1.2) | sigma (1.2) |
|---|---|---|---|---|
| 44 | 189 | 6.68 ± 0.25 | 189 | 6.68 ± 0.25 |
| 70 | 23 | 6.59 ± 0.70 | 24 | 7.18 ± 0.74 |
| 111 | 13 | 5.73 ± 0.81 | 16 | 10.15 ± 1.28 |
| 197 | 12 | 6.52 ± 0.96 | 21 | 12.81 ± 1.40 |
| 377 | 49 | 8.94 ± 0.64 | 77 | 12.58 ± 0.72 |

At the narrow window the profile is **flat** at 6-7 km/s to 200 pc; at the wide one it rises
to 12.8. That factor of two is the selection, not the cluster. Anyone quoting a single
periphery PM dispersion without stating the window is quoting their own cut.

### The honest probe: line-of-sight velocities, no PM selection

| r (pc) | N | sigma_los (km/s) |
|---|---|---|
| 67 | 116 | 6.12 ± 0.48 |
| 83 | 36 | 8.36 ± 1.07 |
| 178 | 5 | 5.22 ± 1.80 |

### Answer

The cluster does **not** stop at the edge of our data, and the profile does not keep falling.
Our outermost Gaia point is 0.231 mas/yr = 5.96 km/s at 57 pc. Beyond it the dispersion
**flattens at 6-8 km/s** out to at least 200 pc, on both the narrow-window proper motions and
the spectroscopy independently. The flattening sets in around 50-100 pc, which is where
r_J(peri) = 90 pc sits.

A plateau at the Jacobi radius is exactly the ambiguity this project exists to resolve:
**potential escapers and unbound debris produce one, and so does a bound dark halo.** Neither
can be preferred from these numbers alone, because the periphery samples have selection
functions we have not modelled (WP6, the Kuzma 2025 footprint) and the outer bins hold 5 to
50 stars.

Code: `kinematics/periphery.py`, plot `plots/periphery_where_the_cluster_ends.png`, tests
`tests/test_periphery.py`. No fits were run.

### Follow-up: is there actually an overdensity out there?

User: "have we tried detecting an overdensity of cluster-like stars at radii corresponding to
the green points (spectroscopic)?" We had not. The project measured star counts against the
MGE out to 42 arcmin (JOURNAL, light profile validation) and stopped; there was never a
background-subtracted excess test at the radii the spectroscopy reaches.

**Method** (`selection/periphery_density.py`). Pristine periphery catalogue, G0 < 16, uniform
coverage: the all-star surface density varies by only 6 per cent between 0.67 and 5.1 deg
(1785 to 2016 per square degree), so counts divided by annulus area need no footprint
correction. Cluster-like means both a proper motion within 0.8 mas/yr of systemic **and**
[Fe/H] < -1.2, the field median being -0.26 against -1.49 for the cluster's inner stars. The
background is measured, not assumed: the identical selection is repeated in **eight control
windows** of the same radius placed at the same distance from the field's own proper-motion
centroid.

| annulus | r (pc) | N | density | background | excess | significance |
|---|---|---|---|---|---|---|
| 0.50-0.70° | 57 | 59 | 78.25 | 0.17 | 78.1 ± 10.2 | **7.7σ** |
| 0.70-0.90° | 76 | 16 | 15.92 | 0.00 | 15.9 ± 4.0 | **4.0σ** |
| 0.90-1.20° | 100 | 5 | 2.53 | 0.00 | 2.5 ± 1.1 | **2.2σ** |
| 1.20-1.60° | 133 | 3 | 0.85 | 0.25 | 0.60 ± 0.51 | 1.2σ |
| 1.60-2.10° | 175 | 3 | 0.52 | 0.30 | 0.22 ± 0.31 | 0.7σ |
| 2.10-2.80° | 232 | 5 | 0.46 | 0.23 | 0.23 ± 0.23 | 1.0σ |
| 2.80-3.60° | 303 | 6 | 0.37 | 0.19 | 0.19 ± 0.16 | 1.2σ |
| 3.60-4.40° | 379 | 4 | 0.20 | 0.32 | -0.12 ± 0.13 | -0.9σ |
| 4.40-5.10° | 450 | 5 | 0.24 | 0.20 | 0.04 ± 0.11 | 0.3σ |

**Answer.** Yes, and it lands where it matters. The two spectroscopic points that carry the
kinematic plateau sit on a 7.7σ excess at 57 pc and a 4.0σ excess at 76 pc, with 2.2σ still
present at 100 pc. Those stars are omega Cen's. The third spectroscopic point (178 pc, 5
stars) has **no** support in a circular average.

The density falls as r^-6 between 57 and 100 pc, so the stellar component **truncates**, and
it truncates right where r_J(peri) = 90 pc sits.

**But a circular average dilutes a tail.** Beyond 1.6 deg the 23 surviving candidates are not
isotropic: an axial Rayleigh test gives a preferred axis at position angle 148° with
p = 0.016, and the along-axis excess is 0.34 ± 0.16 per square degree against -0.02 ± 0.10
across it (2.1σ). The spectroscopic members beyond 1.2 deg lie on that same axis, at position
angles -32°, -24°, -19°, 149°, 149° -- a bipolar arrangement. That is what tidal tails look
like. The axis was chosen after seeing the data and the result is 2σ, so it is a lead, not a
detection, and the proper test needs the Kuzma selection function (WP6).

**Consequence for the dark-matter question.** The stellar tracer truncates at the Jacobi
radius while the velocity dispersion flattens at 6-8 km/s across the same range. A truncated
tracer with a flat dispersion is precisely the configuration in which the mass inside is
degenerate with what the tracer is doing, and it is also what unbound debris produces. The
two spectroscopic points with real density support (57 and 76 pc) are the only ones that can
carry weight in a bound model.

Code `selection/periphery_density.py`, plot `plots/periphery_overdensity.png`, tests
`tests/test_periphery_density.py`. No fits run.

## 2026-09-19 -- Pristine through our own mixture model, and the profile extended outwards

User: "use pristine and our background + cluster models with gaia stars, extend the profiles
outwards and compare." Done, in `kinematics/pristine_profile.py`.

### What was built

The same cluster-plus-field decomposition we use for Gaia, applied to the Pristine periphery
catalogue with **no proper-motion selection at all**, which removes the circularity that made
the earlier window-based numbers meaningless:

* sample: every Pristine star with [Fe/H] < -1.2, a chemical cut independent of kinematics;
* cluster: a two-dimensional Gaussian with its axes along each star's radial direction,
  convolved with that star's full error covariance. The correlation coefficient comes from
  `kuzma2025_periphery_gaia_covariance`, which is **not in the same row order** as the main
  table and has to be joined on `source_id` (a positional join silently scrambles it; there
  is a test for this);
* field: an empirical two-dimensional proper-motion density built from the same catalogue
  under the same metallicity cut, beyond 3.6 deg where the star-count excess is consistent
  with zero;
* fitted per annulus: the cluster fraction, both dispersions and both mean motions.

The depth term is omitted: at most 0.003 mas/yr, an order below the statistical errors here,
and the tracer model it needs is only calibrated inside 66 pc.

### The comparison, on identical annuli

Two catalogues, two selections, two field models, the same nine bins:

| r (pc) | Gaia | Pristine | Pristine / Gaia | cluster stars | field fraction |
|---|---|---|---|---|---|
| 9.1 | 0.4789 ± 0.0373 | 0.5048 ± 0.0286 | 1.054 ± 0.101 | 86 | 0.000 |
| 11.2 | 0.4497 ± 0.0191 | 0.4415 ± 0.0187 | 0.982 ± 0.059 | 181 | 0.000 |
| 13.9 | 0.4350 ± 0.0064 | 0.4465 ± 0.0126 | 1.026 ± 0.033 | 381 | 0.000 |
| 17.3 | 0.4016 ± 0.0060 | 0.4013 ± 0.0113 | 0.999 ± 0.032 | 424 | 0.002 |
| 21.7 | 0.3614 ± 0.0054 | 0.3569 ± 0.0101 | 0.988 ± 0.032 | 443 | 0.007 |
| 27.0 | 0.3204 ± 0.0047 | 0.3142 ± 0.0089 | 0.981 ± 0.031 | 337 | 0.023 |
| 34.4 | 0.2934 ± 0.0084 | 0.2870 ± 0.0122 | 0.978 ± 0.050 | 181 | 0.072 |
| 43.7 | 0.2523 ± 0.0089 | 0.2488 ± 0.0141 | 0.986 ± 0.066 | 105 | 0.161 |
| 56.7 | 0.2314 ± 0.0138 | 0.2240 ± 0.0206 | 0.968 ± 0.106 | 38 | 0.310 |

**Weighted mean ratio 0.995 ± 0.014.** Every bin is within 5 per cent. This is the strongest
external check the outer Gaia measurement has had: nothing is shared between the two except
the sky.

### Extended outwards

| r (pc) | N | field fraction | cluster stars | sigma (km/s) | usable |
|---|---|---|---|---|---|
| 34 | 396 | 0.08 | 364 | 7.25 ± 0.21 | yes |
| 55 | 87 | 0.32 | 59 | 6.00 ± 0.47 | yes |
| 74 | 66 | 0.76 | 16 | 5.90 ± 0.81 | yes |
| 99 | 115 | 0.96 | 4.8 | 3.80 ± 1.39 | no |
| 136 | 214 | 0.99 | 2.0 | 0.11 ± 0.00 | no |
| 176 | 325 | 0.99 | 2.8 | 3.32 ± 0.98 | no |
| 235 | 551 | 0.99 | 4.2 | 4.48 ± 1.20 | no |

**The profile runs out of cluster, not out of signal-to-noise.** Beyond about 100 pc the
fitted cluster fraction falls below 5 per cent and the mixture has two to five effective
stars to work with; the 136 pc bin is degenerate and collapses to exactly zero. Those bins
are flagged `reliable = False` and are drawn as open symbols, never as measurements. This is
the same truncation the star counts found, seen from the kinematic side.

### What this changes

1. The fixed-window reading of the same catalogue rises from 6.6 to 12.6 km/s over
   70-380 pc. The mixture, on the same stars, does not. **The rise was the window.**
2. Where the mixture still has cluster stars, at 74 pc, it gives 5.90 ± 0.81 km/s and the
   spectroscopy gives 6.12 ± 0.48 at 67 pc. The plateau at 6 km/s survives both an
   independent catalogue and an independent tracer.
3. The 8.36 ± 1.07 km/s spectroscopic point at 83 pc is **not** confirmed by the proper
   motions, which give 5.9 ± 0.8 just outside it. With 36 stars against 16 that is a
   1.6 sigma difference, not a conflict, but it should not be leaned on.

Plot `plots/profile_extended_pristine.png`, tests `tests/test_pristine_profile.py`. No fits.

### Master figure: everything the likelihood is shown

`plots/master_datasets.png`, from `plot_master_datasets()`, which reads `DEFAULT_DATASETS`
rather than a hand-written list, so it cannot drift from what is actually fitted. A test
asserts every default key has a drawing style and that the drawn point count equals the
likelihood's.

Top panel: all five datasets converted to km/s at D = 5.43 kpc so proper motions and
line-of-sight velocities share one axis. Behind them, faint and explicitly labelled "not
fitted": the published EDR3 spline, the Pristine mixture points and the periphery
spectroscopy. Bottom panel: each dataset's radial span and point count.

| dataset | kind | N | radial span |
|---|---|---|---|
| HST, radial PM | pmr | 40 | 1.8-311 arcsec (0.05-8.2 pc) |
| HST, tangential PM | pmt | 40 | 1.8-311 arcsec |
| MUSE, line of sight | los | 29 | 5.7-290 arcsec (0.15-7.6 pc) |
| Gaia DR2 PM (published) | pmc | 9 | 180-1747 arcsec (4.7-46 pc) |
| Gaia EDR3 PM (ours) | pmc | 9 | 356-2148 arcsec (9.4-57 pc) |

**127 points total.** Three things the figure makes plain:

1. Inside 300 arcsec (7.9 pc) there is **no usable Gaia at all**; HST and MUSE carry that
   range alone, and they do so with 109 of the 127 points.
2. The whole Gaia contribution is 18 points over 4.7 to 57 pc, and that is the only range
   where a dark halo can show itself against the stars.
3. **Gaia DR2 still contributes 4 points inside 380 arcsec**, the exact range whose EDR3
   equivalent was removed as extrapolation. This is flagged on the figure. It remains an
   open inconsistency in the dataset selection and is the obvious next thing to decide.

All the new figures are now produced by `ocen-dm plot-constraints`, so they regenerate
together.

### Two corrections to the master figure

**1. There ARE usable EDR3 data inside 380 arcsec, and we already use them.** The earlier
phrasing ("the range whose EDR3 equivalent was removed") was wrong and is withdrawn.

| annulus | passing the quality flag | also err < 0.4 sigma | median error | err / sigma |
|---|---|---|---|---|
| 0-100" | 0 | 0 | -- | -- |
| 100-200" | 0 | 0 | -- | -- |
| 200-250" | 2 | 2 | 0.022 | 0.04 |
| 250-300" | 3 | 3 | 0.026 | 0.05 |
| 300-340" | 17 | 17 | 0.026 | 0.05 |
| 340-380" | 36 | 35 | 0.025 | 0.05 |

The true empty zone is **inside 200 arcsec**, where not one star passes. Between 200 and 300
there are 5 stars, too few for a bin at `min_stars = 25`. Between 300 and 380 there are 53,
and **our innermost EDR3 point already uses them**: r = 356 arcsec, 52 stars,
sigma = 0.4789 +- 0.0373, errors only 5 per cent of the signal.

So what was removed from the published spline was its points at 24, 48, 96 and 168 arcsec,
which sit where the flag passes nothing. Its 336 arcsec point sat on real data and we
replaced it with our own measurement of the same stars, not with nothing.

**2. We measure the radial and tangential components and then throw them away.** HST enters
as two datasets, 40 points each. Our Gaia measurement is built from `sigma_pmr` and
`sigma_pmt` and handed to the likelihood only as the combination (`kind = "pmc"`). There is
no reason for that asymmetry, and the discarded information is the interesting part:

| r (arcsec) | sigma_R | sigma_T | T/R |
|---|---|---|---|
| 356 | 0.4885 | 0.4691 | 0.960 |
| 437 | 0.4806 | 0.4165 | 0.867 |
| 533 | 0.4578 | 0.4109 | 0.898 |
| 662 | 0.4242 | 0.3777 | 0.891 |
| 830 | 0.3700 | 0.3525 | 0.953 |
| 1043 | 0.3224 | 0.3184 | 0.988 |
| 1328 | 0.2815 | 0.3049 | 1.083 |
| 1690 | 0.2429 | 0.2614 | 1.076 |
| 2148 | 0.2126 | 0.2489 | 1.171 |

HST's outer points run T/R = 0.892 to 0.843 over 201-311 arcsec, so the radial bias is
continuous across the instrument boundary, and the turn to tangential happens near 1100
arcsec (29 pc). The figure now carries an anisotropy panel showing HST, our Gaia and
Pristine together.

**Proposal, not implemented.** Feed Gaia EDR3 as two datasets, `gaia_edr3_ours_radial` and
`gaia_edr3_ours_tangential`, exactly as HST is fed, taking the likelihood from 127 to 136
points. That is what would let beta(r) be constrained where the dark halo is supposed to
live, instead of being fixed by HST inside 8 pc and extrapolated outwards. It needs the
covariance between the two components per bin, which the 2-D fit can report but currently
does not. Awaiting agreement.

Also withdrawn: the claim that Gaia DR2's four points inside 380 arcsec sit in a region with
no data. They sit in a region where **EDR3** has almost none, but DR2 has its own, brighter
selection and we hold no star list for it, so nothing can be said about its usability from
here. The annotation was removed from the figure.

### One symbol per dataset, fill for the component

User: "use consistent symbols and colors for the datasets, radial and tangential measurements
should use the same symbol." Implemented as a single table in `plotting/style.py` rather than
per-figure choices, so the figures cannot drift apart.

| dataset | colour | symbol |
|---|---|---|
| HST (oMEGACat) | orange, SERIES[1] | circle |
| MUSE, line of sight | amber, SERIES_EXTRA | square |
| Gaia DR2 (published) | ink | plus |
| Gaia EDR3 (our measurement) | blue, SERIES[0] | diamond |
| Pristine periphery | green, SERIES[2] | down triangle |
| periphery spectroscopy | secondary ink | star |

The **component is the marker fill, never a different symbol**: `full` for a scalar quantity
(a combined proper-motion dispersion, or a line-of-sight one), `left` for radial, `right` for
tangential, via matplotlib's half-filled markers. So HST radial and HST tangential are two
halves of the same orange circle, and the Gaia radial and tangential components are two
halves of the same blue diamond. Data shown for context keep their symbol and colour and are
lightened, so nothing that is merely unfitted can be mistaken for a different instrument.

`dataset_style(key, component, fitted=)` returns the marker keywords and
`dataset_label(...)` the matching legend text; `DATASET_KEY_MAP` translates a likelihood
dataset key straight into the pair. Five tests pin the invariants: components share symbol
and colour and differ only in fill, every dataset's (colour, symbol) pair is unique, an
unfitted mark keeps its identity, and every key in `DATASETS` and in `DEFAULT_DATASETS` has
an agreed symbol.

Applied to `master_datasets.png`, `profile_extended_pristine.png` and
`periphery_where_the_cluster_ends.png`. The older figures still use their own colours and
will be migrated when next touched.

## 2026-09-19 -- the HST/Gaia overlap, and why our HST measurement was wrong

User, twice: "we need an overlap between HST and Gaia EDR3 measurements!!!! you keep showing
me that they dont overlap, why???" The answer is that there was never a gap in the data,
only in the product we were reading.

### There was no gap

| | reaches |
|---|---|
| published oMEGACat PM profile | 300 arcsec |
| oMEGACat astrometry quality flag | ~340 arcsec (66 stars at 340-380, **0** beyond 380) |
| **the oMEGACat catalogue itself** | **466 arcsec**, 23730 stars at 340-380 and 7725 at 380-420 |

We had been feeding the likelihood the published profile, so HST stopped at 300 arcsec while
our Gaia measurement started at 356. Measuring HST ourselves in Gaia's own annuli
(`kinematics/hst_profile.py`) produces a genuine **two-bin overlap at 300-460 arcsec**.

Taking the unflagged stars raw does not work: their error-deconvolved dispersion *rises*
outwards (0.80, 0.88, 1.78 mas/yr at 340-380, 380-420, 420-466), which no cluster does. HST
proper motions are relative to the cluster, so members sit at the origin and field stars sit
7.5 mas/yr away; only 1-2 per cent of stars lie beyond 3 mas/yr but in a raw variance they
dominate. The same cluster-plus-field mixture used everywhere else, with the same empirical
Gaia DR3 field template scored in absolute proper motion, handles it.

### The user caught the remaining error

"the published omegacat catalog agrees almost perfectly with our measured Gaia profile, why?
maybe we need to rethink our HST measurements?" Correct. Checked directly:

| our HST measurement | vs the published oMEGACat profile |
|---|---|
| flagged stars only | 1.1 to 2.3 per cent high |
| all stars, uncorrected | 4.7 to 6.1 per cent high |

So our method reproduces theirs when given the same stars, and the excess comes entirely from
the stars oMEGACat rejects. Running the same mixture on each subset where both exist:

| annulus | sigma flagged | sigma unflagged | unflagged / flagged |
|---|---|---|---|
| 150-200" | 0.6268 ± 0.0089 | 0.6846 ± 0.0097 | 1.092 ± 0.022 |
| 200-250" | 0.5852 ± 0.0083 | 0.6367 ± 0.0090 | 1.088 ± 0.022 |
| 250-300" | 0.5511 ± 0.0078 | 0.5945 ± 0.0084 | 1.079 ± 0.022 |
| 300-340" | 0.5257 ± 0.0074 | 0.5531 ± 0.0078 | 1.052 ± 0.021 |

**Weighted mean 1.077 ± 0.011.** The unflagged stars' errors are underestimated by about
8 per cent -- the same pathology as unflagged Gaia, at a fifth of the amplitude. Outside 340
arcsec every HST star is unflagged, so the raw outer points carry the full bias.
`hst_profile` now divides each bin by `1 + f_unflagged * (ratio - 1)` and carries the
ratio's uncertainty as a per-bin systematic.

### The overlap, after the correction

| annulus | HST (r) | Gaia (r) | Gaia/HST before | after, radius-matched |
|---|---|---|---|---|
| 300-380" | 0.5042 (323") | 0.4789 (356") | 0.950 | 0.966 ± 0.077 |
| 380-460" | 0.4548 (397") | 0.4497 (437") | 0.989 | 1.013 ± 0.047 |

**Weighted mean Gaia/HST = 1.000 ± 0.040.** The instruments agree exactly. (The two samples
sit at different radii inside one annulus -- HST's coverage falls outwards, Gaia's flag
passes more stars outwards -- so both are moved to the annulus midpoint on the local slope
before the ratio is formed.)

### What this revises

Every earlier statement in this journal that HST and Gaia disagree by 7-9 per cent was
measuring HST's rejected stars, not the instruments. The correct statement is that the two
agree to 0 ± 4 per cent once each catalogue's own quality flag is respected, and that both
published profiles are reproduced by our pipeline to 1-3 per cent.

A residual 1-2 per cent offset between our flagged measurement and the published oMEGACat
profile remains, from method (our per-annulus mixture against their Voronoi-binned MCMC). It
is a floor on any instrument comparison and is not worth chasing further.

Code `kinematics/hst_profile.py` (with `unflagged_bias()` to remeasure the correction),
figure `plots/hst_gaia_overlap.png`, seven tests in `tests/test_hst_profile.py`. Also fixed
`tests/test_master_plot.py`, which still referenced a style table removed when the master
figure moved to the shared `DATASET_STYLE`.

**Not implemented, awaiting agreement:** extending the fitted HST dataset from 300 to 460
arcsec using this measurement, which would give the likelihood a real instrument overlap
instead of two profiles meeting at a point.

### What the HST quality flag actually is, and whether the correction was justified

`selection_hq_astrometry` is oMEGACat's own published selection (Häberle et al. 2025, ApJ,
https://arxiv.org/abs/2503.04903, Section 2.1), not anything we defined. It keeps 48 per cent
of the stars with proper motions and requires all of:

* temporal baseline longer than **10 years**;
* **N_used / N_found > 0.75** -- the fraction of individual measurements surviving clipping,
  a low value meaning unreliable astrometry;
* **reduced chi-square < 5** in both proper-motion components;
* the proper-motion error inside the **lower 95 per cent of the error distribution in its own
  0.5-magnitude bin**;
* reliable photometry in **both F625W and F814W**: unsaturated (which excludes everything
  brighter than F625W = 13.8), point-spread-function fit quality above the 85th percentile of
  its magnitude bin, and neighbouring flux inside the fit aperture below half the star's own;
* **F625W < 24** outright, because there the error limit reaches 0.3 mas/yr, "similar to half
  of the typical velocity dispersion in the outer regions", and "including stars with errors
  similar to the actual velocity dispersion would ... make it quite sensitive to the modeling
  of the proper motion errors."

That last sentence is the principle, and it is the same one this project reached
independently for Gaia with the `err < 0.4 sigma(R)` cut. The two-filter photometric
requirement is also not cosmetic: F625W and F814W span 2002 to 2022, so requiring both is how
they verify the astrometry holds across the whole baseline.

The `_and_membership` version adds a colour-magnitude cut to the red-giant branch or main
sequence plus a total proper motion below 4.5 mas/yr.

Empirically, flagged against unflagged among stars with proper motions: median F625W 20.86
against 22.77, median PM error 0.06 against 0.12 mas/yr, pass rate falling from 82 per cent
at F625W 17-19 to 38 per cent at 23-24 and zero beyond 24.

**Was the 1.077 correction justified for the outer bins?** Beyond 380 arcsec not one star has
F625W or F814W photometry, so the two-filter requirement rejects them automatically and their
proper-motion errors are in fact decent, 0.08-0.10 mas/yr. If they failed only on
photometry their astrometry might be sound. Tested on inner unflagged stars split by whether
they have photometry:

| annulus | flagged | unflagged, has photometry | unflagged, no photometry |
|---|---|---|---|
| 150-250" | 0.6059 | 0.6606 (1.090) | 0.7114 (1.174) |
| 250-300" | 0.5511 | 0.5950 (1.080) | 0.5873 (1.066) |
| 300-340" | 0.5257 | 0.5605 (1.066) | 0.5486 (1.044) |

Missing photometry is **not** a free pass: those stars are inflated as much as the rest, and
at 150-250 arcsec more. That is what the paper's own reasoning predicts, since the photometry
requirement exists partly to certify the astrometry. The correction stands. The caveat is
that the no-photometry ratio nearest the boundary is 1.044 rather than 1.077, so the outer
correction may be up to 3 per cent too large, which would move Gaia/HST from 1.000 to about
0.97 and remain consistent with agreement.

### The clean overlap: both sides quality-selected, nothing corrected

Answering "can we stay within 340 arcsec, and is there enough Gaia there?" -- yes to both,
with the precision stated honestly. HST restricted to its own flagged stars, Gaia to
quality-flagged low-noise stars, no unflagged-star correction applied to either side:

| window | HST flagged | Gaia flagged | Gaia/HST | precision |
|---|---|---|---|---|
| 300-340" | 16978 stars at r = 311" , 0.5257 ± 0.0074 | 17 stars at r = 318" , 0.5385 ± 0.0685 | **1.028 ± 0.132** | 13 % |
| 300-380" | 16978 stars at r = 311" , 0.5256 ± 0.0074 | 52 stars at r = 356" , 0.4791 ± 0.0373 | **0.933 ± 0.074** | 7 % |

The 300-340 window is the cleaner comparison and the reason is the effective radii:
**311 and 318 arcsec**, essentially the same place, so the slope correction is negligible and
nothing is being extrapolated. Its cost is that Gaia has only 17 usable stars, giving 13 per
cent. Widening to 380 triples the Gaia sample and halves the error but moves Gaia's effective
radius to 356 arcsec against HST's 311, so a 45-arcsec slope correction does part of the
work.

Both windows agree with unity, and both agree with the corrected all-star comparison
(1.000 ± 0.040). Three routes, three different treatments of the unflagged stars, one answer.

**Decision to propose:** fit HST out to 340 arcsec rather than 460, using flagged stars only,
which adds one bin beyond where the published profile stops at 300 and needs no correction
of any kind. Gaia carries 300 arcsec outwards. The two genuinely overlap at 300-340, where
both instruments have quality-selected stars at the same effective radius.

### Flagged stars only, and the fraction dropped

User: "if the unflagged stars should not be used we just used flagged stars and compute their
dispersion." Right, and `hst_profile` now defaults to `require_flag=True`,
`correct_unflagged=False`.

For the record, what the fraction was and why it existed. The dispersion formula never
contained it. It was a rescaling applied after the fit: a bin of flagged stars is unbiased, a
bin of unflagged stars reads 1.077 times high because their quoted errors are too small, and
a bin that is fraction *f* unflagged was assumed to sit linearly in between, so sigma was
divided by `1 + f (1.077 - 1)`.

The fraction only ever did work in bins straddling 340 arcsec, where the flag runs out
mid-bin, and those bins only existed because the HST annuli had been chosen to match Gaia's
edges rather than the flag's own boundary at 340. The goal behind all of it was to push HST
to 460 arcsec to meet Gaia. That goal was unnecessary: **Gaia starts at 300 arcsec, so
300-340 is already an overlap using flagged stars on both sides.**

The measurement, flagged stars only, nothing corrected:

| annulus | N | sigma (mas/yr) |
|---|---|---|
| 150-200" | 142615 | 0.6268 ± 0.0089 |
| 200-250" | 150092 | 0.5852 ± 0.0083 |
| 250-300" | 119576 | 0.5511 ± 0.0078 |
| 300-340" | 16912 | 0.5257 ± 0.0074 |

and the three comparison routes, which agree:

| route | HST | Gaia | Gaia/HST |
|---|---|---|---|
| 300-340", both flagged | 16912 at 311" | 17 at 318" | 1.028 ± 0.132 |
| 300-380", both flagged | 16978 at 311" | 52 at 356" | 0.933 ± 0.074 |
| 300-460", HST corrected | 92568 at 327" | 328 at 430" | 0.966 ± 0.038 |

`plots/hst_gaia_overlap.png` redrawn: HST from flagged stars to 340 arcsec as the
measurement, the corrected extension to 466 drawn faintly as a cross-check only. Eight tests,
one of which pins the default so the correction cannot creep back in.

### The periphery figure was still drawing the published HST profile

It stopped at 300 arcsec while Gaia's first point sits at 356, so the two appeared to meet
without overlapping. Replaced with our own flagged-star measurement, which reaches 340
arcsec, and the Gaia point for the *same* annulus (300-340, 17 stars, 0.5385 ± 0.0685 mas/yr)
is drawn alongside HST's last point (311 arcsec, 0.5257 ± 0.0074). The two now sit on top of
each other at 8.2 and 8.4 pc, with the annulus shaded.

The master figure still shows the published HST profile, correctly, because that is what the
likelihood is currently given. Extending the fitted HST dataset to 340 arcsec is the pending
decision.

### Our HST measurement runs the full radial range, and a correction to an earlier claim

It stopped at 150 arcsec in the figures only because I had run it on the four annuli needed
for the overlap test. There is no obstacle: the catalogue has 1.4 million proper motions from
1.8 arcsec outwards. Run over 3 to 340 arcsec against the published profile:

| r (arcsec) | N | ours | published | ratio |
|---|---|---|---|---|
| 3.8 | 146 | 0.7708 ± 0.0382 | 0.8048 | 0.958 |
| 8.3 | 595 | 0.7944 ± 0.0225 | 0.7854 | 1.012 |
| 18.0 | 3224 | 0.7654 ± 0.0108 | 0.7776 | 0.984 |
| 39.4 | 15574 | 0.7493 ± 0.0106 | 0.7440 | 1.007 |
| 85.9 | 66878 | 0.7070 ± 0.0100 | 0.7079 | 0.999 |
| 176 | 142615 | 0.6268 ± 0.0089 | 0.6236 | 1.005 |
| 270 | 119576 | 0.5511 ± 0.0078 | 0.5504 | 1.001 |
| 311 | 16912 | 0.5257 ± 0.0074 | 0.5199 | 1.011 |

**Median ratio 1.0042, scatter 2.5 per cent**, the scatter coming from the innermost bins
where only a few hundred stars survive the flag.

**Correction.** I reported earlier that our flagged measurement sat 1.1 to 2.3 per cent above
the published profile and called that a method floor. That was an interpolation error on my
side: I interpolated the published profile against its `r_lower` column instead of
`r_median`, which shifts the comparison outwards by half a bin on a falling profile. There is
no 1-2 per cent method floor. Our pipeline reproduces oMEGACat's own numbers to a few tenths
of a per cent wherever the statistics allow.

### Can the overlap be widened? Limits on both sides

| | HST flagged | Gaia usable |
|---|---|---|
| 250-300" | 119576 | 3 |
| 300-340" | 16912 | 17 |
| 340-360" | **66** | ~18 |
| 360-380" | **0** | ~17 |
| 380-466" | 0 | 276 |

The overlap is therefore **300-360 arcsec**, slightly wider than the 300-340 quoted before:
HST's flag survives to 360, on 66 stars. That bin gives HST 0.4912 ± 0.0347 at 341 arcsec
against Gaia 0.4375 ± 0.0402 at 363, a ratio of 0.891 ± 0.103, consistent with the other
routes.

Neither edge can be pushed further within these catalogues.

* **Outward (HST).** The flag requires F625W *and* F814W photometry, and beyond 360 arcsec
  no star in the catalogue has either. It is not a quality threshold that could be relaxed;
  the measurement simply does not exist. Using the unflagged stars there is the 8 per cent
  correction already rejected as a fitting input.
* **Inward (Gaia).** The quality flag passes 5 stars between 250 and 300 arcsec and none
  inside 250. Relaxing it is not an option either: the star-by-star HST comparison showed
  unflagged Gaia carries 0.59 to 0.90 mas/yr of undeclared scatter inside 300 arcsec, which
  is larger than the dispersion being measured.

The one untested route is **oMEGACat II** (arXiv:2404.03722), the standalone HST astrometric
catalogue, which we have not ingested. It is probably the same footprint, but its outer
photometric coverage has not been checked.

## 2026-09-19 -- master figure and the data-preparation write-up

`plots/master_datasets.png` regenerated with our own HST measurement drawn to 360 arcsec
alongside the published profile that the likelihood currently receives, and the 300-360
arcsec overlap shaded.

`docs/data_analysis.tex` written and compiled to `docs/data_analysis.pdf`, 10 pages, 6
figures, with a bibliography of verified links. Sections: the governing principle (a star
whose error approaches the dispersion measures the error model, not the cluster); the
catalogues; the estimator (2-D mixture, empirical field, exact perspective, depth, rotation);
HST (flag definition, validation against the published profile, the 360 arcsec limit and why
it cannot move); Gaia EDR3 (why the published spline is not used, the raw-error problem, the
failure of the density-only inflation, our low-noise measurement); cross-validation
(HST/Gaia overlap by three routes, Pristine at 0.995 +- 0.014); supporting datasets; the
Jacobi radius and the periphery overdensity; the 127-point dataset; and the three open
decisions.

LaTeX build artefacts are gitignored; the .tex and .pdf are tracked.

## 2026-09-19 -- Codex review: two real bugs in the estimator

Ran the Codex review on the data-analysis modules and the write-up (effort high, 419 s, diff
80 kB). Twelve findings. Two are genuine defects in the estimator, verified independently and
fixed; several are valid criticisms of the write-up, now corrected; the rest are recorded as
limitations. The review is kept at `docs/codex_review_2026-09-19.md`.

### Bug 1: wrong sign on the covariance cross-terms in the mean update

`outer_profile.py`, the weighted-least-squares step inside the EM loop. With
$J=[[\cos,-\sin],[\sin,\cos]]$ and $\Sigma^{-1}=(1/\det)[[c_{dd},-c_{ad}],[-c_{ad},c_{aa}]]$,
the normal equations $(\sum r J^\top \Sigma^{-1} J)u = \sum r J^\top \Sigma^{-1} x$ give a
**minus** sign on the `cad` terms in `A12` and `b2`. Both were positive.

Verified three ways before touching anything: hand derivation; a synthetic anisotropic sample
with streaming, where the code returned mean_T = -0.125 for an injected -0.200 while the
corrected algebra returned -0.205; and direct numerical maximisation of the same likelihood,
which returned -0.2048.

`cad` vanishes when the cluster is isotropic and the errors uncorrelated, which is why every
existing mock missed it. Effect on the real data:

* the **combined dispersion changes by at most 0.17 per cent** -- the quantity the likelihood
  currently receives is essentially unaffected;
* the **anisotropy** moves materially: T/R at 437 arcsec goes 0.867 to 0.823, at 533 arcsec
  0.898 to 0.865;
* the **fitted rotation** moves by up to 17 per cent.

The last one is an independent confirmation that the fix is right: our fitted mean tangential
motion against Vasiliev & Baumgardt's published rotation curve, which the estimator never
sees, had a median ratio of **0.90 before the fix and 0.977 after**.

### Bug 2: a 1.4 per cent floor on every quoted uncertainty

The profile-likelihood interval stepped by a fixed $0.02\sigma$ and returned the first point
past the $\Delta\ln L = 0.5$ crossing, with no interpolation. That is a floor of 2 per cent
per component, 1.4 per cent combined, independent of sample size. Every HST bin in the
write-up sat exactly on it: $0.0089/0.6268 = 1.42$ per cent, $0.0083/0.5852 = 1.42$, and so
on for $10^5$ stars.

Now bracketed with a geometrically growing step and interpolated on the likelihood drop.
Recovery on clean synthetic samples: reported error / analytic expectation = 0.98, 1.00, 0.96
at N = 2000, 20000, 100000. The HST uncertainties fall from 1.42 per cent to **0.12-0.39 per
cent**; the Gaia ones from 1.4-6 to 1.1-7 per cent (they were never at the floor).

### Findings accepted into the write-up

* **Pristine is not an independent astrometric check.** Its proper motions are Gaia's, and
  the samples share stars: 33 of the 52 selected Gaia stars in the innermost annulus, 11-37
  per cent further out (verified). It remains a valid test of selection and estimator. The
  ratio, recomputed after the fixes, is 1.005 +- 0.012.
* **The abstract overstated the overlap.** The clean, quality-selected comparison is
  1.03 +- 0.13; the 4 per cent figure belonged to the route using rejected HST stars, an
  empirical correction and a radius adjustment. Fixed, and the three routes are now stated to
  be nested rather than independent.
* **HST covariances are diagonal**, because the catalogue supplies no correlation column. The
  text claimed full covariances.
* **Position angle convention.** Our 148 deg is measured from east through north; the
  astronomical convention gives **122 deg**.
* **The systemic-motion sentence was wrong.** The 5.3 mas/yr artefact is what removing
  nothing produces. Subtracting the constant vector removes the leading term; the exact
  projection removes a residual of order 0.06 mas/yr at 2400 arcsec.
* **1.077 is a ratio of dispersions, not of errors.** Intrinsic and measurement variances add.
* **The 2 per cent bound is approximate**, since the threshold is on the mean of the two
  component errors and against a pilot profile.
* **The midpoint prescription is a sensitivity, not a calibrated systematic.**

New section *Known limitations* records the four that are not fixable today: the field
template does not follow the science selection (restricting it to the same error ceiling
moves the outermost bin by ~1.5 per cent); adding a constant does not restore HST's frame;
the error cut uses the published profile as a pilot and equipartition makes brightness a
kinematic selection; per-bin errors are treated as independent.

`tests/test_mixture_recovery.py` added: seven tests injecting streaming and anisotropy, a
comparison against direct numerical maximisation, and an $N^{-1/2}$ scaling check. The old
mocks all had zero streaming and isotropic dispersions, which is why nothing caught either bug.

### Two scientific conclusions weakened by the error-bar fix

Removing the 1.4 per cent numerical floor made the Gaia uncertainties honest, and two claims
in the journal that rested on the inflated bars no longer stand as stated.

**1. The outer residual is no longer "flat".** Over 460-1500 arcsec the residual against the
K1 model has a weighted mean of $+4.26$ per cent and $\chi^2 = 15.5$ for 6 degrees of
freedom, $p = 0.017$. With the old floor-limited errors it was consistent with a constant.
Statistically there is now 2-sigma-level structure. **It is not established**, because the
field template's selection mismatch moves the outer bins by about 1.5 per cent, which exceeds
the per-point statistical errors of 0.75-1.45 per cent. The test now pins the offset and
requires no single point to run away, instead of asserting flatness.

**2. The measured anisotropy is noisier at 434 arcsec than the old value suggested.** The
sign fix moved that bin from T/R = 0.867 to 0.813, and with 394 stars it is the noisiest
point in the range. The claim that the model and the measurement agree where HST anchors the
fit now holds over 500-800 arcsec (within 0.03), not from 400. The outer disagreement is
unaffected and if anything stronger: model 0.89 against measured 1.03-1.13 beyond 1150
arcsec.

Neither change touches the combined dispersion, which moved by at most 0.17 per cent.

### Before/after figure for the two estimator fixes

`plots/estimator_before_after.png`, built from `kinematics/estimator_audit.py`, which reruns
the real measurement with either defect restored through a `LEGACY` switch in
`outer_profile`. The comparison is therefore recomputed, not quoted from notes. Four tests
pin that the switch reproduces each defect and never leaks outside its block.

| panel | what it shows |
|---|---|
| dispersion, after/before | the quantity the likelihood receives moved by at most **0.17 %** |
| uncertainty vs N | the old 1.41 % floor, and the corrected values following $1/\sqrt{2N}$ |
| anisotropy | more radial at 400-700 arcsec: T/R at 437 arcsec from 0.867 to 0.823 |
| rotation / published curve | median over 356-1328 arcsec from **0.93 to 0.98** |

The last panel is the independent evidence: the published rotation curve is never used by the
estimator, and correcting the algebra moves our fitted rotation onto it.

Full suite clean at **348 passed** after the fixes and the test updates.

### Why HST carries no rotation of its own

Checked after the user asked why the rotation swap leaves HST untouched.

oMEGACat proper motions are **locally corrected**: an HST field is too small to hold enough
extragalactic anchors, so zero motion is defined as the mean motion of the cluster stars in
each patch. Within every patch the mean is then zero by construction, and the cluster's
streaming has been spent on setting the frame. Measured directly with our own mixture:

| annulus | N | our fitted \|mean_T\| | Vasiliev curve | ratio |
|---|---|---|---|---|
| 40-80" | 61850 | 0.0024 | 0.0543 | 0.04 |
| 80-130" | 110161 | 0.0064 | 0.0899 | 0.07 |
| 130-200" | 189791 | 0.0074 | 0.1385 | 0.05 |
| 200-260" | 180903 | 0.0076 | 0.1815 | 0.04 |
| 260-340" | 105677 | 0.0060 | 0.2086 | 0.03 |

**96 per cent of the rotation is gone**, as the procedure implies. So HST's rotation term
must come from outside, and `likelihood.py` already takes it from the Vasiliev & Baumgardt
curve. Our own curve cannot replace it: we fit rotation only from 356 arcsec outwards, while
HST ends at 360.

Consequence for the two comparison figures: `master_datasets_rotation_published.png` is
self-consistent, and `master_datasets_rotation_ours.png` is a hybrid -- our rotation on the
Gaia points, the published one on HST. Relabel before use.

Not a bug, checked: HST's stored `streaming2` is $v_{\rm rot}^2$ while Gaia's is
$v_{\rm rot}^2/2$. HST's dataset is the tangential component alone, Gaia's is the combined
dispersion, and the ratio between those conventions is exactly $\sqrt2$ -- which is why the
comparison came out at 1.41 at every radius.

## 2026-09-19 -- the dataset rebuilt: field template, HST to 360 arcsec, Gaia split

Three changes, all agreed beforehand.

**1. Selection-matched field template.** `field_density_2d` now accepts `err_max` and
`g_range`, and `build_edr3_profile` passes each annulus its own ceiling and magnitude range,
so the template describes the same population as the stars being fitted rather than a
fainter, worse-measured one. Eight of the nine bins get a matched template of 12000-45000
stars; the innermost (300-380 arcsec) falls back to the full template because only 783 stars
survive its tight cuts, and with a field fraction of 0.039 there it barely matters.

Effect: nothing inside 1000 arcsec, then +0.48 per cent at 1690 and **+1.66 per cent
(0.32 sigma) at 2148 arcsec** -- the size the Codex review predicted. That bin has 95 per
cent field, so it is where the template does the most work.

**2. HST extended to 360 arcsec.** `build_hst_profile` writes
`ocen_pm_dispersion_hst_ours.ecsv`: flagged stars only, no correction of any kind, 22 bins
per component from 3 to 341 arcsec (the outermost covering 340-360). This replaces the
published oMEGACat profile in the likelihood, which stopped at 300 and so never overlapped
Gaia. Its uncertainties run 0.1 to 4.5 per cent against the published profile's uniform 1.4.

**3. Gaia split into components.** Checked first whether that is legitimate: over 60
independent synthetic realisations with uniform position angles the correlation between the
fitted sigma_R and sigma_T is **-0.076**, so treating them as two datasets costs about 0.6
per cent in the joint chi-squared. The product now carries `sigma_pmr_err` and
`sigma_pmt_err`, each the midpoint of the two error models with half their separation added.

`DEFAULT_DATASETS` is now
`hst_pm_radial_ours, hst_pm_tangential_ours, muse_los_dispersion, gaia_dr2_pm,
gaia_edr3_ours_radial, gaia_edr3_ours_tangential`, **100 points**, down from 127 but each far
better measured. Every incompatible combination is forbidden in `_FORBIDDEN_TOGETHER`.

The write-up gained a **Rotation** section and the new dataset table, and the master figure
now draws the published HST profile and Gaia spline for context with the fitted data on top.

No fits have been run.

## 2026-09-19 -- the cluster centre was wrong by 10.7 arcsec

Found by the second Codex review, verified three ways, fixed.

### What it was

Four modules hard-coded `(201.696833, -47.476583)` with the comment "Baumgardt catalogue
centre". **`data/raw/baumgardt_gc_catalogue/` is empty** -- that catalogue was never
downloaded. The number came from nowhere: written from memory on 2026-09-17, given a
provenance it did not have, and copied into four files rather than shared.

Both catalogues we use carry their own centre, and neither was consulted:

| source | centre | offset from ours |
|---|---|---|
| oMEGACat, from its pixel grid at (15000, 15000), 272 stars | 201.696833, **-47.479569** | 10.75" |
| Vasiliev & Baumgardt, from the member table's `x, y` | 201.696838, **-47.479339** | 9.92" |

Those two agree with each other to **0.83 arcsec**. A fifth module,
`plotting/data_overview.py`, had `-47.4795` -- the right value -- so the codebase held two
different centres at once.

The symptom Codex reported and I reproduced exactly: the 61 stars labelled as the 2.55-3.24
arcsec bin were really at **7.80-13.90 arcsec** from the true centre.

### The fix

New module `ocen_dm/cluster.py` holds the centre once and `centre()` derives it from either
catalogue. All five modules import from it. Adopted the oMEGACat value: it is set by the
densest, best-measured field in the core, and the 0.83 arcsec residual is far below any bin
width. Six tests in `tests/test_cluster_centre.py`, one of which greps the source tree so the
literal cannot be re-introduced.

### Where it mattered

**Gaia: nowhere.** Every bin moved by 0.5 sigma or less; a 10.7 arcsec shift at r > 300
arcsec is at most 3.6 per cent in radius on a shallow gradient.

**HST: in the inner bins, and less than the first look suggested.** Matched on bin edges the
shifts are at most 2.6 sigma and mostly under 1.5:

| bin | N before -> after | sigma_R before/after | sigma_T before/after | shift |
|---|---|---|---|---|
| 4.1-5.3" | 152 -> 166 | 0.7917 / 0.8265 | 0.7405 / 0.8628 | +0.8 / +2.6 sigma |
| 5.3-6.7" | 226 -> 275 | 0.7786 / 0.8426 | 0.8478 / 0.8271 | +1.8 / -0.6 |
| 13.8-17.6" | 1532 -> 1522 | 0.7611 / 0.7909 | 0.7631 / 0.7688 | +2.1 / +0.4 |
| 22.4-28.5" | 4191 -> 4155 | 0.7660 / 0.7580 | 0.7765 / 0.7560 | -1.0 / -2.5 |
| 340-360" | 66 -> 175 | 0.5109 / 0.5309 | 0.4706 / 0.4648 | +0.7 / -0.2 |

My first comparison table, matching row by row, showed shifts up to 54 sigma. That was
wrong: the two products have 22 and 21 bins, so the rows do not correspond. Matching on bin
edges is the right comparison.

The instrument overlap improves slightly: **1.011 +- 0.123** at 300-340 arcsec (was 1.028 +-
0.125), 0.965 +- 0.069 at 300-380, 0.962 +- 0.029 for the corrected route.

All products, the match and every figure regenerated.

## 2026-09-20 -- the six review findings, applied

### 1. Field template: parallax cut removed, template re-queried

The query required `parallax_over_error < 5`; the science selection does not, so the template
described a different population than the stars it was classifying. Cut removed, field
re-queried: **192 054 -> 231 269 stars**.

Effect on the dispersions: at most **0.33 sigma** (outermost tangential bin), nothing above
0.15 sigma inside 1300 arcsec. Smaller than feared. My earlier inference about the direction
of the bias was also wrong: I argued from the *science* sample split by parallax that the
template overstated the field density at the cluster's proper motion by 2.4x, but measured on
the template's own annulus the fraction within 1 mas/yr of systemic is **0.53 per cent either
way**. The two populations differ in the science annulus, not in the template's.

### 2. Rotation uncertainty now carried

The curve is published with percentiles and we used only the median. An error `dv` moves the
dispersion the model must predict by `v dv / sigma`; that is now added in quadrature.

| dataset | bin | stat error | total after |
|---|---|---|---|
| HST tangential | 223.5" | 0.00088 | 0.00434 (x4.9) |
| HST tangential | 270.5" | 0.00105 | 0.00604 (x5.7) |
| HST tangential | 310.8" | 0.00244 | 0.00769 (x3.1) |
| Gaia tangential | 436" | 0.01832 | 0.02486 (x1.36) |

This was invisible while every error sat on the 1.4 per cent floor. It is a **floor** on the
right correction: the rotation curve moves coherently between bins and is treated here as
independent. A shared nuisance parameter would be correct and is a model change.

### 3. Model averaged over the stars actually measured

The likelihood averaged over complete annuli with `Sigma(R) R` weighting while the data
average a selection whose coverage changes within the bin. Both products now store eight
equal-count quantiles of the selected stars' radii, and `BinnedProfile.r_nodes` makes the
likelihood average over those. HST's 300-340 arcsec bin has a star-weighted mean radius of
**313.2 arcsec against the annulus midpoint of 320**. Effect on the prediction: up to
**1.17 sigma** in the outer HST bins, as the review estimated.

### 4. Shared calibration -- documented, not fixed

The error-model systematic responds to one calibration choice and so moves every bin
coherently, but it is added per bin in quadrature and consumed by a diagonal likelihood. The
correct treatment is a shared nuisance parameter, which is a model change. Recorded in the
write-up's limitations.

### 5. Rotation figures corrected

They subtracted streaming from data that were already dispersions about fitted means, and
used the **combined** term for both components. Now they add the streaming back, per
component, to show the second moment the model predicts. The radial term is ~0.002 and the
tangential ~0.03-0.07 mas^2/yr^2, so the old figures moved the radial points by half the
tangential streaming for no reason.

### 6. Evidence comparisons bound to the observations

`comparison_table` matched dataset *names*, so it would print a Bayes factor between a real
run and a mock. It now fingerprints the datasets, the point count, the input hashes from
`run.yaml` and the mock provenance. Verified: real vs mock now prints "n/a (different data)"
plus an explicit warning; real vs real still gives Delta ln Z = -51.01.

### Also

* the failing point-count test fixed: **98**, not 100 -- the centre fix pushed two HST bins
  below the minimum star count;
* the tautological centre test replaced by one that matches the loader's radii star by star
  against radii computed independently from the catalogue, and checks the old centre would
  have failed it;
* the centre literals carried to 8 decimals so the loader matches the derived value to a
  milliarcsecond;
* the three products predating the centre fix quarantined under
  `data/processed/kinematics/_stale_pre_centre_fix/` with a README, rather than deleted, so
  the older figures stay reproducible.

## 2026-09-20 -- three decisions taken by the user

### 1. Error inflation is an option, not a default

The product averaged the raw-error and eta-inflated fits and carried half their separation as
a per-bin systematic. That hid the choice and treated one coherent calibration decision as
independent noise. Now **the default is `raw`** -- the catalogue's own uncertainties,
unmodified -- and both fits are stored per component (`sigma_pmr`, `sigma_pmr_eta`, ...).
`sigma_sys` is still recorded but **no longer added to the quoted error**.

### 2. The published rotation curve is the default

Gaia's streaming term was our own fitted mean. It is now the **published curve with its
percentiles propagated**, matching how HST has always been treated, with our own means
available as the alternative. The radial component carries no published rotation (their
radial mean is fixed to perspective expansion), so its streaming is zero by default and
non-zero under `rotation="ours"`.

### 3. Gaia DR2 dropped

Four of its nine points sit inside 380 arcsec, the range whose EDR3 equivalent we removed as
extrapolation, and we hold no star list for DR2 so the same check is impossible. Out of the
default; still loadable by key.

### How the variants are run

Both choices are switches rather than assumptions:

```
ocen-dm fit --family K2-cored                                  # raw errors, published rotation
ocen-dm fit --family K2-cored --gaia-errors eta                # inflated errors
ocen-dm fit --family K2-cored --gaia-rotation ours             # our rotation curve
```

The default dataset is now **89 points**: HST 21+21, MUSE 29, Gaia EDR3 9+9. Five tests in
`tests/test_gaia_options.py` pin the defaults, that the product carries both error models,
that each switch moves what it should and nothing else, and that DR2 is out of the default
but still loadable.

The write-up gained a *Choices left explicit* section and is 13 pages.

## 2026-09-20 -- mass-modelling write-up

`docs/mass_modelling.tex` -> `mass_modelling.pdf`, 5 pages, written from the code rather than
from memory (every parameter, prior and equation checked against `mass_models/`,
`kinematics/jeans.py`, `kinematics/anisotropy.py` and `kinematics/fit.py`).

Sections: how the question is posed (K1 nested in K2, two halo slopes run separately rather
than fitting gamma); the four mass components and why the remnants are a separate component
rather than an inflated mass-to-light ratio; the anisotropy family and the An & Evans
constraint; the Jeans solver with its projections, the `r = R cosh u` substitution and the
reverse-accumulation fix; the split-normal likelihood with radial averaging, streaming and
instrument scales; the prior table; nested sampling and the evidence-comparability rule;
the four cross-checks (JamPy, an AGAMA DF, literature presets, injection); five stated
omissions; and the table of runs.

The omissions are stated plainly because each is a route to a spurious halo: sphericity,
rotation removed rather than modelled, equilibrium at 0.7 r_J(peri), the diagonal likelihood
against coherent systematics, and a single tracer population under mass segregation.

## 2026-09-20 -- mass-modelling review: the four fixes agreed

### 2. The likelihood was not a normalised density

Each side of the split normal carried its own `-ln(err)`, so the density **jumped by
`err_hi/err_lo` as the model crossed the datum** -- a factor of 4 for errors of 0.5 and 2,
reproduced exactly. It is now the standard two-piece form with the shared normalisation
`sqrt(2/pi)/(err_lo + err_hi)`: continuous, integrates to 1.000000, and reduces to a plain
normal when the errors are equal. MUSE has asymmetric errors, so this was live.

### 3. The evidence fingerprint missed the switches I had just added

I wrote the fingerprint in the morning and added `--gaia-errors` and `--gaia-rotation` in the
afternoon without going back. Neither was recorded in `run.yaml`, so two runs differing only
in the error model fingerprinted identically. `run_nested` now records a `dataset_options`
block (both switches, the tracer, the backend and the dataset list) and the fingerprint
includes it, plus the mock seed, family and generating parameters -- the previous version
read a `mock_from` field that does not exist.

### 5. The AGAMA prior ran outside what AGAMA accepts

Its Cuddeford DF is undefined below `beta_0 = -0.5` and raised an uncaught `RuntimeError`
mid-run. The prior for that backend is now bounded at `AGAMA_BETA0_MIN = -0.5`; the Jeans
backend keeps `[-1, 0]`.

### 8, 9. Two things the write-up got wrong

The default tracer is **not** a surface-brightness MGE: it is HST star counts inside 25
arcsec spliced onto Trager light outside, and the same shape supplies both the tracer density
and the stellar mass density, which assumes a radius-independent M/L. And the radial
averaging and rotation propagation apply only to the four profiles we measured: **MUSE has
neither**. Its rotation and dispersion were fitted jointly by the survey, so it needs their
covariance rather than an independent error; the first-order rotation contribution in its
innermost bin is 1.22 km/s against quoted errors of 1.90/2.52. Recorded as outstanding.

Ten tests in `tests/test_likelihood_normalisation.py`.

Still open, and deliberately not touched: the An & Evans constraint with a central point mass
(1), the monotonic anisotropy family (4), the AGAMA positivity claim (6) and the coherent
rotation uncertainty (7).

## 2026-09-20 -- the modelling plan: a ladder, not a bigger model

The user, after three review rounds each of which proposed more parameters: "I am now worried
that we are overcomplicating the model, this seems to be the case of premature optimisation
that I would like to avoid at all costs. Instead, we need a simple and robust set of initial
models to test." Correct, and the state of play makes the case: the data are finished, the
machinery has been reviewed three times, and **no fit has been run on the current dataset at
all**.

`docs/MODELLING_PLAN.md` is the plan. Three rungs, same two families on each, climbed only
when the rung below says to:

| rung | anisotropy | scales | K1 | K2 | points/param (K2) |
|---|---|---|---|---|---|
| 1 | constant beta | none | 6 | 8 | 11.1 |
| 2 | beta(r), 3 params | none | 8 | 10 | 8.9 |
| 3 | beta(r) | per instrument | 11 | 13 | 6.8 |

The decision to climb rests on the residuals of the rung below, not on a reviewer's
suggestion. A Delta ln Z that is **stable** across rungs is itself the finding: it says the
answer is not about orbital freedom.

### Taken now, because they remove freedom

* **`beta_0` narrowed to `U[-1, -0.5]`** in both families. An & Evans give `beta <= gamma/2`
  for a self-gravitating cusp, hence `beta_0 <= 0` for a cored tracer, but `gamma >= beta +
  1/2` in a point-mass-dominated potential, hence `beta_0 <= -1/2`. My counter-argument --
  that the sphere of influence is unresolved -- **fails**, and I checked it: `r_infl` reaches
  the innermost datum at 0.11 pc once `M_bh > 7400 Msun`, which is 46 per cent of the
  log-uniform prior.
* **The AGAMA positivity claim dropped.** AGAMA clips negative DF values to zero, so
  construction does not certify positivity. Worse, now that the ceiling is -0.5 and AGAMA's
  floor is also -0.5, **the two ranges meet at a single point**: the AGAMA backend cannot
  represent a physically admissible cored-tracer-plus-black-hole model at all. It keeps its
  own range, is a diagnostic only, and its evidence is never compared.

### Deferred, with the trigger written down

The anisotropy turnover (one parameter, revisit if a rung leaves an anisotropy-shaped
residual) and the shared rotation nuisance (revisit if the outer tangential bins drive the
result). Neither is wrong; neither has been shown to be needed.

### Also

`--constant-beta` was not exposed on the CLI, so rung 1 was not runnable as written. It is
now, and both it and `--no-scales` are recorded in the run provenance. Three tests updated
for the narrowed prior: the unit-cube midpoint now maps to -0.75, and an injection test used
a truth outside the new prior.

### Review of the plan: rung 1 was broken by my own prior change

User: "Does any of these apply to the first mass modelling experiment?" Of the four review
items, only the anisotropy prior touches the first round, and I had applied it wrongly.

An & Evans is a condition **at r -> 0**. With `beta(r)`, `beta_0` is the central value and
the bound belongs there. With a **constant** `beta`, `beta_0` is `beta` everywhere, so the
ceiling of -0.5 I set in the morning had become a global constraint: tangential dispersion
at least 22 per cent above radial at every radius. The data never show that -- projected
`sigma_T/sigma_R` runs 0.84 to 1.18, i.e. `beta` roughly -0.4 to +0.3. Rung 1 as I left it
forbade every anisotropy the data display.

Fixed: the bound applies only to the radially varying family; constant `beta` gets
`CONSTANT_BETA_RANGE = (-1, 0.5)`, wide enough to hold what the data imply with margin.

### The first round, defined

Added **rung 0, isotropic**: `beta = 0` fixed, no scales, K1 with 5 parameters against
K2-cored with 7, plus K2-NFW to bracket the halo shape. It is the model everyone
understands, its Bayes factor has the simplest meaning, and its residuals answer directly
whether anisotropy is needed at all. `--isotropic` added to the CLI and recorded in the
provenance.

| rung | anisotropy | `beta` prior | K1 | K2 | points/param (K2) |
|---|---|---|---|---|---|
| 0 | isotropic | fixed 0 | 5 | 7 | 12.7 |
| 1 | constant | `U[-1, 0.5]` | 6 | 8 | 11.1 |
| 2 | `beta(r)` | `beta_0 ~ U[-1, -0.5]` | 8 | 10 | 8.9 |
| 3 | `beta(r)` + scales | as 2 | 11 | 13 | 6.8 |

Smoke-tested: every rung's two families build and return a finite likelihood in **4 ms**
per evaluation at three points of the unit cube. `tests/test_modelling_ladder.py` pins the
parameter counts, the nesting of K1 in K2, that the central bound touches only the radially
varying family, and that the CLI exposes every rung. 34 fit-side tests pass.

**No fit has been run.** The first round is three commands away:

```
ocen-dm fit --family K1       --isotropic --no-scales
ocen-dm fit --family K2-cored --isotropic --no-scales
ocen-dm fit --family K2-nfw   --isotropic --no-scales
```

### Instrument scales removed from the ladder

User: "this sounds like a fudge to me." Agreed, and the case is stronger now than when the
scales were introduced. They are degenerate with mass by construction, and the data
preparation has made them unnecessary: HST/Gaia agree at 1.01 +- 0.12 in the overlap and
Pristine/Gaia at 1.005 +- 0.012 with no scale anywhere. A remaining disagreement would be a
finding about the data, not something to absorb. The MUSE offset is likely physical
(equipartition, different tracer) and deserves a model, not a factor.

Rung 3 deleted. Scales survive only as a **one-off diagnostic** after the first round: fit
once with them free and check they land at 1.00; a pull-away is a red flag about that
instrument. That run's evidence is never compared.

## 2026-09-20 -- a second route to the DM question: simulate the host's disruption

User's idea, recorded while the rung-0 fits run. Constrain omega Cen's dark-matter content
**from the other end**: instead of inferring it from the present-day kinematics, simulate the
tidal disruption of the host dwarf and see how much DM the surviving nucleus is left with.
Before any simulation, build a small set of distinct progenitor orbit histories. Three
classes:

1. **Current orbit, static axisymmetric potential.** The dwarf was always on today's orbit.
   In McMillan (2017): r_peri = 1.59 kpc, r_apo = 6.99 kpc, e = 0.63, from the systemic
   PM (-3.257, -6.730), v_los = 232.7 km/s, D = 5.43 kpc (computed 2026-09-19 for the Jacobi
   radius).
2. **Same, with dynamical friction integrated back in time.** The dwarf came in on a wider
   orbit and sank; the tidal history is gentler early on.
3. **Deposited by GSE, then migrated inward by the bar** -- Dillamore, Zhang & Belokurov
   2026, arXiv:2606.12516, "Bar-induced migration of omega Centauri away from Gaia
   Sausage-Enceladus". Their result: omega Cen *can* be traced back to the GSE phase-space
   region under a decelerating bar, but only for a present-day pattern speed
   Omega_b <~ 26 km/s/kpc, well below most current estimates.

Discussion in the reply of the same date. Not started.

## 2026-09-20 -- progenitor orbits: assumptions investigated, 3 + 3 + 3 orbits picked

User: "investigate the plausible assumptions about the initial conditions for the dwarf
(important for class 2 and 3) and come up with a small set of most plausible orbits in each
class (say 2 or 3)". Full note: `docs/PROGENITOR_ORBITS.md`; figure `plots/progenitor_orbits.png`;
table `results/tails/progenitor_orbits_summary.ecsv`; code `src/ocen_dm/tails/progenitor_orbits.py`,
`python -m ocen_dm.tails.products`; tests `tests/test_progenitor_orbits.py` (6 pass).

### Bugs found on the way

* galpy's `solarmotion` is the *peculiar* solar motion; I had added v_0 = 233 km/s to V, which
  made the orbit prograde with apo 10.6 kpc. Caught by comparing with an astropy transform
  (now a test). Correct: L_z = -529 kpc km/s, retrograde, E = -1.85e5 km^2/s^2 (McMillan 2017).
* galpy checkout (~/Work/src/galpy, 2024-03) needs two import shims: astroquery's broken
  version metadata and scipy's removed `vectorize1`. Both live in `_import_galpy_safely`.
* galpy's friction force is Python-level: ~30 min per orbit. Replaced by a leapfrog with the
  AGAMA McMillan 2017 host and the same Chandrasekhar formula (~1 s); frictionless run
  reproduces class 1 (1.63/7.04 vs 1.57/7.04 kpc).

### Class 1
McMillan17 1.57/7.04 (e 0.63), MWPotential2014 1.97/6.76 (0.55), Irrgang13I 1.25/7.20 (0.70).
All three kept: the spread is the potential spread.

### Class 2 -- the key assumption is the mass history, not the mass
Constant-mass backward friction is wrong in the important direction (a constant 1e10 Msun
satellite is pumped to 40-140 kpc). Adopted: exponential stripping since infall 10 Gyr ago,
M(t) = M_inf exp(-(t_inf - t)/tau), floored at the nucleus mass, r_h = 1 kpc (M_inf/1e10)^(1/3),
sigma = v_c/sqrt(2). Selection rule: the orbit must be at 50-150 kpc (virial radius of the
young MW) at infall. Scan over M_inf = 1e10 ... 2e11 and tau = 0.5 ... 3 Gyr shows a 1e10 halo
needs slow stripping (tau >= 1.5 Gyr) to have sunk from the virial radius; a GSE-mass halo
needs fast stripping (tau <= 1 Gyr). Picked: (1e10, 2.0) r_inf 65 kpc; (3e10, 1.25) 85 kpc;
(1e11, 0.75) 111 kpc. In all three the orbit is essentially today's for the last 3-5 Gyr and
apo 10-20 kpc at 5 Gyr ago.

Literature for M_inf: NSC-host relation M* ~ 1e9 (Pfeffer+2021); oMEGACat X 4.5e9 "GE-like"
dwarf mass (Souza+2026, arXiv:2603.23589), favouring Sequoia/Thamnos debris; GSE 2e11 halo
(Naidu+2021).

### Class 3 -- GSE-debris orbit before bar migration
Dillamore+2026: retrograde 1:1 resonance of a decelerating bar (eta = 0.003 from ~45 km/s/kpc)
drags omega Cen to lower E and more retrograde L_z; works only for Omega_b <~ 26 today. The
bar is not modelled here; the class-3 initial condition is the debris orbit: apocentres from
the Belokurov+2023 chevrons (11.5, 15.5, 21 kpc), L_z = -300 (less retrograde than today),
inclination 60 deg -> peri 0.6 kpc, e 0.89-0.95, z_max 1.5-2.5 kpc. Pericentre is set by
L_z in the flattened potential and hardly depends on the inclination. Two readings of the
class recorded (omega Cen = GSE's nucleus vs a smaller host in group infall with GSE).

### Still running
The galpy `ChandrasekharDynamicalFrictionForce` cross-check of the fast integrator (M = 3e9
constant, 5 Gyr; fast result peri 1.84 / apo 47.1 kpc) -> `results/tails_galpy_friction_check.log`.
Rung-0 fits (started 03:03, 3 processes) are still in ultranest's sampling phase.

## 2026-09-20 -- class 3 done properly: the Dillamore+2026 bar migration reproduced

User: "yes i want class 3 done properly, reproduce their setup". The static GSE-debris guesses
(L_z = -300) are superseded.

**Key discovery:** `~/Work/Code/oCen_bar` is a clone of Dillamore's own repository
(github.com/adllmr/oCen_bar; authors Dillamore, Belokurov, Zhang) with the scripted pipeline:
Hunter+2024 AGAMA potentials, slowing-bar construction, 10^3 omega Cen samples, backward
integration over an Omega_b,0 grid, Belokurov+2023 GSE contours (`artifacts/gaiadr3_gse_elz.fits`).
The Hunter potential files and contours are referenced from that clone (`OCEN_BAR_DIR`), not
copied. New module `src/ocen_dm/tails/bar_migration.py`; products `python -m ocen_dm.tails.bar_migration`;
tests `tests/test_bar_migration.py` (4) + 1 in `test_progenitor_orbits.py`; note rewritten in
`docs/PROGENITOR_ORBITS.md`.

Set-up (paper + code): bar amplitude Dehnen switch-on 0-1 Gyr, deceleration onset 1-2 Gyr, then
eta = 0.003 to t_f = 8 Gyr; length scales as Omega_b,0/Omega_b(t); Omega_1 = 45.1 for Omega_b,0 = 24
(paper's "~45" -- confirms t1 = 1, t2 = 2 rather than the code defaults t1 = 2, t2 = 3 in
`make_slowing_bar_potential`); astropy default solar frame; bar angle 28 deg.

Results (fraction of the 1000 samples inside the outer GSE contour 8 Gyr ago):
20: 0.44, 22: 0.58, 22.5: 0.85, 23: 0.01, 24: 0.89, 25: 0.79, 26: 0.18, >= 27: 0. Matches the
authors' stored grid (0.87 / 0.69 / 0.08 at 24 / 25 / 26, dip at 23 too) and is seed-stable: the
jagged Omega dependence is resonant-phase sensitivity. Paper's "Omega_b,0 <~ 26" reproduced.
Migration is late: E, L_z sit at GSE values from 8 to ~2.5 Gyr ago and move only in the last
~2.5 Gyr. Picked (Omega_b,0 = 24, E0 quantiles 0.16/0.5/0.84 among successful samples):
pre-migration peri/apo 0.40/10.9, 0.42/11.6, 0.45/12.2 kpc, e 0.93, L_z(t=0) = +69, +163,
+229 (slightly prograde, not retrograde as guessed). Pericentres need the 0.25-Myr
re-integration (`refine_trajectory`): the 20-Myr grid cadence gave 0.5-0.9 kpc, up to 2x too
large, because a pericentre passage at 0.4 kpc lasts ~1 Myr.

Lesson re-learnt: `pytest ... | tail -1` hid a failing test and the previous commit (d7c49be)
went in with 1 failure; fixed here and the chain now uses `set -o pipefail`.

Side result relevant to every class: in the barred Hunter24 potential at Omega_b,0 = 24
(corotation ~9.5 kpc) omega Cen's present orbit has peri/apo 0.81/9.3 kpc over the last Gyr,
against 1.56/7.02 axisymmetric. A slow bar changes "today's orbit".

Also closed: galpy `ChandrasekharDynamicalFrictionForce` cross-check of the fast leapfrog
(M = 3e9 constant, 5 Gyr): galpy peri/apo 1.79/45.95 vs fast 1.84/47.1 kpc (3%), 482 s vs 1 s.

## 2026-09-20 -- LaTeX write-up of the progenitor orbits

User: "produce a latex write-up explaining the model construction for each Class, the
assumptions, analyze the resulting models and compare between them, have sufficient number of
informative plots illustrating the orbits". -> `docs/progenitor_orbits.tex/.pdf` (12 pages,
9 figures, 6 tables). New figure module `src/ocen_dm/tails/writeup_figures.py`
(`plots/po_*.png`, `results/tails/po_comparison.ecsv`, `po_class2_scan.ecsv`).

Fixes made while producing the figures: the fast class-2 integrator stored v_T in the
right-handed (astropy) sense, opposite to galpy's prograde-positive convention used by class 1
-> flipped, so all OrbitSummary objects are prograde-positive; the class-3 L_z was needlessly
sign-flipped in the E-Lz figure (the oCen_bar frame is already prograde-positive).

Comparison numbers (per-Gyr elements, McMillan17 enclosed mass for all): last Gyr, classes 1
and 2 identical (peri 1.6, apo 7.0, T = GM/r^3 = 1.4e4 (km/s/kpc)^2); class 3 peri 0.8, apo
9.3, T = 4e4 because the slow bar is part of its potential. 8-9 Gyr ago: class 2 peri 5-6 kpc,
T ~ 1.5e3 (ten times weaker than today); class 3 peri 0.4 kpc, T ~ 1e5 (seven times stronger).
L_z early: class 2 -> -2000..-4000 (more retrograde), class 3 -> 0..+230 (less).

## 2026-09-20 -- does a GC + DM remnant migrate through the bar like a particle? (friction test)

User asked whether a nucleus with a bound DM remnant is affected by the bar resonance exactly
like any particle; answer: the centre of mass does, to O((r/R)^2), except for mass-dependent
friction, non-coherent debris, and phase-locked tidal shocks. User: "lets code and run the cheap
decisive test". `src/ocen_dm/tails/friction_test.py` (vectorised leapfrog in the time-dependent
AGAMA barred potential + Chandrasekhar friction, constant bound mass = upper bound);
tests `tests/test_friction_test.py` (2). Calibration: nucleus mass gives dLz/dt = 11 kpc km/s/Gyr,
dE/dt = 7e2 (paper: 9.6, 6.3e2); frictionless leapfrog reproduces agama.orbit.

Result (fraction inside GSE 8 Gyr ago, Omega_b,0 = 24): particle 0.87; 1e7 Msun 1.00 (endpoint
shifted, dE +0.13e5, dLz +240); 1e8 Msun 0.00 -- never captured by the resonance, L_z drifts
more retrograde; 1e9 Msun 0.00, apo 74 kpc. Omega grid 20-27: 1e7 keeps the <~26 window; 1e8
never exceeds 0.33 and is 0 for 24-27. Threshold 1e7-1e8 Msun (higher for a stripping system).
Conclusion for the DM question: a class-3 nucleus that migrated from GSE debris did so
essentially without a dark envelope; the class-3 initial condition is only self-consistent if
the envelope was gone before ~2.5 Gyr ago. Written into docs/progenitor_orbits.tex (Sec. 5.4,
Table 6, Fig. 7) and docs/PROGENITOR_ORBITS.md.

Also: stale point count (98 -> 89) in tests/test_master_plot.py fixed; full suite otherwise
392 passed.

## 2026-09-20 -- class 2 in a barred potential

User: "For Class 2, can we investigate if the orbits change at all if the integration is done
in a potential with a rotating bar?" -> `src/ocen_dm/tails/class2_barred.py`, figure
`plots/class2_barred.png`, table `results/tails/class2_barred.ecsv`. The three class-2
histories re-integrated (20 error samples, backward leapfrog with friction) in Hunter+2024:
axisymmetric, and with the Dillamore+2026 growing/decelerating bar ending at 37.5 (Hunter's
value) and at 24 km/s/kpc. The bar is 8 Gyr old; earlier the potential is axisymmetric (checked
that the AGAMA scale modifier holds A = 0 before t = 0).

Result: with the mainstream bar (37.5) nothing changes beyond the sample scatter -- last-Gyr
peri 1.9 vs 1.7 kpc, apo 7.0-7.3 in both, 5-6 Gyr elements within 10%, apocentre before infall
80-99 vs 73-94 kpc (16-84% ranges overlap). With the slow bar (24) the orbit is resonantly
perturbed as in class 3: last-Gyr peri 0.85-1.07 and apo 10-12 kpc, 5-6 Gyr peri 0.9-2.1 vs
1.9-3.4, infall apocentre 95-112 vs 73-94 (still inside the 50-150 selection band). So the
class-2 histories are robust to a fast bar and change only under the slow-bar hypothesis, in
which case the last ~2.5 Gyr look like class 3 (plunging pericentres). The Hunter axisymmetric
baseline itself reproduces the McMillan17 class-2 picks (infall 67/94/79 vs 65/85/111 kpc).

## 2026-09-20 -- resonances of the present orbit (added to the write-up)

`src/ocen_dm/tails/resonances.py`: Omega_r 69.7, Omega_z 63.1, Omega_phi -44.3 km/s/kpc for
today's orbit (Hunter24 axi; McMillan17 70.2, 62.8, -44.0). Omega_b = Omega_phi +
(l Omega_r + n Omega_z)/m: corotation and OLR-type negative (retrograde orbit); retrograde 1:1
at 25.4 (the Dillamore resonance); next physical ones 42.8 (5:4), 48.6 (4:3), 60.3 (3:2). At
33-41 km/s/kpc the present orbit is in no low-order resonance. Written into
docs/progenitor_orbits.tex (Sec. 5.5) together with the class-2-in-bar result (Sec. 5.6,
Fig. 8); test added.

## 2026-09-20 -- which N-body code for the disruption simulations

User: "lets establish how as in what code to use ... I have ../satellite_experiment ... Report
back with your findings". Full assessment: `docs/NBODY_CODE_ASSESSMENT.md`; environment helper
`bin/nemo_env.sh`. Summary: satellite_experiment is a July scaffold (Plummer satellite, analytic
host, gyrfalcON driver) with pipeline-test values; the NEMO/gyrfalcON build behind it works after
three environment fixes (dylib paths, SIP stripping DYLD_* through /usr/bin/time, a stale
dlerror() from LLVM libomp that NEMO's loadobj misreads -> stub dylib). Benchmarks: 1e5
particles 0.08 s/step, 1e6 0.7 s/step (serial). The AGAMA plug-in runs with the rotating
Hunter+2024 bar, so the decelerating-bar potential of class 3 can be used directly.
Recommendation: gyrfalcON + AGAMA plug-in for round one (~5-8 h per 10-Gyr run, 1e6
particles), IC generation with AGAMA self-consistent models, pyfalcon as fallback for
per-particle friction control; Gadget-4 only if serial speed becomes limiting. A symlink
$NEMOOBJ/acc/agama.so -> venv agama.so was created; the stub lives in ~/.cache/ocen_dm.

## 2026-09-20 -- Codex review of the N-body plan: a real bug and several corrections

User: "use ask codex skill to discuss (discuss means not accept but evaluated, agree,
disagreed, reject too)". Review saved as `docs/codex_review_nbody_2026-09-20.md`; my verdicts in
the reply of the same time. The big one: **the fast leapfrogs divided the acceleration by
1.02271 instead of multiplying** ((km/s)^2/kpc = 1.0227 (km/s)/Gyr), kicks 4.4% too weak,
pericentre 4% too large -- exactly the "AGAMA vs galpy" difference I had explained away.
Verified: the corrected integrator reproduces agama.orbit to four figures (1.5730/7.0367 vs
1.5724/7.0368). All AGAMA leapfrogs now run in natural units (kpc/(km/s) = 0.977792 Gyr), which
also removes the Gyr-vs-natural-unit clock mismatch Codex flagged in friction_test/class2_barred;
bar_migration times are documented as the paper's kpc/(km/s) units (t_f = 8 = 7.82 Gyr).

Regenerated: class-2 scan (values shift 10-20%, picks still satisfy 50-150 kpc: 65/86/77 kpc);
friction test (1e7: window intact, 0.8-0.98; 1e8: marginal 0.3-0.5 at Omega 22-24, 0 at >=25 --
previously 0 everywhere; 1e9: 0, apo 79 kpc); class2_barred (same conclusions; slow-bar
last-Gyr peri now 1.3-1.7 rather than 0.85-1.07); products, write-up figures, LaTeX tables and
prose, PROGENITOR_ORBITS.md, NBODY_IC_PROPOSAL.md. Galpy cross-check of the corrected
integrator: peri 1.786 = galpy 1.786, apo 44.3 vs 46.0.

Other Codex points confirmed: falcON kernel default P1 and pair-softening mean -> rigid nucleus
must be an explicit AGAMA potential; Boldrini+2020 use live GCs (my citation wrong; Meadows+2020
use softened point masses, eps = 13 pc); nucleus density at 10 pc is 155 not 2500 Msun/pc^3 ->
DM particles near the nucleus <~ 1e4 Msun; a translating frame does not supply friction ->
prescribed UniformAcceleration a_df(t) in the inertial frame; Jacobi radius at 0.4 kpc is 33-59
pc, not 50-90; two class-2 histories end today above the nucleus mass (6.7e7, 1.0e7) -- model
change, proposed not applied. Disagreed: "live nucleus wrong by construction" is overstated
but the conclusion holds (t_relax(2e4) ~ 0.25 Gyr even with ln Lambda reduced 4x).

## 2026-09-20 -- tidal tracks (Penarrubia+2010, Errani & Navarro 2021) walked along the nine orbits

User wanted to explore the tidal-track idea. `src/ocen_dm/tails/tidal_tracks.py`: EN21 model
(track eq. 5, remnant profile eqs 7-9, time evolution eqs 4, 10-16) with a passage-by-passage
walker using the local host mean density inside each pericentre; tests
`tests/test_tidal_tracks.py` (3). Products: `results/tails/tidal_tracks_budget.ecsv`,
`plots/tidal_tracks_budget.png`; write-up Sec. 8 of `docs/progenitor_orbits.tex`.

Finding 1: our orbits are far outside EN21's calibration (T_mx0/T_peri = 8-65 vs 0.2-2), so
the walk extrapolates to meaningless remnants (1e-12); but the density criterion is robust: a
self-bound remnant needs ~20 x host mean density inside the pericentre = 1.7e10 (classes 1-2)
to 1.5e11 Msun/kpc^3 (class 3); an NFW of 1e10 Msun, c=5.5 has that only inside 4.6 pc
(6.7e3 Msun) / 0.5 pc (1e2 Msun); c=12: 14 pc (1.7e5) / 1.6 pc (2.4e3). A bare DM halo is
stripped below the nucleus mass on all nine orbits.
Finding 2: any DM around omega Cen today is held by the nucleus, so the budget is the initial
DM inside the nucleus's Jacobi radius at the smallest pericentre: class 3 (r_J 33-36 pc)
2-4e5 (1e10, c=5.5), up to 1.0-1.2e6 (c=12); classes 1-2 (r_J 40-70 pc) 0.5-1.5e6, up to
4.4e6. Upper envelopes: no adiabatic contraction (raises), no pericentric shocks inside r_J
(lowers) -- the N-body question.
Finding 3: this shrinks the N-body problem to the inner ~100 pc (2e6 Msun of DM in a 1e10
halo -> 1e5-1e6 particles of ~10 Msun), with the outer halo analytic and stripped per the
tracks; the multi-mass worries largely disappear.

## 2026-09-20 -- predicted DM density inside omega Cen if the DM within r_J is retained

User asked for plots of the predicted DM density for different r_J, initial masses and
concentrations at 2-3 radii inside the cluster including the deep WD field (20 pc, Scalco+2024).
`src/ocen_dm/tails/dm_density_plots.py`: model rho(r) = rho_NFW(r; M200, c) exp(-r/r_J)
(EN21-style truncation at the nucleus's Jacobi radius), grid M200 = 1e9/1e10/1e11, c = 5/10/15
(z = 2), r_J 10-200 pc; stellar density from the MGE light model scaled to the rung-0 K1
M_star = 2.88e6. Figures `plots/dm_density_in_rj_profiles.png` (profiles vs r, three M200
panels, r_J = 35 and 70 pc) and `plots/dm_density_in_rj_vs_rj.png` (rho at 3, 10, 20 pc vs r_J);
table `results/tails/dm_density_grid.ecsv`.

Reading: inside r_J the density is the initial cusp, rho ~ rho_s r_s / r, so r_J matters only
through the truncation factor (<= 2x between class 3 and classes 1-2 at 20 pc); the leverage is
M200^(1/3) and c (factor ~5 each across the grid). At 20 pc the grid spans 0.6-17 Msun/pc^3
(23-640 GeV/cm^3) against 9 Msun/pc^3 in stars: DM comparable to or above the stars there for
M200 >= 1e10 with c >= 10, or any 1e11. At 3 pc stars dominate by 10-200x. The 0.1-3
Msun/pc^3 band assumed in docs/dm_capture_constraints.tex covers only the low-mass/low-c
corner of the grid. Caveats: no adiabatic contraction (raises), no shocks inside r_J (lowers).

## 2026-09-20 -- dm_capture_constraints.tex extended

Added (user request, with the links/explanations from their Google Doc): a Motivation section
(McCullough & Fairbairn 2010; Hooper+2010; the LZ 248 keV event arXiv:2609.02823; inelastic
interpretation McCabe 2609.04181; Higgsino exclusion by solar neutrinos 2609.02775, 2609.07807;
why WDs reach mass splittings up to ~5 MeV; the four open items) and a section "Predicted DM
density inside omega Cen from the progenitor-orbit analysis" with the frozen-cusp-inside-r_J
figure (plots/dm_density_in_rj_profiles.png), the density table at 3/10/20 pc, and the caveats
(shocks inside r_J, contraction/cores, two-body heating). PDF now 15 pages.

## 2026-09-20 -- equilibrium of tidally truncated DM around the nucleus (new section)

User: "would not the DM have to relax into the new state inside r_J?" -> yes: tides truncate in
energy, not radius. New section in docs/dm_capture_constraints.tex with the analytic estimate,
references and a test programme; code `src/ocen_dm/tails/truncated_equilibrium.py` (tests 2).
Analytic result: a rho ~ r^-gamma tracer in the nucleus's Kepler potential has f(E) ~ (-E)^(gamma-3/2);
truncating at Phi(r_J) leaves the fraction F_gamma(x) = 1 - I_x(gamma-1/2, 3/2), x = r/r_J.
For gamma = 1: F(20 pc) = 0.35 (r_J 70) and 0.14 (r_J 35); gamma = 3/2: 0.60 / 0.28 -- much
more than the e^-x stand-in (0.75 / 0.56) used in the profile figure. So the frozen-cusp grid at
20 pc (0.6-17 Msun/pc^3) becomes 0.2-6 (classes 1-2) and 0.1-2.4 (class 3) before shocks.
Shock heating per passage (impulsive + GHO99 adiabatic correction (1+(omega tau)^2)^-1.5):
classes 1-2 (r_p 1.57 kpc, v_p 387): 20 pc protected (dE/E 1e-5 per passage), 50-70 pc eroded
in 30-170 passages; class 3 (r_p 0.42, v_p 507): 35 pc eroded in ~44 passages, 20 pc reaches
unity in ~600 (vs 66 passages available) -> class-3 remnant set by shocks. Experiments E1-E6
listed (AGAMA Eddington truncation; isolated relaxation; static tide at the two pericentres;
class-1 and class-3 orbits with shocks; contracted/cored cusps; stellar heating estimate).
Figure plots/truncated_equilibrium.png. PDF 20 pages.

## 2026-09-20/21 -- rung 0: K1 and K2-cored finished, K2-NFW still sampling

K1 (isotropic, no scales, 5 params): ln Z = -190.02 +- 0.48, chi2_ml = 771.5/89 (4.0 h).
K2-cored (7 params): ln Z = -188.94 +- 0.36, chi2_ml = 761.1/89 (10.0 h). Delta ln Z = +1.1 +- 0.6
in favour of the cored halo: no evidence either way. M_DM(<100 pc) = 1.3e6 [0.19, 2.5] x1e6,
95% upper limit 3.1e6; r_s 261 [46, 663] pc (unconstrained); the DM trades against M_rem
(3.1e5 -> 7.7e5) and M_star (2.88e6 -> 2.3e6). Distance 5.31 +- 0.02 in both (prior 5.43 +- 0.05).
Per-dataset chi2 (K1 / K2): HST radial 160/152 (21 pts), HST tangential 477/484 (21), MUSE 77/76
(29), Gaia radial 19/25 (9), Gaia tangential 38/24 (9). The isotropic model fails on the HST
tangential profile beyond ~100 arcsec (data below model by up to 10 sigma at 150-300 arcsec):
the outer HST field is radially anisotropic (sigma_T < sigma_R) or carries a tangential
systematic -- this is exactly what rung 1 (constant beta) is meant to absorb. JamPy cross-check
of the best samples agrees with our solver to ~1 in ln L (chi2 shifts of ~18 between HST radial
and tangential, the two engines' known difference).
K2-NFW: 11.8 h, 7.0e6 calls, iteration 7641, remainder fraction still 82% -- ultranest is
crawling along the r_s-M_DM ridge (1 accepted draw per ~4000 calls). Not stopped; decision
pending (let it run vs restart with the step sampler).
Report: results/fits/comparison.md; figures plots/fit_posterior_profiles.png,
plots/fit_rung0_K1_posterior_profiles.png, plots/fit_rung0_K2_cored_posterior_profiles.png.

Two tooling fixes: report._family_for now honours the ladder switches (isotropic / constant
beta / no scales) from run.yaml -- it crashed with KeyError beta_0 on rung-0 runs; run.yaml
input hashes are now taken at launch, not at the end (the Gaia EDR3 profile was rewritten at
17:22 during the test-suite run with identical numbers, which made the comparison table flag
the two runs as "different observations" although they fitted the same data). Which code path
rewrote the product is not identified (no test calls build_edr3_profile directly).

## 2026-09-21 -- rung 0 complete (all three), and Section 11 experiment N1 done

K2-NFW finished after 20.6 h / 1.43e7 calls: ln Z = -190.09 +- 0.35, chi2_ml 763.1/89,
M_DM(<100 pc) = 5.3e5 [0.23, 11.9] x1e5 (95% upper limit 1.64e6), r_s 123 [7, 544] pc.

Three-way (89 points, isotropic, no instrument scales):
  K1 no DM     ln Z = -190.02 +- 0.48   max ln L -169.7   chi2 772
  K2 cored     ln Z = -188.94 +- 0.36   max ln L -164.4   chi2 761
  K2 NFW       ln Z = -190.09 +- 0.35   max ln L -165.5   chi2 763
Delta ln Z: cored - K1 = +1.08 +- 0.60; NFW - K1 = -0.07 +- 0.59; cored - NFW = +1.15 +- 0.50.
No evidence for dark matter at rung 0; the NFW form is indistinguishable from no halo. 95%
upper limits on M_DM(<100 pc): 3.1e6 (cored), 1.6e6 (NFW).

Comparability of the evidences verified directly: recomputing K1's chi2 at its stored ML point
with today's files reproduces all five per-dataset values to 5 significant figures, so the
Gaia-product hash change of 20 Sep 17:22 was cosmetic (metadata) and the three runs fitted
identical data. The report's warning is therefore a false positive on these runs.

All three fits are poor in the same place: HST tangential chi2 477-484 for 21 points, data
below model by up to 10 sigma beyond ~100 arcsec, while HST radial is 150-160 and MUSE 76-77.
The DM cannot fix it (it acts outward, the misfit is at 150-300 arcsec where it would need
sigma_T < sigma_R). This is the anisotropy signal rung 1 (constant beta) is meant to absorb;
until then the DM limits are conditional on isotropy.

Section 11 experiment N1 (verification of F_gamma(x) = 1 - I_x(gamma-1/2, 3/2)): passed.
Quadrature of the physical Eddington integral matches the closed form to 1e-15 for gamma = 1
and 1.5; Monte Carlo with a numerical inverse-CDF agrees within 1.1 sigma at x = 0.1, 0.3, 0.6
(now a test). Two sampler bugs found and fixed on the way -- a rejection envelope that is not an
upper bound when gamma < 1.5 (biased the retained fraction 0.60 -> 0.87), and trapezoid CDF
bias at the integrable endpoint singularity (0.6%); both would have poisoned N2-N4, which reuse
this sampling machinery.

## 2026-09-21 -- where the rung-0 chi2 comes from

Diagnosis at the K2-cored ML point (per-point residuals, chi2 split by radius):

  dataset                 R<50"        50-150"        >150"
  HST radial          12.9 (11 pts)  49.5 (5)     90.1 (5)
  HST tangential      13.6 (11)     165.8 (5)    304.4 (5)
  MUSE LOS             9.5 (13)      27.5 (10)    38.9 (6)

Inside 50 arcsec the fit is excellent for all three datasets (1.2, 1.2, 0.7 per point). Every
bit of the bad chi2 is beyond 50 arcsec, and two thirds of the total comes from the HST
tangential component alone.

The *physical* size of the misfit is small: HST radial data are +1.2% above the model at
138-224", HST tangential -3.7 to -10% below at 108-311", MUSE -9.4% at 290". These become 6-11
sigma only because the HST statistical errors there are 0.16-0.35%. A chi2 of 761 corresponds
to a model wrong by a few per cent.

The signature is a *ratio* error, not a mass error. sigma_T/sigma_R in the data falls from
0.99 at 41" to 0.85 at 311"; the model (which already subtracts the published rotation from the
tangential component) only falls to 0.92. The geometric mean of the two components is matched
to 2-3%, i.e. the enclosed mass is fine; what is missing is either radial anisotropy
(beta ~ 0.15-0.25 at 4-8 pc) or more tangential streaming (an extra 5.8 km/s in quadrature at
223", on top of the 4.45 km/s the published curve already supplies -- i.e. a 64% larger rotation
amplitude). The MUSE LOS deficit at the same radius favours anisotropy, since sigma_LOS at large
projected R samples the tangential direction and is measured by an independent instrument.

Error inflation alone cannot rescue it: adding a systematic floor in quadrature gives
chi2(89) = 761 (raw), 428 (0.5%), 271 (1%), 144 (2%), and HST tangential remains the largest
contributor at every floor. The shape is wrong, not just the error bars.

Consequence: all three rung-0 models are misspecified in the same way, the evidence comparison
is between three models that lack a needed degree of freedom, and the DM limits (3.1e6 cored,
1.6e6 NFW within 100 pc) are conditional on isotropy. Rung 1 (constant beta) is required before
any DM statement; the rotation alternative is testable with the --gaia-rotation ours variant and
with the equivalent HST switch.

## 2026-09-21 -- Section 11 test programme written up

`docs/SECTION11_TESTS.md`: ranks the uncertainties in Section 11 (r_J definition 1.8x > gamma_ad
factor 2 > non-Kepler corrections > initial anisotropy > the algebra, which is now verified),
then specifies N0-N6 with set-up, what each verifies, pass criteria and cost. N1 (the closed
form) is done and passed. Minimum useful set N0+N2+N3 = one working day; N4 delivers
rho_DM(20 pc) per orbital class for the WD comparison.

## 2026-09-21 -- rung-0 analysis written up

`docs/RUNG0_ANALYSIS.md`. Beyond the chi2 diagnosis already journalled, one new finding: the
sigma_T/sigma_R ratio, corrected for the rotation the model already subtracts, is -8% at
5.8-8 pc and +12 to +18% at 34-55 pc. Within Gaia alone the data ratio rises monotonically
0.86 (11 pc) -> 1.18 (55 pc), crossing unity at ~27 pc. That is the classical anisotropy
profile of a tidally limited cluster: isotropic core, radial at 4-20 pc (beta ~ +0.25),
tangential beyond ~30 pc (beta ~ -0.2 to -0.4). A single rotation-curve normalisation error
cannot produce the sign reversal, so anisotropy is the natural reading; the rotation variants
remain the cheap competing test.

Consequence: a constant beta (rung 1) cannot fit a sign reversal -- the deferred beta(r)
turnover parameter is now required, and rung 2 is promoted from optional to necessary.

Second new finding: the DM "evidence" decomposes as Delta chi2 = -10 total, of which Gaia
tangential -14 and Gaia radial +6, both at 34-55 pc where the data demand beta < 0. The cored
halo is acting as a proxy for tangential anisotropy; the +1.1 in ln Z is not a statement about
dark matter. Rung-0 M_DM limits must not be quoted.

Also recorded: kinematic distance 5.311-5.313 +- 0.021 kpc in all three runs, 2.4 prior sigma
below the 5.43 +- 0.05 prior and 2.4x tighter -- to be re-checked at rung 1 since anisotropy
changes the PM-to-LOS ratio; M_BH = 4.3-4.7e4 +- 8% in all three; total central mass stable at
3.1-3.2e6 while the star/remnant split moves by 6e5 when a halo is added; and a note that the
NFW family needs the ultranest step sampler at rungs 1-2 (20.6 h and 1.4e7 calls at rung 0,
most of it on the M_DM-r_s ridge).

## 2026-09-21 -- code/analysis audit and documentation corrections

Read the data-to-likelihood pipeline, model families, fit/report machinery and the
orbital/DM-survival calculations against their saved products. The findings and
reproduction recipes are in `docs/CODE_ANALYSIS_AUDIT.md`. Updated the README,
modelling and validation plans, current analyses, source docstrings/comments, and
the four companion LaTeX/PDF notes. Earlier entries remain a chronological record;
the qualifications below supersede their stronger interpretations.

The maximum log likelihoods of all three 89-point rung-0 fits reproduce exactly
with the current default inputs. From the unrounded stored evidences, Delta ln Z
is +1.08 +/- 0.60 (cored minus K1), -0.08 +/- 0.59 (NFW minus K1), and
+1.16 +/- 0.50 (cored minus NFW). All remain poor fits. Reproduction supports the
metadata-only explanation of the Gaia hash difference, but cannot prove that all
original inputs were identical: complete input snapshots were not saved.

The PM ratios motivate anisotropy tests; they do not recover intrinsic beta(r),
establish enclosed mass to 2-3%, or rule out a refitted error model. The floor
experiment used a fixed parameter vector. A central moment fit alone does not
establish a robust IMBH detection. The positive-DF condition for a cored tracer
in a Kepler potential applies to constant beta as well as varying beta. The
current varying-beta law is monotonic; the proposed turnover is not implemented.
There is a Gaia rotation switch but no equivalent HST switch in the present CLI.

Corrected two material survival calculations: at 20 pc, the local class-3 shock
estimate is 66 x 0.001701 = 0.112, not order unity; and the energy-cut retention
factor multiplies the initial density, replacing the exponential taper rather
than multiplying it again. The circular spherical tidal coefficient is
3 - d ln M_host/d ln R. The old 59-pc alternative has not been reproduced.
The N2 baryonic models are alternatives, avoiding the previous double counting
of the cluster mass. Initial aperture budgets and extrapolated tidal tracks are
not validated remnant masses or universal survival bounds.

Confirmed software issues are documented for a separate implementation pass:
report replay ignores saved Gaia variants; mock comparison reads obsolete
metadata keys; mock generation does not inherit the generating model's rung;
and provenance/family reconstruction does not capture all resolved inputs and
settings. A bare fit also retains the legacy anisotropy/scale defaults, so the
adopted recipes now show explicit flags and unique run labels.

Validation: 104 focused tests passed; no full-suite or long-grid rerun. All four
PDFs rebuilt and were visually inspected. The eleven changed Python files have
identical executable syntax trees after removing docstrings, and all 39
snapshotted processed-data and rung-0 result files retain their original hashes.
No numerical behaviour, observations, posterior samples or orbit products changed.

## 2026-09-21 -- reporting, mock handling and run preservation repaired

The user approved this implementation pass and clarified the interpretation of rung 0:
it was designed to establish a baseline and guide subsequent modelling. The documentation
now leads with that purpose. Scientific extensions remain the planned research sequence.

New runs snapshot the actual binned likelihood arrays and resolved model/MGE before
sampling, with parameter order, priors, fixed values, distance, halo taper, engine,
sampler settings, source hashes, git state and available package versions. Reports
replay those inputs, verify checksums, and check the stored maximum likelihood before
plotting. Legacy runs retain their original files and use recorded options with a
warning; an inconsistent reconstruction raises an error.

Mock generation from a saved run now uses that source's model independently of the
fitted family/rung. Bare vectors require an explicit generating family. Mock evidence
comparisons use the actual metadata schema; new runs compare saved observation arrays.
The existing Gaussian-average-error/clipped mock generator is unchanged and named in
provenance. The CLI and driver refuse existing run directories, and interrupted runs
keep their launch snapshots.

Validation: 78 focused tests passed, including 36 new regression cases; all 40 replay/CLI
tests passed again after strengthening two checks. The three archived rung-0 maximum
likelihoods reproduce exactly and all 39 snapshotted data/result files are unchanged.
The data and mass-modelling PDFs were rebuilt and inspected. No production fit was run.
Details, recipes and remaining historical limitations: `docs/CODE_ANALYSIS_AUDIT.md`
and `docs/MODELLING_PLAN.md`.

## 2026-09-21 -- rung-1 production fits launched

User authorised the rung-1 fits while discussion of the two independent dynamical
programmes continues: survival of an existing DM remnant, and disruption of a full
nucleated progenitor. No dynamical simulation was launched in this step.

At 10:46:47 UTC, started `rung1_K1`, `rung1_K2_cored` and `rung1_K2_nfw` in fresh
directories under `results/fits/`, with PIDs 14568, 14570 and 14572 respectively.
All three have one free constant anisotropy parameter, beta uniform on [-1, 0.5],
no instrument scales, the composite tracer and Jeans backend, and the same 89
measurements and remaining priors as rung 0. Gaia uses raw errors and published
rotation. Each run uses 400 live points, dlogz = 0.5, seed 42 and no call limit;
NFW uses the slice step sampler. Each process is restricted to one numerical
library thread and has a job-lifetime idle-sleep inhibitor.

Preflight: beta = 0 reproduces each corresponding rung-0 maximum log likelihood
exactly. Launch metadata confirm 6/8/8 free parameters and identical saved data
fingerprints: `37dfc7fa807b7f2f9862433d9bca9bba1224b5a2871c748a6529baff12aaad99`.
All three processes were verified active, with UltraNest sampling the initial
400 live points. These are launch checks, not convergence or fit results.

Exact commands, environment overrides, PIDs, logs and a source archive are in
`results/fits/rung1_launch_20260921T104646Z/`; each fit also has its own resolved
model, data and source-hash snapshots in `run.yaml` and `data_snapshot.json`.
The archived rung-0 run directories were not reused.

## 2026-09-21 -- controlled dynamical experiment pipeline implemented

Added `ocen_dm.dynamics` and the `ocen nbody prepare/run/analyse/plot` commands for
the two independent programmes: survival of an existing compact DM remnant and
disruption of a nucleated progenitor containing stars and DM. Usage, assumptions
and reproducible commands are in `docs/DYNAMICAL_EXPERIMENTS.md`; four pilot
templates and a short software-check configuration are in `configs/dynamics/`.

Preparation builds spherical constant-beta DFs in the joint nucleus+DM+stellar
potential, independently checks the signed inversion and verifies recovered
densities and anisotropy. The Plummer nucleus is counted once. Optional importance
sampling by orbital pericentre provides variable particle masses with inverse
selection weights. It protects low-pericentre orbits rather than selecting by
instantaneous radius; effective counts and later heavy-particle contamination
are recorded. The templates still require particle-count/softening/time-step
convergence before scientific use.

The compact core/cusp pair is isotropic and matched in M_DM(<20 pc). The broad
progenitor core fails the signed DF check when isotropic in this particular
combined potential. Its tangential beta_DM=-0.3 variant passes; both progenitor
templates use beta_DM=-0.3 to keep their anisotropy matched, with isotropic stars.
They share total masses and scales, so their central DM masses differ.

Live evolution uses gyrfalcON particle self-gravity and the analytic rigid nucleus;
frozen controls retain the complete initial satellite potential. Isolation,
McMillan17 class-1 and slowing-bar class-3 tracks preserve all Cartesian coordinates
and the host's absolute clock. Force comparisons validate the saved host assets.
Both modes use prescribed nucleus trajectories without dynamical friction or
back-reaction on the nucleus. Those are explicit future extensions for historical
reconstruction. No long production dynamical experiment was launched.

Validation exposed and corrected three interface traps: AGAMA's generic export
does not retain analytic Plummer parameters; gyrfalcON's output-key option alone
does not read particle keys; NEMO's standalone time header rounds to six significant
digits. Analytic potentials are now serialized explicitly, a run-local no-op
manipulator requests input keys, and full-precision time columns are parsed.

Each preparation and evolution has fresh-directory protection, saved configurations,
input checksums, source archives, code metadata and status records. Live runs also
record executable hashes and commands. Diagnostics include component mass and
shell-density profiles, effective counts, anisotropy, mass-weighted energy changes,
a spherical iterative bound-mass proxy and an explicitly defined tidal-tensor scale.
All particles remain in the simulation.

Validation: 124 focused tests passed in 45.70 seconds, including analytic Plummer
inversion, circular point-mass and flat-rotation tidal scales, phase-space refinement,
both progenitor templates, frozen energy conservation, and direct NEMO tests of
external gravity, binary self-gravity and arbitrary particle IDs. A subsequent
manifest/config-consistency guard passed its focused regression. Ten short
end-to-end evolutions (five preparations, live and frozen) completed and were replayed:
isolation, class 1, class 3 starting at 2000 Myr in the bar clock, and cored/cusped
progenitors on class-1 orbits. Each spans about 0.1 Myr with 512 remnant particles or
1000 particles per progenitor component. These validate software and integration,
not long-term survival or numerical convergence. Products and the inspected comparison
figure are in `results/dynamics/validation_20260921_v2/`.

The three rung-1 fit processes were verified active while this work completed.

## 2026-09-21 -- first isolation pilots launched

Prepared all four dynamical templates at their configured particle counts and durations.
Every model passed the signed-DF, density/beta recovery, saved-potential and input-checksum
checks. Compact remnants have 100,000 DM particles and run for 20 Myr; progenitors have
200,000 DM plus 100,000 stellar particles and run for 100 Myr. Each has a live evolution
and a frozen control from the same saved particles. These runs test initial equilibrium
and numerical drift before tidal passages; they do not establish resolution convergence.

Two independent queues started at 11:32:26 UTC under
`results/dynamics/isolation_20260921T113007Z/`: remnant worker PID 17802 and progenitor
worker PID 17804. Each queue runs cusp live/frozen, then core live/frozen, and writes
comparison figures. At most two evolutions run simultaneously, each using one numerical
library thread and nice level 5. Job-lifetime sleep inhibitors are attached to the workers.
The three rung-1 fits remain active at their existing priority.

The full preparations expose a resolution limitation: the cored progenitor has only
16 DM particles within 20 pc (the cusped progenitor has 539). The cored model's pilot
can assess broader equilibrium behaviour, but its central density needs finer sampling.
This limitation and all aperture effective counts are saved in the launch record.

Both first live evolutions started gyrfalcON successfully, with the expected particle
counts and external nucleus. Actual maximum/minimum steps are 0.0149199253 and
0.00186499066 Myr. Commands, configurations, input/source hashes, logs, process IDs and
queue statuses are preserved in the batch directory. A first startup attempt was blocked
by the sandbox's process-priority restriction before any evolution began; its records
are retained, and the authorised retry runs outside that sandbox. No tidal run was queued.

## 2026-09-21 -- remnant single-passage pilots and numerical controls launched

The user approved advancing remnant survival independently of progenitor isolation.
At 11:53:05 UTC, launched three low-priority, single-threaded queues in
`results/dynamics/remnant_passage_20260921T114939Z/`: cusp passage (worker 18978),
core passage (18980), and numerical controls (18982). Both tidal live runs and
the first timestep control were verified advancing in gyrfalcON. The existing
progenitor queue and all three rung-1 fits remained active.

The class-1 McMillan17 track starts at the nearest past apocentre of the nominal
present orbit, refined by zero radial velocity. It spans one pericentre at
1.99817 kpc, reached after 45.11369 Myr, and ends at the next apocentre after
88.09158 Myr. Binary-step alignment makes the actual integration 88.10216 Myr.
The initial/final apocentre radii are 6.95887/6.91760 kpc. Both profiles have
byte-identical orbit and host assets; maximum sampled interpolation acceleration
error is 4.01e-6. The selected phase-space state and orbit search are preserved.

Each tidal model has 100,000 DM particles, 0.2 pc softening, and matched live/frozen
evolutions. Both also have live/frozen isolation references spanning the same full
88.10216 Myr, to distinguish tides from longer numerical drift. The initial models
are spherical isolated equilibria placed at apocentre, without prior tidal relaxation
or adiabatic growth of the host field. The nucleus remains rigid and prescribed.

For each core/cusp model, three 20.00762 Myr isolation controls change one numerical
setting: N=200,000 at fixed softening/timestep; softening=0.1 pc at fixed N/timestep;
or half the maximum and minimum timesteps at fixed N/softening. The doubled-N models
have their own frozen controls. Softening and timestep controls reproduce the original
particle positions, velocities, masses and IDs exactly and reuse the original frozen
reference. All 20 Myr controls share the completed baseline's exact endpoint.

Ten preparations passed DF, density/anisotropy, input-checksum, particle-matching and
orbit checks. The bounded batch contains 16 evolutions plus comparison figures, with
at most three new evolutions at once. Commands, operational scripts, source hashes,
resolved configurations, preflight records, process IDs and logs are preserved in the
batch directory. These are sensitivity checks and first-passage pilots; convergence
of tidal results still needs assessment from completed outputs and, if warranted,
further resolution tests.

## 2026-09-21 -- first tidal results, NFW rung-1 completion, and control failure

Both remnant class-1 single-passage live/frozen pairs completed 88.10216 Myr.
Diagnostics checksums were verified and the comparison figures inspected. Relative
to their initial samples, live M_DM(<20 pc) changes by -1.205% (cusp) and -3.415%
(core); the corresponding frozen changes are -1.343% and -2.521%. The live spherical
bound-mass proxies decrease by 6.876% and 8.579%, respectively. These are preliminary
single-passage results: the matched-duration isolation references are still running,
and numerical sensitivity checks are incomplete. No escape fraction or converged
survival limit is inferred from the spherical binding proxy.

Both half-timestep isolation runs completed. At the common 20.00762 Myr endpoint,
their M_DM(<20 pc) differs from the original live run by about +0.116% (cusp) and
-0.062% (core). Maximum logged fractional energy drift fell to 6.17e-5 and 2.66e-5,
compared with 2.93e-4 and 1.68e-4 in the original runs.

The cusp 0.1-pc-softening control aborted in gyrfalcON with a NaN position for body
96346, stopping the numerical-control queue. An unchanged rerun in a fresh
`evolution_live_repeat` directory reproduced the same error. The final saved snapshot
at 1.50691 Myr is finite; body 96346 belongs to its closest particle pair, separated
by 0.01112 pc. This observation does not yet identify the numerical cause. Both
failed runs and all logs remain preserved; this softening setting is not validated.

The already-approved, unstarted doubled-N controls (both profiles, live/frozen),
followed by the cored softening control, were moved to a separate bounded continuation
queue in `results/dynamics/remnant_controls_continue_20260921T122459Z/`.
Worker 29353 was verified running the cusp doubled-N case. Original queue records
remain unchanged. The failure does not prevent these independent controls from running.

The cusped progenitor completed 100.00826 Myr in live and frozen isolation; its bound
DM and stellar mass proxies remain constant. Central shell estimates are noisy at
roughly 200 effective DM particles. The cored progenitor remains in live isolation.

The rung-1 NFW fit completed and its saved maximum likelihood replay passed. It has
beta=0.10846 [0.10391, 0.11316] (median and 16th/84th percentiles), chi2=316.01 for
89 measurements, and lnL_max=58.1051411. The corresponding rung-0 NFW chi2 was 763.12.
Stored ln Z is 30.53163 +/- 0.40937. These are conditional results for the current
constant-beta NFW model. K1 and cored rung-1 fits remain active but sample inefficiently;
the three-family comparison is not yet available.

## 2026-09-21 -- P1 overflow repaired and slow fits replaced with slice sampling

The user authorised a bounded repair of the NaN failure and sampler inefficiency.
Instrumenting the installed scalar P1 Taylor kernel identified the cause: at order 3,
X=5.43901e7 and the preceding raw derivative=1.3316e30 produce an intermediate
7*X*D3 about 5.07e38, exceeding float's range. Multiplication by eps^2/2 would bring
the correction back into range, but the overflow had already occurred. The resulting
NaN is numerical; neither the physical softening nor the timestep needed changing.

The new coefficient helper evaluates the same P1 derivatives using a dimensionless
correction, without forming the unnecessary next derivative. New live runs build a
run-local falcON library that replaces only kernel.o, checks the upstream source
block, and records compiler, source/object and library hashes. Unsupported source
versions fail explicitly. Missing or altered runtime libraries also fail instead
of silently falling back. The installed NEMO source and library remain unchanged,
verified by hashes. The helper and builder are included in source archives.

A 64-particle close-cell reproducer aborts with the original library and completes
with the repaired library. Its maximum force error relative to direct P1 summation
is 7.25e-5. The focused regression suite passed 95 tests, including real NEMO external
gravity, binary self-gravity and the new force comparison. After strengthening the
runtime integrity check, all five kernel/checksum tests passed.

The original 100,000-particle, 0.1-pc-softening cusp configuration completed the full
20.00762 Myr in `evolution_live_p1_fixed`, with 41 analysed snapshots and validated
IDs, masses, endpoint and diagnostic checksum. M_DM(<20 pc) changes from 100128.84
to 100435.89 Msun; the bound-mass proxy remains constant. Maximum logged fractional
energy drift is 1.57e-4. At the final valid pre-crash snapshot, the repaired and old
runs differ by at most 7.27e-5 pc in position and 2.47e-4 km/s in velocity. Both failed
attempts remain preserved. The comparison figure was inspected. The doubled-N,
cored-softening and full-duration remnant isolation controls have also completed.

K1 and cored rung-1 replacements were launched as `rung1_K1_slice` (PID 30585) and
`rung1_K2_cored_slice` (30587), using the existing SliceSampler with 2*ndim steps.
Their exact saved observations, model, priors, seed and convergence settings match
the old runs; beta=0 again reproduces the corresponding rung-0 likelihoods. After
verifying replacement startup, the old jobs were stopped with SIGINT. Their readable
HDF5 checkpoints were copied, hashed and retained, and their manifests now record
the interruption and replacement labels. No existing project run was resumed.

Early comparisons over matched likelihood ranges show about three times fewer
likelihood evaluations per accepted iteration: 48.9 versus 145.7 for K1 and 53.1
versus 163.6 for the cored model. These are stage-specific efficiency measurements,
not completed fits or predicted total speedups. Both replacement fits remain active.
Evidence, scripts, checkpoints and repair verification are under
`results/repairs/numerics_20260921T125417Z/`. No new scientific grid was launched.

## 2026-09-21 -- completed dynamical batch analysed

Analysed all 24 completed evolutions from 14 preparations, including the repaired
cusp softening control; excluded and preserved the three failed attempts. Verified
prepared-input and diagnostic hashes, manifest links, durations, mass invariance,
nested apertures and matched initial conditions. Independently reread 12 final
snapshots, checked IDs/masses/epochs and reproduced the saved aperture masses.
Extended the progenitor endpoint profiles to 10 kpc. The script is
`bin/analyse_dynamical_batch.py`; results, source copy, six PNG/PDF figures and
checksums are in `results/dynamics/batch_analysis_20260921_v2/`. The earlier analysis
directory is retained. The full interpretation is in
`docs/DYNAMICAL_BATCH_ANALYSIS.md`.

At 88.10216 Myr, live M_DM(<20 pc) is 1.226% lower than its isolation control for
the cusp and 3.386% lower for the core. Final-10-Myr median deficits are 1.406%
and 3.309%; frozen medians are 1.430% and 3.028%. These temporal summaries are not
confidence intervals. Live M_DM(<70 pc) endpoint deficits are 5.351% and 7.120%,
and spherical bound-mass proxy declines are 6.876% and 8.579%. Final-10-Myr shell
density depletion is about 17% at 70 pc; beta there is -0.298 and -0.285, versus
near zero in isolation. The comparison fixes M_DM(<20 pc), not total DM mass.

The saved orbit crosses the disc at 24.42066 and 53.27006 Myr, at cylindrical
radii 5.09361 and 3.00873 kpc. Pericentre occurs at 45.11369 Myr. The largest
compressive host-force-gradient pulse accompanies the second crossing, followed
by the strongest outer depletion. The frozen 70-pc cohort gains much more
reference energy than the 20-pc cohort. This suggests a disc-shocking contribution;
the batch does not isolate it from the preceding pericentre response.

Halving timestep or softening changes the 20-Myr remnant endpoint M20 by at most
0.116%. Doubled-N live/frozen differences using their own matched samples are
+0.011% and -0.155%. Isolation energy drift decreases with half timesteps; the
88-Myr baseline runs reach 0.135% and 0.0728%. Recommend extending timestep and
particle-number controls to the full tidal interval before repeated passages.

Both progenitor live/frozen pairs completed 100.00826 Myr. Bound component masses
remain constant; no particles above ten times the minimum component mass enter
the 100-pc protection aperture at saved epochs. At sampled radii from 100 pc to
10 kpc, live/frozen endpoint masses differ by less than 1% in both components.
The core has only 16 initial DM particles inside 20 pc, ending with 18 live and
23 frozen; its large central fluctuations are unresolved sampling behaviour.
The existing models support an extended-body disruption pilot, while central
remnant measurements need improved sampling and its isolation validation. These
programmes remain independent. No new evolution or fit was launched.

## 2026-09-21 -- LaTeX report of the first dynamical batch

Created `docs/dynamical_batch.tex` and `docs/dynamical_batch.pdf`: an
11-page write-up with six vector figures, six tables, the initial-condition
and diagnostic equations, complete batch specifications, numerical controls,
force-kernel repair, interpretation and next experiments. Regenerated the six
plots at page width using `bin/plot_dynamical_writeup.py`; their input and
figure checksums are in `docs/dynamical_batch_figures/provenance.json`.
The timestep setting `fea=0.2` was checked against the installed gyrfalcON
manual and described as the acceleration-based timestep factor.

The PDF compiled without warnings, unresolved references or overfull/underfull
boxes. Rendered pages were inspected for table, equation and figure layout.
The portable archive `output/source/dynamical_batch_latex.zip` contains the
source, all six figures, provenance, checksums and build instructions. Its
contents were verified and compiled independently. No simulation products or
physical settings changed.

## 2026-09-21 -- Add the 20-pc density summary to the write-up

Added Section 8 and Figure 7 on page 11 of `docs/dynamical_batch.pdf`, with
all 24 completed evolutions and their distinct comparison references. The
text separates local shell density from enclosed mass, gives the tidal
medians, and discusses positive isolation offsets and the unresolved core.
`bin/plot_density_changes_20pc_paper.py` lays out the audited summary values
at the report's text width. The updated report has 12 pages and seven figures.
The PDF remains beside its LaTeX source. The portable LaTeX archive was updated
and compiled independently; its PDF text matches the report. No simulations
or fit settings changed.

## 2026-09-21 -- Lifetime exposure and the surviving-cluster constraint

The user requested multi-Gyr density evolution and stressed that Omega Cen
must survive to the present. Recorded `docs/LIFETIME_EVOLUTION.md`, linked
from the experiment guide and the earlier N4 plan. It separates stellar age,
MW accretion, envelope loss and present-orbit exposure; specifies density,
stellar-structure and joint-survival plots; and requires long numerical controls.
The existing fixed nucleus cannot test stellar survival. A responsive nucleus
and an assessment of physical relaxation are required for that inference.
The class-1 10-Gyr experiment is a controlled baseline; the implemented bar
history spans 7.82234 Gyr. Historical comparisons must end at the observed epoch.
Published exposure/migration references were checked against ADS and arXiv.
No long-duration run was launched and no lifetime density result was inferred
from the first-passage percentages.

## 2026-09-21 -- Responsive nuclei, staged lifetime queue, completed rung 1

Implemented a live Plummer nucleus sampled in the combined potential, with its
analytic gravity removed during self-gravitating evolution. Its signed DF and
density/moment recovery are checked like the DM. An isotropic Plummer nucleus
fails the cusp positivity check; β = −0.1 passes for both compact halos and is
the common controlled choice. Stellar diagnostics follow a shrinking-sphere
centre and include bound mass, three projected half-mass radii, 3D half-mass
radius, velocity dispersions and displacement from the reference orbit.
Initial/final kinematic residuals use archived observations and flag sparse
bins; they are diagnostic Cartesian projections, not a survival classification.

Added present-endpoint orbital rewinding, preserving the bar's absolute clock.
At accuracy 1e-12 the 10-Gyr static reference returns within 0.0165 pc and
0.000933 km/s; the 7.82234-Gyr bar track returns within 0.422 pc and 0.0161 km/s.
Their measured pericentre/disc-crossing counts are 114/209 and 67/135. The
local physical DM-heating estimate uses the initial DF and an explicit stellar
mass-spectrum bracket, without treating simulation particles as real stars or
turning the timescale into a depletion law. Orbital averaging remains needed.

Launched `results/dynamics/lifetime_20260921T162000Z` with controller PID 44098
and four workers. It has 12 full-passage controls, 16 live-nucleus 500.004-Myr
experiments, then eight conditional 7.822/10-Gyr frozen-potential comparisons.
Failed runs or numerical checks prevent stage advancement. The controller
executes preserved, hashed source; existing experiments are unchanged. Initial
verification found all four first-stage runs active with growing snapshots.
No new long-duration density or stellar-survival result is claimed yet.

Verified 54 focused tests in total: the existing dynamics/kernel checks and
12 lifetime cases, including actual NEMO integration in isolation and a host,
DF rejection, centring, units, observational bins, numerical-gate failure and
checksum rejection. The new compiled isolation test checks total momentum
rather than imposing an arbitrary tolerance on the noisy low-N density centre.
Plots and current run state are in `docs/LIFETIME_BATCH_REPORT.md`; refresh it
with `bin/report_lifetime_batch.py` and the batch path.

All three rung-1 models are complete and their maximum likelihoods replay
exactly on the same 89 data points. χ² is 342.51 (K1), 316.01 (NFW), 315.29
(cored), compared with 771.53, 763.12, 761.10 at rung 0. Cored versus NFW gives
Δln Z = 0.93 ± 0.54. Added `docs/RUNG1_COMPARISON.md`, profile figures and a
rung-0/rung-1 χ² comparison. The results support the intended next step in the
modelling ladder; constant anisotropy still leaves structured outer residuals.

## 2026-09-21 -- PNG plot convention and the Gaia anisotropy tradeoff

Moved analysis figures to `plots/` and adopted PNG as the sole plot format,
including the older diagnostic figures stored in result directories. Converted
the paper figures at 300 dpi and removed PDF plot exports after verifying their
PNG replacements. Numerical simulation inputs and fit products were preserved.
Updated generators, Markdown/LaTeX references, README and project AGENTS.md;
`results/maintenance/plot_location_migration.json` records old locations and checksums.
The 12-page dynamical document remains beside its LaTeX source. Rebuilt and
visually checked all pages, then verified the portable archive compiles with
its seven PNGs and reproduces the document text.

Investigated the rung-1 Gaia tangential deterioration by replaying the fits.
Cored χ² changes from 23.56 at rung 0 to 58.50 at rung 1. Holding rung-1 masses
and distance fixed and changing only β to zero gives χ² = 14.93 for Gaia
tangential, while HST radial/tangential worsen to 666.41/362.91. The equivalent
NFW check gives 17.78 instead of 62.23. The fixed-mass check is not a new fit.
Saved the bin-level diagnostics in
`results/fits/gaia_rung_anisotropy_check_20260921.json` and expanded the rung-1
comparison. The worst cored Gaia tangential residual is 4.31σ at 662 arcsec.

Confirmed that rung 2 introduces β(r), but the existing form is monotonic and
uses β(0) <= -0.5. It cannot combine a tangential centre, radial intermediate
region and an outer decline. Documented the distinction between a turnover
extension and a broad-central-prior monotonic diagnostic. No rung-2 fit was
launched while answering this modelling question.

## 2026-09-21 -- Analytical predictions versus the completed dynamical pilot

Added `bin/compare_dynamical_analytics.py` and three PNG figures under `plots/`.
The actual initial-DF energy cut predicts 20-pc density losses of 10.18%/14.71%
for cusp/core, compared with live losses of 2.50%/4.66%. It over-removes phase
space relative to the finite first passage. The scalar pericentre heating
formula, evaluated on the actual 1.998-kpc orbit, instead underpredicts the
full-orbit energy increment. Initial-radius cohorts contain extended internal
orbits; their contribution, and continued work on escaping particles, separates
the mean from the median. A large mean/prediction ratio is not a density-loss
fraction or a calibration of the adiabatic index.

The analysis checks saved-particle IDs, masses, epochs and input hashes,
reproduces the existing energy medians, subtracts matched isolation energy
drift, and verifies the independent signed-DF density integration with a finer
grid. Three tests pass for exact Kepler limits and finite-shell averaging.
Expanded `docs/dynamical_batch.tex` with equations, two tables and three figures;
updated the companion Markdown analysis and the interpretation of the proposed
N3 criterion. No simulation or rung-2 fit was launched for this analysis.

## 2026-09-21 -- Flexible rung-2 fits launched

Implemented the authorised two-transition stellar anisotropy profile with
three independent beta control levels and two ordered radial scales. Convex
weights prevent overshoot; the exact integrating factor supports the existing
spherical Jeans solver. The profile allows peaks, troughs and monotonic shapes.
The central prior retains beta_0 in [-1, -0.5] for the current cored MGE plus
point mass; this necessary condition is not a global DF-positivity certificate.
Legacy profiles and replay remain supported, and unsupported JAM/AGAMA or
constant-beta combinations fail explicitly. The full prior specification is
in `docs/MODELLING_PLAN.md`.

Verified 105 focused tests, including 27 new tests for the profile, numerical
Jeans integration, special-case projections, validation and saved-model replay.
Rechecked the 10 ladder tests after adopting the new rung-2 parameter counts.
The launcher independently replayed all three completed rung-1 likelihoods
and checked 258 finite prior points per rung-2 family. It preserves the same
89 measurements, MGE, mass and distance priors, rotation/error options, and
excludes instrument scales. The central-beta prior is explicitly different
from the diagnostic rung-1 prior.

Launched `rung2_K1_turnover` (PID 48144), `rung2_K2_nfw_turnover` (48150)
and `rung2_K2_cored_turnover` (48152), with 400 live points, dlogz=0.5,
seed 42 and 2*ndim slice steps. All three use preserved source and one numerical
library thread; the separate dynamical queue is unchanged. Manifests and logs
are in `results/launches/rung2_turnover_20260921/`. Initial checks found all
three samplers advancing past several thousand likelihood evaluations with
the exact common data fingerprint. No result or completion-time estimate is
inferred from this early progress.

## 2026-09-21 -- Simple single-transition companion fits launched

Following the agreed intermediate step, launched the existing smooth profile
`beta(r) = beta_0 + (beta_inf-beta_0) r²/(r²+r_beta²)` with both endpoints
uniform on [-1, 0.5] and `r_beta` log-uniform on [0.5, 100] pc. This retains
rung 1's diagnostic central assumptions; the central DF restriction is a
separate control. The five-parameter turnover runs retain their own priors.

Extended `bin/run_rung2_fits.py` with `--profile simple`, preserving its default
turnover recipe. The driver copies each saved rung-1 mass/distance model and
89-point data snapshot, stores the full resolved priors, and checks the
constant-beta maximum-likelihood limit at transition radii 0.5, 10 and 100 pc.
Global solver and constructor defaults are unchanged. Custom parent mass
priors survive reconstruction, and existing batch/output directories are
refused before any worker starts.

All 115 focused tests passed, including 10 new tests covering endpoint priors,
saved-model replay, all three projected velocity moments in the constant-beta
limit, preservation of custom parent settings, the turnover recipe, and run
preservation. Each production family passed 258 finite prior likelihood
checks and reproduced its parent likelihood in the constant limit. The
launched configurations and data fingerprints were checked independently.

The new labels are `rung2_K1_simple_beta` (PID 52088, 8 parameters),
`rung2_K2_nfw_simple_beta` (52090, 10), and `rung2_K2_cored_simple_beta`
(52092, 10). All use 400 live points, dlogz=0.5, seed 42, 2*ndim slice steps,
archived source and one numerical-library thread. Records, logs and test
verification are in `results/launches/rung2_simple_beta_20260921/`. Initial
checks found all three sampling beyond 3,000 likelihood evaluations, with no
tracebacks. The turnover K1 fit had already completed in 5,999.7 seconds;
the two turnover halo fits and the dynamical controller remain active.
The existing turnover source/configuration hashes were verified unchanged.

## 2026-09-21 -- Keep the plots directory image-only

Moved the nine JSON/ECSV plot-support files to `results/plot_data/` and the
earlier location-migration log and Finder metadata to `results/maintenance/`.
Updated all seven figure/report generators, current document references,
README and AGENTS output conventions. The dynamical report's portable LaTeX
bundle uses the same separation: PNGs in `plots/`, JSON in `results/plot_data/`.
The compiled report remains beside its LaTeX source in `docs/`.

All 112 PNG files and both numerical ECSV tables remain byte-for-byte
unchanged. Three generators were exercised against the relocated inputs;
all seven updated scripts parse. Migration checksums and the distinction
between regenerated metadata and historical generator hashes are recorded in
`results/maintenance/plots_support_migration_20260921.json`. Fit and simulation
run directories, including archived source used by active jobs, were untouched.

## 2026-09-22 -- Independent and conditional DM-density posterior checks

Launched two flexible-core fits at 06:52:53 UTC in batch
`results/fits/densitycheck_20260922T065245Z`: a free-density repeat (PID 77776)
and a rho20 = 2 solar-mass-per-cubic-parsec conditional posterior (PID 77778).
Both use 800 live points, dlogz 0.25, four times ndim slice steps, independent
seeds 1729/1730, archived source/data and one numerical-library thread.
Original model priors, including the million-solar-mass stellar floor, remain.

Added an explicitly serialised conditional model. Conditioning the original
log-uniform halo-mass prior leaves the scale-radius prior unchanged here;
full mass support is checked before sampling. The conditional/free evidence
ratio will estimate the marginal posterior PDF at density 2. It is not a
probability atom or a replacement for DM/no-DM model averaging.

All 71 targeted tests and a limited conditional UltraNest smoke test passed.
Both production workers were verified alive and advancing beyond their prior
initialisation. Original fits and the four active dynamical simulations were
preserved. See `docs/DM_DENSITY_POSTERIOR_CHECKS.md` for the prior derivation,
file locations, checks and interpretation planned after completion.

## 2026-09-22 -- Spatial density and stellar DF audit

Added `df_consistency.py`, the all-run driver `bin/audit_df_consistency.py`,
and `docs/DF_CONSISTENCY_AUDIT.md`. Inventoried 29 run records, including two
interrupted attempts and two still-active posterior checks. Checked 358,746
saved posterior rows (103,468 unique vectors), 800 checkpoint live points,
and 30 saved density-profile/repeat-start optima. All spatial densities are
non-negative. Reconstructed all 7,202 archived mass-profile draws exactly;
legacy likelihood replay differences are recorded separately.

Twenty of the 23 completed fits outside the invalid old solver archive fail
the cored-tracer/central-black-hole necessary condition beta0 <= -1/2. All
three turnover best samples pass that central bound but fail P'' >= 0 near
0.085--0.18 pc under a separable augmented-density assumption. Approximately
87--90% of their posterior samples fail the same screen. This excludes that
DF construction, not every nonseparable DF. No passing necessary-condition
screen is labelled a DF-positivity certificate. The checks concern the stellar
tracer; dark-matter/remnant DFs are unspecified by the Jeans fits.

All 71 focused tests passed (DF consistency, turnover anisotropy and run
replay). Best-sample classifications agree after doubling the radial grid.
Repeating all turnover posterior checks at 1,024 and 2,048 points changes
failure fractions by at most 0.11 percentage points relative to 512 points;
narrow negative intervals account for the borderline changes. The PNG figure
was inspected, and all numerical products live outside plots/. Existing fit
samples, priors and active workers were not changed.

## 2026-09-22 -- Positive stellar action-DF branch in AGAMA

Implemented a separate stellar DF model and joint photometric/kinematic pilot
fitter. Analytic AGAMA DoublePowerLaw DFs and their positive mixtures generate
both stellar density and anisotropy, with self-consistent stellar gravity.
The gravitational model retains an exact point-mass BH and optional static
Plummer remnants and truncated cored/cusped DM. Positivity is established for
the stellar DF; static remnants/DM do not yet have specified DFs.

An independent velocity integral exposed a central density discrepancy in
the standard AGAMA self-consistency helper. The new branch uses explicit
spherical velocity quadrature and AGAMA actions/potential construction,
checking off-grid density closure at doubled velocity order and total mass.
No negative DF or density values are clipped. The saved streaming subtraction
remains an explicitly documented moment approximation, with impossible values
rejected. Photometric errors must be supplied explicitly; the diagnostic pilot
adopts 0.1 mag / sqrt(relative weight), with one profiled zero-point.

All 94 focused tests passed (23 new controls plus 71 previous consistency,
anisotropy and replay tests). Four real-data evaluations -- no DM, cored DM,
cusped DM and a two-component stellar mixture -- passed independent projection
and complete resolution-doubling checks. The largest dispersion change was
0.011 observational errors. The mixture generates intermediate radial and
outer tangential anisotropy while keeping its DF positive.

A 96-call no-DM local fitting pilot completed in 294 seconds. Kinematic chi2
fell from 15881.77 to 761.52, with photometric chi2 662.30. The optimizer hit
its evaluation limit, so these are workflow-validation results, not converged
fits or new DM/BH constraints. The pilot's numerical checks passed. Added
`docs/AGAMA_DF_MODELS.md`, example configs, snapshots and PNG figures under
the established output conventions. Previous Jeans fits and dynamical runs
were not modified or restarted.

Mixture weights can also be fitted through log mass ratios followed by a
softmax, keeping every DF coefficient positive and the fractions normalized.
An independent fresh-process replay of the pilot reproduces its joint
objective and all 89 predictions exactly from the saved inputs.

## 2026-09-22 -- DF infrastructure methods note

Added `docs/agama_df_models.tex` and its 11-page compiled PDF beside the
source. The note derives the positive action-DF family, mixture weights,
self-consistent gravity, velocity integrals, projection and joint data
objective. It documents the numerical controls, the 94-test verification
record, the bounded unconverged no-DM pilot and current model scope.

Four numbered figures show example density/anisotropy profiles, numerical
accuracy, optimizer progress and the pilot against the data. The six-panel
data comparison occupies a landscape page. Added `bin/plot_df_writeup.py`
to regenerate the two new diagnostic figures from saved results, with
source checksums and plotted arrays in `results/plot_data/`. All figures
are PNGs under `plots/`. The LaTeX build has no unresolved references or
box warnings; all pages were rendered and inspected. No fits were launched
or modified for this documentation task.

## 2026-09-22 -- Independent fixed-potential DF flexibility challenge

Implemented an independent positive energy-based mock generator. Its
Osipkov--Merritt lowered-isothermal populations solve Poisson's equation jointly,
including luminous stars, concentrated dark remnants and an optional extended
dark component. The three cases comprise a single isotropic population and
two stellar mixtures with and without a halo. Their exact projected expectations
are evaluated at the saved 89 kinematic bins and 82 photometric radii; error
scales are adopted from the existing data/pilot, with no random noise or rotation.

Added a fixed-potential action-DF fitter with cached actions and analytic
derivatives. One, two and three positive DoublePowerLaw shapes were fitted from
three starts each. Wider-bound controls repeated the one- and two-DF fits.
All six batches finished: 45 starts, 6,605 evaluations, 44 optimizer terminations
and one capped alternative start whose solution was superseded. Each selected
fit passed resolution, outer-domain, independent AGAMA projection, true-potential
and mock-projection controls. The 104 focused tests passed.

Three DFs meet the predefined accuracy target in all mocks. At doubled
resolution the RMS/worst-bin residuals, in adopted observational errors, are
0.0346/0.1735 (single isotropic truth), 0.0098/0.0452 (mixed stars/remnants) and
0.0231/0.1135 (mixture plus halo). Two DFs approach this accuracy but retain
worst-bin residuals of 0.45--0.63 errors even with wider search bounds. One DF
performs substantially worse. Some parameters still contact the wider bounds,
so this is an attained-capacity comparison, not an impossibility theorem.

Added `docs/DF_CAPACITY_CHALLENGE.md`, five inspected PNGs in `plots/`, combined
plot data/checksums in `results/plot_data/df_capacity_20260922.json`, and preserved
run/code/input histories under `results/df/capacity*_20260922_*`. Components of
the fitted light DF are basis functions, not identified physical populations.
The fixed true potential already contains all gravity; these runs do not infer
DM or remnant masses. Noisy mass recovery, independent mass/light weighting,
instrument population selection and real-data fitting remain subsequent stages.

## 2026-09-22 -- Free-potential DF mass recovery launched

Implemented independent positive stellar mass and light mixtures in
`df_mass_recovery.py`. Each trial solves stellar self-consistency from the
mass DF and projects the light DF. Added explicit, separately serialized
dark-template masses and radial dilations. The first recovery test gives the
fitter matched dark density families while freeing their amplitudes/scales;
stellar action shapes are independent of the positive energy-based mock truth.
Dark spatial density is positive, but fitted dark DFs are not certified.

The 111 focused tests passed; a further seven-test rerun checked the added
shared-M/L path. Both mock starting models passed resolution, native projection
and cold replay controls. A four-evaluation, fixed-shape optimizer smoke test
reduced Q from 13.647 to 0.0901 and passed the same numerical checks. That smoke
test was deliberately capped and is not a completed mass-recovery inference.

Launched the preserved batch `results/df/mass_recovery_20260922T151408Z/` with
two single-threaded workers. Six noiseless fits (two injected clusters, three
starts each) fit weights, then action scales, then all 27 coordinates. For each
case passing a converged, numerically checked noiseless accuracy gate, the batch
adds one shared-M/L control and eight noise realizations with two starts each.
The full planned batch has 40 jobs. The noise ensemble is a pilot study of
recovery scatter, not a posterior-coverage certification or a real-data fit.

All workers import frozen scientific source and saved data/initial-shape inputs.
The new reporter writes PNGs under `plots/`, numerical/provenance records under
`results/plot_data/`, and `docs/DF_MASS_RECOVERY_STATUS.md`. Methods and limitations
are documented in `docs/DF_MASS_RECOVERY.md`. Previous fits and simulations were
not restarted or stopped.

## 2026-09-23 -- Compact DF matched-family recovery infrastructure

Added `compact_recovery.py` and `bin/run_compact_df_recovery.py` for the next
recovery stage. Two refined regularized-DF equilibria inject either zero DM or
a cored halo with rho(20 pc)=2 solar masses per cubic parsec; both contain
concentrated remnants. The planned six fits use two starts per experiment:
free halo on each mock, plus a separate no-halo control on the no-DM mock.
All five stellar controls and both remnant coordinates are free. Free-halo
fits also vary rho20 (including exactly zero) and scale radius.

The driver preserves source and mock snapshots, checkpoints each improved
solution atomically, and resumes interrupted fits from their saved parameters.
Two single-threaded workers and per-attempt evaluation budgets bound the run.
Numerical convergence, observable residuals, and physical recovery have
separate gates; an accurate observable fit cannot certify a correct mass split.
The acceptance criteria and scope are in `docs/COMPACT_DF_RECOVERY.md`.

The focused verification suite passed 63 tests, including five new checks of
mock construction, asymmetric residuals, full resolution refinement, and the
distinction between observable fit quality and physical recovery. End-to-end
mock preparation and bounded fitting are the next launch checks. Independent
stellar families, noisy mocks, and real-data inference remain later stages.

## 2026-09-23 -- Recovery smoke test exposed a polar-mapper singularity

Both v1 mock truths passed refinement checks (maximum shifts 0.00703 and
0.00686 adopted errors). The two-evaluation optimizer smoke test saved its
checkpoints but failed while refining the selected perturbed model. The
native action mapper had returned infinite velocities at cylindrical R=0
for a polar representative, causing the downstream frequency integration to
raise `integrateGL: order is too high (not implemented)`.

Preserved the diagnostic scripts, model checkpoint, potential, and input orbit
in version control. The first diagnostic was incorrectly placed in `/tmp`;
it has been moved into the repository and AGENTS.md now records the requirement
that scientific/debugging code remain reviewable, together with regular commits
and pushes. No experiment code will rely on that removed temporary script.

The frequency evaluator now uses an equatorial orbit with the same Jr and L.
A new regression catches the angular-momentum error in the old polar mapping.
All 64 focused tests pass, and the previously failing complete refined-model
rebuild succeeds. V1 inputs and the failed smoke-test record remain preserved.
The corrected source will be frozen for a fresh v2 recovery batch.

## 2026-09-23 -- Compact recovery v2 launched

Regenerated both mock truths from corrected commit `01992f0`; nominal-resolution
shifts remain 0.00703 and 0.00686 adopted errors. A separate two-evaluation smoke
test completed its checkpoint, budget-stop, cold-replay, and refined validation
workflow. Its numerical checks passed: maximum cold shift 0.0228 errors,
refinement shift 0.0000662 errors, and native projection discrepancy 0.0326%.
It is deliberately unconverged and is not a successful mass-recovery result.

Launched `results/df/compact_recovery_20260923_v2/` at 00:02:56 UTC with two
single-threaded workers, six jobs, and 180 trial equilibrium evaluations per
job. The first jobs fit a free halo to the cored-DM and no-DM mocks. Results,
checkpoints, and source/input snapshots remain in the batch directory; the
progress PNG and its data follow the project output convention. Scientific
recovery gates remain pending until the full fits and validations finish.

## 2026-09-23 -- Compact recovery v2 finished: all six fits fail the gates

All six jobs in `results/df/compact_recovery_20260923_v2/` exited cleanly by
01:34 UTC (`report.json`). None terminated within the 180-evaluation budget, so
every fit is recorded as unconverged and fails the combined gate.

| Fit | Numerical | RMS / max (err) | M_star error | rho20 (true) |
|---|---|---|---|---|
| cored, free halo, start0 | pass | 0.006 / 0.025 | 4.1% | 2.53 (2.0) |
| cored, free halo, start1 | fail | 0.006 / 0.022 | 8.5% | 3.08 (2.0) |
| no DM, free halo, start0 | pass | 0.015 / 0.073 | 14.9% | 1.42 (0) |
| no DM, free halo, start1 | pass | 0.024 / 0.131 | 24.4% | 2.34 (0) |
| no DM, no halo, start0 | pass | 0.076 / 0.356 | 0.7% | fixed 0 |
| no DM, no halo, start1 | pass | 0.037 / 0.147 | 0.2% | fixed 0 |

Columns: numerical = cold replay, refinement, and projection checks; RMS/max
= refined noiseless-mock residuals in adopted errors; rho20 = halo density at
20 pc in solar masses per cubic parsec.

Free-halo fits match the noiseless observables to 0.006-0.024 errors RMS, and
total enclosed mass stays within 0.2-4.8% at every diagnostic radius. The
star/halo split is not recovered, however. In the no-DM mock a spurious halo
(rho20 1.4-2.3) replaces 15-24% of the stellar mass. In the cored mock the two
starts give rho20 = 2.53 and 3.08. Because no fit converged, part of this may
be optimizer budget; the residual levels nonetheless show a nearly flat
objective along the stellar-mass/halo-density direction. The no-halo start1
objective was still falling at the cap (0.255 to 0.238 over the last 10
evaluations). The start0 no-halo fit misses the maximum-residual gate
(0.356 > 0.3).

Cored start1 fails numerical validation: independent projected moments
differ by 0.81% in pmr and 0.46% in los, against the 0.5% threshold. Its cold
and refinement shifts pass (0.0021 and 0.000089 errors).

The DF models are therefore not ready for DM inference. Candidate follow-ups,
not yet launched: a profile scan in fixed rho20 to map the degeneracy
directly, an explicitly recorded extension of the capped fits, and a check of
projection accuracy at the cored start1 solution. DF-focused unit tests
(`test_regularized_df.py`, `test_compact_recovery.py`) pass: 24 tests.

## 2026-09-23 -- Degeneracy figure for compact recovery v2

Added `bin/plot_compact_degeneracy.py`, which reads saved per-trial logs and
refined profiles only (no model evaluations). Output:
`plots/compact_recovery_20260923_v2_degeneracy.png`, with data and provenance in
`results/plot_data/compact_recovery_20260923_v2_degeneracy.json`. The top row
shows optimizer trials in (M_star, rho20) coloured by objective; these are
optimizer paths, not a profile likelihood. The bottom row shows refined
stellar and halo enclosed mass relative to the injected total.

Every free-halo path reaches the low-objective region (below about 1) after a
few steps, then drifts slowly along a direction of rising rho20 and falling
M_star. It stops at the budget cap well away from the truth. The refined
profiles show the halo absorbing the missing stellar mass beyond about 5 pc.

The objective at the injected truth, evaluated at fitting resolution, is
2.3e-4 in both mocks (RMS 0.0012 errors over 171 points). The best cored fits
reach 6.3e-3 and 6.8e-3, and the no-DM free-halo fits 0.041 and 0.10. The truth
is therefore still a strictly better minimum and the fits are unconverged, not
exact degeneracies. However, these objective differences are far below the
chi-square scatter expected with real noise (sqrt(2 x 171) = 18.5). The
stellar-mass/halo-density split is thus practically unconstrained by these
observables. A fixed-rho20 profile scan would measure that constraint directly.

## 2026-09-23 -- Compact recovery v2 mocks compared with the observed data

Added `bin/plot_compact_mocks.py` (reads saved snapshots only). Output:
`plots/compact_recovery_20260923_v2_mocks.png`, with data in
`results/plot_data/compact_recovery_20260923_v2_mocks.json`. It shows both
noiseless mocks with their adopted errors, the observed values they replace,
and (cored - no DM)/error per bin. The data span 0.05-68 pc projected, at
0.0263 pc per arcsec for 5.43 kpc.

The injected truths are not Omega Cen-like. Central mock dispersions are
12.8 km/s (LOS) and 0.49 mas/yr (PM) for no DM, against about 18-21 km/s and
0.76-0.86 mas/yr observed, roughly 35-40% low. The mock surface-brightness
profile has a flat core out to about 10 pc and falls later than observed.
The observed profile declines from about 3 pc. The 3e6 Msun stellar mass
with J0 = 300 pc km/s therefore gives a cluster too light and extended for
the real data. The recovery test is a self-consistency check of the
machinery; its degeneracy amplitude need not transfer to the observed regime.

The two mocks are strongly separable bin by bin. The HST PM dispersions differ
by up to 43 (radial) and 20 (tangential) adopted errors near 3-6 pc, and the
Gaia PM, MUSE LOS, and outer photometry by about 5. The fitted degeneracy is
therefore not a failure to see the halo's effect on the data. The stellar DF
and remnant parameters absorb that effect instead.

## 2026-09-23 -- Why the v2 fits hit the budget

With 9 free coordinates, each forward-difference Jacobian costs 10 equilibrium
evaluations, so 180 evaluations allow only about 18 Gauss-Newton steps (about
28 min per fit). The best objective was still falling at the cap in every fit:
between evaluations 150 and 180, cored start0 went from 0.016 to 0.0038 and
cored start1 from 0.12 to 0.0072.

The scipy tolerances (ftol=1e-5 relative, xtol=2e-4, gtol=2e-3) are tighter
than the evaluation noise floor. Warm-started and cold-replayed objectives at
the same parameters differ by up to 0.0025 (cored start0: 0.00377 warm,
0.00629 cold), with maximum residual shifts of 0.002-0.02 errors. Over the
0.002 finite-difference step, that noise enters the Jacobian at order 1-10
per unit coordinate. Near the minimum the relative-ftol test cannot be
satisfied reliably, so a larger budget alone may still end at the cap.

Proposal, not implemented: an absolute stopping rule (e.g. objective change
below 0.01 over two successive accepted steps, well below a meaningful
delta-chi2 of 1), a larger cap (about 600 evaluations), and cold or tighter
equilibrium solves inside the Jacobian so that its noise stays below the
tolerance.

## 2026-09-23 -- Fitter convergence controls and observed-data DF fits

Jacobian noise probe (`bin/probe_compact_jacobian_noise.py`; record in
`results/maintenance/jacobian_noise_compact_recovery_20260923_v2_cored_dm_free_halo_start0.json`)
at the v2 cored start0 optimum. With iteration tolerance 5e-4, repeated
base-point evaluations with different seeds shift residuals by up to 0.0060
errors. Forward-difference columns differ from cold-solve columns by 0.1-7%
(chained seeds) or 1.5-2% (base-point seeds). With tolerance 5e-5 these fall
to 0.00059 errors and 0.01-0.65%, at 16.7 s instead of 13.5 s per evaluation.
Jacobian noise was therefore modest; the 180-evaluation cap was the main cause
of the v2 non-convergence.

Driver changes in `bin/run_compact_df_recovery.py` (v2 keeps its frozen code):
- `AbsoluteStop` (in `compact_recovery.py`, 2 new tests) ends a fit when the
  objective at accepted iterates improves by less than `--stop-delta` over
  `--stop-window` steps; the result is recorded as optimizer-terminated.
- `--jacobian-seed base` seeds every probe from the base point's stars.
- `--iteration-tolerance` sets the batch fitting tolerance.
- `prepare-real` fits the observed pilot snapshot (89 kinematic bins,
  82 photometric radii) with free-halo and no-halo branches, two starts each.
  Starting objectives are 2604 and 8943.

Observed-data acceptance criteria, fixed before fitting and judging fit
quality only: chi2_kin/N < 1.3 (89 bins), chi2/n < 2 in every dataset,
photometric RMS < 0.05 mag (errors are adopted, not measured), optimizer
terminated by the rule, and numerical validation passed. These do not test the
mass decomposition.

Smoke batch `results/df/observed_df_smoke_20260923/` (12 evaluations per job)
completed the budget-stop, validation, gate, and report path. All numerical
checks passed (cold shift at most 0.0057 errors, refinement at most 0.0014,
projection at most 0.09%). The capped fits reached chi2_kin/N = 8.5-11.6, as
expected at 12 evaluations.
