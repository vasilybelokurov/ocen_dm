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
