# Known gaps and inconsistencies

Found in the 2026-10-13 audit. Each item is open unless marked otherwise.


1. **Baseline counts use the old (schema-1) products.** *In progress.*
   - `dftwo_counts_20260923`, `rho20scan_*` and `mpsscan_*` used schema-1 counts (star-weighted bin average, biased; see `counts.py`). Their frozen inputs still replay exactly; the old products are kept as `*_schema1_20260923.*`.
   - The canonical `ocen_counts_*.json` are schema 2 since commit `1af8c33`. The base model was refitted with them on 2026-10-13 (`results/df/dftwo_counts_area_20261013`); the scans have not been rerun.
2. **Docstring and metadata drift.**
   - *Fixed 2026-10-13:* the `bin/build_count_profiles.py` docstring described the old quantile nodes.
   - `selection/field_template.py:107` writes `"parallax_over_error_max": 5` into the product's `cuts` metadata, but no parallax cut is applied (`parallax_cut: none`). The metadata is misleading.
   - `likelihood.py:251 _edr3_ours` admits that its note string describes an obsolete construction.
3. **Stale provenance records.**
   - `provenance/manifest.json` and `checksums.txt` list `vasiliev2021_gc_members/clusters.zip` at `data/raw/vasiliev2021_gc_members/clusters.zip`. That directory is empty, and the copy actually used is the external `<code>/apogee_halo_rotation/.../vasiliev_baumgardt_2021_clusters.zip`.
   - The manifest also keeps an obsolete `kuzma2026_members.*` row with status `manual` and `exists=false`, next to the real `kuzma2026_members.txt` entry.
   - `data/raw/baumgardt_gc_catalogue/` is empty because the file is referenced from `CAT`, which is intended.
   - The T3 product's `ocen_raw.manifest` says `unregistered`, because it was built before the manifest entry existed.
4. **Registry incomplete.**
   - `provenance/datasets.yaml` `omegacat_vi_kinematics` has `files: []`, although 22 files are downloaded and md5-verified in the manifest.
   - `datasets.yaml` omits the following:
     - STREAMFINDER 2021 and 2024, HST LP (all in `configs/data.yaml` only).
     - Tango Zenodo 4300977, `gc_catalog_full.fits`, Gaia XP/Andrae pulls.
     - Hunter, Portail, McMillan and MWPotential2014_mod potentials; DB98; GSE contours.
     - MIST/minimint BC tables; the K8 WSDB field template; T2 (covered by product metadata only).
     - Baumgardt+2018 RVs, planned in `docs/DATA_GAPS_PLAN.md` but never ingested.
5. **Products without builders.**
   - S2 `stream5455_gaia_errors.npz`, which every stream fit depends on.
   - S5 XP pulls.
   - The S1 MRT→FITS conversion, which has no checksum for the FITS.
6. **External path dependencies outside the repository and CAT.**
   - The VB21 zip under `<code>/apogee_halo_rotation`.
   - The oCen_bar potentials and GSE contours under `OB`, reached through two different path conventions.
   - McMillan17 from the installed AGAMA package, so results depend on the AGAMA version.
7. **Inconsistent systemic values across modules.** *Resolved 2026-10-16.*
   - `cluster.py` is now the single reference for centre, systemic PM, v_los and distance, with published sources (see [cluster data](cluster.md)); about 110 typed-in copies in 47 files import it, and `tests/test_cluster_reference.py` bans new ones.
   - Products built before that date keep the values they were built with (e.g. the K8 template frame, (−3.2514, −6.7444) from K6, v_los 232.6).
8. **Selection inconsistencies within the light and tracer chain.**
   - L2 (Jeans tracer and K4 depth term) uses oMEGACat `selection_hq_f625w` star counts inside 25″.
   - C1 (DF baseline) deliberately avoids that flag as 12–15 % radially incomplete.
   - K10 `hst_quality` uses a flag that contains a CMD cut (project journal, 2026-10-08; flagged and not changed).
9. **Unverified references** (per `datasets.yaml` or the project journal):
   - Trager 1995 ADS link (`verified: false`).
   - Baumgardt+2019 DOI (`verified: false`).
   - Baumgardt GC database version.
   - Kuzma+2021 (registered, no data, `verified: false`).
   - Hunter+2024 paper.
   - Origins of the Portail17 and MWPotential2014_mod .ini files.
   - GSE contours.
   - Andrae+2023 reference.
   - The meaning of the VB21 profile `undocumented_col12`.
10. **Registered or present but unused.**
    - The O1 oMEGACat files: alternative binnings, equipartition, `.dat` twins, notebook.
    - S3 STREAMFINDER 2021.
    - The unused Hunter24 variants, `pot_bar_mod`, PriceWhelan19/22, MWPotential2014.
    - The large XP catalogues in `CAT`.
    - `datasets.yaml` entries with no data: `noyola2006_sbp` (no NGC 5139 rows), `omegacat_ii_astrometry` (no files, "not used until validated"), `kuzma2021_edr3_tails`, `youakim2023_debris`.
11. **Quarantined stale products.**
    - `proc/kinematics/_stale_pre_centre_fix/` holds three products built before the 2026-09-19 centre correction: `ocen_outer_fit_quality`, `ocen_pm_dispersion_ours` and `ocen_pm_dispersion_mixture`. Its README says they feed only historical figures.
    - `plotting/constraints.py:437` writes a fresh `ocen_outer_fit_quality.ecsv` to the non-stale path when `save=True`. None currently exists there.
    - No code reads the other two stale files.
12. **Rebuild date versus snapshot.**
    - `ocen_pm_dispersion_edr3_ours.ecsv` was rebuilt on 2026-10-11, after the baseline. The values were verified identical today, but the loaders auto-build missing products, so a fresh checkout can silently produce data that differ from the frozen snapshots (README:162–165).
    - Run directories freeze their inputs (`data_snapshot.json`, count JSON with sha256), which protects the recorded runs.
13. **Heavy and uncommitted items.**
    - `results/streams/m5_gaia_knee.npz` (242 MB) is a one-off diagnostic cache.
    - `bin/streams/fit_bar_ocen.py` has uncommitted edits, and `check_stream_coordinate.py` (a consumer of S4) is untracked.
