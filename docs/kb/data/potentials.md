# Potentials and auxiliary inputs

Shorthands (`OB` = oCen_bar clone, `CAT`) are defined on the [Data products](index.md) page.

## P1–P7. Potentials and auxiliary model inputs

**P1. Hunter+2024 MW potentials.**
- **Location.** `OB/agama_potentials/MWPotentialHunter24_*.ini`, from the oCen_bar clone (github.com/adllmr/oCen_bar per JOURNAL.md:3767–3772).
- **Path conventions.** `streams/restricted.py:37 HUNTER24_DIR` is a relative path; `tails/bar_migration.py:25` uses the `OCEN_BAR_DIR` env variable with default `$OCEN_BAR_DIR`.
- **Files used:**
  - `axi` (73,032 B): `restricted.py`, `lmc.py:63`, `mw_variants.py:45`, `bar_migration.py:100` and others.
  - `full` (227,828 B): `restricted.py:92`.
  - `baryon_full` and `baryon_axi`: bar = A × (baryon_full − baryon_axi).
  - `spiral` (735,071 B): `streams/spiral.py:82`.
- **Unused:** `nonaxi`, `rotating`, `rotspiral`, `pot_bar_mod`, `PriceWhelan19/22`, `MWPotential2014.ini`.
- **Citation.** arXiv:2403.18000 (JOURNAL.md:6896). JOURNAL.md:7442 says the paper itself was not retrieved: **UNVERIFIED**.

**P2. Portail17.ini and Portail17_axi.ini** (241,764 B and 47,001 B).
- Used through opt-in `bar_model='portail17'` (`restricted.py:72,98`).
- Paper: Portail+2017, arXiv:1608.07954 (JOURNAL). The provenance of the .ini files is **UNVERIFIED**.

**P3. MWPotential2014_mod.ini** (469 B).
- Used by `tails/bar_migration.py:155`. Source **UNVERIFIED**.

**P4. McMillan17.ini.**
- The code uses the AGAMA package copy (`site-packages/agama/data/McMillan17.ini`, AGAMA 1.0.152) via glob: `tails/progenitor_orbits.py:222`, `restricted.py:114`, `bin/streams/scope_jacobi.py:39`, `bin/nbody/run_nbody.py:52`.
- A second copy exists at `OB/agama_potentials/McMillan17.ini` (797 B).

**P5. DB98_Model1.ini.**
- File: `configs/potentials/DB98_Model1.ini` (1,284 B). This is the only file in `configs/potentials`.
- Model: Dehnen & Binney 1998, MNRAS 294, 429, arXiv:astro-ph/9612059, Model 1 (the Ibata+2019 host), built from DB98 Tables 1 and 3.
- Checks at JOURNAL.md:6258–6260: v_c(R0) = 222.0 km/s; Σ(|z| < 1.1 kpc) = 67.7 M☉ pc⁻².
- Consumers: `restricted.py`, `backtrack_fimbulthul.py:34`, `sgr_*`, `run_db98_ibata_frame.sh`, `run_rotation_tests.sh`.

**P6. Vasiliev, Belokurov & Erkal 2021 "Tango for three" archive.**
- **Location.** `CAT/vasiliev21_tango_zenodo4300977/`.
- **PROVENANCE.txt records:**
  - MNRAS 501, 2279; arXiv:2009.10726
  - https://doi.org/10.5281/zenodo.4300977 (v2; the paper cites v1, 4038141)
  - md5 `151950f0…`; downloaded 2026-10-12
  - units 232,500 M☉, G = 1, time unit 0.9778 Gyr
- **Files used:**
  - `potentials_{axisym,triax}/trajlmc.txt`
  - `accel.txt`
  - `potentials_axisym/potential_evolving.ini`
  - `mw00.pot`
- **Consumers.** `bin/streams/lmc_tango_check.py`, `lmc_vs_tango.py`. The LMC parameters in `streams/lmc.py` follow Tango Section 3.1.
- **Pitfall** (JOURNAL.md:7497). AGAMA reorders the components of `potential_evolving.ini`, so select them by type.
- **Registration.** Not registered anywhere in `provenance/` or `configs/`.

**P7. GSE E–Lz contours.**
- **File.** `OB/artifacts/gaiadr3_gse_elz.fits` (17,280 B), labelled "Belokurov+2023" in JOURNAL.md:3770 with no reference: **UNVERIFIED**.
- **Consumer.** `tails/bar_migration.py:27,157`.

**Internal tails outputs, not external data.**
- `results/tails/{progenitor_orbits_summary.ecsv, progenitor_orbit_tracks.npz}`, from `tails/products.py:48`.
- `bar_migration_*`, `friction_test.ecsv`, `class2_barred.ecsv`, `resonances.ecsv`, `dm_density_grid.ecsv`, `po_*.ecsv`.
- The N-body scripts read only `results/nbody/<model>/{ics.npz, orbit/snap_today.npz}`.
