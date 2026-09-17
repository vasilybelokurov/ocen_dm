# Plan: closing the data gaps before Milestone 3

*Written 2026-09-16 after inspecting the real VizieR table structures and confirming the
Gaia join on WSDB. Nothing below is assumed from memory; unverified items are marked.*

## Goal

Give the kinematic fit (Milestone 3) and the orbit/stream models (Milestone 6) every input
they need as provenance-tracked processed products: a light profile for the MGE, the
cluster's systemic phase space with covariance, Gaia PM correlations for the tail
catalogue, an outer kinematic profile, and the tail survey's selection footprint.

## The finding that reorders the priorities

The specification says the tails are sensitive to extra mass at **50–200 pc** (section 1). At
5.43 kpc that is **32′–127′**. The HST kinematics we hold stop at 346″ ≈ 6′ ≈ 9 pc. So the
radii that matter most for the dark-matter question are covered by *none* of the internal
kinematics ingested so far. They are covered by:

- **Baumgardt+ 2019, table4** — Gaia DR2 proper-motion dispersion profiles of 102 clusters
  (`Name, NPM, r, sigma, E_sigma, e_sigma`), reaching the cluster outskirts;
- **Baumgardt+ 2018, tabled** — individual radial velocities with membership probability
  (`Cluster, RAJ2000, DEJ2000, RV, e_RV, Dcen, PMemb`), an outer σ_LOS dataset;
- our own **Kuzma 2025 members** (1.7′–300′, but only G₀ < 16, i.e. bright giants).

What was going to be "optional" (WP5) is therefore essential. Without it the DM constraint
would come entirely from the tails, with the internal kinematics constraining only the
inner 9 pc — which is exactly the degeneracy the project is about, but with no inner-outer
lever arm on the internal side.

## Assumptions and unknowns

| Item | Status |
|---|---|
| Trager+ 1995 `tables`: `logr` is log10 of radius in arcsec, `muV` is V-band SB in mag/arcsec², `Weight` and `DataSet` flag heterogeneous sources | FIELD units to be read from the VOTable on ingest; **UNVERIFIED** until then |
| Noyola & Gebhardt 2006 `photdata`/`smdata`: HST core SB, calibrated to V | filter/zero-point **UNVERIFIED**; read the ReadMe |
| Only the *shape* of the light profile enters the MGE — the stellar M/L is a free nuisance parameter (spec section 3.1) — so the SB zero-point and the extinction correction cancel | Verified by construction; removes the largest risk in WP2 |
| Vasiliev & Baumgardt 2021 `tablea1.corr` is the pmRA–pmDE correlation of the *systemic* PM | consistent with the column set; confirm from the ReadMe |
| Distance: oMEGACat VI kinematic 5.43 kpc vs Baumgardt & Vasiliev catalogues (`Rsun`, `Dist`) | a real tension to record; **keep D as a parameter**, not a constant |
| WSDB `gaia_dr3.gaia_source` carries `pmra_pmdec_corr`, `ruwe`, `parallax`; a Kuzma 2025 `source_id` resolves | **Verified** 2026-09-16 |
| Kuzma 2025 footprint is a ~5.1° circle with a G₀ < 16 limit | to be reconstructed from the catalogue, not assumed |

## Status update, same day

Two of the gaps closed from the user's canonical catalogue store `~/data/catalogues/`
(retrieval kind `local`, verified in place, never copied):

- **WP3 done** via `gc_catalog_updated.fits` (Baumgardt database compilation, UNVERIFIED against
  upstream): systemic PM, D, RV, mass, r_h. The distance tension assumed above is absent —
  both sources give 5.43 kpc.
- **WP5 largely done** via `gc_members_gaia_vasiliev.fits` (Vasiliev & Baumgardt 2021): 228,055
  stars to G = 21 within 63 pc with PM covariance, plus the authors' σ_PM and v_rot profiles.
  Still to fetch: the Baumgardt 2019 DR2 profile for radii beyond 40′.
- **New constraint on WP2/Milestone 3:** EDR3 σ_PM is 20–25% below HST σ_PM at 1–3′,
  converging by 5′. Candidate causes are crowding systematics or energy equipartition
  (bright giants vs faint stars). The kinematic likelihood must model σ(m) or use a common
  magnitude range; the two profiles must not be stacked.
- **WP4 verified feasible:** one `sqlutilpy.local_join` returns all 157,481 Kuzma 2025
  source_ids in 162 s with 100% match; median |pmra_pmdec_corr| = 0.11.

**2026-09-17:** WP1 done for Trager 1995 and Baumgardt 2019 (Noyola 2006 checked: ω Cen absent,
not applicable). **WP2 done** — `light_model.py`, 11-Gaussian MGE, rms 0.18 mag, R_h = 280″ =
7.37 pc. **WP5 done** for the DR2 profile (4.7–46 pc). Open: WP4, WP6, star-count cross-check.

## Work packages

### WP1 — Registry, fetch, loaders for five VizieR catalogues (~1 h)

Extend `provenance/datasets.yaml` (kind `vizier`, already supported) with the exact tables,
each filtered to NGC 5139 at retrieval:

| Dataset | Table | Rows we need |
|---|---|---|
| `trager1995_sbp` | `J/AJ/109/218/tables` | `Name` = NGC 5139 |
| `noyola2006_sbp` | `J/AJ/132/447/photdata`, `smdata` | `NGC` = 5139 |
| `vasiliev2021_gc_pm` | `J/MNRAS/505/5978/tablea1` | one row |
| `baumgardt2019_gc` | `J/MNRAS/482/5138/table1`, `table4` | one row; the σ_PM(r) profile |
| `baumgardt2018_gc` | `J/MNRAS/478/1520/table2`, `tabled` | one row; individual RVs |

Loaders follow the existing pattern (roles, declared units, `build_product`), with the
VizieR byte-reproducibility caveat recorded as for Fimbulthul. Tests on synthetic fixtures
with the real column names.

**Check:** row counts printed; the VB21 row's `pmRA, pmDE` against the Kuzma 2025 member
centroid we already measured, (−3.2, −6.7) mas/yr.

### WP2 — Light profile → MGE (`src/ocen_dm/light_model.py`) (~2 h)

1. Merge Trager (ground, wide) and Noyola (HST, core) V-band profiles on a common
   log-radius grid, honouring Trager's `Weight`; overlap region used to check relative
   zero-points.
2. Convert μ → I ∝ 10^(−0.4 μ) (shape only, see assumptions).
3. Fit a 1-D MGE in projection: Gaussians with widths fixed on a log grid, amplitudes by
   non-negative least squares, zero components pruned. No external MGE package needed.
4. Deproject analytically (spherical): Σ_j ∝ exp(−R²/2σ_j²) ↔ ρ_j ∝ exp(−r²/2σ_j²) with the
   same σ_j — which is exactly our `MGE` class.
5. Angular → physical scale through the distance *parameter*.

**Checks:** reproject and compare to the data (rms in mag, residual figure to `plots/`);
half-light radius against Harris/Baumgardt; total luminosity against the catalogue M_V
for an assumed M/L. **Decision needed:** V-band only (Trager+Noyola) for v1, or also the
star-count profile from our own catalogues as a cross-check (recommended as check, not
input).

### WP3 — Systemic phase space (`configs/ocen_phase_space.yaml` + loader) (~0.5 h)

RA, Dec (SIMBAD/Harris), μ_α*, μ_δ, corr (VB21), v_sys (Baumgardt 2018 ⟨RV⟩; oMEGACat VI),
distance (three values, each with source), all with uncertainties and citations.

**Checks:** VB21 PM vs our member centroid; literature ⟨RV⟩ vs our Kuzma 2026 member mean
(234.8) and oMEGACat median (232.6); the three distances tabulated side by side.

### WP4 — Gaia DR3 covariance join for Kuzma 2025 (`src/ocen_dm/selection/gaia_covariance.py`) (~0.5 h)

`sqlutilpy.local_join` of the 157,481 `source_id`s against `gaia_dr3.gaia_source` for
`pmra, pmdec, pmra_error, pmdec_error, pmra_pmdec_corr, parallax, parallax_error, ruwe`.
Product: `kuzma2025_periphery_gaia.ecsv`.

**Checks:** 100% join completeness; catalogue `pmRA/pmDE` equal to Gaia's to 1e-6 (they
are DR3 values); `corr` ∈ [−1, 1]; RUWE distribution reported.

### WP5 — Outer kinematics (~1 h)

From WP1's Baumgardt tables: the NGC 5139 σ_PM(r) profile (DR2) and the individual RVs
with `PMemb`. Products in `data/processed/kinematics/` with the same asymmetric-error schema
structure as the oMEGACat profiles, plus a **shared-stars / independence note**: the DR2
profile and our Kuzma 2025 members may overlap; the Baumgardt RVs and the oMEGACat MUSE
sample overlap in the core.

**Checks:** in the 3′–6′ overlap, Baumgardt σ_PM vs oMEGACat σ_PM after converting with a
common D; a combined figure σ(r) from 0.1 pc to the outskirts.

### WP6 — Selection footprint for Kuzma 2025 (`configs/selection_kuzma2025.yaml`) (~0.5 h)

Reconstruct from the catalogue: angular footprint (circle radius, centre), magnitude limits
(edge of the G₀ distribution), any Pristine coverage holes (map of stars with finite
`CaHK_0`). Store as a declarative selection; implement `selection/kuzma2025.py` that applies
it identically to real and mock catalogues.

**Check (spec section 13, item 12):** the same function applied to the real catalogue and to
a mock produces identical bookkeeping columns.

## Order and dependencies

```
WP1 (fetch) ──► WP2 (MGE) ──► Milestone 3
     │
     ├──────► WP3 (phase space) ──► Milestone 6
     │
     └──────► WP5 (outer kinematics) ──► Milestone 3
WP4 (Gaia join, WSDB)  ── independent ──► Milestone 8
WP6 (footprint)        ── independent ──► Milestone 8
```

Total ≈ 5–6 h. Each WP ends in a commit and a journal entry; figures for WP2 and WP5.

## Outputs

- `provenance/datasets.yaml`: +5 datasets; `data/processed/literature/*.ecsv`
- `src/ocen_dm/light_model.py`, `plots/light_profile_mge.png`
- `configs/ocen_phase_space.yaml`, `src/ocen_dm/data/phase_space.py`
- `src/ocen_dm/selection/{gaia_covariance,kuzma2025}.py`, `configs/selection_kuzma2025.yaml`
- `data/processed/kinematics/{baumgardt2019_sigma_pm,baumgardt2018_rv}.ecsv`, `plots/kinematics_all_radii.png`
- tests for every loader and for the MGE fit (reprojection error, deprojection identity)

## Decisions requested

1. Distance: treat as a **free parameter** with the three literature values as priors/cross-checks (recommended), or fix to 5.43 kpc?
2. Light profile: V-band Trager+Noyola only for v1 (recommended), with our star counts as a cross-check?
3. Include the Baumgardt 2018 individual RVs as an outer σ_LOS dataset (different instrument systematics from MUSE) — yes, flagged as a separate likelihood term (recommended)?
