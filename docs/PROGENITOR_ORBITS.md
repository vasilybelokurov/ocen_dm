# Progenitor orbits for the host-disruption simulations

*2026-09-20. Code: `src/ocen_dm/tails/progenitor_orbits.py`; products: `python -m ocen_dm.tails.products`
→ `results/tails/progenitor_orbits_summary.ecsv`, `results/tails/progenitor_orbit_tracks.npz`,
`plots/progenitor_orbits.png`; tests: `tests/test_progenitor_orbits.py`.*

## Purpose

Constrain ω Cen's dark-matter content from the other end: simulate the tidal disruption of the
host dwarf and see what the surviving nucleus keeps. Before any simulation we need a small set
of distinct orbital histories. Three classes were defined by the user; this note fixes the
initial-condition assumptions and picks 2–3 orbits per class.

## Present-day phase-space point (common to all classes)

Observables: centre from `cluster.py`, D = 5.43 kpc, μ = (−3.257, −6.730) mas/yr,
v_los = 232.7 km/s. Solar parameters: R₀ = 8.21 kpc, v₀ = 233.1 km/s, z☉ = 20.8 pc,
(U,V,W)☉ = (11.1, 12.24, 7.25) km/s (McMillan 2017 / Schönrich 2010).

| quantity | value |
|---|---|
| galactocentric (x, y, z) | (−4.90, −4.07, 1.41) kpc, r = 6.52 kpc |
| (v_x, v_y, v_z) | (95.2, −28.9, −87.1) km/s (right-handed, disc L_z < 0) |
| L_z | −529 kpc km/s (retrograde; prograde-positive convention) |
| E (McMillan 2017) | −1.85 × 10⁵ km² s⁻² |

Checked: galpy and astropy agree to 0.01 kpc / 0.1 km/s once galpy's `solarmotion` is given the
peculiar motion only (galpy adds v₀ itself; including it double-counts 233 km/s and makes the
orbit prograde — the first version of the code had this bug).

## Class 1 — today's orbit in a static axisymmetric potential

No assumptions beyond the potential. Three standard potentials bracket the uncertainty:

| potential | r_peri | r_apo | e | z_max |
|---|---|---|---|---|
| McMillan 2017 | 1.57 | 7.04 | 0.63 | 2.88 |
| MWPotential2014 (Bovy 2015) | 1.97 | 6.76 | 0.55 | 2.72 |
| Irrgang 2013 model I | 1.25 | 7.20 | 0.70 | 3.06 |

(kpc; backward integration 10 Gyr, galpy `dop853_c`.) McMillan 2017 reproduces the 2026-09-19
AGAMA orbit (1.59 / 6.99) used for the Jacobi radius; the AGAMA and galpy implementations of
McMillan 2017 differ by 4% in pericentre. **Pick: all three** — the class-1 spread is the
potential spread, 1.25–1.97 kpc in pericentre.

## Class 2 — today's orbit with dynamical friction integrated backwards

### What has to be assumed

Chandrasekhar friction scales with the satellite's bound mass, so the backward orbit depends on
**M(t)**, not on a single mass. Integrating back with a constant mass is wrong in the
direction that matters: a constant 10¹⁰ M☉ satellite is pumped to 40–140 kpc within 10 Gyr,
i.e. it "arrives" from outside the virial radius although the real object was already stripped
to a fraction of that mass for most of its orbit.

Assumptions adopted:

1. **Infall mass M_inf** from the progenitor literature: the NSC–host relation gives
   M★ ≈ 10⁹ M☉ (Pfeffer et al. 2021, MNRAS 500, 2514,
   https://doi.org/10.1093/mnras/staa3407); oMEGACat X gives a "GE-like" dwarf mass
   4.5 (−2.4, +5.7) × 10⁹ M☉ (Souza et al. 2026, arXiv:2603.23589) and favours Sequoia/Thamnos
   over GSE as debris. Halo masses M₂₀₀ for M★ = 10⁸–10⁹ at z ~ 2 are 10¹⁰–10¹¹ M☉; GSE proper
   is 2 × 10¹¹ M☉ (Naidu et al. 2021, arXiv:2103.03251). Scanned: 10¹⁰, 3 × 10¹⁰, 10¹¹, 2 × 10¹¹.
2. **Infall time** t_inf = 10 Gyr ago (z ≈ 2; GE merger 8–11 Gyr ago). Also scanned 8 Gyr.
3. **Stripping**: bound mass M(t) = M_inf exp[−(t_inf − t_lb)/τ], floored at the nucleus mass
   3.55 × 10⁶ M☉; τ scanned 0.5–3 Gyr. Half-mass radius r_h = 1 kpc (M_inf/10¹⁰)^{1/3}, fixed.
4. **Friction**: galpy's Chandrasekhar formula (σ = v_c/√2, ln Λ = ln[r / max(r_h, GM/v²)]),
   host = McMillan 2017 through AGAMA, leapfrog with 0.25 Myr steps (~1 s per orbit). The
   frictionless run reproduces class 1 (1.63 / 7.04 kpc).
5. **Consistency requirement** that selects (M_inf, τ): the backward orbit must be at the
   virial radius of the young Milky Way, **50–150 kpc, at t_inf**. Histories that do not reach
   50 kpc cannot have sunk to today's orbit by friction alone; histories beyond 150 kpc need
   more friction than the mass allows.

### Scan result: apocentre in the last Gyr before t_inf (kpc), t_inf = 10 Gyr

| M_inf \ τ [Gyr] | 0.5 | 0.75 | 1.0 | 1.25 | 1.5 | 2.0 | 2.5 | 3.0 |
|---|---|---|---|---|---|---|---|---|
| 10¹⁰ | 23 | 29 | 36 | 47 | 58 | 65 | 72 | 86 |
| 3 × 10¹⁰ | 46 | 46 | 57 | 85 | 87 | 99 | 120 | 166 |
| 10¹¹ | 78 | 111 | 127 | 148 | 139 | 154 | 307 | 367 |
| 2 × 10¹¹ | 53 | 93 | 125 | 226 | 233 | 229 | 228 | 497 |

(The table is noisy at the 20% level because the phase of the last apocentre matters.)
Consequences: a **Sequoia-scale 10¹⁰ M☉ halo reaches the virial radius only if stripping is slow
(τ ≳ 1.5 Gyr)**; a GSE-scale halo must be stripped fast (τ ≲ 1 Gyr) or it would have come from
far outside the virial radius. In every viable history the satellite's bound mass is below
10⁹ M☉ by 5 Gyr ago, and the orbit has essentially today's shape (apo 7–12 kpc) for the last
3–5 Gyr: **the present orbit is a good proxy for the last half of the history in all cases**.

### Picked orbits

| label | M_inf | τ | r at t_inf | 5–6 Gyr ago (peri–apo) | 2–3 Gyr ago |
|---|---|---|---|---|---|
| Sequoia-like | 10¹⁰ | 2.0 | 65 kpc | 2.9–17.9 | 2.0–9.2 |
| intermediate | 3 × 10¹⁰ | 1.25 | 85 kpc | 2.2–11.6 | 1.7–7.4 |
| GSE-like | 10¹¹ | 0.75 | 111 kpc | 1.7–7.5 | 1.6–7.1 |

## Class 3 — deposited by GSE, migrated inward by the bar (done properly)

*Code: `src/ocen_dm/tails/bar_migration.py` (`python -m ocen_dm.tails.bar_migration`); figures
`plots/bar_migration_elz.png`, `plots/bar_migration_class3_picks.png`; tables
`results/tails/bar_migration_grid.ecsv`, `results/tails/bar_migration_class3_picks.ecsv`.*

### Set-up reproduced

Dillamore, Zhang & Belokurov 2026 (arXiv:2606.12516), following the authors' own pipeline in
`github.com/adllmr/oCen_bar` (local clone `~/Work/Code/oCen_bar`, referenced not copied):

* **Potential**: Hunter et al. 2024 (A&A 692, A216) Milky Way with the Sormani et al. 2022 bar
  (M_bar = 1.83 × 10¹⁰ M☉), as AGAMA files. Split into the axisymmetric part plus the baryonic
  bar part, scaled and rotated, minus the axisymmetrised baryonic part (so the m = 0 mass is
  unchanged while the bar grows).
* **Bar history**: amplitude switch-on (Dehnen 2000 eq. 4) from t = 0 to t₁ = 1 Gyr; pattern
  speed constant at Ω₁ until t₁, smooth onset of deceleration to t₂ = 2 Gyr, then constant
  η = −Ω̇/Ω² = 0.003 until t_f = 8 Gyr (bar age). Ω₁ follows from the present-day Ω_b,0:
  Ω₁ = 45.1 for Ω_b,0 = 24 (the paper's "≈ 45"), 72 for 30, 110 for 35. Bar length scales as
  S = Ω_b,0/Ω_b(t) (≈ corotation), amplitude fixed.
* **ω Cen**: 10³ samples of (D = 5.43 ± 0.05 kpc, μ = (−3.257, −6.730) ± 0.025 mas/yr,
  v_los = 232.7 ± 0.21 km/s), astropy default frame (R₀ = 8.122, v☉ = (12.9, 245.6, 7.78)),
  bar angle 28°, integrated from t_f back to 0 with `agama.orbit`.
* **Test**: (E, L_z) at t = 0 in the axisymmetric potential against the Belokurov et al. 2023
  GSE contours (energy zero-point matched at R = 8.2 kpc). Success metric here: fraction of
  samples inside the outermost contour.

### Results

| Ω_b,0 [km/s/kpc] | 20 | 21 | 22 | 22.5 | 23 | 23.5 | 24 | 24.5 | 25 | 25.5 | 26 | 26.5 | ≥ 27 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fraction inside GSE at t = 0 | 0.44 | 0.38 | 0.58 | 0.85 | 0.01 | 0.34 | **0.89** | 0.81 | 0.79 | 0.54 | 0.18 | 0.02 | 0.00 |

* Reproduces the paper: overlap only for Ω_b,0 ≲ 26; nothing for the mainstream 30–40.
* Cross-check against the authors' stored grid (`oCen_bar/artifacts/03_orbit_grid.npz`, which
  used t₁ = 2, t₂ = 3): 0.87 (24), 0.69 (25), 0.08 (26), 0.09 (23), 0.80 (22.5) — same
  pattern, same jagged Ω dependence. The dip at 23 is robust to the random seed (0.013 / 0.028)
  and to the bar timings: it is resonant-phase sensitivity, not noise.
* **Migration happens late**: E and L_z of the successful samples are constant at GSE-like
  values from 8 to ~2.5 Gyr ago and move to today's (E, L_z) only in the last ~2.5 Gyr, as
  Ω_b falls from ~30 to 24 and the retrograde 1:1 resonance sweeps over ω Cen.
* **The bar changes today's orbit too.** In the barred Hunter24 potential at Ω_b,0 = 24
  (corotation ≈ 9.5 kpc, so ω Cen at 6.5 kpc is inside corotation) the present orbit has
  peri/apo = 0.81/9.3 kpc over the last Gyr, versus 1.56/7.02 in the axisymmetric version and
  1.57/7.04 in McMillan 2017. Class 1's "today's orbit" therefore depends on whether a slow
  bar is present — relevant to the disruption simulations of all classes.

### Picked orbits (Ω_b,0 = 24; samples at the 16/50/84th percentiles of E(t = 0) among those ending inside GSE)

| sample | E(t=0) [10⁵] | L_z(t=0) | pre-migration (7–8 Gyr ago) peri / apo / e / z_max | last Gyr peri / apo |
|---|---|---|---|---|
| 986 | −1.301 | +69 | 0.87 / 10.9 / 0.85 / 4.8 | 0.84 / 9.0 |
| 133 | −1.265 | +163 | 0.60 / 11.6 / 0.90 / 4.9 | 0.83 / 9.3 |
| 240 | −1.233 | +229 | 0.50 / 12.2 / 0.92 / 5.2 | 0.84 / 9.1 |

(L_z prograde-positive.) The pre-migration orbit is a **most-bound GSE-debris orbit**: apocentre
11–12 kpc, pericentre 0.5–0.9 kpc, slightly prograde to zero L_z — not retrograde. The earlier
static guess (L_z = −300) was wrong in sign and is kept only as `class3_gse_debris_orbit`.

### What is still not modelled

* Anything before the bar formed (> 8 Gyr ago): GSE's own infall and disruption (Naidu et al.
  2021 fiducial: M★ = 5 × 10⁸, M_DM = 2 × 10¹¹, from R_vir at z ≈ 2, circularity 0.5,
  retrograde) and where ω Cen's host sat in it. The two readings of "deposited by GSE" (ω Cen
  = GSE's nucleus; or a smaller host in group infall) both start from the picked debris orbit.
* Dynamical friction on the nucleus (negligible at 3.6 × 10⁶ M☉, as the paper notes) and on the
  host while it still existed (class 2 covers this).
* The bar-history caveat recorded in the oCen_bar journal: with fixed η = 0.003 a present-day
  Ω_b,0 = 35 implies an initial Ω₁ = 110 km/s/kpc, so the scan over Ω_b,0 is also a scan over
  implausible early histories; the low-Ω_b,0 requirement is tied to this slowdown law.

## What differs between the classes (for the simulations)

| | early orbit (≥ 5 Gyr ago) | late orbit | pericentre |
|---|---|---|---|
| class 1 | today's | today's | 1.3–2.0 kpc |
| class 2 | wider, apo 10–20 kpc at 5 Gyr, infall from R_vir | today's | 1.6–2.9 kpc early |
| class 3 | most-bound GSE debris, apo 11–12 kpc, L_z ≈ 0 to +230, until ~2.5 Gyr ago | today's, reached by bar migration in the last ~2.5 Gyr | 0.5–0.9 kpc |

## Caveats

* Chandrasekhar friction with σ = v_c/√2 and a fixed r_h is a 20–30% level approximation; the
  mass-history dependence is far larger than that.
* Only McMillan 2017 for class 2 (AGAMA host); class 3 uses Hunter et al. 2024 (the only
  potential with the fitted bar); class 1 carries the potential spread.
* The galpy `ChandrasekharDynamicalFrictionForce` (Python force, ~30 min per orbit) is kept as
  `class2_friction_backwards` for cross-checks; the fast leapfrog is the working tool.
* galpy checkout at `~/Work/src/galpy` is from 2024-03 and needs two import shims (astroquery
  metadata, scipy `vectorize1`), both in `progenitor_orbits._import_galpy_safely`.
