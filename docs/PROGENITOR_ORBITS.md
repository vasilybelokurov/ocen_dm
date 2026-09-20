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

## Class 3 — deposited by GSE, migrated inward by the bar

### What the literature gives

* Dillamore, Zhang & Belokurov 2026 (arXiv:2606.12516): modified Hunter et al. 2024 barred
  potential; bar grows 0–1 Gyr, decelerates with η = −Ω̇/Ω² = 0.003 from ≈ 45 km/s/kpc;
  ω Cen is caught by the **retrograde 1:1 resonance** (Ω_φ − Ω_b + Ω_r = 0), which scatters
  orbits to **lower E and more retrograde L_z**; ω Cen overlaps the GSE debris only if today's
  Ω_b ≲ 26 km/s/kpc. Back-integration of 10³ samples of the observed phase-space point.
  Friction is irrelevant at the nucleus mass (L̇_z ≈ 9.6 kpc km/s/Gyr).
* GSE debris (Belokurov et al. 2023, MNRAS 518, 6200, arXiv:2208.11135): |L_z| < 700 kpc km/s,
  chevrons with apocentres 11.5, 15.5, 21, 23, 25 kpc; most-bound debris E ≈ −1.4 × 10⁵ in
  their potential. In McMillan 2017 a near-radial orbit with these apocentres has
  E = Φ(r_apo) = −1.65, −1.50, −1.34, −1.26 × 10⁵ km² s⁻² (11.5, 15.5, 21, 25 kpc).
* GSE's own infall (Naidu et al. 2021 fiducial): M★ = 5 × 10⁸, M_DM = 2 × 10¹¹, c = 4,
  from R_vir at z ≈ 2, circularity 0.5, inclination 15°, retrograde; MW at z = 2:
  M₂₀₀ = 5 × 10¹¹, disc 6 × 10⁹, bulge 1.4 × 10¹⁰.

### Two readings of "deposited by GSE"

(a) **ω Cen is GSE's nuclear cluster** (Pfeffer et al. 2021). Then the disruption simulation is
GSE's own (Naidu et al. 2021 set-up above) and the dwarf's initial orbit is GSE's infall orbit;
the nucleus ends on a GSE-debris orbit, which the bar then drags inward.
(b) **ω Cen's host was a smaller dwarf that fell in with GSE** (group infall; the oMEGACat X
reading with Sequoia/Thamnos). Then the host disrupts on a GSE-debris orbit.

In both readings the nucleus spends the time between GSE's disruption (~9–10 Gyr ago) and the
bar's deceleration to Ω_b ≈ 25 (late) on a GSE-debris orbit. That orbit is the class-3 initial
condition for the host's tidal history. **Bar migration itself is not modelled here**: it needs
the Hunter et al. 2024 barred potential and a low pattern speed, and it moves the nucleus, not
the (already dispersed) host.

### Picked orbits (static McMillan 2017; start at apocentre, inclination 60°)

| label | apo | L_z | peri | e | z_max | E [10⁵ km²/s²] |
|---|---|---|---|---|---|---|
| GSE debris, most bound | 11.5 | −300 | 0.66 | 0.89 | 1.5 | −1.646 |
| GSE debris, middle | 15.5 | −300 | 0.62 | 0.92 | 1.9 | −1.492 |
| GSE debris, outer | 21 | −300 | 0.58 | 0.95 | 2.5 | −1.342 |

L_z = −300 (less retrograde than today's −529, as the resonance makes L_z more retrograde).
The pericentre is set by L_z in the flattened potential and hardly depends on the inclination
(0.6–0.8 kpc for 45–70°); these orbits are much more plunging than class 1 or 2, which is the
physical point of class 3: **stronger tides at pericentre, weaker at apocentre**.

## What differs between the classes (for the simulations)

| | early orbit (≥ 5 Gyr ago) | late orbit | pericentre |
|---|---|---|---|
| class 1 | today's | today's | 1.3–2.0 kpc |
| class 2 | wider, apo 10–20 kpc at 5 Gyr, infall from R_vir | today's | 1.6–2.9 kpc early |
| class 3 | GSE debris, apo 11–21 kpc | today's after bar migration | 0.6 kpc |

## Caveats

* Chandrasekhar friction with σ = v_c/√2 and a fixed r_h is a 20–30% level approximation; the
  mass-history dependence is far larger than that.
* Only McMillan 2017 for classes 2 and 3 (AGAMA host); class 1 carries the potential spread.
* The galpy `ChandrasekharDynamicalFrictionForce` (Python force, ~30 min per orbit) is kept as
  `class2_friction_backwards` for cross-checks; the fast leapfrog is the working tool.
* galpy checkout at `~/Work/src/galpy` is from 2024-03 and needs two import shims (astroquery
  metadata, scipy `vectorize1`), both in `progenitor_orbits._import_galpy_safely`.
