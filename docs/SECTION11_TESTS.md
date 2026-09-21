# Numerical tests of Section 11 (tidally truncated DM around a heavy nucleus)

*2026-09-21. Target: Section 11 of `dm_capture_constraints.tex` ("The equilibrium of tidally
truncated dark matter around a heavy nucleus"). Code: `src/ocen_dm/tails/truncated_equilibrium.py`;
tests `tests/test_truncated_equilibrium.py`.*

## What Section 11 claims, and what is actually uncertain

The section makes three quantitative claims:

1. **Energy truncation.** A tracer ρ ∝ r^−γ in the nucleus's Kepler potential has the isotropic
   DF f(E) ∝ (−E)^(γ−3/2); truncating it at Φ(r_J) leaves the fraction
   F_γ(x) = 1 − I_x(γ−½, 3/2) of the initial density at r = x r_J. At 20 pc:
   0.35 (γ = 1, r_J = 70 pc), 0.14 (γ = 1, r_J = 35 pc), 0.60 / 0.28 for γ = 3/2.
2. **Re-adjustment is small**, because the nucleus (3.55 × 10⁶ M☉ inside ~7 pc) dominates the
   potential over the ≲ 10⁶ M☉ of DM inside 70 pc, and phase mixing takes ~10⁸ yr.
3. **Shock heating** per pericentre is ΔE/|E| ≈ (Δv/σ)² (1 + ω²τ²)^−γ_ad with γ_ad ≈ 1.5, giving
   20 pc adiabatically protected for classes 1–2 but eroded for class 3.

Ranked by how much they can move the answer, the uncertainties are:

| rank | uncertainty | size of the effect on ρ_DM(20 pc) |
|---|---|---|
| 1 | **the definition of r_J** — the simple r_peri[M/3M_host]^{1/3} gives 33 pc at the class-3 pericentre, the form with the d ln M/d ln r term gives 59 pc | ~3× (through F_γ(r/r_J)) |
| 2 | **γ_ad and the impulsive approximation** at ωτ ~ 1 | factor ~2 in the number of passages to erosion |
| 3 | **departures from Kepler** (Plummer core inside ~10 pc, stars, DM self-gravity near r_J) | tens of per cent |
| 4 | **anisotropy of the initial cusp** (radial loses more, tangential less) | tens of per cent |
| 5 | the closed-form algebra itself | **none — verified, see N1** |

The experiments below are ordered by cost. N1 is done; N0 + N2 + N3 is one working day and
covers ranks 1–3; N4 produces the number the white-dwarf comparison needs.

---

## N1 — the truncation formula (DONE, minutes)

**Set-up.** Two independent checks with G = M = r_J = 1.
(a) Direct quadrature of the physical Eddington integral
ρ_t(r)/ρ(r) = ∫_{Φ(r)}^{Φ(r_J)} (−E)^(γ−3/2) √(E−Φ(r)) dE ÷ ∫_{Φ(r)}^{0} (…) dE,
without the substitution that produces the beta function.
(b) Monte Carlo: sample p(v) ∝ v² (−E)^(γ−3/2) at fixed r by numerical inverse-CDF, cut
E > Φ(r_J), measure the surviving fraction. Uses no beta-function identity.

**Result — passed.** Quadrature agrees with 1 − I_x(γ−½, 3/2) to 2 × 10⁻¹⁶ (γ = 1) and
6 × 10⁻¹⁷ (γ = 3/2) at x = 0.1, 0.3, 0.6. Monte Carlo with 4 × 10⁵ samples agrees within
1.1 σ at all six points. Now a unit test (`test_truncation_formula_against_monte_carlo`).

**Two sampler bugs found on the way**, both of which would have propagated into N2–N4 because
they reuse this machinery:
* a rejection-sampling envelope that is not an upper bound when γ < 3/2 (the weight
  v²(−E)^(γ−3/2) diverges as v → v_esc): it under-sampled exactly the loosely bound particles
  that truncation removes, and inflated the retained fraction from 0.60 to 0.87 — a 24 σ error
  that looked like a failure of the formula;
* trapezoidal CDF bias at the integrable endpoint singularity (0.6 %, 5 σ at 2 × 10⁵ samples),
  removed by sampling in θ with v = v_esc sin θ.

---

## N0 — what is r_J, really? (minutes, no dynamics)

**Set-up.** Along each of the nine orbits, evaluate the tidal tensor T_ij = −∂²Φ/∂x_i∂x_j from
`agama.Potential.forceDeriv` in the frame co-rotating with the orbit, add the centrifugal term,
diagonalise, and solve λ_max r³ = G M_nuc for the tidal radius. Compare with
(i) r_peri [M_nuc/3M_host(<r_peri)]^{1/3} and (ii) the same with the (2 − d ln M/d ln r) factor.
Report the value at pericentre and its variation around the orbit; do it in McMillan 2017,
Hunter et al. 2024 axisymmetrised, and the barred Hunter potential at Ω_b = 24 (where the
non-axisymmetric term contributes directly).

**Verifies.** Rank-1 uncertainty. It is the input to every F_γ number in the section, and the two
conventions currently in use differ by 1.8× at the class-3 pericentre (33 vs 59 pc), which is
~3× in ρ_DM(20 pc).

**Pass criterion.** A single stated convention with the spread between conventions quoted as a
systematic; the barred-potential value within 20 % of the axisymmetric one (if not, class 3
needs its own r_J).

**Cost.** Minutes. No new machinery beyond `forceDeriv`.

---

## N2 — truncated equilibrium in the true potential (minutes, AGAMA, no dynamics)

**Set-up.** Build the total potential: nucleus Plummer (3.55 × 10⁶ M☉, a = 5.4 pc if 7 pc is the
3-D half-mass radius), stars (the project's MGE light model scaled to the rung-0 K1
M★ = 2.9 × 10⁶ M☉), DM cusp (NFW 10⁹/10¹⁰/10¹¹ M☉ at z = 2 concentrations, plus c = 12 and a
contracted γ = 3/2 and a cored γ = 0 variant). Obtain the DM distribution function by Eddington
inversion **in the total potential** (`agama.DistributionFunction(type='QuasiSpherical')`), set
f(E) = 0 above Φ(r_J) for the r_J values from N0, recompute ρ_DM(r), and iterate the potential
once so the DM's own (small) contribution is consistent.

**Variants.** (a) apocentre criterion r_apo(E, L) > r_J instead of the energy cut — these differ
because a high-L orbit with E above Φ(r_J) never reaches r_J; the two bracket the truth.
(b) radially and tangentially anisotropic initial DFs (β = ±0.3), which is rank-4.

**Verifies.** Rank 3 and 4: how much the Plummer core (inside ~10 pc), the stellar mass and the
DM self-gravity break the Kepler idealisation behind F_γ.

**Pass criterion / deliverable.** ρ_t(r)/ρ(r) at 3, 10, 20 pc tabulated against F_γ; the
difference quoted as the "non-Kepler correction". **This output replaces the e^{−r/r_J} stand-in
in the density figure of Section 9** — the one number in the note that is currently a placeholder.

**Cost.** Minutes per model, ~1 h for the whole grid.

---

## N3 — one pericentre passage, test particles (≈ 1 h)

**Set-up.** 10⁵ massless particles sampled from the N2 equilibrium, orbiting an analytic nucleus
whose centre follows the precomputed class-1 (r_peri = 1.57 kpc, v_peri = 387 km/s) or class-3
(0.42 kpc, 507 km/s) orbit; integrate with `agama.orbit` in the inertial frame, the potential
being host + Plummer nucleus with a time-dependent `center` read from the orbit track. Run from
apocentre to apocentre across one pericentre. Measure ΔE in the nucleus frame per particle,
binned by initial radius (10, 20, 35, 50, 70 pc).

**Verifies.** Rank 2 — the entire shock table of Section 11, including the adiabatic correction,
and it **fits** γ_ad instead of assuming 1.5. No self-gravity, no softening, no relaxation, so
the result is clean.

**Pass criterion.** The impulsive scaling ΔE ∝ r² recovered where ωτ ≪ 1; the measured
suppression at ωτ ≳ 1 within a factor 2 of (1 + ω²τ²)^−γ_ad; γ_ad fitted to ±0.3. A measured
ΔE/|E| more than 3× the analytic value at 20 pc would overturn the "20 pc is protected"
conclusion for classes 1–2.

**Cost.** ~1 h including set-up. This is the highest value-per-hour experiment in the list.

---

## N4 — many passages: the number the WD comparison needs (hours)

**Set-up.** As N3 but integrating the full histories: the class-1 orbit for 10 Gyr (~114
pericentres) and the class-3 debris orbit for 8 Gyr (~66), the latter also once in the
time-dependent barred potential of the companion note. Track the bound mass (energy in the
nucleus frame < 0) and ρ_DM at 3, 10, 20, 35, 50, 70 pc as functions of time, with output every
~50 Myr and densely around pericentre.

**Verifies.** The combination of truncation and shocks, i.e. the central claim of the section:
whether the density at 20 pc is **truncation-limited** (classes 1–2 — the prediction is that it
settles at F_γ × the initial cusp and then stays put) or **shock-limited** (class 3 — the
prediction is continued erosion, with the effective boundary migrating inward from r_J).

**Pass criterion.** ρ_DM(20 pc, t_final) within a factor 2 of F_γ × ρ_NFW for classes 1–2; a
measurable decline for class 3. Convergence: two particle numbers (10⁴, 10⁵), two time steps.

**Deliverable.** The corrected ρ_DM(20 pc) per orbital class, which is what enters the white-dwarf
heating comparison and the local DM mass fraction f_DM.

**Cost.** Hours (test particles in an analytic potential; no tree gravity).

---

## N5 — is the test-particle approximation valid? (one run)

**Set-up.** Repeat one N4 case with a live DM component (10⁵ particles, pyfalcon or gyrfalcON,
softening ≤ 1 pc) around the rigid nucleus potential.

**Verifies.** That neglecting DM self-gravity is safe — justified a priori because the DM inside
70 pc (≲ 10⁶ M☉) is subdominant to the nucleus (3.55 × 10⁶ M☉), but the re-adjustment claim
(item 2 of Section 11) has not been tested.

**Pass criterion.** Bound mass and ρ(20 pc) within 20 % of the test-particle run. If it passes,
N3–N4 can be used for the full grid of halo masses and concentrations at negligible cost.

---

## N6 — stellar two-body heating (estimate only, no simulation)

Diffusion coefficients from the Spitzer formulae with the K1 stellar density and remnant
population, integrated over 10 Gyr, to check that the DM inside ~5 pc is not heated out by the
stars. An N-body treatment is out of scope; the estimate decides whether it needs to be in
scope.

---

## Summary of the programme

| | what it settles | cost | status |
|---|---|---|---|
| N1 | the closed-form truncation factor | minutes | **done, passed** |
| N0 | the definition of r_J (largest single uncertainty) | minutes | to do |
| N2 | non-Kepler corrections; replaces the e^{−r/r_J} model in the figure | ~1 h | to do |
| N3 | the shock table; calibrates γ_ad | ~1 h | to do |
| N4 | ρ_DM(20 pc) per orbital class — the number the WD comparison needs | hours | to do |
| N5 | validity of test particles | one run | to do |
| N6 | stellar heating inside 5 pc | estimate | to do |

Minimum useful set: **N0 + N2 + N3** (one working day), which fixes the r_J convention, replaces
the placeholder density model with a computed one, and calibrates the shock formula. N4 then
produces the final numbers.
