# Numerical tests of Section 11 (tidally truncated DM around a heavy nucleus)

*2026-09-21. Target: Section 11 of `dm_capture_constraints.tex` ("The equilibrium of tidally
truncated dark matter around a heavy nucleus"). Code: `src/ocen_dm/tails/truncated_equilibrium.py`;
tests `tests/test_truncated_equilibrium.py`.*

## Status after the 2026-09-21 audit

**Code update:** [DYNAMICAL_EXPERIMENTS.md](DYNAMICAL_EXPERIMENTS.md) describes the new
initial-condition and live/frozen evolution drivers. The instantaneous tidal diagnostic
has analytic circular-limit tests, and smoothly tapered joint equilibria have signed-DF
and density-recovery checks. These provide machinery for the experiments below; the
short software validation runs do not complete N2's hard energy-cut iteration or the
long tidal experiments N3–N5.

N1 is the only completed experiment in the full grid specified below. The energy-cut formula is verified for its stated
Kepler, isotropic, power-law assumptions; it is not a validated tidal-loss prescription.
At 20 pc the class-3 shock estimate accumulates to **0.112 over 66 identical passages**,
not order unity. At 35 pc it reaches 1.52. Whether outer losses subsequently deplete the
20-pc density is a question for N3/N4. The previous N2 setup also counted the cluster's
stellar mass twice; the revised component definitions below remove that ambiguity.

**Completed pilot comparison (21 September 2026).** Section 6 of
[the dynamical write-up](dynamical_batch.pdf) now evaluates a fixed-potential
energy cut in the actual cusp/core initial DFs and compares the scalar shock
estimate with the completed 88-Myr live/frozen pilot. This supplies partial
N2/N3/N5 evidence, without completing their proposed parameter grids or N2's
self-consistent iteration. At 20 pc the instantaneous energy cut predicts
10.18%/14.71% density loss, while the live simulations give 2.50%/4.66%.
The scalar pericentre energy estimate is much smaller than the full-orbit
mean; extended internal orbits, disc crossings and continued work on escapers
prevent interpreting that discrepancy as a fitted adiabatic correction.
See `results/plot_data/dynamical_analytics.json` and `bin/compare_dynamical_analytics.py`.

## What Section 11 estimates, and what is actually uncertain

The section uses three approximations:

1. **Energy truncation.** A tracer ρ ∝ r^−γ in the nucleus's Kepler potential has the isotropic
   DF f(E) ∝ (−E)^(γ−3/2); truncating it at Φ(r_J) leaves the fraction
   F_γ(x) = 1 − I_x(γ−½, 3/2) of the initial density at r = x r_J. At 20 pc:
   0.35 (γ = 1, r_J = 70 pc), 0.14 (γ = 1, r_J = 35 pc), 0.60 / 0.28 for γ = 3/2.
2. **Re-adjustment is small**, because a nucleus with total mass 3.55 × 10⁶ M☉ and approximate 3-D half-mass radius 7 pc
   can dominate low-DM models. Some grid choices have comparable DM mass; the small-response
   approximation must be checked model by model. The phase-mixing estimate is ~10⁸ yr.
3. **Shock heating** per pericentre is ΔE/|E| ≈ (Δv/σ)² (1 + ω²τ²)^−γ_ad with γ_ad ≈ 1.5, giving
   much stronger local heating in class 3. The current arithmetic does not establish
   order-unity heating at 20 pc within that history.

Ranked by how much they can move the answer, the uncertainties are:

| rank | uncertainty | size of the effect on ρ_DM(20 pc) |
|---|---|---|
| 1 | **the definition of r_J** — the simple r_peri[M/3M_host]^{1/3} gives 33 pc at the class-3 pericentre, a historical slope-corrected estimate was 59 pc, but has not been reproduced consistently | ~3× (through F_γ(r/r_J)) |
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

**Result — passed in the recorded N1 experiment.** Quadrature agrees with 1 − I_x(γ−½, 3/2) to 2 × 10⁻¹⁶ (γ = 1) and
6 × 10⁻¹⁷ (γ = 3/2) at x = 0.1, 0.3, 0.6. Monte Carlo with 4 × 10⁵ samples agrees within
1.1 σ at all six points. The committed unit test (`test_truncation_formula_against_monte_carlo`) uses 4 × 10⁴
samples per case and a four-sigma tolerance; it is a smaller regression check than the
recorded 4 × 10⁵-sample experiment. The sampler is local to the test, not yet a reusable
phase-space initial-condition generator.

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

**Set-up.** First recover the circular, spherical reference
`r_J³ = G M_sat / (Omega² - d²Phi/dR²)`, with `Omega² = G M_host(<R)/R³`.
Thus the enclosed-mass denominator is **3 − d ln M_host/d ln R**, not 2 minus that slope.
For a point-mass host the coefficient is 3; for a flat circular-speed curve it is 2.
This follows by differentiating `dPhi/dR = GM_host/R²`; see the
[galpy tidal-radius definition](https://docs.galpy.org/en/v1.5.0/reference/potentialrtide.html).

Then evaluate the tidal tensor `T_ij = -d²Phi/dx_i dx_j` along each orbit, add the
instantaneous centrifugal matrix for an explicitly defined rotating frame, and solve
`lambda_max r³ = G M_sat(<r)` where a disruptive eigenvalue is positive. Compare with the
existing coefficient-3 proxy and the circular spherical reference. An eccentric orbit in
a slowing bar has no conserved Jacobi integral: this is an instantaneous scale, not an
exact escape surface. Record the angular velocity and treatment of time-dependent terms.
Report the value at pericentre and its variation around the orbit; do it in McMillan 2017,
Hunter et al. 2024 axisymmetrised, and the barred Hunter potential at Ω_b = 24 (where the
non-axisymmetric term contributes directly).

**Verifies.** Rank-1 uncertainty. It is the input to every F_γ number in the section, and the historical estimates differ by 1.8× (33 vs 59 pc), enough to change the retention
factor substantially. Recalculate both under explicit conventions rather than treating
the unverified 59-pc estimate as a validated reference.

**Validation.** Recover the point-mass and flat-rotation circular limits and converge the
force derivatives. Report the convention spread as a systematic. A barred/axisymmetric
difference above 20% is a physical result to investigate, not a numerical test failure.

**Cost.** Minutes. No new machinery beyond `forceDeriv`.

---

## N2 — truncated equilibrium in the true potential (minutes, AGAMA, no dynamics)

**Set-up.** Define the baryonic mass once. Use either (a) a rigid Plummer approximation
with total nucleus mass 3.55 × 10⁶ M☉ and a = 5.4 pc for a 7-pc 3-D half-mass radius, or
(b) the full K1 stellar MGE + remnant + point-mass model, reconstructed from one saved
parameter vector and its distance. **Do not add a 2.9 × 10⁶ M☉ cluster MGE on top of a
3.55 × 10⁶ M☉ Plummer representing that same cluster.** An extra dwarf stellar envelope,
if included, needs its own mass and profile.

Add the DM model (NFW 10⁹/10¹⁰/10¹¹ M☉ with stated z = 2 concentrations, plus contracted
and cored variants). Compute the DM DF in the combined potential, apply the energy cut
at Phi(r_J), and iterate the DM potential to a stated convergence tolerance. Check
recovered density and moments before and after sampling. AGAMA clips negative DF values,
so a successful `QuasiSpherical` construction is not a positivity test. The gamma = 0
variant cannot use the gamma > 1/2 Kepler power-law formula.

Measure `M_DM(<r)/M_baryon(<r)` throughout the domain. If DM is not subdominant, bring
self-gravity into the experiment at this stage; do not assume all grid points qualify
for N3/N4's test-particle treatment and postpone the check until N5.

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
ΔE/|E| more than 3× the analytic value at 20 pc would reject the amplitude
of that local approximation in a matched shock experiment. It would not
by itself establish destruction of the 20-pc density. The completed pilot
shows why: a small fraction of extended or escaping internal orbits can
dominate the full-orbit cohort mean while most inner density remains.
Isolated-shock calibration must match the forcing history, internal orbital
frequencies, statistic and energy normalization before applying this criterion.

**Cost.** ~1 h including set-up. This is the highest value-per-hour experiment in the list.

---

## N4 — many passages: the number the WD comparison needs (hours)

**Updated scope (21 September 2026).** The
[long-term evolution plan](LIFETIME_EVOLUTION.md) specifies the density histories,
stellar-survival comparisons and convergence controls now required. The frozen,
rigid-nucleus experiments below establish conditional DM evolution. A present-day
survival constraint additionally requires a responsive nucleus. The implemented
class-3 history spans 8 natural units = 7.82234 Gyr; its original clock must be
preserved. The costs below are planning estimates, not measured production runtimes.

**Set-up.** As N3 but integrating the full histories: the class-1 orbit for 10 Gyr (~114
pericentres) and the class-3 debris orbit for 8 Gyr (~66), the latter also once in the
time-dependent barred potential of the companion note. Track the bound mass (energy in the
nucleus frame < 0) and ρ_DM at 3, 10, 20, 35, 50, 70 pc as functions of time, with output every
~50 Myr and densely around pericentre.

**Tests.** The combination of truncation and shocks:
whether the density at 20 pc is **truncation-limited** (classes 1–2 — the prediction is that it
settles at F_γ × the initial cusp and then stays put) or continues to decline as outer orbits are heated (a class-3 hypothesis). The local
20-pc heating estimate alone does not predict rapid erosion there.

**Scientific comparison.** Test whether the final class-1/2 density lies within a factor
2 of F_γ × ρ_initial, and measure rather than require a class-3 decline. Disagreement can
falsify the approximation. **Numerical validation:** converge with two particle numbers
(10⁴, 10⁵) and two time steps; separate numerical error from the physical discrepancy.

**Deliverable.** The corrected ρ_DM(20 pc) per orbital class, which is what enters the white-dwarf
heating comparison and the local DM mass fraction f_DM.

**Cost.** Hours (test particles in an analytic potential; no tree gravity).

---

## N5 — is the test-particle approximation valid? (one run)

**Set-up.** Repeat one N4 case with a live DM component (10⁵ particles, pyfalcon or gyrfalcON,
softening ≤ 1 pc) around the rigid nucleus potential.

**Verifies.** That neglecting DM self-gravity is safe — checked first from N2's enclosed-mass ratios. Some initial profiles exceed 10⁶ M☉
inside 70 pc, so this approximation is not valid across the grid by assumption.

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
