# Long-term DM evolution subject to stellar-cluster survival

Research requirement recorded on 21 September 2026: measure density evolution
over the cluster's exposure to the Milky Way, and use the observed survival of
Omega Cen to constrain the allowed histories. This is the next-stage design;
the completed batch contains one 88.1-Myr tidal passage and isolation controls,
not measurements over the cluster's lifetime. The implementation and staged queue
are now active; see the [batch report](LIFETIME_BATCH_REPORT.md) for current scope,
tests, physical heating estimates and run states.

## The inference we want

Determine the present DM density, especially at 20 pc, among histories that
leave a stellar remnant consistent with today's cluster. A surviving compact
stellar core can coexist with substantial loss of its extended halo or dwarf
envelope. Survival therefore constrains the joint initial model and orbital
history; it does not require unchanged stellar mass or imply unchanged DM.

The few-per-cent first-passage changes establish an initial response. Neither
linear accumulation nor a fixed fractional loss per passage predicts the
long-term density: successive shocks act on a distribution already altered by
stripping, heating and re-equilibration. Measure whether the inner density
continues to decline, approaches a plateau, or changes regime.

## Duration and orbital history

Distinguish stellar age, time since accretion into the Milky Way, time since
loss of the dwarf envelope, and time on the present orbit. They are separate
quantities. Proposed reporting epochs are 0.1, 0.5, 1, 3, 5, approximately 8,
and 10 Gyr where a chosen history covers them; they are not measured ages.

| History | Proposed coverage | Interpretation |
|---|---|---|
| Static class 1 | Up to 10 Gyr | Controlled repeated forcing in the present host; a baseline for accumulated tides |
| Growing/slowing bar, class 3 | Full stored history, 7.82234 Gyr | Time-dependent orbital comparison with the original absolute bar clock |
| Accreted massive progenitor, class 2 | Bracket accretion/exposure histories, initially 4-10 Gyr | Requires forward orbital evolution coupled to changing bound mass and host assumptions |

The class-3 endpoint is **8 natural time units**, each 0.97779222 Gyr, rather
than exactly 8 Gyr. Do not extend the saved bar history beyond its endpoint.
For historical comparisons, rewind from the observed present phase-space
point and integrate forward to that epoch. Extending the current class-1
pilot from its nearest-apocentre start does not automatically end at today's
orbital phase. Actual shock and pericentre counts must be measured from each
track, rather than repeating the first 88-Myr interval.

The literature motivates multi-Gyr exposure without selecting one unique
history. Stream models give accretion lower limits of 7 or 4 Gyr depending on
the assumed initial stellar mass ([Youakim, Lind & Kushniruk 2023,
MNRAS 524, 2630-2650](https://doi.org/10.1093/mnras/stad1952)). Bar migration
offers an alternative to a fixed past orbit, conditional on the adopted bar
parameters ([Dillamore, Zhang & Belokurov 2026, MNRAS 551,
stag1449](https://doi.org/10.1093/mnras/stag1449)). These published records
were verified against ADS and arXiv. Neither supplies a unique exposure time
for the compact remnant. A historical interpretation over the full interval
must assess host growth and the progenitor's changing protection of its nucleus.

## Use the surviving stars as an observational constraint

The first-batch nucleus was an imposed Plummer potential with fixed mass and size.
It cannot expand, lose stars or disrupt. The live `stars` component in the
progenitor runs represents the extended stellar body, not a responsive compact
nucleus. Current runs therefore cannot pass or fail a stellar-cluster survival
test, regardless of their duration.

The new pipeline supports a live nucleus, sampled in equilibrium with the other
components, and removes its analytic contribution when its particles supply the
gravity. Diagnostics follow a measured stellar centre; the rigid calculations
remain controls. Both live nuclei use β = −0.1, which passes the signed-DF checks;
an isotropic Plummer nucleus is inadmissible in the cusped compact-halo model.
Vary initial stellar mass and scale as well as halo parameters; do not require
the initial nucleus to equal the present cluster or retain an arbitrary fixed
fraction of its stars.

At the observed endpoint, compare the surviving remnant with the stellar
mass/luminosity, projected surface-density profile, half-light radius, and
line-of-sight and proper-motion kinematics. Include rotation and shape as the
model permits, and compare extra-tidal material using the observational
selection. Determine tolerances from the data and model uncertainties rather
than imposing an arbitrary survival threshold. The current rigid-nucleus
parameters are simulation inputs, not the observational acceptance criteria.
If kinematic-fit posteriors guide these comparisons, avoid counting their
underlying observations twice in a joint likelihood.

## Numerical and physical checks over long durations

1. Extend the repaired-kernel timestep and particle-number comparisons through
   the full 88-Myr forcing, then through a 0.5-1-Gyr bridge with matched isolation
   and at least a second seed. Record the runtime and storage before fixing the
   full production grid. A short isolation check does not establish Gyr stability.
2. Carry long-duration isolation references alongside tidal runs. Compare live
   and frozen satellite potentials for selected cases, reassessing the frozen
   approximation as the mass distribution evolves. Plot numerical differences
   alongside the inferred tidal effect; do not linearly extrapolate energy drift.
3. Check central shell counts, softening and unequal particle masses. In live
   nucleus runs, demonstrate that massive simulation DM particles do not
   artificially heat the stars. The current progenitor core's 15 effective
   particles in the 20-pc shell are inadequate for this measurement.
4. Estimate physical stellar relaxation and gravitational heating of DM from
   the stellar/remnant mass spectrum over the adopted history, extending N6.
   A softened collisionless tree run does not automatically reproduce those
   processes. If important, add a calibrated collisional treatment and relevant
   stellar mass evolution. Gravitational heating of cluster DM by stars is a
   separate effect from Galactic shocks ([Bertone & Fairbairn 2008,
   Section III](https://doi.org/10.1103/PhysRevD.77.043515)); its size must be
   evaluated for Omega Cen rather than imported from M4.

## Required analysis and plots

| Product | What it will show |
|---|---|
| Density histories | Absolute DM shell density and ratios to both initial conditions and matched isolation at 3, 10, 20, 35, 50 and 70 pc; median trends over orbital phase, resolved shock excursions, and effective counts |
| Evolution with radius | Density profiles at the reporting epochs, with enclosed DM masses reported separately; all-particle and bound-proxy measurements distinguished |
| Stellar survival histories | Bound stellar mass, half-light radius, central profile and velocity dispersion through time, with present observations shown at the endpoint |
| Depletion versus stellar response | Present rho_DM(20 pc) versus final stellar mass and size; identify which initial models and histories reproduce the observed remnant |
| Conditional density summary | DM density at 20 pc across observationally compatible histories, keeping orbit/initial-model spread separate from seed and numerical uncertainty |
| Forcing and convergence | Pericentre/disc events, host tidal gradients, isolation drift, timestep/particle-number comparisons and the estimated collisional contribution |

Sparse saved snapshots should accompany denser scalar diagnostics and event
sampling; the present uniform full-snapshot output must be costed before
multi-Gyr production. No lifetime density curve should be labelled a result
until it is measured from an evolved run.

The immediate computational step is the full-passage convergence comparison
and the 0.5-1-Gyr bridge. Development of a responsive nucleus and the physical
relaxation estimate can proceed alongside that work. Fixed-nucleus lifetime
runs remain useful conditional DM experiments; the survival-selected result
requires the stellar response. The launched batch contains 12 full-passage
controls, 16 responsive-nucleus 500-Myr cases, and eight conditional long
frozen-potential comparisons. Four workers run concurrently. Each later stage
waits for the previous numerical checks; a failure stops advancement. A fully
live, observationally selected multi-Gyr grid awaits the measured stellar
response and assessment of physical relaxation.
