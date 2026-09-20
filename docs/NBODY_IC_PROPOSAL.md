# Proposal: particle-mass scheme and initial conditions for the host-disruption runs (2026-09-20, draft for discussion)

Context: `docs/NBODY_CODE_ASSESSMENT.md` (code choice: gyrfalcON + AGAMA plug-in),
`docs/PROGENITOR_ORBITS.md` and `docs/progenitor_orbits.pdf` (the nine orbital histories).
Goal of the simulations: how much dark matter remains bound to omega Cen's nucleus after the
host dwarf is tidally disrupted, for each orbital class.

## The two timescales that dictate the scheme

1. Self-relaxation of a live N-body nucleus. omega Cen: N ~ 1e7 stars, t_relax ~ 10 Gyr. A live
   nucleus of 2e4 particles has t_relax ~ 0.1 N / ln N x t_cross ~ 0.1 x 2e4 / 10 x 0.3 Myr ~ 60 Myr:
   it core-collapses and evaporates spuriously within a Gyr; 1e6 particles still give only ~3 Gyr.
   Softening barely helps. A live nucleus over 10 Gyr is therefore wrong by construction.
2. Heating of a light species by a heavy one (Spitzer): rate ~ m_heavy rho_heavy, to be compared
   with self-relaxation ~ m_light rho_light. Inside the nucleus rho_star ~ 2500 Msun/pc^3 while an
   NFW halo of 1e10 Msun has rho_DM ~ 2 Msun/pc^3 at 10 pc; with m_DM = 1e3, m_star = 180 Msun the
   ratio is ~0.004 at the centre; but rho_star(10 pc) = 155 Msun/pc^3 for a Plummer a = 7 pc, so
   the ratio there is 0.07 for m_DM = 1e3 and 7 for m_DM = 1e5 (Codex correction): DM particles in
   the nucleus region must be <~ 1e4 Msun. Heating of the dwarf's stellar body (DM-dominated) is
   the "Ludlow" case, set by m_DM and rho_DM, not by the stellar particle mass.

## Proposed scheme (round one)

* Nucleus: a rigid Plummer potential (M = 3.6e6 Msun, a = 5.4 pc if 7 pc is the 3D half-mass
  radius, 7 pc if projected) that moves with a tracer particle -- implemented as an AGAMA
  potential with a time-dependent centre, NOT as a softened falcON particle: falcON's default
  kernel is P1 (not Plummer) and pair softening is (eps_i + eps_j)/2, so a 7 pc particle does not
  produce a 7 pc Plummer field (Codex). Optionally 1e4 massless tracers from the nucleus DF.
  Literature: Meadows et al. 2020 (arXiv:1910.11887) use single softened point masses (eps = 13 pc)
  for Fornax GCs; Boldrini, Mohayaee & Silk 2020 (arXiv:1909.07404) instead use live 1e4-particle
  GCs with equal particle masses -- the earlier draft cited them wrongly.
* Dwarf stars: 1e5 equal-mass particles (1e3-1e4 Msun for M* = 1e8-1e9), eps ~ 20-50 pc.
* Dark matter: radially graded particle masses (Zemp et al. 2008 multimass): ~1e3 Msun inside
  ~0.3 kpc rising to 1e5-1e6 Msun outside; ~1e6 particles total; softening graded with mass;
  implemented by sampling the AGAMA DF, thinning outer particles with probability p and
  reweighting by 1/p. Keep mass ratios between species sharing a region <~ 10-30 and heavy-particle
  softening >= the light species' inter-particle spacing.
* Equilibrium: AGAMA self-consistent model of the three components (nucleus potential included as a
  fixed Plummer), sampled into particles; short isolated run (~10 dynamical times) to verify.
* Host: AGAMA INI files already used in oCen_dm.tails (class 1: McMillan17, MWPotential2014 and
  Irrgang13 as in the orbit set; class 2: McMillan17; class 3: Hunter24 with the decelerating bar).
  Class-2 friction: a translating frame does NOT supply friction (the residual -a_df at the origin
  makes the satellite drift off the prescribed path; Codex). Instead run in the inertial frame and
  apply the precomputed Chandrasekhar acceleration a_df(t) as a `UniformAcceleration` time series to
  the satellite particles; it is also applied to escaped debris (small) and is not self-consistent
  with a changing bound mass -> iterate the orbit with the measured bound mass, or couple the drag
  online to the remnant's bulk motion.
* Convergence checks: inner DM particle mass x3; nucleus softening x2; one short live-nucleus run
  (N ~ 1e5, ~1 Gyr) to confirm the point-mass approximation for the bound-DM diagnostic.
* Diagnostic: mass of DM (and stars) bound to the nucleus within its instantaneous Jacobi radius
  vs time; snapshots every ~50 Myr.

## Points I am unsure about (for review)

* Is a rigid nucleus adequate when the question is DM bound *to* it? The Jacobi radius of the
  nucleus at pericentre is 33-59 pc at 0.4 kpc (class 3) and 69-107 pc at 1.6 kpc (classes 1-2)
  in Hunter24, i.e. 5-15 x the Plummer scale; the Plummer force differs from a point mass by 45% at
  10 pc, 3% at 50 pc.
* Is the multimass thinning safe for DM that later falls into the nucleus's region (heavy outer
  particles migrating inward through tidal shocks)?
* Resolution target: >= 1e3 effective survivors (N_eff = (sum m)^2 / sum m^2) for an integrated
  bound mass, 1e4 for a profile; a 1e5 Msun remnant then needs inner DM particles of ~100 Msun.
  Radius-only thinning is insufficient (heavy particles near apocentre plunge through the centre;
  Zemp et al. 2008 refine on pericentre with mass ratios ~2 per shell and a protection radius >> the
  inner shell) -> refine on pericentre/energy and monitor heavy-particle contamination.
* Alternatives not chosen: full live nucleus with equal masses (relaxation), two-stage resimulation
  with a live nucleus in the recorded tidal field (round two), pyfalcon/Python integrator for
  per-particle friction control.
