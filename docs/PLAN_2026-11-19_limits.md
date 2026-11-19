# Plan 2026-11-19: inner-beta fix, null calibration, DM upper limit (after adversarial review 4)

Goal: like-for-like spherical vs JAM models that obey the central BH anisotropy limit, then a calibrated 95% upper
limit on the dark mass inside the data radius. Authorised by the user to run without further approval.

## Steps, tests, decision rules (fixed before running)
1. Inner node. beta(2") = beta0_max = -0.5 fixed; b10 prior [-1, 0]; other nodes unchanged ([-1, 0.95]).
   Spherical: family option `inner_node_arcsec`; JAM: `--binner 2`. JAM A (turnover beta, beta_0 <= -0.5, r_beta
   ~2.5") already satisfies it and is not refitted.
   Tests: (a) inner_node_arcsec=None reproduces the old likelihood bit-for-bit at the saved ML point; (b) ln g exact
   vs quadrature with the inner node; (c) prior bounds; (d) config round trip.
   Runs: spherical local refits (K1, NFW), UltraNest K1 + NFW (3 ranks each); JAM C6 K1 + NFW, 18 starts, 6 workers.
   DF check: An & Evans / Ciotti-Morganti necessary condition gamma_tracer(r) >= 2 beta(r) (and gamma >= beta + 1/2
   where the BH dominates the enclosed mass) on a radial grid at each optimum.
   Rule: spherical is co-equal with JAM only if (DF condition holds) and chi2 <= chi2(JAM best) + 15 and
   |D - 5.43| < 2 sigma (0.10 kpc). Otherwise spherical stays a cross-check.
2. Bookkeeping: report data lnL, distance-prior term and lnpost separately; polish nested ML points locally.
3. Null calibration (spherical): 200 parametric mocks from the polished K1 optimum (Gaussian noise with the data
   errors, split-normal), each fitted locally with K1 and NFW (warm starts from the truth / the real-data optima);
   q95, q99 of 2 dlnL. Injections: NFW (r_s 300 pc) with M_dm(<100 pc) = 1e6 and 2e6 Msun, 20 mocks each.
   Rule: a halo is claimed only if the real 2 dlnL > q99 and injections are recovered.
4. Profile limit on M_dm(<50 pc) (50 pc ~ 1950" at 5.3 kpc, the outer data edge): spherical NFW on r_s
   {30, 100, 300, 1000} pc x M_dm grid, local optimisation of all other parameters incl. D; threshold = calibrated
   q95 of the null (one-sided; nominal 2.71). JAM C6 NFW on the same r_s grid with warm starts (fewer M grid points).
   Report the largest limit over r_s.
5. Rotation-error covariance: replace the independent per-bin rotation-error term by the covariance of the
   rotation-subtracted T (and R) profiles over the rotation-field draws; refit spherical K1/NFW locally. Report both
   if the halo gain moves by > 0.5 lnL or the limit by > 20%.
6. Sensitivity (local refits): M_rem/M_star HalfNormal(0.10) cap 0.3; M_star floor raised to a literature
   population-synthesis M/L_V of living stars; inclination 50 deg on JAM C6. N-body: tides change outer sigma_r < 1%.

Resources: <= 10 cores, < 24 GB; OMP/VECLIB threads 1.
