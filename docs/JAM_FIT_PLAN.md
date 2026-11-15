# Flattened JAM fits of omega Cen: plan after review 2 (2026-11-15)

Context: bounded least-squares JAM fits (bin/fit_jam_halo.py --method lsq) with q' free. Real data: JAM K1 -lnpost
-178.62, JAM NFW -180.14 (Delta lnL 1.73, Delta chi2 3.46; halo at r_s = 1000 pc, the box edge). Review 2 (Codex, Fable)
found: the spherical-vs-JAM comparison was not on identical data; the halo gain is not significant by any calibrated
test; the two outermost Gaia radial bins (field fractions 0.85, 0.95) drive the residual misfit; several code issues.

## Done
- 2026-11-15: all rotation terms from ONE rotation field (the fitted field, fit_rotation_field.py, cylindrical MAP).
  Before, only HST used it; MUSE used the published v_rot^2/2 and Gaia-T the VB21 curve
  (tune_jam_to_data.Setup._corrected_data). Gaia-R unchanged (the field's mean radial motion is 0).

## Steps (revised 2026-11-15 after the critic/defender debate and discussion with the user)
1. DONE: fitter fixes (screened/separated starts, restarts, failure counts, safe outputs); flattened enclosed mass.
2. DONE: data fixes -- rotation-term error from the spread of the subtracted quantity (sqrt<vbar^2> for HST/MUSE,
   <vbar> for Gaia); Gaia ring terms weighted by the tracer surface density. fit_jam_halo --rmax (fit bins R <= rmax).
3. Base fits on identical data, one global q' (no shape constraint): JAM K1, JAM NFW, spherical K1, spherical NFW
   (all stars per bin for the rotation terms); plus JAM K1/NFW fitted to R <= 1100" only (like-for-like base for the
   shape-constrained fits). Report: data/model/residual plots, per-bin halo gain radial vs tangential per dataset.
4. Shape measurement from star counts (method agreed 2026-11-15, revised after Codex + Fable reviews):
   a. Refit the circular count MGE with a denser outer sigma grid (ratio ~1.4 beyond ~300", >= 1 Gaussian beyond
      772"); Gaia beyond ~850" from PM + parallax members (position-independent cut); check the outer slope
      (members ~ -5.6 at 1550-2200" vs current MGE -6.7).
   b. Ring table: m=1, m=2 (eps, PA), m=3 for HST (F625W < 19, < 18) and Gaia (all, PM members), with plots.
   c. Unbinned Poisson point-process fit of the star positions: q' a smooth function of log sigma (3-5 nodes), NOT
      limited to q' > cos i (consistency with i = 44 checked afterwards); Gaussian amplitudes refitted (flattening at
      fixed sigma changes the ring average, +7.4% at R = 3 sigma for q' 0.85); PA fixed (104.2) and free; centre and
      m=1 nuisance per dataset; Gaia field uniform + gradient (variant: quadratic). HST 60-250"; Gaia 300-2400".
      Beyond 2400" (Kuzma 2025 members: 12-55 per ring): PA / tail direction only.
   d. Mocks (>= 20 each; pass = within 1 sigma, mean bias < 0.01): round; flattened; a step in q'(R); wrong outer
      slope; 5" centre offset; magnitude-dependent m=1; curved field; PM-selection injection (rotation, perspective).
      Completeness maps from artificial stars only if bright cuts or mocks show a problem.
   e. Plot eps(R), PA(R) with literature (van de Ven et al. 2006 Table 2, verified from the paper: q' 1.000, 0.934,
      0.877, 0.848, 0.850, 0.836, 0.866, 0.926 at sigma' 0.153, 1.477, 2.525, 3.691, 5.219, 7.534, 11.07, 17.55
      arcmin; q'_j = observed flattening (Geyer et al. 1983) at R = sigma_j); decide the reliable radius.
   f. JAM inputs: q'_k per Gaussian; remnants q' from the same q'(log sigma) or round; fitter takes q' arrays.
   Noted: periphery code angle convention atan2(y,x) is not PA east of north; N-body orbit_today stopped at 1955.6 Myr
   but is labelled 2000 Myr.
5. JAM with the measured shape: measured q'_k at 60-1100"; inner q' free; outer Gaussians bracketed (inner value vs
   round); remnants: both variants -- same q' as the stars, and round; PA variant as in step 4 (JAM's axis is the
   kinematic one; a different photometric PA is a consistency check). Fit R <= 1100" (main) and all bins (check only).
   Quote dark mass inside the fitted radius (~30 pc).
6. Later, if still needed: halo-mass profile; null mocks (K1 truth; no-DM N-body snapshot for the outer bins);
   power test with injected halos (CSD3 for large batches; needs approval).
Dropped: Gaia mean radial motion subtraction (double counting); dropping outer bins as a halo test (superseded).
Debate findings (2026-11-15): the old halo gain came only from Gaia radial (dlnL +1.79), Gaia-T worse (-0.87); outer
T/R misfit (K1 1.22/1.42 vs data 1.10/1.18) is outside JAM's beta control at fixed q'; HST rotation-draw spread of the
gain s.d. 0.2 (+0.5 for the spherical-law field).

## Reporting rules
- Quote -lnpost differences to +-0.2 until restarts bring the run-to-run scatter below 0.05.
- No significance from chi2 p-values (halo shape undefined at zero mass; r_s at the box edge).
