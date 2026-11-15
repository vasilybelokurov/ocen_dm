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
4. In parallel: shape measurement from star counts -- per-Gaussian q' (completeness, field background, trusted radii).
   Method to be agreed with the user before running. Known first pass (2026-11-13): eps 0.12-0.15 at 114-250",
   0.15-0.17 at 300-850", 0.10 at 850-1100", 0.03-0.06 beyond (uncertain: tides, field, crowding).
5. JAM with the measured shape: q' per Gaussian fixed inside ~1100"; outer Gaussians bracketed (inner value vs round);
   fit R <= 1100" (main) and all bins (check only). Quote dark mass inside the fitted radius (~30 pc).
6. Later, if still needed: halo-mass profile; null mocks (K1 truth; no-DM N-body snapshot for the outer bins);
   power test with injected halos (CSD3 for large batches; needs approval).
Dropped: Gaia mean radial motion subtraction (double counting); dropping outer bins as a halo test (superseded).
Debate findings (2026-11-15): the old halo gain came only from Gaia radial (dlnL +1.79), Gaia-T worse (-0.87); outer
T/R misfit (K1 1.22/1.42 vs data 1.10/1.18) is outside JAM's beta control at fixed q'; HST rotation-draw spread of the
gain s.d. 0.2 (+0.5 for the spherical-law field).

## Reporting rules
- Quote -lnpost differences to +-0.2 until restarts bring the run-to-run scatter below 0.05.
- No significance from chi2 p-values (halo shape undefined at zero mass; r_s at the box edge).
