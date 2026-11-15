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
4. In parallel: shape measurement from star counts (method agreed with the user 2026-11-15).
   Known first pass (2026-11-13): eps 0.12-0.15 at 114-250", 0.15-0.17 at 300-850", 0.10 at 850-1100", 0.03-0.06
   beyond (uncertain: tides, field, crowding). Method:
   - unbinned Poisson likelihood of the star positions over each survey's footprint; count-MGE Gaussian widths fixed,
     one q'_k per Gaussian free (q'_k > cos i = 0.72 for i = 44 deg); no mge_fit_sectors (fails on count images);
   - PA: both variants -- fixed at the kinematic 104.2 deg, and free;
   - data: HST F625W < 19 at 60-250", Gaia G < 17 from 300" out to the edge of the member sample (catalogue reaches
     ~2400"), no quality-flag or membership cuts; Gaia field = uniform + linear gradient. Beyond ~2400": Gaia DR3 from
     WSDB in a wider field (PM + CMD selected, selection must not depend on position angle) and/or the Kuzma et al.
     2025 Pristine periphery members (data/processed/tails/kuzma2025_periphery.ecsv, to 5.1 deg; check that its
     membership has no circular spatial prior). Measure as far out as the
     members go; decide AFTER seeing the results (null/recovery mocks, bright cuts, field-model variants, PA(R))
     out to which radius the shape is reliable and equilibrium-like (tides) (user, 2026-11-15);
   - centre (< ~60"): excluded, no shape constraint (no m=2 information there; centroid/crowding problems) -- the
     inner Gaussians' q' is left free in the JAM fit (one shared value), not measured;
   - outer Gaussians: beyond the radius judged reliable (to be set from the measurement; first-pass guess ~1100") --
     bracketed: last reliable value vs round;
   - crowding check: bright cuts (F625W < 18, G < 16) must give the same q'; else adopt the bright cuts;
   - tests: round mocks through the real footprints/field/completeness return round; flattened mocks recover q' to
     0.01-0.02; literature cross-check (van de Ven et al. 2006, https://ui.adsabs.harvard.edu/abs/2006A&A...445..513V;
     values read from the paper).
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
