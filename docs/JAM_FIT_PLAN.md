# Flattened JAM fits of omega Cen: plan after review 2 (2026-11-15)

Context: bounded least-squares JAM fits (bin/fit_jam_halo.py --method lsq) with q' free. Real data: JAM K1 -lnpost
-178.62, JAM NFW -180.14 (Delta lnL 1.73, Delta chi2 3.46; halo at r_s = 1000 pc, the box edge). Review 2 (Codex, Fable)
found: the spherical-vs-JAM comparison was not on identical data; the halo gain is not significant by any calibrated
test; the two outermost Gaia radial bins (field fractions 0.85, 0.95) drive the residual misfit; several code issues.

## Done
- 2026-11-15: all rotation terms from ONE rotation field (the fitted field, fit_rotation_field.py, cylindrical MAP).
  Before, only HST used it; MUSE used the published v_rot^2/2 and Gaia-T the VB21 curve
  (tune_jam_to_data.Setup._corrected_data). Gaia-R unchanged (the field's mean radial motion is 0).

## Steps (in order)
1. Code fixes in fit_jam_halo.py / report_jamfit_lsq.py (no long runs):
   - failed models: every start (standard, given, screened) must give a valid model; failures counted per run;
   - start variety: starts picked best-first with a minimum separation in the normalised parameter box;
   - restarts: each least-squares run restarted from its end point until the gain is < 0.02 (fresh Jacobian),
     to cut the 0.05-0.2 run-to-run scatter and stalls;
   - outputs: no overwrite without --overwrite; save arguments, data values/errors/rotation terms, git commit;
   - mock report: label "at or below truth + 0.1"; enclosed mass of the flattened components computed properly.
   Test: short run on the q'0.90 + halo mock; must reach the truth score.
2. Gaia mean radial motion check: our mixture fit's mean_pmr in the outer rings; if non-zero, remove it like rotation.
3. Refit JAM K1 and NFW on the corrected data (~40 min); refit spherical K1/NFW on the same data (all stars per bin).
4. Drop the Gaia bins beyond 1500" (K1 and NFW, ~40 min): does the halo gain survive?
5. Halo-mass profile: ~8 fits with M_dm_100 fixed, incl. exactly zero (~2-3 h) -> upper limit.
6. Null distribution: ~30 noisy no-halo mocks (truth = K1 best fit, noise from the data errors), each fitted with K1 and
   NFW; Delta chi2 distribution. Real gain beyond the 95th percentile = hint, beyond the 99th = detection.
   ~8-10 h on the Mac (6 starts per fit) or on CSD3 CPUs (to decide).
Deferred, needs approval (model change): q' varying with radius (rounder outer model, measured ellipticity
0.03-0.06 beyond 1100"), only if steps 2-4 do not explain the outer radial misfit.

## Reporting rules
- Quote -lnpost differences to +-0.2 until restarts bring the run-to-run scatter below 0.05.
- No significance from chi2 p-values (halo shape undefined at zero mass; r_s at the box edge).
