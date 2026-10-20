# Fable review: omega Cen centre (oMEGACat VI kinematic centre; rigorous design) -- 2026-10-20

Same brief as docs/reviews/codex_centre_20261020.md. Fable ran short read-only checks on the oMEGACat catalogue.
Claude verified the depth asymmetry: faint fraction N(22-24)/N(<24) by 45-deg PA sector is 0.398-0.419 at 20-60" and
0.409-0.428 at 60-100" (symmetric), but 0.420 (PA 270-315, W/NW, shallowest) to 0.452 (PA 90-135) at 120-200", where the
bright (F625W<19) density is highest at PA 270-315 (0.644 vs 0.547-0.607 per arcsec^2).

## Findings (Fable)
- Reproduced the oMEGACat VI centre on their 2434 PM Voronoi bins: 2-D Gaussian + floor (-1.45+-0.74, +1.39+-0.68)". The
  offset comes from R > 100" (R < 100": +0.1+-1.8, +2.2+-1.7); without the floor the centre jumps to +4.3" N.
- Magnitude split in the same bins: m < 20.5 (-0.6+-1.0, +1.1+-1.0), 20.5-22 (-1.8+-1.4, +1.1+-1.3), m > 22 (-3.5+-1.4,
  +2.1+-1.2): the westward offset is carried by faint, large-error stars. Shape mismatch, ellipticity and equipartition bias
  the Gaussian centre by <= 0.25" (synthetic tests). A 13% differential PM-error scale between shallow and deep halves of the
  field is enough for a ~2" shift in the faint half. Net: statistical 0.7" correct; systematic ~2" (W).
- Our count offset is not independent: depth-asymmetric footprint beyond ~100" (SE deep, W/NW shallow). The odd (m=1) wedge
  asymmetry at 80-200" is anti-correlated between the depth-dominated 23-25 slice and every brighter slice (r = -0.54 to
  -0.88) and grows with brightness (apparent W shift -3" at 21-23, -12" at m < 17.5): a catalogue/exposure pattern, not a
  physical m=1 mode. Our W offset grows with aperture exactly where the depth asymmetry starts.
- Fisher precision per axis: counts (all < 23, R < 200") 0.28"; PM dispersion gradient 0.45" (R < 200"); vlos rotation 1.5-2.5"
  statistical and much worse systematically (a 1 km/s MUSE pointing zero point at 100" mimics ~25"). Core R < 20" alone:
  10" (counts), 12" (dispersion); the central sigma rise ~3"; fast stars ~1" (UNVERIFIED).
## Design (Fable)
Joint unbinned likelihood, shared centre: counts with a per-cell completeness map C(x) (from N(22-24)/N(19-22) in 20" cells)
and an m=1 term with a radial law (= a drifting centre), F625W 19-22 only; per-star dispersion with per-cell error scales
f_k and a magnitude term; vlos only as a sanity check. Mocks: real positions/errors, injected offsets 0-5" and depth-correlated
error-scale errors 1.0-1.3; require bias < 0.3" and correct coverage. Expected: ~0.4" statistical, 0.5-1" after
completeness marginalisation (global); the core cannot define its own centre.
## Which centre (Fable)
Dip tests: AvdM10, marginalise over a 1.5" prior. Inner Jeans/DF: AvdM10 fixed (1.5" moves sigma(10") by ~0.5% < 2-3% bin
errors), +-2" robustness rows. IMBH: centre as a free nuisance (AvdM10 +- 1.5" or the fast-star position) for M <~ 2e4.
Do not adopt the oMEGACat VI centre. Counter-argument: a mass-dependent real m=1 mode cannot be excluded; the joint fit with
the depth map settles it at ~0.5".
