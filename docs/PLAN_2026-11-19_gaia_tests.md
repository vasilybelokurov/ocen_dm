# Plan: end-to-end validation of the Gaia EDR3 outer dispersion profile (2026-11-19)

Product under test: `ocen_pm_dispersion_edr3_ours` (src/ocen_dm/kinematics/outer_gaia.py): 9 rings, 300-2400",
per-ring 2D mixture (cluster Gaussian in R/T + empirical 2D field template in absolute PM), raw catalogue errors,
perspective and depth corrections, membership-weighted model nodes.

Motivation: the outermost Gaia-R point (0.212 +- 0.016 mas/yr, 95.5% field) sits ~2 sigma above every model and
carries most of the residual halo gain. The ring's mixture passes a posterior-predictive check
(plots/gaia_ring8_ppc.png), so the question is whether the estimator is unbiased and its error honest.

Already done (not repeated): unit recoveries with synthetic Gaussian fields (tests/test_outer_profile.py); sigma_R/T
correlation -0.08 over 60 synthetic realisations; selection-matched template (+1.7% at 2148"); parallax cut removed
from template (<= 0.33 sigma); raw vs eta error model (VB21 replication, 2-6%); Gaia completeness vs HST LP fields;
node weighting fix.

## Tests, ordered by expected impact on the outer rings

### T1. End-to-end injection with the REAL field (bias and error calibration)  [main test]
- Per ring: draw field stars from the template region (r > 1 deg) with the ring's field count, matched in G and
  per-star errors, placed at the ring's positions; add N_cluster synthetic members (the fitted count) with known
  sigma_R, sigma_T, rotation and the members' real error distribution; perturb by the errors; run the production
  `dispersion_2d` with the production template (which then also contains those field stars: test both with and
  without leave-out).
- Variants: truth = fitted values; truth with a sigma gradient across the ring (+-10%) and the measured rotation
  gradient; non-Gaussian cluster (Student-t nu = 8 or a radial-anisotropy mixture).
- 200 realisations per ring per variant.
- Pass: |bias| < 0.3 sigma_stat; pull s.d. in [0.9, 1.1]; coverage of the 68% interval 0.63-0.73.

### T2. Gaia error-model calibration (the largest known systematic: raw vs eta moves the outer rings 3-8%)
- Internal: in each ring split members by per-star error (quartiles) and by G; fit sigma per split. Correct
  errors give the same sigma in all splits; error underestimation by k makes sigma rise with e. Fit k(G) jointly.
- External: HST-Gaia star-by-star overlap (300-360", plots/hst_gaia_star_by_star.png): distribution of
  dmu / sqrt(e_G^2 + e_H^2) per G bin gives k directly.
- Pass: k consistent with 1 (raw adopted) or a measured k(G) replaces the VB21 eta (needs your approval: it changes
  the errors).

### T3. Field-model sensitivity
- Template variants: region (1-1.5 deg, 1.5-2 deg, four quadrants, avoiding the tidal-tail directions), KDE
  bandwidth 0.10 / 0.15 / 0.25 mas/yr, with and without G/error matching.
- Independent field model: a free 2D Gaussian mixture (2-3 components) fitted in the ring itself, no template.
- Simple cross-check: stars within |dmu| < 1 mas/yr, flat local field from the 1-2 mas/yr annulus.
- Posterior-predictive check (bin/gaia_ring_ppc.py) for every ring, with Poisson deviance.
- Report the spread of sigma per ring against sigma_stat.

### T4. Sample and astrometric systematics
- Bright-star frame offset: split G < 13 / >= 13 and apply the Cantat-Gaudin & Brandt 2021 correction
  (2021A&A...649A.124C, https://ui.adsabs.harvard.edu/abs/2021A%26A...649A.124C); a PM offset between the two
  groups inflates a combined sigma.
- Small-scale correlated systematics (Lindegren et al. 2021, 2021A&A...649A...2L,
  https://ui.adsabs.harvard.edu/abs/2021A%26A...649A...2L): estimated variance contribution at the ring scale.
- Quality cuts: RUWE < 1.2 / 1.4, ipd_frac_multi_peak, astrometric_excess_noise (unresolved binaries/blends).
- Pass: shifts < 0.5 sigma_stat.

### T5. Geometry and physics corrections
- Position-angle jackknife (4 quadrants; along vs across the tidal-tail axis): escapers/tails in the 50-63 pc ring.
- Ring edges: halve the outer rings, shift edges by half a ring; check sigma vs radius is smooth.
- Rotation: per-star mean from the rotation-field model vs one mean per ring (rotation gradient inflates sigma_T).
- Perspective and depth: vary distance (5.24-5.43 kpc) and systemic PM within errors.

### T6. External comparison
- Extend the VB21 replication to the current pipeline (we are -5% vs VB21 at 2122" with eta).

## Outputs
- results/gaia_tests/*.json, plots per test, a summary table: ring x {bias, pull width, error-model shift, field
  spread, sample shift, PA jackknife}.
- Decision rule: any shift > 0.5 sigma_stat in a ring is reported with its sign; the product is changed only for a
  bias (T1) or a calibrated error model (T2), and only with your approval. No added error floors.

## Cost
- T1: ~5 s per fit; 9 rings x 3 variants x 200 = 5400 fits ~ 7.5 core-h -> ~1 h on 8 cores (~1 GB each).
- T2-T6: minutes each. Needs one WSDB query for RUWE/ipd columns (T4).
- My time: roughly 2-3 h including plots.
