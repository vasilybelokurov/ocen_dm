ask-codex: target=default base=ff92e1d effort=high diff=79723 bytes repo=/Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm git=1
ask-codex: review saved to /var/folders/jb/7f0bq_sx0435s1wvlt2yfkcr0000gn/T/tmp.jwZTFPyzk5/codex-review.md (419s, effort high)

## Findings

1. **Severity: high — Incorrect covariance signs corrupt the fitted streaming and anisotropy.**  
   **Location:** [outer_profile.py:542](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/outer_profile.py:542>) and line 544.  
   **Problem:** The `cad` terms in `A12` and `b2` have the wrong signs for the inverse covariance used by the likelihood. The mean update therefore does not maximise that likelihood.  
   **Failure scenario:** A 3,000-star synthetic sample with `(σR, σT)=(0.50, 0.25)` and nonzero streaming returned `(0.453, 0.314)`; correcting those signs in memory returned `(0.499, 0.253)`. In the actual Gaia 380–460″ bin, the correction changed the tangential mean from −0.221 to −0.258 mas/yr, affecting the streaming correction used for Jeans inference.  
   **Confidence:** high — derived the normal equations and reproduced both synthetic and real-data effects.

2. **Severity: high — Pristine is not an independent astrometric validation.**  
   **Location:** [docs/data_analysis.tex:263](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:263>).  
   **Problem:** Pristine supplies a different selection, but its proper motions are Gaia measurements; DR3 reuses EDR3 astrometry. The ratio uncertainty also omits covariance between the overlapping samples. [Pristine paper](https://arxiv.org/html/2502.01135v1), [Gaia release documentation](https://www.cosmos.esa.int/web/gaia/dr3).  
   **Failure scenario:** A shared astrometric or estimator bias can leave the ratio near unity while both profiles are wrong. I found shared source IDs in every comparison annulus, including 33 of the 52 selected Gaia stars in the innermost bin.  
   **Confidence:** high — checked catalogue provenance, joined source IDs, and inspected the independent-error propagation at `constraints.py:1395`. Positive covariance need not make that ratio error too small, but the quoted uncertainty is not a calibrated independent-validation precision.

3. **Severity: high — The fixed field template does not follow the science selection or its measurement errors.**  
   **Location:** [outer_gaia.py:74](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/outer_gaia.py:74>); [hst_profile.py:162](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/hst_profile.py:162>).  
   **Problem:** Both fits use an observed Gaia KDE with different magnitude, quality, parallax and error selections; only its normalisation can change. Neither selection matching nor a correction for the template’s different measurement-error convolution is performed.  
   **Failure scenario:** In the outermost Gaia bin, approximately 95% of selected stars are field stars, with an error ceiling near 0.092 mas/yr, whereas the template’s median error is 0.282 mas/yr. Restricting the template to that error ceiling changed the combined fitted dispersion by about 1.5%—well above the quoted approximately 0.2% error-model systematic.  
   **Confidence:** high for the mismatch and measured sensitivity; the actual bias requires a validated field model. Both diagnostic fits used the corrected mean-update signs.

4. **Severity: medium — Adding the systemic motion does not restore HST’s absolute proper motions.**  
   **Location:** [hst_profile.py:116](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/hst_profile.py:116>).  
   **Problem:** The locally corrected HST motions have lost position-dependent cluster streaming, including rotation. Adding only a constant systemic vector leaves those corrections missing when evaluating the absolute Gaia field template. The source paper explicitly states that the relative PMs contain no rotation signal. [oMEGACat VI](https://arxiv.org/html/2503.04903v1).  
   **Failure scenario:** Where the removed local rotation is approximately 0.25 mas/yr, an HST field star is scored against the Gaia KDE at a position displaced by that amount; membership responsibilities and dispersion can consequently change.  
   **Confidence:** high for the frame mismatch — checked the loader, `absolute_pm`, and catalogue description; its numerical impact remains unmeasured.

5. **Severity: medium — The reported statistical errors contain a numerical floor.**  
   **Location:** [outer_profile.py:584](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/outer_profile.py:584>).  
   **Problem:** Confidence limits advance in fixed `0.02 × sigma` steps and return the first crossing without interpolation or root-finding. This imposes a minimum approximately 2% component uncertainty, or 1.4% combined uncertainty near isotropy.  
   **Failure scenario:** Increasing a clean Gaussian sample to 100,000 stars should give approximately 0.16% combined statistical precision; this implementation still reports approximately 1.4%. Several large HST bins in the write-up sit at precisely that numerical floor.  
   **Confidence:** high — checked the interval algorithm and its propagation into the tables.

6. **Severity: medium — The error cut does not establish the claimed 2% bound.**  
   **Location:** [docs/data_analysis.tex:196](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:196>); [outer_gaia.py:71](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/outer_gaia.py:71>).  
   **Problem:** The code thresholds the arithmetic mean of the two errors against the **published** dispersion, whereas the bound requires control of error variance relative to the dispersion actually inferred. Neither individual component nor the mean squared error has the stated bound.  
   **Failure scenario:** For `σpublished=0.4`, errors `(0.02, 0.29)` pass because their mean is 0.155. Rescaling these errors by 1.13 can reduce a raw combined dispersion of 0.4 by approximately 3.7%, exceeding the claimed bound.  
   **Confidence:** high — checked the implemented selection and calculated the variance change.

7. **Severity: medium — The midpoint prescription is presented as a calibrated systematic without establishing a bracket.**  
   **Location:** [docs/data_analysis.tex:188](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:188>) and line 199.  
   **Problem:** The magnitude comparison does not establish that the true dispersion lies between the raw and inflated fits: it uses PM-derived membership cuts, unmatched magnitude populations and broad radial bins. Moreover, its faint sample has `G>19`, while the final selected Gaia sample contains no stars that faint.  
   **Failure scenario:** Membership selection suppresses the faint sample’s velocity wings, producing a lower dispersion even with correctly calibrated errors; this is then misidentified as over-inflation and used to justify the midpoint. A shared error-model defect can bias both endpoints in the same direction.  
   **Confidence:** high — inspected `magnitude_consistency()` and the selected catalogue. Half-separation is a sensitivity measure; its coverage as a systematic uncertainty has not been demonstrated.

8. **Severity: medium — The abstract upgrades a 13% comparison to 4% precision.**  
   **Location:** [docs/data_analysis.tex:24](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:24>).  
   **Problem:** The quality-selected overlap gives `Gaia/HST=1.028±0.132`, not `0±4%`. The approximately 4% uncertainty belongs to the route using rejected HST stars, an empirical correction, and substantial radius adjustment.  
   **Failure scenario:** A 10% instrument or population offset remains compatible with the clean overlap, yet the abstract would encourage treating it as excluded. The three routes also share stars and cannot be counted as independent confirmations.  
   **Confidence:** high — checked the overlap table and `hst_gaia_overlap_tests()`.

9. **Severity: medium — The full-covariance description is false for HST.**  
   **Location:** [docs/data_analysis.tex:76](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:76>).  
   **Problem:** HST explicitly sets every error correlation to zero; the supplied FITS table has no correlation column. The analysis therefore assumes diagonal HST errors rather than using measured full covariances.  
   **Failure scenario:** With nonzero correlations and uneven angular coverage, omitted terms `±2ρ εα εδ sinφ cosφ` are attributed to intrinsic radial/tangential variance, potentially creating or suppressing anisotropy.  
   **Confidence:** high for the documentation/implementation discrepancy — inspected the FITS schema and loaded arrays. Whether actual HST correlations produce a material bias is unverified.

10. **Severity: medium — The quoted tail position angle uses the wrong convention.**  
    **Location:** [docs/data_analysis.tex:300](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:300>).  
    **Problem:** `load_periphery()` measures angles from east towards north using `atan2(y,x)`, whereas astronomical position angle is measured from north through east. The computed 147.81° corresponds to astronomical PA **122.19°**, modulo 180°.  
    **Failure scenario:** Comparing the reported 148° with an orbit or published stream PA gives a spurious approximately 26° misalignment.  
    **Confidence:** high — recomputed the statistic and coordinate conversion from the 23 candidates.

11. **Severity: low — Constant systemic subtraction is confused with no subtraction.**  
    **Location:** [docs/data_analysis.tex:95](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:95>).  
    **Problem:** The approximately 5.3 mas/yr artefact arises from projecting the unsubtracted systemic motion. Subtracting its constant equatorial vector removes that leading term and leaves much smaller perspective/basis residuals.  
    **Failure scenario:** For a rigid cluster at 2400″, the repository’s projection gives approximately 0.057 mas/yr residual component RMS after constant subtraction, not 5.3.  
    **Confidence:** high — evaluated the repository’s exact projection on a synthetic annulus.

12. **Severity: low — A dispersion excess is incorrectly translated into the same percentage error underestimate.**  
    **Location:** [docs/data_analysis.tex:158](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/data_analysis.tex:158>).  
    **Problem:** A dispersion ratio of 1.077 does not imply that measurement errors are underestimated by approximately 8%; intrinsic and measurement **variances** add. Population differences and contamination also need exclusion before attributing the excess to errors.  
    **Failure scenario:** With intrinsic dispersion 0.55 and quoted error 0.10 mas/yr, that dispersion excess would require true errors near 0.242 mas/yr if error miscalibration alone caused it.  
    **Confidence:** high — checked the variance identity and the empirical correction’s interpretation.

## 1. Unverified assumptions the change relies on

- **Selection bias:** An error-only cut does not inherently truncate velocities or force the answer towards the published profile. Its unbiasedness requires the error/quality selection to be independent of kinematics at fixed radius and tracer population. Selecting brighter, more massive stars can lower the measured dispersion through partial equipartition, which is measured in ω Cen. [oMEGACat VI](https://arxiv.org/html/2503.04903v1).
- **Jacobi radii:** I could not reproduce 90/185 pc from repository code; they are hard-coded in `periphery.py`. The calculation needs the potential, orbital state, tidal derivative/tensor and angular-frequency convention. The separate `CompositeMassModel.jacobi_radius()` uses the point-mass factor `3Mgal`, which is not generally appropriate for an extended, flattened Galactic potential. [Tidal-radius formulation](https://arxiv.org/abs/1111.5013).
- **Rayleigh null:** Doubling the angles is appropriate for axial data, but uniform angular selection after the metallicity and PM cuts is unestablished. Uniform *radial* counts do not establish this.
- **Spectroscopy:** Applying no additional PM cut does not establish that target selection or the supplied membership flag was PM-independent.
- **Joint inference:** Shared Gaia stars/observations, shared field templates, spatial astrometric correlations and common error calibration are not represented by independent per-bin error bars.

## 2. Missing or weak tests

- Five targeted existing tests passed, but the mixture mocks use zero streaming and uniform angular coverage, missing the demonstrated sign error.
- No confidence-interval coverage or \(N^{-1/2}\) scaling test detects the numerical uncertainty floor.
- No end-to-end recovery tests vary field selection, HST reference-frame corrections, HST correlations, or template uncertainty.
- No selection robustness test varies the pilot dispersion, error threshold and magnitude/mass distribution independently.
- Several tests require the dispersion to decline or the surveys to agree. Those assertions can reject a real outer rise while accepting shared biases.
- No reproducible test establishes the Jacobi numbers or the astronomical PA convention.

## 3. What the author should verify before releasing

- Correct the likelihood algebra and interval calculation, then regenerate dispersions, anisotropy, streaming corrections and comparison uncertainties.
- Validate selection-matched field models and quantify HST frame/covariance sensitivity.
- Recalibrate error inflation on the population actually used; propagate shared nuisance parameters across bins.
- Recompute Pristine/Gaia comparisons with overlap-aware covariance and describe them as selection checks.
- Publish the Jacobi calculation and uncertainty; distinguish a density non-detection from evidence for a physical truncation.
- Calibrate the tail statistic against the survey selection and any searched radial/chemical/PM thresholds. The Rayleigh resultant already accounts for fitting the axis; the subsequent along-axis wedge comparison is selected using the same data.