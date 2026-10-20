# Codex review: omega Cen centre (oMEGACat VI kinematic centre; rigorous design) -- 2026-10-20

Claude verified the code findings (symmetry prefilter at 60", b at bound 5 for 120" fits, no convergence check) and that the brief misstated the symmetry results (21-23: +4.1, -3.8).

## Brief
# Determining the centre of omega Cen (NGC 5139): critique a published kinematic centre and design a rigorous measurement

## Context
We model omega Cen's inner kinematics (Jeans and action-DF models, IMBH question). The adopted centre is Anderson & van der
Marel 2010 (AvdM10, 2010ApJ...710.1032A): RA 13:26:47.24, Dec -47:28:46.45, from HST ACS star counts (isodensity contours and
"pie slices"), confirmed within 1-2" by an HST PM kinematical centre and 2MASS unresolved light; quoted accuracy ~1".
Goldsbury+2010 (2010AJ....140.1830G, ACS GC survey ellipse fits within 2'): +0.41" E, +0.35" N of AvdM10 (quoted 0.1").
Haeberle+2024 (Nature) IMBH position from fast-star accelerations: 0.77" NE of AvdM10.

## The published kinematic centre (oMEGACat VI, Haeberle+2025, 2025ApJ...983...95H, sec. 4.2), text quoted
"we fit a 2D Gaussian function to velocity dispersion values determined in the Voronoi bins [vorbin, target 250 stars per
bin, 2,434 bins, median uncertainty ~0.02 mas/yr]. A single 2D Gaussian provides a decent fit with a reduced chi2 = 0.99 and
independently recovers the cluster center (dRA = -1.31 +- 0.72", dDec = 1.44 +- 0.66") with respect to the photometric
center determined in AvdM10. The best-fit position angle is 101.4 +- 2.4 deg ... flattened ... 1 - sigma_y/sigma_x = 0.09
... the residuals show that there is a significant rise of the velocity dispersion within r < 10". ... Allowing a second
Gaussian component (with the same center as the first component) ... reduce[s] the reduced chi2 to 0.91. The kinematic
center, the position angle, and the flattening of the outer component are still successfully recovered."
Their PM catalogue (oMEGACat II, 2024ApJ...970..192H): 1.4 M stars to the half-light radius (~5'), F625W to ~25, PMs from
~20-yr HST baselines, LOCALLY CORRECTED: "this local approach to measure proper motions also removes any signature of
rotation" (mean PM field ~0 everywhere; dispersions preserved). Partial energy equipartition is detected (sigma depends on
stellar mass). The HST field is a mosaic of many programmes with varying depth/exposure.
Separation from AvdM10: 1.95" (chi2 = 8.1 for 2 dof, p = 0.018); the authors call it "in agreement".

## Our own star-count centres (oMEGACat catalogue, all stars with F625W, x = dRA* east+, y = dDec, arcsec from AvdM10)
HST coverage: every 3" cell occupied to R = 250" (C1 coverage >= 0.956 inside 250").
- Inside 60" the projected density changes by only ~15% (core radius ~140"), so the core alone fixes the centre to 1-5".
- Symmetry (opposite wedges, 8 pairs, annulus 5-55"): consistent with 0 within +-2".
- Unbinned Poisson fits of an elliptical cored power law Sigma ~ (1+m^2/a^2)^-b in circular apertures about AvdM10:
  R<120": F625W<19 (-3.1+-1.1, +1.9+-1.3); 19-21 (-0.8+-1.0, +1.0+-0.8); 21-23 (-1.8+-0.8, +0.9+-0.7); all<23 (-1.7+-0.6, +1.1+-0.4)
  R<200": F625W<19 (-4.0+-0.7, +1.0+-0.7); 19-21 (-2.4+-0.5, +1.0+-0.4); 21-23 (-1.8+-0.5, +0.3+-0.5); all<23 (-2.3+-0.4, +0.6+-0.2)
  (bootstrap errors, 20 draws; q and PA poorly constrained, b and a vary with aperture).
- A PM rotation-field centre is impossible with these locally corrected PMs (mean tangential PM 0.001-0.008 mas/yr to 200").
- Motivation: about AvdM10, the bright stars (F625W<19) show a ~1.8-sigma central dip inside ~3" and a bump at 5-8"; fainter
  stars are flat (slope 0.00 +- 0.07 over 2-8.5"); about the oMEGACat VI centre nothing is seen.
Data available locally: oMEGACat VI catalogue (RA, DEC, x, y, pmra/pmdec_corrected + errors, vlos + err for ~300k MUSE stars,
f625w, f814w, quality/membership selections), Gaia DR3/EDR3 members (absolute PMs, sparse in the core), JWST NIRCam mosaics of
the centre (GO 8322, 4343; bright stars saturated), HST raw data not reduced by us (no artificial stars).

## Questions
1. What are the assumptions behind the oMEGACat VI kinematic centre and how robust is it? Identify concrete failure modes
   (model shape vs the true sigma(R), asymmetric field footprint / depth, magnitude mix + equipartition, PM error
   under-estimation varying by field, Voronoi binning, local correction, rotation-projection plumes, counter-rotating core),
   and estimate how large a bias each could produce in arcsec.
2. Is our star-count offset (~2" W, 1" N at R<200") independent evidence, or could it share systematics with theirs (same HST
   footprint and depth pattern)? What would distinguish a real large-scale centre offset (lopsidedness / m=1 mode, offset
   core vs envelope) from a systematic?
3. Design a better, more rigorous centre measurement with this data: what likelihood (counts, dispersion, LOS rotation,
   combined?), what selection (magnitude-matched, completeness-limited, symmetric footprint), what nuisance terms (m=1
   mode, radial-dependent centre, ellipticity), how to calibrate with mocks (injected offsets, realistic depth maps), and
   what precision is realistically achievable for (a) the core centre (R < 20") and (b) the global centre. Be critical and
   specific; give the expected precision with a short calculation.
4. What is the right centre to use for (i) inner profile/dip tests, (ii) Jeans/DF modelling of the inner 10", (iii) the
   IMBH question? Should the centre be a nuisance parameter in the dynamical fits?
Separate what you verified from what you assumed. Argue against your own recommendation once.

## Answer
The saved results contradict the packet’s claim that the wedge centres are all consistent with zero within ±2″. The faint-star result is **(+4.09″, −3.81″)**, with a quoted x error of 0.74″, near the search boundary. That result needs investigation before using the symmetry test as corroboration. (`results/plot_data/hst_density_centre.json:129`, `bin/hst_density_centre.py:43`)

**My recommendation is to retain AvdM10 as the reporting origin, treat the physical centre as uncertain, and fit centre coordinates jointly with the inner dynamics. Neither the published dispersion centre nor the present count fits establish a subarcsecond core centre.**

1. **How robust is the published kinematic centre?**

   It estimates the centre of an empirical, approximately elliptical **dispersion field**. Interpreting that as the potential centre assumes that tracer populations, selection and measurement errors do not introduce an asymmetric dispersion pattern. The two-Gaussian fit imposes a common centre, so its agreement does not independently demonstrate that the inner excess and outer cluster are concentric. Reduced χ² near one tests average residual size, not bias in two particular parameters. [oMEGACat VI, §4.2](https://arxiv.org/html/2503.04903v2#S4.SS2)

   The packet’s p=0.018 treats AvdM10 as exact. Adding an **assumed independent 1″ uncertainty per coordinate** for AvdM10 changes χ² from 8.07 to 2.57, giving p≈0.28. Thus “in agreement” is defensible; neither calculation establishes which centre is correct. The actual reference covariance and meaning of its quoted accuracy are needed.

   The following are **conditional sensitivity estimates, not measured biases**. For illustration, take a dispersion gradient \(g_\sigma=10^{-3}\ {\rm mas\,yr^{-1}\,arcsec^{-1}}\). An unmodelled dipole of amplitude \(A_\sigma\) can mimic approximately \(A_\sigma/g_\sigma\) arcseconds. Actual responses require the fit’s spatial weights and nuisance-parameter covariance.

   | Failure mode | Concrete mechanism and conditional scale |
   |---|---|
   | Wrong radial profile | A radial model error alone does not select a direction under symmetric sampling. If asymmetric weights convert 10–30% of a 0.01 mas/yr radial residual into a dipole, the shift is **1–3″**. |
   | Footprint/depth | Missing area alone need not bias a correctly specified conditional velocity likelihood. It changes the weighting of model errors and population mixtures. For counts, a **1–3% completeness dipole** divided by a density gradient of 0.0065 arcsec⁻¹ mimics **1.5–4.6″**. |
   | Magnitude mix/equipartition | With \(\sigma\propto m^{-\eta}\), assumed \(\eta=0.09,\sigma=0.65\), a dipole of 0.05–0.2 in effective \(\ln m\) gives **3–12″**. This is a stress test; the local mass mixture must be measured. |
   | Underestimated PM errors | \(\Delta\sigma\simeq(e_{\rm true}^2-e_{\rm quoted}^2)/(2\sigma)\). Spatially asymmetric 10–20% error underestimation with quoted errors 0.05–0.15 mas/yr gives roughly **0.4–8″**. Uniform underestimation changes the profile, not necessarily its centre. |
   | Voronoi bins | Compare the model with the constituent stars, rather than just one bin coordinate. For illustrative bin widths 10–30″ and curvature scales 50–150″, the curvature error corresponds to approximately **0.03–0.8″**; misplaced representative coordinates can contribute more. These bin sizes are assumptions. |
   | Local corrections | Removing mean rotation does not guarantee independent, unchanged dispersions. The reduction uses overlapping reference sets and adds correction uncertainty. An *unaccounted* reference-noise difference between 50 and 150 stars gives about **4″** in a simplified mean-subtraction calculation. Correct accounting removes that particular bias; correlations can still enlarge uncertainty. [oMEGACat II, §5.3](https://researchonline.ljmu.ac.uk/id/eprint/24270/1/oMEGACat.%20II.%20Photometry%20and%20Proper%20Motions%20for%201.4%20Million%20Stars%20in%20Omega%20Centauri%20and%20Its%20Rotation%20in%20the%20Plane%20of%20the%20Sky.pdf) |
   | Rotation-projection plumes | Symmetric paired plumes produce **zero first-order translation bias** with symmetric sampling. Ten-percent odd leakage of projected streaming variance for assumed \(V=3–5\) km/s and \(\sigma=20\) km/s corresponds to approximately **1–3″**. Rotating radial/tangential axes alone preserves the covariance trace. |
   | Counter-rotating core | A concentric reversal alone gives **zero translation bias**. Asymmetric unresolved streaming can produce the preceding scale; an actually displaced component shifts a common-centre fit according to its statistical weight. Counter-rotation should be allowed in the LOS mean model. [Pechetti et al.](https://arxiv.org/abs/2401.15149) |

   These effects overlap; adding those ranges in quadrature would be unjustified.

2. **Are the count offsets independent evidence?**

   They provide a different observable, but **not an independent control of HST selection systematics**. Your fit selects finite F625W and positions, then normalizes the density over an unmasked circular aperture. It contains no spatial completeness or depth term. (`bin/hst_density_centre.py:65`, `bin/hst_density_centre.py:75`, `bin/hst_density_centre.py:93`)

   Occupancy is insufficient evidence for completeness: the footprint estimator marks cells occupied by *any* catalogue star and explicitly assigns full area inside 30″. A cell can contain many detections while missing a few percent of the relevant population. (`bin/build_count_profiles.py:85`, `bin/build_count_profiles.py:99`)

   There are additional reasons to distrust the quoted count uncertainties:

   - Twenty individual-star bootstraps hold field selection errors fixed. Their estimated standard deviations themselves have roughly 16% sampling uncertainty. (`bin/hst_density_centre.py:106`)
   - The optimizer returns parameters without checking or saving convergence status. (`bin/hst_density_centre.py:84`)
   - Both faint and combined 120″ fits reach \(b\simeq5\), the imposed upper bound. That makes their errors conditional on a restrictive profile family. (`bin/hst_density_centre.py:76`; `results/plot_data/hst_density_centre.json:164`, `:225`)
   - The wedge search prefilters stars to 60″, although a 55″ annulus displaced by more than 5″ needs stars beyond 60″. The saved faint solution has a 5.6″ displacement, so a crescent of its trial annulus is unavailable. (`bin/hst_density_centre.py:42`, `:52`; `results/plot_data/hst_density_centre.json:129`)

   A real offset should persist in common magnitude/colour populations, separate programmes, independent instruments, and symmetric masks. A systematic should track exposure boundaries, crowding, error distributions or photometric thresholds. Repeat fits in disjoint radial ranges to distinguish a common translation from an inner/outer displacement.

   Crucially, a small translation produces
   \[
   \delta\ln\Sigma=-\boldsymbol{\delta c}\cdot\nabla\ln\Sigma.
   \]
   This is itself an \(m=1\) pattern. An arbitrary radial dipole and an arbitrary centre are therefore non-identifiable without an explicit convention or external constraint.

3. **A more rigorous measurement**

   I would use a **marked spatial point-process likelihood with conditional stellar velocities**:

   \[
   \log L=\sum_i\log\lambda(\mathbf x_i,m_i,k_i)
          -\int_W\lambda\,dA\,dm
          +\log p(\{\mathbf v_i\}\mid\{\mathbf x_i,m_i,k_i\}),
   \]
   where \(k_i\) identifies population or observing configuration, and
   \(\lambda=C(\mathbf x,m,k)\Sigma_k(\mathbf x)+\lambda_{\rm field}\).

   Concrete choices:

   - **Selection:** use fixed sky apertures and a fixed mask containing diametrically opposite coverage. Fit narrow magnitude/colour slices with common support. For kinematics, retain quality cuts and model their selection; matching magnitude distributions alone does not match stellar masses across populations.
   - **Counts:** use a flexible monotonic radial profile, ellipticity and position angle; integrate selection over the actual window. Without artificial-star tests, constrain relative completeness through overlaps and independent imaging, then propagate uncertainty rather than declaring completeness from occupancy.
   - **PMs:** fit Cartesian components with their measurement covariance, mass dependence and spatial error-scale terms. Forward-model local correction effects and shared-reference covariance where possible.
   - **LOS velocities:** fit the mean and dispersion simultaneously, allowing rotation reversal, perspective effects, programme zero points and binaries. A solid-body mean field alone cannot determine both centre coordinates: one translation is degenerate with systemic velocity and the perpendicular coordinate is unconstrained.
   - **Geometry:** first fit one centre. Then test separately centred inner and outer components and constrained dipole terms. Define or regularize the dipole so it cannot freely reproduce a translation.
   - **Combination:** count positions and velocities conditional on those positions can be combined. Repeated binnings or overlapping magnitude samples are not additional independent observations.
   - **Mocks:** inject 0, 0.5, 1, 2 and 5″ offsets into catalogues with realistic depth, crowding, mass mixtures, error fields, rotation and correction correlations. Include radial profiles outside the fitted family. Reapply selection and binning. Measure bias, interval coverage and false dip detections; use programme/spatial-block resampling as a complementary check.

   **Precision:** distinguish a core measured locally from a centre extrapolated inward.

   For a smooth core \(\Sigma\simeq\Sigma_0e^{-\alpha R^2}\), the stated 15% decline by 60″ gives
   \[
   \alpha=-\ln(0.85)/60^2=4.5\times10^{-5}\ {\rm arcsec^{-2}},\qquad
   {\rm SE}(c_x)\simeq\frac{1}{\alpha R_{\max}\sqrt N}.
   \]

   The saved 60″ samples contain 86,075 combined and 13,481 bright stars. Assuming approximately constant central density gives \(N_{20}\approx9,600\) and 1,500. Thus **R<20″ counts alone give about 11″ and 29″ per coordinate**, even before completeness and profile uncertainty. These are illustrative Fisher limits, not confidence intervals. (`results/plot_data/hst_density_centre.json:213`, `:30`)

   For the **global** counts, the saved 200″ models imply an approximate RMS logarithmic gradient of 0.0065–0.0070 arcsec⁻¹. Using
   \[
   {\rm SE}(c_x)\simeq\frac{\sqrt2}{\sqrt N\,g_{\rm rms}}
   \]
   gives **0.28″ combined and 0.66″ bright**, conditional on the fitted circularized profile and known selection. The underlying model parameters and sample sizes are saved at `results/plot_data/hst_density_centre.json:239` and `:56`.

   A reasonable *design target* is **1–2″ total global uncertainty after calibration**. Current results do not establish that accuracy. A distinct central dispersion excess could provide **1–3″** localization if, illustratively, \(N_{\rm eff}\sim1,000\) and \(|d\ln\sigma/dR|\sim0.01\) arcsec⁻¹; the two-component Gaussian Fisher estimate is about 2.2″ per coordinate. That locates the excess, not automatically the smooth stellar core.

4. **Which centre should enter the science?**

   - **Dip tests:** retain a predeclared reporting origin and marginalize over an independently estimated centre distribution. Do not select the centre that maximizes or removes the dip. The existing scan changes observed counts while retaining fixed model expectations, so it is a sensitivity diagnostic rather than a fitted-centre significance calculation. (`bin/check_count_centre.py:51`)
   - **Inner Jeans/DF fits:** make centre coordinates nuisance parameters and recompute stellar radii, component projections and selection integrals. The current HST loader fixes these using `OCEN_RA/DEC`; the DF objective consumes already reduced radial data and count nodes. A new centre cannot be represented by simply shifting those radial coordinates. (`src/ocen_dm/kinematics/hst_profile.py:125`; `src/ocen_dm/kinematics/df_fit.py:145`, `:164`)
   - **IMBH inference:** distinguish the cluster/potential centre from the compact object’s position. Fit their relationship explicitly and propagate astrometric registration uncertainty. The 0.77″ NE position is model-dependent: the 2024 paper reports accelerations individually consistent with zero within 3σ. Your plot additionally substitutes a 45° direction for “NE”; that marker is not a measured coordinate pair. [Häberle et al.](https://arxiv.org/html/2405.06015v2); `bin/hst_density_centre.py:19`, `:36`

**Verified versus assumed.** I verified the selection, likelihood normalization, search limits, bootstrap procedure, optimizer handling, saved results, fixed-centre projection and radial DF inputs cited above. The existing centre tests check catalogue-origin consistency and propagation, not recovery of a physical centre under asymmetric completeness; the count tests exercise amplitude recovery and integration. (`tests/test_cluster_centre.py:25`, `:47`; `tests/test_counts.py:18`, `:91`, `:105`)

The bias table, central-count extrapolation and achievable total precision are assumptions or analytical sensitivity calculations. I did not reproduce the catalogue fits or inspect the authors’ fitting implementation. Test execution was unavailable because the active Python lacks pytest; it also lacks NumPy/Astropy. Nothing was modified.

**Strongest argument against my recommendation.** A real common displacement remains plausible: all three disjoint 200″ magnitude slices lie west and north of the origin, and the combined result resembles the published dispersion offset. (`results/plot_data/hst_density_centre.json:47`, `:108`, `:169`, `:230`) Excessively flexible completeness and dipole terms could absorb that physical signal, while an AvdM10-centred prior could pull the answer back artificially. Report prior sensitivity and independently identifiable centre definitions.

**What I would check before trusting the result:**

- Full published centre covariance, reference-centre uncertainty and RA-versus-RA* convention.
- Actual bin averaging, fit weighting and centre stability under binning/profile changes.
- Expanded wedge searches with sufficient input coverage; optimization convergence, slope-bound and quadrature checks.
- Magnitude-dependent completeness and PM-error maps by programme, including local-correction covariance.
- Independent HST/JWST/Gaia registration and population-matched comparisons.
- End-to-end injected-offset coverage tests, including centre selection in the dip statistic and joint centre–IMBH inference.