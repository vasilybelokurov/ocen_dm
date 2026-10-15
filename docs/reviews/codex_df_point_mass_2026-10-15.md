The scratch patch leaves incorrect finite-centre assumptions in place: the frequency ratio is forced to 2 at \(L=0\), and the spline enforces the corresponding boundary derivative. With a point mass, the \(L\to0^+\) frequency ratio is 1. Bypassing the guard therefore does **not** establish a correct BH extension. ([regularized_df.py:95–128](<src/ocen_dm/kinematics/regularized_df.py#L95>), [regularized_df.py:147–158](<src/ocen_dm/kinematics/regularized_df.py#L147>))

**My recommendation: proceed as a conditional compatibility experiment, but change the depletion coordinate, fix the radial boundary, and validate observables before fitting. This is more than a 30-line production change.**

1. **What \(D(q)\) actually does.** It suppresses small *total action*, rather than all low-angular-momentum orbits. A radial orbit with large \(J_r\) survives. Away from a Kepler potential, \(q\) is not an energy coordinate, so multiplication changes the angular weighting of the existing DF; its original contour construction does not prevent that. The implementation explicitly distinguishes \(q\) from its contour label \(L_c\). ([regularized_df.py:225–233](<src/ocen_dm/kinematics/regularized_df.py#L225>), [regularized_df.py:260–283](<src/ocen_dm/kinematics/regularized_df.py#L260>))

   For example, in a harmonic core,
   \[
   E-\Phi(0)=\Omega(2J_r+L),\qquad
   q=\frac{E-\Phi(0)}{2\Omega}+\frac L2.
   \]
   At fixed energy, increasing \(D(q)\) favours larger \(L\). Thus it can introduce tangential bias and a central density depression even without a BH. It does not necessarily empty the centre: finite-action stars still pass through it.

   The small-action slope also needs correction. In an isotropic Kepler control, writing positive binding energy as \(\mathcal E\),
   \[
   q=\frac{GM_\bullet}{\sqrt{2\mathcal E}},\qquad
   f\sim q^p\sim\mathcal E^{-p/2},
   \]
   and
   \[
   \rho=4\pi\sqrt2\int_0^\psi f(\mathcal E)
                 (\psi-\mathcal E)^{1/2}\,d\mathcal E.
   \]
   Consequently,
   \[
   \rho(r)\sim
   \begin{cases}
   r^{-(3-p)/2},&0\le p<2,\\
   r^{-1/2}\log(r_0/r),&p=2,\\
   r^{-1/2},&p>2.
   \end{cases}
   \]
   The proposed quadratic depletion is the **logarithmic borderline**, not a pure \(r^{-1/2}\) cusp. Every nonzero \(J_c\) has that same limiting behaviour; \(J_c\) changes where the transition occurs. These are analytical deductions, consistent with the isotropic bound in [An & Evans](https://arxiv.org/abs/astro-ph/0511686).

   Renormalization is mandatory: include the same depletion factor in both the callable DF and its action integral, then recompute the amplitude during every potential rebuild. Currently those expressions contain only the exponential. Merely multiplying the returned DF changes the requested stellar mass. ([regularized_df.py:217–223](<src/ocen_dm/kinematics/regularized_df.py#L217>), [regularized_df.py:276–305](<src/ocen_dm/kinematics/regularized_df.py#L276>))

2. **Form I would choose.** Use depletion in **\(L_c\)**. With zero anisotropy bias and correct frequencies, that preserves the energy contours; multiplying by \(D(q)\) generally does not. A function of \(L_c\) also preserves the existing contour derivative ratios wherever its derivative is nonzero. This follows from the characteristic equation used to construct \(L_c\). ([regularized_df.py:177–179](<src/ocen_dm/kinematics/regularized_df.py#L177>), [regularized_df.py:225–233](<src/ocen_dm/kinematics/regularized_df.py#L225>))

   For an initial one-parameter compatibility screen, I would use
   \[
   D(L_c)=\frac{(L_c/J_c)^4}{1+(L_c/J_c)^4},
   \]
   with an explicit \(J_c=0\) branch returning 1, and compare against the undepleted model. Call \(J_c\) a depletion scale.

   If the objective is to infer the **asymptotic slope**, free \(p\), with \(\gamma=(3-p)/2\) for \(0\le p<2\). Freeing both \(p\) and the transition scale requires two parameters. No single scale parameter independently controls both quantities. Neither scalar depletion represents a loss cone or permits an intrinsically flat isotropic stellar core around a BH.

3. **Numerical risks.**
   
   - Fix both the radial frequency boundary and the spline boundary. Replacing only the frequency assignment leaves the spline enforcing the wrong condition. ([regularized_df.py:128](<src/ocen_dm/kinematics/regularized_df.py#L128>), [regularized_df.py:147–158](<src/ocen_dm/kinematics/regularized_df.py#L147>))
   - The dangerous velocity regime is the thin shell near escape speed: after strong high-binding depletion, stars with finite binding energy dominate the central density. The fixed Gauss–Legendre speed rule can miss that shell. Density closure alone does not certify pressure accuracy: the closure check extracts only density, while the observable table subsequently uses the nominal velocity order. ([positive_df.py:157–175](<src/ocen_dm/kinematics/positive_df.py#L157>), [regularized_df.py:447–461](<src/ocen_dm/kinematics/regularized_df.py#L447>), [regularized_df.py:468–471](<src/ocen_dm/kinematics/regularized_df.py#L468>))
   - Outside the action lookup grid, evaluation integrates contours directly; normalization nevertheless covers a finite action interval. Doubling node counts does not enlarge that interval or change `r_min`. Test domains separately from resolution. ([regularized_df.py:260–284](<src/ocen_dm/kinematics/regularized_df.py#L260>), [compact_recovery.py:95–109](<src/ocen_dm/kinematics/compact_recovery.py#L95>))
   - The unresolved stellar cusp is continued using the fitted innermost logarithmic slope. Vary `r_min` and check its enclosed mass and force; an integrable cusp does not automatically imply an accurate continuation. ([positive_df.py:180–199](<src/ocen_dm/kinematics/positive_df.py#L180>))
   - Add explicit finite-value checks to final closure diagnostics: the regularized constructor uses threshold comparisons that do not themselves reject NaNs. ([regularized_df.py:452–462](<src/ocen_dm/kinematics/regularized_df.py#L452>))

4. **Identifiability.** Taking the supplied influence radii and radial coverage as assumptions, \(8200\,M_\odot\) is marginally resolved; \(4\times10^4\,M_\odot\) should affect more bins. That supports testing compatibility, not assuming either mass is identifiable.

   Profile \(M_\bullet\) against depletion scale, anisotropy, stellar mass, remnant mass/scale and distance. Counts constrain tracer shape but do not independently fix stellar mass when their amplitudes are profiled—the saved base fit uses free count amplitudes. ([df_fit.py:153–156](<src/ocen_dm/kinematics/df_fit.py#L153>), [summary.json:183–202](<results/df/dftwo_counts_area_20261013/fits/observed_no_halo_start3/summary.json#L183>))

   Require repeated noisy injections at \(0,\ 8200,\ 40000\,M_\odot\), spanning depletion and remnant alternatives. Measure false exclusions, mass bias and interval coverage. The existing mock generator supplies Gaussian dispersion noise and Poisson counts, but deliberately drops streaming; it is not by itself a test of streaming-systematic robustness. ([compact_recovery.py:31–64](<src/ocen_dm/kinematics/compact_recovery.py#L31>))

5. **Fixed Plummer remnants.** Acceptable as a declared first sensitivity case; insufficient for a robust BH measurement. The component is prescribed gravity with no remnant DF, so stellar self-consistency establishes nothing about remnant equilibrium. ([regularized_df.py:339–343](<src/ocen_dm/kinematics/regularized_df.py#L339>), [regularized_df.py:415–423](<src/ocen_dm/kinematics/regularized_df.py#L415>))

   Missing concentrated remnant mass could be assigned to the BH, although refitting can change the direction of the final bias. Bracket plausible centrally concentrated remnant profiles before quoting a BH constraint. A specific relaxed remnant cusp requires additional physical assumptions.

6. **Engineering beyond the DF expression.**
   
   - Persist and validate \(M_\bullet,J_c\) in configuration. The scratch BH lives in a mutable global outside `cfg`, so its serialized configuration cannot reproduce the potential. ([packet:2546–2566](scratch feasibility script, not in the repository))
   - Add fitting coordinates and starting values explicitly: the driver constructs a fixed coordinate list, rejects unknown bound overrides, and copies only named parameters for `--start-from`. ([run_compact_df_recovery.py:175–181](<bin/run_compact_df_recovery.py#L175>), [run_compact_df_recovery.py:290–317](<bin/run_compact_df_recovery.py#L290>), [run_compact_df_recovery.py:900–905](<bin/run_compact_df_recovery.py#L900>))
   - Handle zero boundaries deliberately; logarithmic coordinates exclude zero. ([df_fit.py:183–201](<src/ocen_dm/kinematics/df_fit.py#L183>))
   - Fix mass accounting: `mass_profiles()["total"]` sums stars, remnants and halo only. Its recovery metrics could therefore validate a “total” that omits the BH. Also extend diagnostics below the current 0.5-pc reporting minimum. ([compact_recovery.py:112–141](<src/ocen_dm/kinematics/compact_recovery.py#L112>), [run_compact_df_recovery.py:489](<bin/run_compact_df_recovery.py#L489>))

7. **Simpler experiment.** First correct point-mass support and refit the original DF at fixed \(M_\bullet=0,\ 8200,\ 40000\,M_\odot\), without depletion. Then add \(D(L_c)\) and repeat with matched nuisance freedom. Finally profile BH mass and calibrate conclusions with injections.

   The baseline already fails its own data-fit gates: \(\chi^2_{\rm kin}=156.39\) for 89 points, HST count deviance \(96.92\) for 20 bins, and `data_fit_passed=false`. Numerical validation passes. Therefore a relative fit improvement alone is insufficient evidence for a calibrated BH constraint. ([summary.json:1199–1241](<results/df/dftwo_counts_area_20261013/fits/observed_no_halo_start3/summary.json#L1199>))

Separately:

1. **Verified versus assumed.** I verified the implementation, configuration/driver paths, cited tests and saved baseline diagnostics by reading them. I did not execute tests or reproduce the scratch results. The regularized test currently requires BH rejection; the existing Kepler moment test exercises the separate DoublePowerLaw branch. ([test_regularized_df.py:91–92](<tests/test_regularized_df.py#L91>), [test_positive_df.py:73–94](<tests/test_positive_df.py#L73>)) The cusp scaling is my derivation under the stated isotropic Kepler assumptions. The supplied influence radii and physical remnant-response picture remain assumptions.

2. **Strongest argument against my recommendation.** The scientific question concerns resolved bins, not the mathematical origin. The radial-boundary defect and \(D(q)\)-induced angular changes might be negligible at observational precision. If corrected and uncorrected implementations give indistinguishable bin predictions, the original quadratic proposal could suffice for a narrow compatibility test. Its positivity is not the problem; interpretation and demonstrated accuracy are.

3. **Checks needed to be sure.** Require pure-Kepler frequency/action round trips across eccentricity; \(L_c=J_r+L\) and isotropy controls including the radial boundary; the derived \(p<2,\ p=2,\ p>2\) density limits; exact old-model recovery at both new parameters zero; independent action-mass normalization; separate radial-domain, velocity, contour and projection convergence; direct AGAMA projected-moment comparisons at actual selection nodes; and noisy BH recovery under alternative remnant/anisotropy models. Use [Baes & Dejonghe](https://arxiv.org/abs/astro-ph/0402527) as an analytic numerical benchmark, supplemented by the shallow-cusp controls above. Also require saved-config replay to reproduce the BH force and reported total mass.