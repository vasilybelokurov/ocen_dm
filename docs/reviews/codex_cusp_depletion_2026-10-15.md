The mass-budget premise needs correcting first: with proper renormalization, depletion redistributes the fixed \(M_\star\); it does **not** remove total stellar mass or require replacement by remnants. The amplitude is explicitly normalized to the requested mass. ([regularized_df.py:519–525](<../../src/ocen_dm/kinematics/regularized_df.py#L519>))

**I recommend fixed \(p=4\), one free \(J_c\), as a conditional compatibility experiment.** Include \(M_\bullet=0\), retain the undepleted branch, and test sensitivity to fixed smaller \(p\). Do not interpret this experiment as measuring a cusp slope, identifying a depletion mechanism, or establishing a black-hole mass.

**A. Fixed \(p=4\) is adequate for the first experiment; free \(p\) is not yet justified.**

For an isotropic Kepler limit, the asymptotic result is

\[
\rho(r)\propto
\begin{cases}
r^{-(3-p)/2},&0\le p<2,\\
r^{-1/2}\log(r_0/r),&p=2,\\
r^{-1/2},&p>2.
\end{cases}
\]

Thus \(\gamma=(3-p)/2\) cannot be extrapolated to \(p=4\). Increasing \(p\) above two changes the transition, not the limiting slope. The isotropic lower bound \(\gamma\ge1/2\) is consistent with [An & Evans](https://arxiv.org/abs/astro-ph/0511686).

A three-dimensional \(r^{-1/2}\) cusp has finite central surface density and can look flat in projection. Whether its slope is sufficiently small **between the actual observed radii** depends on its amplitude, transition scale and surrounding cluster. Flat projected counts therefore neither require a three-dimensional core nor establish \(p=4\).

The alternatives have different interpretations:

- **Energy depletion:** equivalent to action depletion in the unbiased pure-Kepler limit, since \(J_r+L=GM_\bullet/\sqrt{2\mathcal E}\). It offers no immediate identifiability advantage. Away from that limit, prefer \(D(L_c)\) because it preserves the existing contours; those contours depend on the potential and anisotropy prescription. ([regularized_df.py:527–559](<../../src/ocen_dm/kinematics/regularized_df.py#L527>))
- **A positive DF floor:** restores the \(r^{-3/2}\) cusp sufficiently far inward. It can hide that cusp below resolution, but adds an unconstrained crossover.
- **A loss cone:** suppresses low angular momentum at fixed energy, unlike suppression of low \(L_c\). Use it only for a specific destruction/scouring hypothesis and test the predicted tangential anisotropy. The distinction is explicit in the loss-cone treatment of [González Prieto et al.](https://arxiv.org/html/2507.06316v2#S2.SS1).
- **Relaxation or mass segregation:** does not directly motivate \(p=4\). Relaxed light populations generally retain approximately \(r^{-3/2}\) cusps; heavier populations become steeper. A shallow depleted population instead needs an evolutionary explanation and a cusp-regrowth timescale. ([Alexander & Hopman](https://arxiv.org/abs/0808.3150), [Merritt](https://arxiv.org/abs/0909.1318))

Fit fixed \(p=4\), then repeat with, for example, \(p=0.4\) and \(p=1\). These span limiting slopes \(1.3\), \(1\), and \(0.5\) without simultaneously fitting two poorly separated shape parameters. Only free \(p\) if likelihood profiles and injections demonstrate that the observed radial range separates it from \(J_c\).

**B. A depleted tracer can retain a dispersion rise; an upward mass shift is plausible, not compulsory.**

In an isotropic, black-hole-dominated power-law region, Jeans gives

\[
\sigma_r^2=\frac{GM_\bullet}{(\gamma+1)r}.
\]

Changing \(\gamma\) from \(3/2\) to \(1/2\) therefore **increases local variance by \(5/3\)** at fixed mass. Depletion preferentially removes tightly bound orbits; surviving stars passing through the centre can have high velocities.

Projection works in the opposite direction by reducing their contribution relative to foreground/background cluster stars. Under idealized asymptotic assumptions:

- \(\gamma=3/2:\ \sigma_{\rm proj}^2\propto R^{-1}\);
- \(\gamma=1/2\), with finite central \(\Sigma:\ \sigma_{\rm proj}^2\propto R^{-1/2}\).

The rise becomes shallower, but survives. These are analytical limits, not predictions for your finite bins. The implementation projects density and pressures separately, so changing their radial shapes changes both numerator and denominator. ([positive_df.py:322–338](<../../src/ocen_dm/kinematics/positive_df.py#L322>))

Consequently, there is a mass–depletion degeneracy, but no universal monotonic rule \(J_c\uparrow\Rightarrow M_\bullet\uparrow\). Calculate the response in every central PM/LOS bin before refitting, then profile \(M_\bullet\) against \(J_c\). The supplied information cannot determine the mass shift quantitatively.

**C. Keep Plummer remnants for the first comparison, but make the conclusion conditional on that choice.**

A prescribed Plummer component is mathematically usable as an external gravitational source. It does not demonstrate a physically admissible remnant DF or collisional equilibrium; the implementation supplies its potential without constructing either. ([regularized_df.py:648–654](<../../src/ocen_dm/kinematics/regularized_df.py#L648>))

Also, your reported compensation is **not an increase in central remnant density**. From the saved parameters,

\[
\frac{\rho_{\rm rem,new}(0)}{\rho_{\rm rem,base}(0)}
=\frac{503869}{273161}\left(\frac{2.180}{4.045}\right)^3
\simeq0.29.
\]

The larger remnant mass is substantially more extended. ([batch.json:362–363](<../../results/df/dfbh40000_seed_20261015/batch.json#L362>), [summary.json:205–206](<../../results/df/dfbh40000_seed_20261015/fits/observed_no_halo_start0/summary.json#L205>))

Do not add freely adjustable stellar depletion and a freely adjustable remnant cusp simultaneously in the first experiment. Subsequently compare a specified concentrated-remnant alternative. The cited evolutionary models predict stellar and stellar-BH cusps with different slopes, supporting that sensitivity test. ([González Prieto et al., Appendix A.3](https://arxiv.org/html/2507.06316v2#A1.SS3))

If depletion is meant to affect only the bright tracer through mass segregation, applying it to the entire self-gravitating stellar DF is an additional assumption. The current count likelihood uses the same model surface density, with a scalar amplitude per count product; it cannot represent a radial change in tracer-to-stellar-mass ratio. ([df_fit.py:151–156](<../../src/ocen_dm/kinematics/df_fit.py#L151>))

**D. The principal implementation risks are numerical, especially the escape shell.**

- **Normalization:** multiply the same depletion factor into the callable and normalization integrand. Evaluate it stably, handle \(J_c=0\) explicitly, and check integration-domain extension independently of order doubling. Both existing orders use the same finite domain. `totalMass()` repeats the normalization quadrature, so it is not independent spatial mass verification. ([regularized_df.py:578–614](<../../src/ocen_dm/kinematics/regularized_df.py#L578>))
- **Contour boundaries:** retain direct integration outside the grid. Test continuity and derivatives at both boundaries; depletion introduces another action scale that can expose interpolation errors. ([regularized_df.py:562–575](<../../src/ocen_dm/kinematics/regularized_df.py#L562>))
- **Underflow cutoff:** its high-action envelope remains valid for unchanged contours and \(0\le D\le1\), using the **new normalized amplitude**. Nonmonotonicity of the depleted DF does not invalidate that envelope. Never add a hard low-action cutoff. ([regularized_df.py:595–608](<../../src/ocen_dm/kinematics/regularized_df.py#L595>))
- **Action guard:** the claim that an overestimated action can only remove mass becomes false for a nonmonotonic depleted DF. The guard recomputes finite actions beyond 300 pc only when their approximate DF is positive; finite errors inside that radius remain unchecked. Compare against exact actions across the depleted family, including erroneous actions landing near its zero or maximum. ([regularized_df.py:294–318](<../../src/ocen_dm/kinematics/regularized_df.py#L294>))
- **Escape shell:** fixed speed quadrature can miss the finite-binding-energy stars dominating a strongly depleted central population. In a Kepler potential their fractional speed-shell width shrinks roughly as \(r/r_{\rm depletion}\). Agreement between two inadequate orders is insufficient. Use an independent binding-energy/adaptive integration reference. The production rule integrates directly in \(v/v_{\rm esc}\), while equilibrium closure checks density only. ([positive_df.py:162–180](<../../src/ocen_dm/kinematics/positive_df.py#L162>), [regularized_df.py:759–774](<../../src/ocen_dm/kinematics/regularized_df.py#L759>))

Use separate undepleted jobs and positive, preferably logarithmic, \(J_c\) starts. Near \(J_c=0\), \(p=4\) gives a weak response that can defeat finite-difference optimization. The driver uses fixed steps in normalized coordinates and resumes saved coordinate vectors, so use fresh batches when changing the parameterization. ([run_compact_df_recovery.py:503–510](<../../bin/run_compact_df_recovery.py#L503>), [run_compact_df_recovery.py:583–588](<../../bin/run_compact_df_recovery.py#L583>))

**E. Retaining the tension is defensible; cutting at \(3''\) is only a sensitivity test.**

A \(3''\) cut cannot remove a discrepancy extending through \(13''\). Report excluded-bin predictions and compare masses using identical retained datasets. The centre diagnostic covers eight bins from \(2''\) to approximately \(13''\), using fixed model expectations without refitting. ([check_count_centre.py:40–56](<../../bin/check_count_centre.py#L40>))

Flat magnitude ratios also cannot exclude incompleteness shared by all slices; the diagnostic explicitly admits that ambiguity. ([check_hst_count_completeness.py:4–7](<../../bin/check_hst_count_completeness.py#L4>))

1. **Verified versus assumed.**

   I verified the implementation issues above and inspected the saved refits. Both stop on `ftol`, with reported optimality approximately 25 and 55; both pass numerical validation but fail the configured data-fit criteria. These are candidate fits, not established global minima. ([8200 summary:112–115](<../../results/df/dfbh8200_seed_20261015/fits/observed_no_halo_start0/summary.json#L112>), [40000 summary:112–115](<../../results/df/dfbh40000_seed_20261015/fits/observed_no_halo_start0/summary.json#L112>), [8200 gates:1314–1344](<../../results/df/dfbh8200_seed_20261015/fits/observed_no_halo_start0/summary.json#L1314>), [40000 gates:1314–1344](<../../results/df/dfbh40000_seed_20261015/fits/observed_no_halo_start0/summary.json#L1314>))

   I did not rerun tests, refits or catalogue diagnostics. I assume the packet’s empirical completeness and centre summaries are accurate. The asymptotic calculations assume spherical equilibrium and the stated isotropic Kepler limit; their applicability at observed radii remains unverified. Nothing was modified.

2. **Strongest argument against my recommendation.**

   Fixed \(p=4\) imposes the shallowest isotropic asymptote and may manufacture compatibility through an extreme tracer choice. The evolutionary alternative with \(\gamma\simeq1.3\) could yield a materially different mass conclusion. An improved objective could therefore demonstrate model flexibility while weakening physical interpretability. The published \(\gamma\simeq1.3\) result is a stellar **mass-density** fit over \(0.001\)–1 pc, not a direct prediction for your magnitude-selected counts. ([González Prieto et al., Appendix A.3](https://arxiv.org/html/2507.06316v2#A1.SS3))

3. **Minimum checks needed before trusting refits.**

   - Exact legacy recovery at \(J_c=0\), schema/optimizer round-trips, positivity, origin behavior and stable extreme-parameter evaluation. Existing tests assert an undepleted, nonzero origin; that does not validate depletion. ([test_regularized_df.py:25–54](<../../tests/test_regularized_df.py#L25>))
   - Independent Kepler density **and pressure** benchmarks spanning \(p<2\), \(p=2\), \(p>2\), including shell resolution and projected asymptotes. The existing constant-DF benchmark does not exercise this regime. ([test_regularized_df.py:209–220](<../../tests/test_regularized_df.py#L209>))
   - Independent action-space normalization, spatial mass closure, exact-action rebuilds, contour-boundary checks, and convergence of every count expectation and dispersion bin.
   - Multiple starts and numerical/Jacobian sensitivity checks at all three masses; then profiles in \((M_\bullet,J_c)\) and fixed-\(p\)/remnant alternatives.
   - Repeated noisy injections, including undepleted truth, depleted truth and mismatched tracer/remnant families, before interpreting significance or coverage. Existing mocks supply Gaussian dispersion noise and Poisson counts but remove streaming, so they do not validate streaming systematics. ([compact_recovery.py:44–64](<../../src/ocen_dm/kinematics/compact_recovery.py#L44>))