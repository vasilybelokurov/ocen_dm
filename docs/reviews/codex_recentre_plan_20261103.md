**The code contradicts the rotation premise.** `DFJointProblem` clears `streaming2` only in its prediction helper, then subtracts it explicitly during scoring. Rotation remains active; a test checks this behavior ([df_fit.py:133–156](src/ocen_dm/kinematics/df_fit.py#L133), [test_positive_df.py:172–183](tests/test_positive_df.py#L172)).

**My recommendation: proceed as a sensitivity calculation, but do not interpret the proposed rerun as validated marginal inference.** Ranked concerns:

1. **The main problem is changing the data summary, not reusing stars across centres.** Averaging alternative likelihoods does not itself double-count measurement noise. The correct object is \(\int p(D\mid\theta,c)p(c)\,dc\), for one fixed observed catalogue and selection window.

   Here, membership and even the selected sample change with centre: both count builders use displaced radial cuts, and HST dispersions are re-estimated from displaced annuli ([centre_marg.py:99–147](src/ocen_dm/centre_marg.py#L99)). Normalized Poisson and split-normal terms do **not** establish that these different summaries have comparable likelihoods ([counts.py:178–196](src/ocen_dm/kinematics/counts.py#L178), [likelihood.py:504–516](src/ocen_dm/kinematics/likelihood.py#L504)). Compression can discard centre-dependent factors; uncertainty could consequently be understated or overstated. Numerical convergence cannot settle that. Also, integrating a likelihood with amplitudes profiled separately at each centre is a hybrid procedure, not full Bayesian marginalization over those amplitudes ([df_fit.py:171–185](src/ocen_dm/kinematics/df_fit.py#L171)).

2. **Rotation needs a consistent per-centre treatment.** MUSE refits a sinusoidal mean but returns only dispersion; the published correction and errors remain unchanged when the dispersion changes ([centre_marg.py:150–175](src/ocen_dm/centre_marg.py#L150), [centre_marg.py:211–215](src/ocen_dm/centre_marg.py#L211)). Its \(v_{\rm rot}^2/2\) correction assumes an azimuthal average ([likelihood.py:177–185](src/ocen_dm/kinematics/likelihood.py#L177)). Recompute the fitted rotational contribution with the actual selected stars’ angular weights, excluding systemic velocity, or demonstrate that retaining it changes model comparisons negligibly.

   HST already recalculates its external tangential correction at each rebuilt median radius ([likelihood.py:287–305](src/ocen_dm/kinematics/likelihood.py#L287)). Whether that adequately represents locally corrected PMs under displacement requires validation; reading corrected columns cannot establish it.

3. **Scrambled Sobol is reasonable; 200–400 draws is not an acceptance criterion.** Use powers of two—128, 256, 512—with the same centres for every model and several independent scrambles ([SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html)). Require uncertainty in each paired marginal \(\Delta(-2\log L)\), estimated across scrambles, comfortably below 0.5.

   Report normalized contributions
   \[
   \alpha_k=\frac{w_k e^{\ell_k}}{\sum_j w_j e^{\ell_j}},
   \qquad N_{\rm eff}=1/\sum_k\alpha_k^2,
   \]
   maximum contribution, and locations of dominant centres, separately for each model. ESS diagnoses concentration; it does not certify QMC convergence.

   **T4 is not fully explained by bin crossings.** Narrow likelihood peaks and optimizer failures are alternatives. The estimator does not check optimizer success, and interval scans can return an endpoint without finding a crossing ([outer_profile.py:578–605](src/ocen_dm/kinematics/outer_profile.py#L578), [outer_profile.py:618–640](src/ocen_dm/kinematics/outer_profile.py#L618)).

4. **Warm starts are legitimate, but the proposed checks are insufficient.** Compare maximized likelihoods, fitted means, field fractions, errors, and final model comparisons—not only dispersions. Errors directly affect both residual penalties and likelihood normalization ([likelihood.py:514–516](src/ocen_dm/kinematics/likelihood.py#L514)). Audit dominant centres, tails and diagonal offsets; one cold start is not an independent guarantee.

   A cheaper first optimization is caching the field density evaluated at every star: those Cartesian PM inputs remain unchanged while masks and projection angles change ([hst_profile.py:135–153](src/ocen_dm/kinematics/hst_profile.py#L135), [outer_profile.py:514–516](src/ocen_dm/kinematics/outer_profile.py#L514)). Benchmark the full rebuild: removing starts leaves the repeated interval optimizations intact ([outer_profile.py:598–643](src/ocen_dm/kinematics/outer_profile.py#L598)).

5. **Gaia can be omitted only after a likelihood-level check.** One 2″ displacement is insufficient. Test both signs on both axes, diagonals and representative prior-tail offsets; compare complete model likelihood ratios, including changed errors and radial nodes. The existing builder leaves Gaia kinematics unchanged through its fallback branch ([centre_marg.py:203–218](src/ocen_dm/centre_marg.py#L203)). Small displacement divided by radius alone does not bound bin-membership effects.

6. **Defer a comprehensive catalogue model, but pursue a fixed-data middle path now.** Use unbinned positions or fixed sky pixels for counts. For inner PMs, reuse the existing heteroscedastic Cartesian mixture construction with DF moments evaluated per star; keep nuisance parameters shared over fixed sky subsets ([outer_profile.py:524–563](src/ocen_dm/kinematics/outer_profile.py#L524)). A hybrid inside 60″ is defensible if the split is fixed on the sky and every star enters once.

   The current T8 does not isolate binning error: its unbinned window is fixed about AvdM10, while its binned window follows the trial centre ([centre_marginal.py:76–94](bin/centre_marginal.py#L76), [centre_marginal.py:140–162](bin/centre_marginal.py#L140)).

7. **Several acceptance claims need strengthening.** T5 checks an arbitrary total-count difference below 200, rather than the documented boundary-crossing identity; T6 checks one offset rather than the grid ([test_centre_marg.py:95–104](tests/test_centre_marg.py#L95)). T9’s recorded discrepancy reaches 0.62σ against a 0.5σ requirement ([CENTRE_MARGINALISATION.md:49](docs/CENTRE_MARGINALISATION.md#L49), [CENTRE_MARGINALISATION.md:81](docs/CENTRE_MARGINALISATION.md#L81)).

   The reported central-deficit probability averages Poisson CDFs whose **observed count changes with centre**. That is not automatically a calibrated significance measure; simulate the entire procedure or define a fixed observable ([centre_marginal.py:150–163](bin/centre_marginal.py#L150), [centre_marginal.py:183–184](bin/centre_marginal.py#L183)). Fixed saved fits also cannot establish the preference after refitting ([centre_marginal.py:48–55](bin/centre_marginal.py#L48)).

**Verified versus assumed.** I read the implementations, tests, earlier reviews and saved diagnostics. The saved T4 difference is 3.9045 ([output:558–561](results/plot_data/centre_marginal_20261103.json#L558)). I did not reproduce catalogue estimates or timings. I attempted the four synthetic centre tests, but the active Python lacks `pytest`. Assumptions requiring evidence include the 2″ prior’s calibration and independence, estimator-difference transferability, and negligible residual PM correlations.

**Strongest argument against my answer.** A calibrated binned approximation may preserve the scientific comparison while a richer catalogue model introduces larger systematic errors. The existing counts cross-check preserves the sign ([CENTRE_MARGINALISATION.md:80](docs/CENTRE_MARGINALISATION.md#L80)). That supports continuing the diagnostic, although its unequal windows prevent it from validating the approximation.

**Needed to be sure.** End-to-end displaced-centre BH/no-BH mocks; comparison against a fixed-data likelihood; rotation and covariance sensitivity; optimizer and interval audits at dominant centres; replicated integration convergence; and refits before interpreting model preference. Baseline reproduction and a delta-prior identity test implementation consistency, not inferential calibration.