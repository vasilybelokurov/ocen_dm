**The proposed `r_s ≥ 10 pc` cut does not establish the physical separation you claim.** For the implemented cored halo, neglecting the distant taper, the half-density radius is only \(0.260\,r_s\): a 10 pc scale still permits a roughly 2.6 pc core. It therefore does not automatically eliminate remnant-like concentration. [regularized_df.py:325](<src/ocen_dm/kinematics/regularized_df.py:325>)

**My recommendation: correct and validate the count likelihood first, then build a constrained multi-mass pilot using mass-resolved photometry and PMs.** Continue reporting conditional enclosed-mass sensitivity curves meanwhile. Several independently flexible DFs would be premature.

**a) Which next step, and what is minimally useful?**

There is a correctness issue upstream of the proposed model expansion. Count-profile nodes are quantiles of the observed stars’ radii, and predicted density is averaged equally over those nodes. Expected annular counts instead require integrating density over the footprint and completeness function. [build_count_profiles.py:54](<bin/build_count_profiles.py:54>), [counts.py:146](<src/ocen_dm/kinematics/counts.py:146>)

For a complete annulus sampled from the correct density, the current average approaches \(\int\Sigma^2\,dA/\int\Sigma\,dA\), whereas the required average is \(\int\Sigma\,dA/A\). This systematically overweights dense regions. My analytic check gives a 1.96% excess for \(\Sigma\propto R^{-2}\) across the adopted HST logarithmic bin width. That is an illustrative calculation, not a measurement of the actual fit bias. [build_count_profiles.py:122](<bin/build_count_profiles.py:122>)

The mocks repeat the same averaging prescription while retaining the observed nodes, so successful recovery cannot expose this error. [compact_recovery.py:55](<src/ocen_dm/kinematics/compact_recovery.py:55>)

After that correction, use:

- Two or three observed stellar-mass bins, plus an extrapolated low-mass population and a separate WD population; retain a phenomenological NS/BH component initially.
- A shared DF family with mass-dependent action scale and a small number of segregation parameters. Do not give every bin independent shape and anisotropy freedom.
- Joint radial counts by mass bin and PM dispersions by mass and radius, with selection-specific tracer predictions.

Mass-resolved PMs supply useful additional information about segregation; they do **not** independently fix unseen low-mass stars, WD retention, or DM. There is already directly relevant omega Cen evidence: oMEGACat VI measures PMs over approximately 0.29–0.69 \(M_\odot\), finding weak, radially varying partial equipartition and different radial/tangential behaviour. Do not impose full equipartition. [arXiv:2503.04903](https://arxiv.org/html/2503.04903v2)

Verified literature matches:

| Reference | Relevant method or correction |
|---|---|
| Baumgardt & Hilker 2018 — [1804.08359](https://arxiv.org/abs/1804.08359) | N-body grids fitted to density, kinematics and available stellar MFs. |
| Gieles & Zocchi 2015 — [1508.02120](https://arxiv.org/abs/1508.02120) | Self-consistent single/multimass LIMEPY models. |
| Zocchi, Gieles & Hénault-Brunet 2018 — [1806.02157](https://arxiv.org/abs/1806.02157) | Omega Cen: remnant populations and anisotropy can mimic central dark-mass signatures. |
| Hénault-Brunet et al. 2019 — [1810.05167](https://arxiv.org/abs/1810.05167) | **Mock M4 method comparison**, including a minimal three-component model. |
| Hénault-Brunet et al. 2020 — [1908.08538](https://arxiv.org/abs/1908.08538) | 47 Tuc multimass inference with variable MF and remnant content. |
| Watkins et al. 2022 — [2206.05300](https://arxiv.org/abs/2206.05300) | Mass-dependent HST PM dispersions and partial equipartition. |
| Evans, Strigari & Zivick 2022 — [2109.10998](https://arxiv.org/abs/2109.10998) | **Axisymmetric dynamical modelling, not a multimass LIMEPY application.** |

**b) Is `r_s ≥ r_half` defensible?**

As an explicitly conditional “extended halo” hypothesis, yes. As a necessary physical prior, no. Absence of stellar mass segregation in particle DM does not imply this radius inequality. Baryonic contraction and stellar heating can both alter DM concentration; neither supplies a universal cutoff at the stellar half-mass radius. [arXiv:1104.2548](https://arxiv.org/abs/1104.2548)

Nor does “non-segregated” uniquely identify DM. WD masses overlap the stellar mass range, and their distribution need not resemble a compact NS/BH subsystem. The photometric calculation itself assumes a mean WD mass of 0.55 \(M_\odot\). [build_mass_function.py:63](<bin/build_mass_function.py:63>)

Report \(M_{\rm DM}(<r)\) at several radii, alongside total mass and total unseen baryonic-plus-DM mass. This improves parameter interpretation; it cannot resolve an unidentified decomposition. Enclosed component masses are already calculated. [compact_recovery.py:112](<src/ocen_dm/kinematics/compact_recovery.py:112>)

**c) What makes the architecture easy or hard?**

- **Existing foundation:** `MassLightDFModel` already separates gravitating mixture weights from tracer weights and clears cached moments when reweighting. However, it uses the older `DFModelConfig`/DoublePowerLaw components and produces one combined tracer. Its strictly positive light fractions also prohibit a truly invisible component within that mixture. [df_mass_recovery.py:95](<src/ocen_dm/kinematics/df_mass_recovery.py:95>), [df_mass_recovery.py:140](<src/ocen_dm/kinematics/df_mass_recovery.py:140>), [positive_df.py:72](<src/ocen_dm/kinematics/positive_df.py:72>)
- **Required solver change:** the regularized configuration contains one stellar DF. Each potential iteration rebuilds its potential-dependent DF, then recomputes density. A multimass version must rebuild all component DFs in the common trial potential, sum their densities, and verify closure for each component. Freezing their initial maps would change the implemented prescription. [regularized_df.py:353](<src/ocen_dm/kinematics/regularized_df.py:353>), [regularized_df.py:415](<src/ocen_dm/kinematics/regularized_df.py:415>)
- **Main interface limitation:** every count product and kinematic profile currently receives the same projected tracer. Multimass work needs explicit selection-to-component weights; combine surface densities and pressure moments before dividing to obtain dispersions. [df_fit.py:130](<src/ocen_dm/kinematics/df_fit.py:130>), [likelihood.py:453](<src/ocen_dm/kinematics/likelihood.py:453>)
- **Remnants are prescribed gravity:** the regularized branch constructs its Plummer potential directly; changing `RemnantPlummer` alone would not change this branch. Neither that potential nor the prescribed halo establishes a positive equilibrium remnant/DM DF. [regularized_df.py:339](<src/ocen_dm/kinematics/regularized_df.py:339>)

AGAMA supports iterative multicomponent self-consistency. The harder work here is the population/selection likelihood and identifiability, not availability of a potential solver. [arXiv:1802.08239](https://arxiv.org/abs/1802.08239)

**d) Overclaims and cheaper discriminating tests**

“Demonstrably wrong” overstates what the photometric budget establishes. The script extrapolates local observed MFs below the detection range and above the turnoff, then assigns inferred WDs to each zone without retention or migration modelling. Its Poisson MF fit does not include completeness uncertainties. The trend is evidence against constant conversion, not a direct measurement of the true total stellar mass profile. [build_mass_function.py:183](<bin/build_mass_function.py:183>), [build_mass_function.py:219](<bin/build_mass_function.py:219>)

Also check the denominator: the MF catalogue requires both F625W and F814W, whereas the HST count product requires F625W only. The discrepancy matters if F814W availability varies with radius. [build_mass_function.py:67](<bin/build_mass_function.py:67>), [build_count_profiles.py:64](<bin/build_count_profiles.py:64>)

A cheaper scientific test is a two-population tracer/mass separation constrained by corrected radial count ratios, followed by prediction of held-out mass-resolved PM dispersions. That tests segregation using additional observables; tightening the halo bound merely restricts the answer.

Do not calibrate errors by forcing global \(\chi^2/N=1\). The saved fit fails its own data-fit gate, including HST count deviance 98.9 over 20 bins. Separate selection/model errors from measurement errors. [summary.json:1204](<results/df/dftwo_counts_20260923/fits/observed_no_halo_start0/summary.json:1204>)

An odd-\(L_z\) term can fit streaming means while preserving the even DF’s raw second moments. It cannot automatically cure all dispersion residuals. The current moment quadrature explicitly assumes a spherical even DF, so adding rotation requires a compatible moment calculation. [positive_df.py:148](<src/ocen_dm/kinematics/positive_df.py:148>)

**1. Verified versus assumed**

Verified: the implementation findings above; the saved baseline parameters and objective; nine focused count/mock tests passed with bytecode and pytest caching disabled. Those tests generate synthetic counts using the same density-node averaging, so they do not validate annular integration. [best.json:2](<results/df/dftwo_counts_20260923/fits/observed_no_halo_start0/best.json:2>), [tests/test_counts.py:7](<tests/test_counts.py:7>), [tests/test_compact_recovery.py:218](<tests/test_compact_recovery.py:218>)

Assumed: the remaining quoted campaign statistics, completeness claims, and astrophysical interpretation. I did not independently reproduce every fit or mock. Twelve point-estimate recovery experiments are not interval-coverage validation. Nothing was modified.

**2. Strongest argument against my answer**

Even constrained multimass inference may transfer the uncertainty into unobserved low-mass stars, WD retention and the segregation prescription. The measured equipartition is weak. If those uncertainties dominate, corrected single-tracer fits with openly conditional mass-decomposition curves could deliver nearly the same defensible conclusion at substantially lower cost. [arXiv:2503.04903](https://arxiv.org/html/2503.04903v2)

**3. What I would check to be sure**

- Independent footprint/completeness-integrated count predictions and catalogue-generated mocks, rather than reused observed nodes.
- Artificial-star completeness, binary/population effects, and matched photometric/kinematic selections.
- WD bookkeeping: WDs are currently included in the photometric stellar budget; moving them to a separate component must remove that contribution. [build_mass_function.py:234](<bin/build_mass_function.py:234>)
- Identical real/mock parameter bounds: defaults currently allow `r_s ≤ 150 pc` for real fits versus `500 pc` for mocks. [run_compact_df_recovery.py:293](<bin/run_compact_df_recovery.py:293>), [run_compact_df_recovery.py:419](<bin/run_compact_df_recovery.py:419>)
- Multistart identifiability and actual interval coverage using segregated, rotating and alternative-remnant truths.
- Preserve MF inputs and distinct output tags: the script overwrites same-tag outputs and records a plot hash, without recording the raw-catalogue hash. [build_mass_function.py:274](<bin/build_mass_function.py:274>)