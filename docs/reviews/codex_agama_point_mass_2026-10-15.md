The scaled-energy hypothesis is **not supported as stated**. AGAMA explicitly handles \(\Phi(0)=-\infty\): its coordinate
\[
X=\log[1/\Phi(0)-1/E]
\]
then becomes \(\log(-1/E)\). There is no automatic compression of weak binding energies caused by an infinite central potential. The grid is adaptive in log radius. See [`scaleE`, actions_spherical.cpp:42](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L42), [`createInterpolationGrid`, potential_utils.cpp:738](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/potential_utils.cpp#L738), and the [upstream implementation](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp).

Also, your frequency mode named `"exact"` still trusts **finite, bound mapper outputs without checking their action residual**. It inverts \(J_r(E,L)\) only for nonfinite or unbound outputs. Accurate quadrature at the wrong energy remains a possible error. [`regularized_df.py:119`](src/ocen_dm/kinematics/regularized_df.py:119)

**A. AGAMA has a direct route, but neither suggested option guarantees your required accuracy.**

- **`ActionFinder(potential, interp=False)` does not disable spherical interpolation.** The factory selects `ActionFinderSpherical` without consulting `interpolate`; that argument reaches the axisymmetric branch only. [`actions_factory.cpp:35`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_factory.cpp#L35)
- **`agama.actions(potential, points)` selects direct spherical integration automatically.** No method argument is needed. It dispatches to `evalSpherical`, which computes energy from the supplied potential, finds turning points, and integrates \(J_r\). [`interface_python.cpp:3724`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/interface_python.cpp#L3724), [`actions_factory.cpp:20`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_factory.cpp#L20), [`actions_spherical.cpp:585`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L585)
- **Direct does not mean arbitrarily accurate.** Its quadrature order is capped at 40, with a heuristic documented for approximately constant accuracy down to \(r_p/r_a=10^{-3}\). That is inadequate assurance for your extreme eccentricities. [`actions_spherical.cpp:30`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L30), [`math_core.h:560`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/math_core.h#L560)

I checked the installed binary in pure Kepler, \(M=8200\), \(J_r+L=100\):

| \(L/(J_r+L)\) | `agama.actions` relative \(J_r\) error | `ActionFinder` error |
|---:|---:|---:|
| 0.001 | \(+4.49\times10^{-4}\) | \(+3.61\times10^{-4}\) |
| 0.01 | \(-1.42\times10^{-5}\) | \(-1.07\times10^{-3}\) |
| 0.3 | \(+8.78\times10^{-13}\) | \(-4.05\times10^{-7}\) |

`interp=False` reproduced the default finder. Thus the \(4.5\times10^{-4}\) Kepler error also exists **without action interpolation**.

There are nevertheless concrete interpolation failure mechanisms worth investigating:

- The finder computes energy using its **internal interpolated potential**, rather than the original composite. Your comparison therefore mixes potential-interpolation and action-interpolation errors. [`actions_spherical.cpp:649`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L649)
- It interpolates \(J_r/(L_c-L)\), clips coordinates to its table boundaries, and clamps negative results to zero. This provides a mechanism for spurious zero actions, although I have not established which mechanism caused your start-2 failure. [`actions_spherical.cpp:337`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L337), [`actions_spherical.cpp:73`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L73)
- The mapper extrapolates its energy table and applies one first-order correction. That is not a converged action inversion. [`actions_spherical.cpp:679`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L679), [`actions_spherical.cpp:725`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L725)

**B. Establish an all-direct BH reference before adopting a hybrid.**

I would first use converged direct \(J_r\) quadrature throughout BH models. Use that implementation to determine where approximation is safe. Replacing only weakly bound actions prematurely assumes the error depends primarily on energy; the eccentricity dependence above already contradicts that simplification.

Your existing helper needs hardening before becoming a general action backend: it assumes `r0` lies inside the orbit, uses fixed-order quadrature, and clips negative radicands. Require valid brackets, explicit radial/circular limits, and convergence checks rather than treating every returned number as certified. [`regularized_df.py:179`](src/ocen_dm/kinematics/regularized_df.py:179), [`regularized_df.py:213`](src/ocen_dm/kinematics/regularized_df.py:213)

For a production hybrid:

- Validate an error envelope over **energy, circularity, and allowed potential parameters**, including actual velocity quadrature nodes.
- Choose the boundary from errors in density, pressures, projected predictions, and likelihood—not an arbitrary radius.
- A hard switch is acceptable if the disagreement is negligible throughout the overlap. Otherwise moving model parameters can move quadrature nodes across a finite jump. Smooth blending removes that jump but does not remove bias or guarantee monotonic actions.

Your velocity nodes are fixed in escape-speed fraction, so their physical energies change with the potential; a fixed Gauss rule does not eliminate this switching problem. [`positive_df.py:157`](src/ocen_dm/kinematics/positive_df.py:157)

If direct evaluation is too expensive, your own table is reasonable. However, AGAMA already uses Kepler-compatible scaled energy and circularity coordinates. The improvement must be **validated table construction, refinement, and fallback**, rather than merely choosing those coordinates again. [`actions_spherical.cpp:347`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/actions_spherical.cpp#L347)

**C. Small action errors can matter, but their magnitude alone cannot establish fit bias.**

For your implemented exponential DF, at fixed amplitude,
\[
\delta\log f \simeq
-\alpha(L_c/J_0)^\alpha\,\delta L_c/L_c.
\]
Errors can therefore be amplified in the steep tail. A relative error in \(J_r\), especially near \(J_r=0\), is not directly a relative error in dispersion. [`regularized_df.py:438`](src/ocen_dm/kinematics/regularized_df.py:438)

The decisive test is:

1. In the **same frozen potential and DF**, replace only the action calculation; compare density, both pressures, and actual bin predictions.
2. Rebuild self-consistency with each backend and compare the full objective.
3. Refit representative starts and the BH-mass profile; measure parameter shifts and likelihood changes.

Use the joint prediction displacement
\[
D^2=\sum_i(\Delta\sigma_i/s_i)^2+
\sum_j(\Delta\mu_j)^2/\mu_j
\]
as an initial diagnostic, then evaluate the actual likelihood with its nuisance parameters treated consistently. The implemented objective combines kinematic residuals and count deviances; count normalizations and backgrounds may be profiled. [`df_fit.py:130`](src/ocen_dm/kinematics/df_fit.py:130), [`counts.py:155`](src/ocen_dm/kinematics/counts.py:155)

For scale: a coherent **dispersion** shift of \(0.1\%\) across 89 independent bins with \(2\%\) errors gives \(D^2\approx0.22\). That does not imply a \(0.1\%\) action error produces that shift.

**D. The main diagnostic traps are these.**

- **Mapper frequency agreement can hide wrong actions.** In Kepler the frequency ratio is independent of energy. In composites, verify the returned \(J_r(E,L)\) and actual angular momentum before trusting a bound mapper point. Your current fallback does not perform those checks. [`regularized_df.py:127`](src/ocen_dm/kinematics/regularized_df.py:127)
- **Export corruption extends beyond the BH.** This AGAMA version writes analytic component types without their parameters. Plummer reconstruction defaults include `mass=1`, `scaleRadius=1`. Reloading can therefore change remnants as well as remove the point-mass singularity. Reconstruct every static component from configuration and verify potential/force equality before comparison. [`potential_factory.cpp:995`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/potential_factory.cpp#L995), [`potential_factory.cpp:1032`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/potential_factory.cpp#L1032), [`potential_factory.cpp:158`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/potential_factory.cpp#L158)
- **`rmax` is not an orbit boundary.** Multipole evaluates an asymptotic continuation outside its grid. Direct integration there can accurately measure the continuation rather than establish its physical accuracy. Test radial-domain convergence separately. [`potential_multipole.cpp:1478`](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/src/potential_multipole.cpp#L1478)
- **Current validation shares the suspect backend.** Closure uses `self.af`, and `GalaxyModel` receives that same finder; direct projected moments therefore do not independently validate actions. The BH action regression checks `agama.actions` only for \(c>0.1\) and \(q<10^5\), excluding the reported difficult regime. [`regularized_df.py:605`](src/ocen_dm/kinematics/regularized_df.py:605), [`positive_df.py:338`](src/ocen_dm/kinematics/positive_df.py:338), [`test_regularized_df.py:318`](tests/test_regularized_df.py:318)

1. **Verified versus assumed.** I read the cited project and local AGAMA sources, checked that the relevant local AGAMA source files match those shipped with the installed package, and ran read-only numerical comparisons. Your radial quadrature recovered the sampled Kepler actions accurately and exposed the direct-routine eccentricity errors above. I did **not** reproduce the base/start-2 failures, orbit integrations, or full fits. I assume your reported observations and uncertainties are accurate. The available AGAMA environment lacks SciPy, so I did not run the project suite.

2. **Strongest argument against my recommendation.** All-direct evaluation may impose substantial cost for errors that largely cancel in integrated observables. A validated hybrid could deliver indistinguishable inference much faster. The missing evidence is the full prediction-and-refit comparison, not further isolated action accuracy.

3. **Checks needed for confidence.** Confirm the actual fitting binary/version; reproduce start-2 from the live potential; separate internal-potential error from action-table error; establish quadrature convergence in composite potentials; test finite-but-wrong mapper outputs; verify export/reconstruction force equality; and compare rebuilt models, refits, and BH profiles against the direct reference.