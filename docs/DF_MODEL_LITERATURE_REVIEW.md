# DF model choice for a stellar cluster embedded in dark matter

Literature review, 22 September 2026. This note reassesses the model family;
it does not replace the saved fits or change the running recovery batch.

The compact model proposal is now available as a [LaTeX note](compact_df_model.tex)
and [compiled PDF](compact_df_model.pdf), with equations, a conditional parameter
budget, three illustrations and the validation programme. It describes the proposed
model, not an implemented or validated replacement for the current mixture.

The note now distinguishes two formulations. The recommended first fit uses
one regularized exponential stellar DF with an anisotropy transition, plus
parametrized DM and optional remnant **density**
profiles, retaining stellar self-gravity. Explicit DFs for the dark components
are an alternative for stronger equilibrium checks and live N-body models;
they are not required to evaluate the stellar likelihood.

**The target is a self-gravitating stellar cluster with a possible extended DM
component, concentrated stellar remnants, and eventually a central black hole.**
Neither a stars-only GC nor a massless stellar tracer in a dominant halo is
an adequate default. The model must cover the transition between these limits.

The review supports retaining AGAMA and positive DFs, but revising the choice
of stellar family before treating the present 27-parameter mixture as our
inference model. The number of free parameters alone does not establish
overfitting. The stronger concern is that the current unlabelled mass/light
decomposition is weakly identified by the fitted observables.

The references below were checked against NASA ADS. The
[BibTeX file](df_model_literature.bib) contains 22 published records, including
the object-origin reference and the LIMEPY erratum; the
[metadata record](../results/diagnostics/df_literature_20260922/ads_metadata.json)
preserves the searches. The discussion distinguishes published constructions
from recommendations for this project.
The added Bekki & Freeman (2003) record is preserved in the
[object-reference metadata](../results/diagnostics/compact_df_writeup/object_reference_ads.json).
The three observational papers added after review have a separate
[ADS metadata record](../results/diagnostics/compact_df_writeup/review_observational_references_ads.json).

The closest precedents are these:

| Work | What it establishes | Consequence for this project |
|---|---|---|
| [Pascale, Binney, Nipoti & Posti (2019)](https://arxiv.org/abs/1904.08447), sections 3.2–3.3 and 4.1 | A compact stellar action DF can coexist with a halo DF and a central point mass in a jointly solved potential. Its GC examples include Omega Centauri. | This is the closest starting architecture. The Omega Cen example is a surface-brightness fit with illustrative kinematics, not a demonstrated joint fit to our present data or evidence for a halo. |
| [Pascale, Posti, Nipoti & Binney (2018)](https://academic.oup.com/mnras/article/480/1/927/5054060), sections 2–3 | Stars and DM each have an action DF; their densities and common potential are computed together. The observational likelihood uses individual velocities and foreground modelling. | A simple stellar family is compatible with a separate DM component. The density and anisotropy of either component depend on their common potential. |
| [Williams & Evans (2015)](https://arxiv.org/abs/1412.4640), section 3 | Action DFs can be designed around density asymptotes and an anisotropy transition with its own action scale. | A radial change in anisotropy can be built into one DF. Adding stellar basis components is not the only way to gain this freedom. |
| [Amorisco & Evans (2011)](https://arxiv.org/abs/1009.1813) and [Amorisco & Evans (2012)](https://arxiv.org/abs/1106.1062), section 3 of the latter | Lowered-isothermal and Michie–King stellar DFs can describe populations embedded in dark haloes. Distinct observed populations can jointly constrain a common potential. | Energy-based cluster DFs are also candidates for an embedded cluster. Their dwarf-galaxy tracer approximation must be replaced by a calculation including stellar self-gravity here. |
| [Jeffreson, Sanders, Evans et al. (2017)](https://arxiv.org/abs/1704.07833), section 3 | A photometrically motivated stellar family can be extended to flattened, rotating action-based cluster models. | Rotation and flattening deserve their own controlled extension. This paper's stars-only implementation cannot simply be reused unchanged after adding a halo. |
| [Vasiliev (2019)](https://arxiv.org/abs/1802.08239), sections 4–5 | AGAMA supplies analytic DFs, density-to-DF inversion, projected observables and iterative construction of composite equilibria. | The numerical infrastructure does not require our particular mixture. A spherical inversion model is a useful independent control. |
| [Binney & Vasiliev (2023)](https://arxiv.org/abs/2206.03523), sections 3.2 and 5.3 | Component DFs jointly generate gravity; the spheroidal DF construction also addresses problematic low-angular-momentum velocity distributions. | Physical component definitions and the detailed orbital prescription matter alongside positivity. Their Galactic component count is not a prescription for Omega Cen. |
| [Vasiliev, Feldmeier-Krause & Sormani (2026)](https://arxiv.org/html/2603.29502v1), sections II–III | A stellar action DF is fitted with a central BH and other gravitating components, using discrete LOS velocities and PMs. It treats missing coordinates, measurement errors and selection explicitly. | Positivity can be combined with full velocity-distribution information. Conditional velocity likelihoods require stated assumptions about velocity and line-of-sight selection. |
| [Binney (2026)](https://arxiv.org/html/2602.23127v1), sections 3–4 and 8 | Simple action combinations can produce cusps in velocity distributions near zero tangential/azimuthal velocity; derivative conditions and improved constructions address this. | Enforce the spherical derivative ratio of two in the finite-centre limit and use smooth contour labels for the stellar DF. This is separate from positivity and stability; saved fits retain their diagnostic role. |

Both formulations include all mass components in the total potential:

\[
\Phi_{\rm tot}=\Phi_\star+\Phi_{\rm rem}+\Phi_{\rm DM}+\Phi_\bullet,
\qquad
\rho_\star(\boldsymbol x)=\int
 f_\star[\boldsymbol J(\boldsymbol x,\boldsymbol v;\Phi_{\rm tot})]\,d^3v.
\]

The extended components jointly source Poisson's equation. A central black
hole contributes its point-mass potential. In the recommended starting model,
the dark densities are prescribed within each trial. In the all-component DF
alternative, the same density integral applies to DM and remnants; every DF
must be defined and checked in the **combined** potential. Holding a DF fixed in action space
as the potential changes generally changes its spatial density. A present-day
equilibrium constructed this way does not, by itself, specify the system's
assembly or stripping history. This interpretation follows the composite
construction of [Pascale et al. (2018)](https://academic.oup.com/mnras/article/480/1/927/5054060).

Using one effective stellar DF says nothing about whether a halo is present.
It is an assumption about the represented stellar population. A constant
stellar M/L ties stellar mass to stellar light; independent remnant and DM
components still allow the **total** mass-to-light ratio to vary with radius.
Conversely, multiple stellar populations become especially useful when their
membership, density or kinematics are constrained separately, as in
[Amorisco & Evans (2012)](https://arxiv.org/abs/1106.1062).

The primary spherical candidate is now a **regularized exponential action
DF**, adapted from the envelope of [Pascale et al. (2019)](https://arxiv.org/abs/1904.08447)
and the contour construction of [Binney (2026)](https://arxiv.org/html/2602.23127v1):

\[
f_\star=A_\star\exp[-(L_c/J_{0,\star})^\alpha],\qquad
\frac{dL}{dJ_r}=-g(J_r,L).
\]

Integrating to the circular-orbit axis assigns the contour label \(L_c\).
The DF is numerically normalized. Set \(q=J_r+L\),
\(c=L/(L+J_r)\), and

\[
g=g_H\exp[-b(q)\sin(\pi c/2)],\qquad
b(q)=b_{\rm out}\frac{(q/J_a)^2}{1+(q/J_a)^2}.
\]

An exact frequency ratio \(g_H=\Omega_r/\Omega_t\) gives an ergodic
reference; the epicycle interpolation in the LaTeX note is a possible
acceleration to validate. The frequency map and normalization are rebuilt
as the potential changes. This is an explicitly potential-dependent family,
not a fixed-action DF carried through adiabatic evolution. No implementation
or validation of this family is claimed yet.

The prescribed-density inference model has the following initial budget:

| Component | Free coordinates | Count |
|---|---|---:|
| Stars | \(M_\star,J_{0,\star},\alpha,b_{\rm out},J_a\) | 5 |
| DM | \(\rho_{20},r_s\) | 2 |
| Remnants | \(M_{\rm rem},a_{\rm rem}\) | 2 |
| **Spherical model with both dark components** | | **9** |

Light normalization (or stellar M/L) adds one coordinate; other observational
nuisance parameters are additional. The two anisotropy coordinates replace
one constant \(\eta\), increasing the stellar count by one. A two-transition
law adds two more coordinates if matched-tracer data require a turnover.
These action controls are not a prescribed \(\beta(r)\), and they still affect
the stellar density. Photometry and all velocity components must be fitted
jointly.

Compare cored and cusped halo profiles using \(\rho_{20}\) as the sampled
normalization, with the outer cutoff fixed within each run and varied between
sensitivity tests. This makes the prior on the target explicit but does not
remove prior sensitivity: uniform density and uniform log-density are different
choices. The no-DM comparison is separate from a bounded logarithmic prior.
Formal total halo mass is derived and can be dominated by extrapolation.

The DM component needs an independent density and outer truncation
prescription. An orbital prescription is needed if we construct a dark DF
for equilibrium certification or live evolution. It must not be treated as
a stellar mass bin subject to stellar equipartition, or forced to end at the
luminous cutoff. At fixed dark density, dark orbital distributions do not
change the smooth potential or the equilibrium stellar likelihood. AGAMA
explicitly supports density components without DFs
([Vasiliev 2019, section 5.2](https://academic.oup.com/mnras/article/482/2/1525/5114593)).

For the finite-centre no-BH baseline, the smoothness condition is concrete:
\(f_{J_r}/f_L\to2\) as \(L\to0\). The precursor
\(f\propto\exp[-((J_r+\eta L)/J_0)^\alpha]\) has ratio \(1/\eta\),
so only \(\eta=1/2\) meets that limit. It remains an analytic/reproduction
model, not the primary candidate. The regularized contour construction
enforces the ratio through \(g\to2\). Its joint action origin, global velocity
distributions, numerical closure and stability still require checks. A central
point mass changes the limiting analysis. The 2026 paper's suggestion about
radial-orbit instability is not a blanket rejection of every saved radially
biased model.

The required independent stellar control is an embedded generalized
lowered-isothermal/LIMEPY DF. A restricted inversion supplies an additional
diagnostic: a stellar density inferred from photometry can be inverted to a
DF in each trial total potential.
AGAMA's [spherical DF documentation](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/doc/reference.tex)
describes the Cuddeford–Osipkov–Merritt family, with
\(\beta(r)=(\beta_0+r^2/r_a^2)/(1+r^2/r_a^2)\).
The inversion can fail positivity and must then reject the trial; it cannot
represent arbitrary anisotropy turnovers. Both the inversion diagnostic and
the independent energy-based control must retain stellar self-gravity.

[Gieles & Zocchi (2015)](https://arxiv.org/abs/1508.02120), with their
[2018 erratum](https://doi.org/10.1093/mnras/stx3144), supplies a useful family
for stellar mass groups and energy truncation. Its radial-anisotropy restriction
and usual isolated-cluster setup need explicit consideration before applying
it to a stellar-plus-DM equilibrium. The Omega Cen studies by
[Zocchi, Gieles & Hénault-Brunet (2017)](https://arxiv.org/abs/1702.00725) and
[Zocchi et al. (2019)](https://doi.org/10.1093/mnras/sty1508) show why anisotropy
and concentrated stellar remnants deserve separate controls in central-mass
inference. These are methodological precedents, not adopted conclusions about
our data.

For a BH extension, we must let the stellar central structure respond.
[An & Evans (2006)](https://arxiv.org/abs/astro-ph/0511686) derive central
slope–anisotropy restrictions, including \(\gamma\geq\beta+1/2\) in the
Keplerian limit. An arbitrary fixed stellar core cannot simply be carried
unchanged into every BH model. These are necessary central conditions, not
a substitute for constructing the DF throughout phase space.

The recent [Pascale et al. (2026)](https://arxiv.org/html/2604.24855v1)
analysis of Draco and Ursa Minor is a useful counterexample to parameter
counting alone: its model has 30 parameters and uses chemically distinguished
stellar components. It assumes negligible stellar self-gravity, however.
That approximation cannot be imported into our compact stellar cluster.
The relevant lesson is to constrain population freedom with population data.

The observational priorities are supported by [Häberle et al. (2025)](https://arxiv.org/abs/2503.04903),
who measure radial variation in partial equipartition and a flattened dispersion
field, and [Vasiliev & Baumgardt (2021)](https://arxiv.org/abs/2102.09568),
who find an outer anisotropy transition in Gaia EDR3. The samples and radial
ranges differ. [Pechetti et al. (2024)](https://arxiv.org/abs/2401.15149)
report central counter-rotation. The even/odd axisymmetric DF decomposition
is useful, but the even part fixes raw second moments; dispersions about the
mean also change when the odd part changes streaming.

For our existing experiment, the source code uses three double-power-law
stellar action DFs without an explicit action cutoff. The generating mocks
instead have finite energy boundaries. A need for several basis terms may
therefore reflect approximation of the chosen mock shapes. This is a
hypothesis to test, not an established explanation for the residuals.
The [capacity challenge](DF_CAPACITY_CHALLENGE.md) measured approximation
accuracy in a supplied potential; it did not select the statistically preferred
family for the real data. The [mass-recovery experiment](DF_MASS_RECOVERY.md)
also gives the basis terms independent mass and light weights. Its very good
observable fits with different mass decompositions are grounds for testing
identifiability, rather than declaring either that DM is unconstrained by the
real data or that 27 parameters are inherently excessive.

My recommendation is a bounded family comparison before another production
inference run:

1. Compare the regularized exponential action DF and the required embedded
   generalized lowered-isothermal control,
   retaining stellar self-gravity, optional DM and remnants in both.
2. Test observable adequacy and recovery of \(\rho_{\rm DM}(20\,\mathrm{pc})\)
   on independent stars-only and stars-plus-DM truths, including shape mismatch
   and known M/L differences. A correct-data-model control separates optimizer
   failure from family mismatch; fitting the same family alone is insufficient.
3. Include one anisotropy transition from the outset. Add a second transition
   only if matched-tracer profiles require it, and demonstrate the resulting
   physical anisotropy profile rather than infer it from input signs.
4. Audit and match HST/MUSE/Gaia stellar-mass and luminosity selections, or
   require a two-tracer recovery test. Test rotation and flattening before
   interpreting a spherical result as a final DM constraint. Count stellar
   mass only once when assigning catalogue-specific tracer weights.
5. Check numerical closure and stellar velocity distributions. For selected
   fitted dark density profiles, investigate whether positive dark DFs exist
   in the combined potential, and construct velocity distributions for live
   dynamical tests. A failed isotropic inversion does not exclude every
   anisotropic realization. These additional consistency and stability checks
   are distinct from evaluating the initial stellar likelihood; the existing
   fitted dark density templates have not been certified as equilibrium DFs.

The appropriate stopping criterion for this comparison is reliable mass
recovery with adequate predictions at the data's precision. Near-exact
reproduction of a noiseless mock is useful as a numerical diagnostic, but is
not sufficient to choose the inference model.
