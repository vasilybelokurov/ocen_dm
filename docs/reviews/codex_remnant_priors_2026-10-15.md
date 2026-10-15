Research mode used: ordinary Codex web research fallback
Native Deep Research: unavailable from this invocation path
Scientific source rules: yes
Model effort: high   Date: 2026-10-15

---

## 1. Executive summary

- **Total remnant mass and your added Plummer mass are different quantities.** A substantial WD population is already included in a dynamically inferred light-following \(M_\star\). Applying a prior for *all* WD+NS+BH mass directly to your additional component would double-count some remnants.
- **A canonical IMF gives an old population a remnant mass fraction of roughly 25–40% before preferential stellar escape**, depending especially on BH formation and retention. WDs dominate. Dynamical evolution can raise the fraction substantially: published multimass fits include clusters with approximately 50–70% total remnant mass. These are model-dependent population estimates, not universal fractions. [Dickson et al. 2023](https://arxiv.org/abs/2303.01637), [Baumgardt & Makino 2003](https://arxiv.org/abs/astro-ph/0211471).
- **For ω Cen, several independent models support a substantial stellar-mass BH population.** Published estimates put approximately 4.5–6% of the total mass in BHs; July–September 2026 preprints obtain approximately 8% and \(\gtrsim7\%\), respectively. This supports a broad prior, rather than adopting one published posterior as a precise constraint. [Zocchi et al. 2019](https://arxiv.org/abs/1806.02157), [Baumgardt et al. 2019](https://arxiv.org/abs/1907.10845), [Dickson et al. 2024](https://arxiv.org/abs/2308.13037), [Dickson et al. 2026a,b](https://arxiv.org/abs/2607.07010).
- **Your action-DF component, \(2.7\times10^5\,M_\odot\), is physically sensible.** The Jeans component, \(5.4\times10^5\,M_\odot\), is also possible as a mixed non-light-following component, but is more demanding if interpreted entirely as BHs. Its mass is not ruled out by a universal \(5\times10^5\,M_\odot\) ceiling: no such robust ceiling exists.
- **My practical default for your extra component is**
  \[
  q\equiv M_{\rm rem}/M_\star\sim {\rm HalfNormal}(0.10),
  \qquad 0\le q\le0.30,
  \]
  with a broad scale prior centred around \(a/r_{h,\star}\sim0.3\). These are recommended regularizing priors, **not published confidence limits**. A separate, broader analysis should allow a mixed WD/faint-star component and alternative formation histories.

## 2. Detailed findings

### Mass bookkeeping: the essential distinction

Define
\[
M_{\rm live}=\text{mass in surviving non-remnant stars},
\]
including faint main-sequence stars, and
\[
M_R=M_{\rm WD}+M_{\rm NS}+M_{\rm BH}.
\]
For a cluster without a separately counted IMBH or dark-matter component,
\[
M_{\rm cl}=M_{\rm live}+M_R,\qquad
f_R=\frac{M_R}{M_{\rm cl}},\qquad
\frac{M_R}{M_{\rm live}}=\frac{f_R}{1-f_R}.
\]

Consequently, total remnant fractions of 30%, 50% and 70% correspond to remnant-to-living-star mass ratios of 0.43, 1 and 2.33. Here “living” is a population definition; it does not mean only the bright stars entering a kinematic catalogue.

Your decomposition is instead
\[
\rho_{\rm tot}(r)=M_\star\,\widehat{\rho}_{L}(r)
                  +M_{\rm rem}\,\widehat{\rho}_{P}(r)
                  +M_\bullet\delta(\boldsymbol r).
\]
Because \(M_\star\) includes remnants assigned to the light-following distribution,
\[
M_{\rm rem}\ne M_R.
\]

More precisely, the added component measures mass that the adopted light-following profile does not accommodate. Calling it “remnants” is an interpretation requiring population and spatial-distribution assumptions. This distinction is also explicit in the multicomponent Jeans analysis of Bañares-Hernández et al.: their photometric mass component includes remnants traced by that distribution. [Bañares-Hernández et al. 2025, §2](https://arxiv.org/html/2408.00939v4).

The physical connection to the IMF is
\[
M_i(t)=A\int \xi(m)\,
       m_{i,\rm final}(m,Z,Y)\,
       P_{i,\rm bound}(m,t)\,dm,
\]
where the integration domain selects progenitors of remnant type \(i\). The corresponding surviving-star integral fixes \(M_{\rm live}\). The normalization \(A\) is shared: remnant mass cannot be chosen independently of the progenitor population unless the formation and stripping history is also allowed to change.

The retention factor must distinguish natal escape, later dynamical ejection, and tidal stripping. It is generally mass-dependent.

### Canonical IMF: a numerical baseline

There is no unique “12-Gyr canonical remnant fraction” without specifying the IMF mass limits, IFMR, kick prescription and subsequent escape.

Dickson et al. explicitly state that WDs contribute approximately **30% of the total mass at approximately 10 Gyr for a canonical IMF**, before preferential loss of low-mass objects. [Dickson et al. 2023, introduction](https://arxiv.org/html/2303.01637v2).

For a transparent calculation closer to your requested age, I integrated a Kroupa IMF with:

- \(dN/dm\propto m^{-1.3}\) from \(0.1\)–\(0.5\,M_\odot\), and \(m^{-2.3}\) from \(0.5\)–\(100\,M_\odot\);
- turnoff mass \(0.8\,M_\odot\);
- WD progenitors \(0.8\)–\(8\,M_\odot\), with \(m_{\rm WD}=0.109m+0.394\);
- \(1.4\,M_\odot\) NSs from \(8\)–\(20\,M_\odot\);
- an illustrative mean BH mass \(15\,M_\odot\) above \(20\,M_\odot\).

The IMF and empirical WD relation come from [Kroupa 2001](https://arxiv.org/abs/astro-ph/0009005) and [Kalirai et al. 2008](https://arxiv.org/abs/0706.3894). The progenitor boundaries and fixed BH mass are simplifying assumptions, not a modern low-metallicity population-synthesis calculation.

**The following numbers are my calculation, not quotations from those papers.**

| Illustrative population at 12 Gyr | WD/total | NS/total | BH/total | All remnants/total | WD/living | NS/living | BH/living |
|---|---:|---:|---:|---:|---:|---:|---:|
| All remnants retained; no preferential escape | 22.6% | 2.0% | 8.0% | 32.6% | 33.5% | 2.91% | 11.9% |
| 10% NS retention; 50% BH retention | 24.0% | 0.21% | 4.3% | 28.4% | 33.5% | 0.29% | 6.0% |
| Same retained remnants, but half the living-star mass removed | 37.3% | 0.32% | 6.6% | 44.3% | 66.9% | 0.58% | 11.9% |

The final row isolates the denominator effect; actual evolution also removes WDs and BHs.

The useful baseline is therefore:

- **WDs:** approximately 20–30% of total retained mass without strong preferential stellar depletion.
- **NSs:** approximately a few per cent if all were retained; usually much less than 1% after realistic kicks.
- **BHs:** a few per cent, potentially approaching 10% with favourable low-\(Z\) formation and retention.
- **All remnants:** approximately 25–40%, with substantial prescription dependence.

A Chabrier IMF provides another canonical baseline, but its precise normalization and single-star/system convention should be specified rather than treating it as numerically identical to Kroupa. [Chabrier 2003](https://arxiv.org/abs/astro-ph/0304382).

### What dynamical evolution changes

**Preferential stellar loss can raise the remnant fraction without increasing remnant mass.** Conversely, BH binary interactions can reduce the BH inventory while heating and expanding the luminous cluster.

Baumgardt & Makino’s stellar-evolution-only fitting formula gives \(f_{\rm WD+NS}\simeq0.283\) at 12 Gyr when evaluated directly. Their dissolving models can become approximately 90% remnants near destruction. However, their IMF ended at \(15\,M_\odot\), and they retained all NSs: this is evidence for the dynamical trend, not a current BH-inclusive population prescription. [Baumgardt & Makino 2003, equations 2 and 15, §3.4](https://arxiv.org/pdf/astro-ph/0211471).

Present-day mass functions support preferential low-mass-star depletion. Sollima & Baumgardt find correlations between MF slopes, relaxation time and dark-remnant fraction; their companion \(N\)-body analysis finds that Kroupa/Chabrier starting IMFs can reproduce the observed MFs. [Sollima & Baumgardt 2017](https://arxiv.org/abs/1708.09529), [Baumgardt & Sollima 2017](https://arxiv.org/abs/1708.09530).

High BH fractions are also possible without a top-heavy IMF. Gieles et al.’s best Pal 5 model has **22% of its present bound mass in BHs**, following loss of approximately 92% of its initial mass. This is a best-fitting evolutionary model, not a confidence interval or a direct BH census. It demonstrates why a universal 10% upper bound is false; it does not establish that ω Cen followed Pal 5’s history. [Gieles et al. 2021](https://arxiv.org/html/2102.11348).

There is an additional complication: BH heating changes WD escape as well as stellar escape. Rostami-Shirazi et al.’s direct simulations find that reproducing the high inferred remnant fractions requires low BH natal kicks in their model grid. Thus “preferential loss always retains every WD” is too simple. [Rostami-Shirazi et al. 2025; published online 27 November 2024](https://doi.org/10.1093/mnras/stae2644).

### Natal kicks and retention

For isolated remnants initially at rest, a Maxwellian kick distribution gives
\[
P_{\rm ret}(v_{\rm esc},\sigma_k)=
 {\rm erf}\!\left(\frac{v_{\rm esc}}{\sqrt2\sigma_k}\right)
-\sqrt{\frac2\pi}\frac{v_{\rm esc}}{\sigma_k}
 e^{-v_{\rm esc}^2/(2\sigma_k^2)}.
\]

Hobbs et al. infer a Maxwellian with one-dimensional dispersion \(265\ {\rm km\,s^{-1}}\) for young pulsar birth velocities. This is an observational result for that pulsar sample, not a universal prescription for every NS formation channel. [Hobbs et al. 2005](https://arxiv.org/abs/astro-ph/0504584).

**Calculated retention probabilities under the simplified formula:**

| Kick dispersion | \(v_{\rm esc}=60\) km/s | 70 km/s | 120 km/s | 140 km/s |
|---|---:|---:|---:|---:|
| 265 km/s | 0.30% | 0.48% | 2.3% | 3.6% |
| 50 km/s | 30% | 42% | 88% | 95% |
| 20 km/s | 97% | 99% | approximately 100% | approximately 100% |

Real calculations must include progenitor velocities, birth positions, binary disruption and binary centre-of-mass recoil.

For NSs, binary evolution and low-kick channels matter enormously. Pfahl et al. find that conventional assumptions give retention below approximately 5% in massive clusters, motivating alternative channels. Modern CMC models distinguish core-collapse NS kicks, \(\sigma=265\ {\rm km\,s^{-1}}\), from electron-capture channels, \(\sigma=20\ {\rm km\,s^{-1}}\). [Pfahl et al. 2002](https://arxiv.org/abs/astro-ph/0106141), [Kremer et al. 2020, §2](https://arxiv.org/pdf/1911.00018).

For BHs, commonly used prescriptions reduce kicks through fallback:
\[
v_{k,\rm BH}=(1-f_{\rm fallback})\,v_{k,\rm NS}.
\]
Lower metallicity changes winds, final masses and fallback; the BH mass spectrum and retention must therefore be calculated together. The result depends on the explosion prescription, rather than metallicity supplying a unique retention fraction. [Fryer et al. 2012](https://arxiv.org/abs/1110.1726).

For ω Cen, Baumgardt et al. quote a present central escape speed of **63 km/s**. Their population calculation uses an initially more massive, smaller cluster and finds that the approximately 4.6% present BH fraction is compatible with several stellar-evolution/kick combinations, after accounting for dynamical ejection. [Baumgardt et al. 2019, §4.1.3](https://arxiv.org/pdf/1907.10845).

**Prior implication:** 100% NS retention is an extreme assumption; a fixed 10% NS retention is a useful approximation, not an observation. High BH natal retention in ω Cen is plausible, but “natal retention”, “fraction surviving to today”, and “mass fraction in BHs today” must never be interchanged.

### Which models impose an IMF-derived inventory?

**LIMEPY itself does not impose an IMF.** It needs component masses and normalizations; the surrounding population model determines whether these are constrained or free.

| Model/application | How remnant content enters |
|---|---|
| **Hénault-Brunet et al. 2019** | Method comparison using mock **M4**, not a 47 Tuc measurement. Includes flexible remnant-MF and dark-mass parameters; therefore not exclusively an IMF-imposed inventory. The correct citation is **MNRAS 483, 1400**, not 482. [Paper](https://arxiv.org/abs/1810.05167). |
| **Hénault-Brunet et al. 2020: 47 Tuc** | Broken stellar MF with fitted high-mass slope; stellar evolution/IFMR generates remnants. NS retention fixed at 10%; BH retention fitted, with massive BHs removed first to represent dynamical ejection. [Paper, §3.2](https://arxiv.org/html/1908.08538v2). |
| **Dickson et al. 2023/2024** | Fitted MF slopes below 0.5, between 0.5–1, and above \(1\,M_\odot\); metallicity-dependent IFMRs; natal-kick prescription and fitted present BH mass-retention parameter. Spatial segregation is tied to component mass through \(s_j\propto m_j^{-\delta}\). [Paper I](https://arxiv.org/html/2303.01637v2). |
| **Zocchi et al. 2019: ω Cen** | Two-component DF: BHs versus other stars/remnants. Fits BH mass fraction and object-mass ratio; approximately 5% is subsequently checked against population expectations. It is not fixed by an explicit full IMF. [Paper, §2.2](https://arxiv.org/pdf/1806.02157). |
| **Baumgardt et al. 2019: ω Cen** | Evolved \(N\)-body models: IMF and stellar evolution supply remnants; initial BH retention is varied. BH ejection and segregation then occur dynamically. [Paper](https://arxiv.org/pdf/1907.10845). |
| **CMC/Kremer models** | IMF, stellar/binary evolution and kicks generate remnants; encounters, escape and tidal evolution determine survival. These are evolutionary predictions conditional on initial conditions and prescriptions. [Catalogue paper](https://arxiv.org/abs/1911.00018). |
| **Dickson et al. 2026a,b — preprints** | Couples multimass equilibrium models to clusterBH evolution. Removes the earlier independent BH-retention parameter: dynamical BH depletion is coupled to cluster evolution and energy demand. [July preprint](https://arxiv.org/html/2607.07010), [September preprint](https://arxiv.org/html/2609.10955). |

Peuten et al. validate multimass lowered-isothermal models against \(N\)-body populations with different remnant retention. Their success supports the DF machinery; it does not establish a universal remnant abundance. [Peuten et al. 2017](https://arxiv.org/abs/1702.01712).

Weatherford et al. instead infer BH populations from the anticorrelation between observable stellar mass segregation and retained BHs, calibrated against simulations. Their 2020 analysis **excluded ω Cen because of insufficient radial field coverage**. It should not be quoted as an independent ω Cen BH measurement. [Weatherford et al. 2018](https://arxiv.org/abs/1712.03979), [Weatherford et al. 2020, §4](https://arxiv.org/pdf/1911.09125).

### Verified numerical results for 47 Tuc and ω Cen

**47 Tuc is a useful demonstration of WD-rich, BH-poor composition, but is more metal-rich and dynamically different from ω Cen.**

Hénault-Brunet et al. infer
\[
M_{\rm BH}=430^{+386}_{-301}\,M_\odot
\]
for 47 Tuc, with quoted errors spanning the 16th–84th posterior percentiles. They also find approximately \(6\times10^4\,M_\odot\) in WDs individually more massive than \(1\,M_\odot\). Thus a small BH fraction does not mean little centrally important remnant mass. Their BH masses depend strongly on the assumed preferential ejection of massive BHs. [Hénault-Brunet et al. 2020, §§5.1.2–5.1.3](https://arxiv.org/html/1908.08538v2).

For ω Cen:

| Study | Reported result |
|---|---|
| Zocchi et al. 2019 | Anisotropic two-component fit: approximately **4.5%**, \(1.35\times10^5\,M_\odot\), in the heavy component; a fitted model value, not a quoted confidence bound. [§4.2](https://arxiv.org/pdf/1806.02157). |
| Baumgardt et al. 2019 | Best-fitting evolved model: approximately **4.6%**, \(1.65\times10^5\,M_\odot\), in stellar-mass BHs. [§4.1](https://arxiv.org/pdf/1907.10845). |
| Dickson et al. 2024 | \(f_{\rm BH}=5.9^{+0.2}_{-0.2}\%\), \(M_{\rm BH}=182000^{+5000}_{-6000}\,M_\odot\); total remnant fraction \(50.3^{+0.9}_{-0.3}\%\). All are approximately 68% posterior intervals. [Table 1](https://arxiv.org/html/2308.13037v2). |
| Dickson et al. 2026a,b | July preprint: approximately **8.3%** BH mass fraction; September preprint: \(\gtrsim7\%\) under each explored prescription. These are characteristic fitted values, not confidence lower limits. [July, §IV.3](https://arxiv.org/html/2607.07010), [September, Fig. 9](https://arxiv.org/html/2609.10955). |

Dickson et al.’s 2024 model also gives 47 Tuc a total remnant fraction \(46.9\pm0.3\%\), despite only approximately 0.046% in BHs. For ω Cen, subtracting its BH contribution leaves approximately **44% in WDs+NSs**—my arithmetic, not an independently measured WD fraction. The authors caution that their statistical errors underestimate total uncertainty; mock tests suggest BH-fraction errors can need inflation by approximately 2.5. [Dickson et al. 2024, Table 1 and §3](https://arxiv.org/html/2308.13037v2).

The 2026 evolution-coupled models explicitly acknowledge additional systematic uncertainty and cannot include an IMBH self-consistently. They therefore update the plausible BH scale without settling the central mass decomposition. [Dickson et al. 2026a, §IV.3](https://arxiv.org/html/2607.07010).

### White dwarfs: abundant does not mean more concentrated than turnoff stars

The relevant comparison is object mass, not luminosity.

A \(0.55\)–\(0.60\,M_\odot\) WD is:

- heavier than the typical low-mass main-sequence star;
- lighter than a normal approximately \(0.8\,M_\odot\) turnoff star;
- therefore expected, after relaxation, to be more concentrated than low-mass stars but **less concentrated than turnoff/giant tracers**, other things being equal.

This is an inference from mass-dependent relaxation, not a claim that all WD populations have identical profiles. The illustrative IMF integral above gives a number-weighted mean WD mass \(0.59\,M_\odot\); old massive WDs form an important tail.

There is direct observational support for the change in spatial distribution after mass loss: Heyl et al. measure young WDs diffusing outward through 47 Tuc’s core. [Heyl et al. 2015](https://arxiv.org/abs/1502.01890).

ω Cen is additionally complicated by helium-rich populations. Bellini et al. interpret its split bright WD sequence using approximately \(0.55\,M_\odot\) normal-population WDs and approximately \(0.46\,M_\odot\) products associated with helium-rich populations. This is a stellar-evolution interpretation of the sequences, not a global dynamical WD census. [Bellini et al. 2013](https://arxiv.org/pdf/1305.0265).

**For your models:** a large WD budget does not automatically justify placing that entire budget in a compact positive Plummer addition to a turnoff-star light profile. Massive WDs, NSs and BHs can supply a concentrated contribution; ordinary low-mass WDs may follow a similar or more extended distribution.

### Upper limits and unusual formation histories

Remnant mass cannot be arbitrary at fixed IMF, surviving-star normalization and retention history. However, **there is no formation-history-independent maximum fraction well below unity**: preferential stripping can leave remnant-dominated survivors.

Top-heavy IMFs can increase remnant production, but also increase early mass loss, expansion and disruption. Haghi et al. explicitly model this coupling; one cannot increase the progenitor count while keeping the subsequent cluster history unchanged. [Haghi et al. 2020](https://arxiv.org/abs/2012.07095).

The counter-evidence is that Dickson et al.’s multimass fits favour high-mass IMF slopes compatible with canonical values and find no significant metallicity trend. This disfavors using a strongly top-heavy IMF as the default explanation for present GC mass budgets, conditional on their remnant prescriptions. [Dickson et al. 2023, §6.1.1](https://arxiv.org/html/2303.01637v2).

Multiple populations and a nuclear-cluster origin broaden the possible histories, but do not supply a free mass multiplier:

- each population has its own surviving stars and remnants;
- helium changes lifetimes and IFMRs;
- stripping a host galaxy does not deposit all its remnants in the nucleus;
- inward transport, capture or cluster mergers require an explicit dynamical history.

**My assessment for ω Cen:** total remnant fractions around **30–55%** are straightforward to entertain; approximately **60–70%** is a useful generous sensitivity range. Fractions above that are not mathematically forbidden, but deserve explicit extreme-IMF or selective-stripping models rather than ordinary prior weight. These ranges are my synthesis, not published confidence limits.

Evans et al.’s often-repeated approximately \(5\times10^5\,M_\odot\) remnant allowance is a model-based estimate using selected CMC populations, **not a universal upper bound**. It also cannot be transferred directly to your extra component without reconciling the luminous-mass definition. [Evans, Strigari & Zivick 2022, discussion](https://academic.oup.com/mnras/article/511/3/4251/6517696).

### Assessment of your two fits

For a Plummer sphere,
\[
r_{h,\rm rem}=\frac{a}{\sqrt{2^{2/3}-1}}=1.305a.
\]

Using your inputs:

| Quantity | Action-DF | Jeans |
|---|---:|---:|
| \(M_{\rm rem}/M_\star\) | 0.090 | 0.200–0.225 |
| Extra-component fraction of total, including the stated point mass | 8.26% | 16.45–18.10% |
| \(r_{h,\rm rem}\) | 2.87 pc | 6.92 pc |
| \(r_{h,\rm rem}/r_{h,\star}\), taking \(r_{h,\star}=10\) pc | 0.287 | 0.692 |

Thus 16.5% corresponds approximately to the upper end of your Jeans \(M_\star\) range, rather than the whole range.

**Action-DF:** \(2.7\times10^5\,M_\odot\) at \(a=2.2\) pc is a credible extra remnant component. Its mass and scale are close to the extended central component preferred by Bañares-Hernández et al., who report approximately \(2\)–\(3\times10^5\,M_\odot\) and \(a=1.5\)–2.2 pc at their stated \(3\sigma\) level. Agreement is reassuring but not independent validation if observational inputs overlap. [Bañares-Hernández et al. 2025, §6](https://arxiv.org/html/2408.00939v4).

**Jeans:** \(5.4\times10^5\,M_\odot\) is physically possible as a mixed additional component. As **BH mass alone**, a 16–18% fraction is substantially above existing ω Cen estimates and requires unusually favourable production, retention or stellar depletion. Its broader profile can also absorb radial population/MF changes that a constant light-following \(M/L\) misses. That possibility should be investigated before assigning the whole mass to BHs.

The scale \(r_{h,\rm rem}/r_{h,\star}\simeq0.69\) is not by itself unphysical for a massive BH subsystem. Balanced-evolution theory predicts subsystem size to depend on both BH mass fraction and object-mass ratio; large BH populations need not be extremely compact. [Breen & Heggie 2013, equation 4](https://doi.org/10.1093/mnras/stt628).

### Priors I would use

**Preferred physical formulation:** normalize a population model to star counts, sample IMF/IFMR/kick and depletion uncertainties, and predict both the remnant inventory and its spatial distribution. Constrain the light-following \(M/L\) jointly. This prevents an unconstrained exchange between “stellar” and “remnant” mass.

For your present two-component parameterization, I recommend the following **explicitly phenomenological** priors.

**Extra mass:**
\[
q=\frac{M_{\rm rem}}{M_\star}
\sim {\rm HalfNormal}(\sigma_q=0.10),
\quad 0\le q\le0.30,
\qquad M_{\rm rem}=qM_\star.
\]

This gives approximately 95% prior probability below \(q\simeq0.19\), after truncation. At \(M_\star=2.4\)–\(3.0\times10^6\,M_\odot\), that is approximately \(4.6\)–\(5.8\times10^5\,M_\odot\). The upper support is \(7.2\)–\(9.0\times10^5\,M_\odot\).

It favours the literature-scale addition while retaining a tail for the Jeans result. **The cutoff is a modelling choice, not a physical theorem.** For a broader mixed-component sensitivity analysis, use \(\sigma_q=0.15\), \(q_{\max}=0.5\).

**Scale**, writing \(x=a/r_{h,\star}\):
\[
\ln x\sim\mathcal N(\ln0.30,\;0.8^2).
\]

For a BH-dominated interpretation I would truncate approximately to
\[
0.03<x<0.77,
\]
corresponding to \(r_{h,\rm rem}<r_{h,\star}\). This retains both fitted scales. For a mixed WD/faint-star residual, extend the upper support to \(x\simeq1\), while examining degeneracy with \(M_\star\). The numerical form and bounds are my recommendation, informed by segregation and subsystem physics—not calibrated population constraints.

**Total remnant inventory, if separately modelled:** a broad starting prior such as
\[
f_R\sim\mathcal N(0.40,0.12^2),\qquad0.15<f_R<0.70
\]
is defensible for sensitivity work on ω Cen. It must be imposed on **embedded plus additional remnants together**, and cannot substitute for a population model.

A continuous half-normal includes arbitrarily small masses but assigns zero probability to exactly zero. For formal component-existence comparisons, fit a separate \(M_{\rm rem}=0\) model or use a spike-and-slab prior.

## 3. Disagreements, uncertainties and open questions

1. **A large remnant population does not settle the IMBH question.** Häberle et al. report central fast stars and an escape-speed lower bound of approximately \(8200\,M_\odot\); their velocity-plus-acceleration analysis gives \(M_\bullet>21100\,M_\odot\) at 99% confidence under their positional assumptions. [Häberle et al. 2024](https://www.nature.com/articles/s41586-024-07511-z). Bañares-Hernández et al.’s combined Jeans/pulsar analysis instead gives \(M_\bullet<6000\,M_\odot\) at stated \(3\sigma\). These are materially inconsistent conditional results, not a combined allowed interval. [Bañares-Hernández et al. 2025](https://arxiv.org/abs/2408.00939).

2. **Your \(4.3\times10^4\,M_\odot\) point mass is not forbidden by remnant supply.** Evolutionary models starting with IMBH seeds can grow approximately \(5\times10^4\,M_\odot\) objects while matching ω Cen’s broad structure. That establishes a possible formation channel, not a measurement. [González Prieto, Rodriguez & Cabrera 2025](https://arxiv.org/abs/2507.06316). The 2026 timing-only preprint gives a much weaker \(M_\bullet<10^5\,M_\odot\) limit at 90% confidence; it does not reproduce the stringent combined-analysis exclusion. [Colom i Bernadich et al. 2026, preprint](https://arxiv.org/abs/2603.21845).

3. **Individual BH detections establish existence, not the total inventory.** Whitaker et al. report an astrometric companion mass \(4.46^{+1.22}_{-1.01}\,M_\odot\), with 16th–84th percentile errors, assuming a \(0.78\,M_\odot\) visible star. They favour a BH interpretation while discussing rare alternatives. This also shows that low metallicity does not require every retained BH to be massive. [Whitaker et al. 2026](https://arxiv.org/html/2606.18350).

4. **The dominant uncertainties are systematic:** BH IFMR and kicks, initial density, primordial binaries, rotation, helium-rich populations, mass-dependent escape, and the translation from observable tracer density to total stellar density. Agreement among models sharing stellar-evolution prescriptions is not fully independent evidence.

5. **Weight of evidence:** a WD-rich cluster with a substantial retained BH population is well supported as a physically viable ω Cen model. A uniquely determined BH fraction, a universal \(5\times10^5\,M_\odot\) remnant ceiling, or a settled exclusion/detection of the IMBH is not supported. For your fits, the decisive unresolved issue is how much radial \(M/L\) variation the extra Plummer component is representing.

## 4. Source list

All entries below are primary research. Journal papers are peer-reviewed; the three specified 2026 modelling/timing entries remain preprints in the records checked.

| Source | Publication and identifiers |
|---|---|
| Kroupa — *On the variation of the IMF* | 2001, MNRAS 322, 231. [astro-ph/0009005](https://arxiv.org/abs/astro-ph/0009005); DOI **10.1046/j.1365-8711.2001.04022.x** |
| Chabrier — *Galactic Stellar and Substellar IMF* | 2003, PASP. [astro-ph/0304382](https://arxiv.org/abs/astro-ph/0304382) |
| Kalirai et al. — *Initial–Final Mass Relation* | 2008, ApJ 676, 594. [0706.3894](https://arxiv.org/abs/0706.3894); DOI **10.1086/527028** |
| Baumgardt & Makino — *Dynamical evolution…in tidal fields* | 2003, MNRAS 340, 227. [astro-ph/0211471](https://arxiv.org/abs/astro-ph/0211471); DOI **10.1046/j.1365-8711.2003.06286.x** |
| Sollima & Baumgardt — *Global mass functions…I* | 2017, MNRAS. [1708.09529](https://arxiv.org/abs/1708.09529) |
| Baumgardt & Sollima — *Global mass functions…II* | 2017, MNRAS 472, 744. [1708.09530](https://arxiv.org/abs/1708.09530); DOI **10.1093/mnras/stx2036** |
| Hobbs et al. — *233 pulsar proper motions* | 2005, MNRAS 360, 974. [astro-ph/0504584](https://arxiv.org/abs/astro-ph/0504584); DOI **10.1111/j.1365-2966.2005.09087.x** |
| Pfahl, Rappaport & Podsiadlowski — *NS retention* | 2002, ApJ 573, 283. [astro-ph/0106141](https://arxiv.org/abs/astro-ph/0106141) |
| Fryer et al. — *Compact Remnant Mass Function* | 2012, ApJ 749, 91. [1110.1726](https://arxiv.org/abs/1110.1726); DOI **10.1088/0004-637X/749/1/91** |
| Peuten et al. — *Testing lowered isothermal models…II* | 2017, MNRAS 470, 2736. [1702.01712](https://arxiv.org/abs/1702.01712); DOI **10.1093/mnras/stx1311** |
| Hénault-Brunet et al. — *Mass modelling…method comparison* | 2019, MNRAS **483**, 1400; online 2018. [1810.05167](https://arxiv.org/abs/1810.05167); DOI **10.1093/mnras/sty3187** |
| Hénault-Brunet et al. — *BH content and IMF of 47 Tuc* | 2020, MNRAS 491, 113; online 2019. [1908.08538](https://arxiv.org/abs/1908.08538); DOI **10.1093/mnras/stz2995** |
| Dickson et al. — *Multimass modelling…I* | 2023, MNRAS 522, 5320. [2303.01637](https://arxiv.org/abs/2303.01637); DOI **10.1093/mnras/stad1254** |
| Dickson et al. — *Multimass modelling…II* | 2024, MNRAS 529, 331; published 14 February. [2308.13037](https://arxiv.org/abs/2308.13037); DOI **10.1093/mnras/stae470** |
| Zocchi, Gieles & Hénault-Brunet — *Stellar-mass BHs and ω Cen kinematics* | 2019, MNRAS 482, 4713. [1806.02157](https://arxiv.org/abs/1806.02157) |
| Baumgardt et al. — *No evidence for IMBHs in ω Cen and NGC 6624* | 2019, MNRAS 488, 5340. [1907.10845](https://arxiv.org/abs/1907.10845) |
| Weatherford et al. — *Predicting stellar-mass BH populations* | 2018, ApJ. [1712.03979](https://arxiv.org/abs/1712.03979) |
| Weatherford et al. — *Dynamical survey…50 GCs* | 2020, ApJ. [1911.09125](https://arxiv.org/abs/1911.09125); DOI **10.3847/1538-4357/ab9f98** |
| Kremer et al. — *CMC Cluster Catalog* | 2020, ApJS 247, 48. [1911.00018](https://arxiv.org/abs/1911.00018); DOI **10.3847/1538-4365/ab7919** |
| Gieles et al. — *Supra-massive BH population in Pal 5* | 2021, Nature Astronomy 5, 957. [2102.11348](https://arxiv.org/abs/2102.11348); DOI **10.1038/s41550-021-01392-2** |
| Breen & Heggie — *BH subsystems in idealized star clusters* | 2013, MNRAS 432, 2779. [Publisher](https://doi.org/10.1093/mnras/stt628); DOI **10.1093/mnras/stt628** |
| Heyl et al. — *A Measurement of Diffusion in 47 Tucanae* | 2015, ApJ 804, 53. [1502.01890](https://arxiv.org/abs/1502.01890); DOI **10.1088/0004-637X/804/1/53** |
| Bellini et al. — *Double WD Cooling Sequence in ω Cen* | 2013, ApJL 769, L32. [1305.0265](https://arxiv.org/abs/1305.0265); DOI **10.1088/2041-8205/769/2/L32** |
| Haghi et al. — *Lifetimes…Top-heavy IMF* | 2020, ApJ 904, 43. [2012.07095](https://arxiv.org/abs/2012.07095); DOI **10.3847/1538-4357/abbfb0** |
| Rostami-Shirazi et al. — *Origin of high dark remnant fractions* | 2025, MNRAS 536, 1332; online 2024. [Publisher](https://doi.org/10.1093/mnras/stae2644); [2508.12167](https://arxiv.org/abs/2508.12167) |
| Evans, Strigari & Zivick — *Dark and luminous mass components of ω Cen* | 2022, MNRAS 511, 4251. [2109.10998](https://arxiv.org/abs/2109.10998); DOI **10.1093/mnras/stac261** |
| Bañares-Hernández et al. — *Combined stellar kinematics and pulsar timing* | 2025, A&A 693, A104. [2408.00939](https://arxiv.org/abs/2408.00939); DOI **10.1051/0004-6361/202451763** |
| Häberle et al. — *Fast-moving stars…ω Cen* | 2024, Nature 631, 285. [2405.06015](https://arxiv.org/abs/2405.06015); DOI **10.1038/s41586-024-07511-z** |
| González Prieto et al. — *Growing the IMBH in ω Cen* | 2025, ApJL 990, L69. [2507.06316](https://arxiv.org/abs/2507.06316); DOI **10.3847/2041-8213/adfd4a** |
| Whitaker et al. — *Long-period stellar-mass BH binary in ω Cen* | 2026, ApJL 1006, L1. [2606.18350](https://arxiv.org/abs/2606.18350) |
| Dickson et al. — *Fast Dynamical Modelling…I* | **Preprint**, submitted 8 July 2026. [2607.07010](https://arxiv.org/abs/2607.07010) |
| Dickson et al. — *Fast Dynamical Modelling…II* | **Preprint**, submitted 10 September 2026. [2609.10955](https://arxiv.org/abs/2609.10955) |
| Colom i Bernadich et al. — *Joint MeerKAT and Parkes view* | **Preprint**, submitted 23 March 2026; under review in checked record. [2603.21845](https://arxiv.org/abs/2603.21845) |