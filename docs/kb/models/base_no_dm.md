# Base model: two-transition DF, no dark matter

The reference model for the cluster mass distribution. Every dark-matter
variant is compared with it.

!!! info "Status (2026-10-13)"
    Refitted on 2026-10-13 with the corrected (area-weighted) count
    likelihood: `results/df/dftwo_counts_area_20261013`, best job
    `observed_no_halo_start3`. The fit is numerically validated but **fails
    its own data-fit gate**. It supersedes the 2026-09-23 fit
    (`results/df/dftwo_counts_20260923`), whose star-weighted count averages
    were biased; parameters moved by ≤ 4% and the mass per counted star by
    −1.8%.

## Mass components

| Component | Form | Free parameters (range) | Best fit |
|---|---|---|---|
| Stars | Self-consistent spherical action DF (below); both the counted tracer and the mass that follows the light | \(M_\star\) (1–6×10⁶ M☉, log) | 3.00×10⁶ M☉ |
| Remnants | Plummer potential; gravity only, no DF, not a tracer | \(M_{\rm rem}\) (10²–10⁶ M☉, log), \(a_{\rm rem}\) (0.3–20 pc, log) | 2.73×10⁵ M☉, 2.18 pc |
| Dark matter | none (\(\rho_{20}=0\)) | — | — |
| Central black hole | none | — | — |

### Stellar DF

The regularized exponential action DF (contour construction of Binney 2026,
MNRAS 549, stag854), rebuilt at every iteration and iterated to
self-consistency in the total potential (stars + remnants):

\[
f(\mathbf J)=A\,\exp\!\left[-\left(L_c/J_0\right)^{\alpha}\right],
\qquad q=J_r+L .
\]

The contour label \(L_c\) follows the spherical frequency ratio, biased by the
anisotropy function \(b(q)\): \(b>0\) favours radial orbits, \(b<0\)
tangential ones. Two transitions:

\[
b(q)=b_{\rm out}\,s(q;J_a)+\left(b_{\rm outer}-b_{\rm out}\right)s(q;J_{\rm outer}),
\qquad s(q;J)=\frac{q^2}{q^2+J^2}.
\]

The core is isotropic (\(b\to0\) as \(q\to0\)); the model is radially
anisotropic at intermediate radii and tangential outside.

| Parameter | Meaning | Range | Best fit |
|---|---|---|---|
| \(J_0\) | envelope scale (pc km/s) | 30–800, log | 56.0 |
| \(\alpha\) | envelope sharpness | 0.6–3 | 0.80 |
| \(b_{\rm out}\) | intermediate anisotropy | −2 to 2 | +0.77 (radial) |
| \(J_a\) | first transition (pc km/s) | 5–500, log | 65.5 |
| \(b_{\rm outer}\) | outer anisotropy | −4 to 2 | −0.33 (tangential) |
| \(J_{\rm outer}\) | second transition (pc km/s) | 2–100 × \(J_a\) | 312 |

### Remnant component

\[
\rho_{\rm rem}(r)=\frac{3M_{\rm rem}}{4\pi a_{\rm rem}^3}
\left(1+\frac{r^2}{a_{\rm rem}^2}\right)^{-5/2}
\]

It stands for all mass that does not follow the counted light: neutron stars
and stellar-mass black holes, any mass-segregated white dwarfs and heavy
stars, and part of the effect of a central black hole. At the best fit it
holds 8% of the total mass, has a half-mass radius of \(1.305\,a_{\rm rem}=2.8\)
pc, and dominates the mass inside ~1.5 pc (60% inside 1 pc).

### Distance, geometry, rotation

- Distance free in 4.8–6.0 kpc with a Gaussian prior \(N(5.43, 0.05)\) kpc;
  best fit 5.356 kpc (1.5σ low).
- Spherical and non-rotating. The published rotation is removed at the
  moment level, \(\sigma^2_{\rm pred}=\langle v^2\rangle_{\rm DF}-v^2_{\rm rot}\);
  a negative result is rejected, never floored.

## Data and likelihood

Objective \(=-2\ln L\) (up to a constant):

| Term | Data | Likelihood |
|---|---|---|
| Kinematics, 89 bins | HST PM radial and tangential (21 + 21), MUSE LOS (29), Gaia EDR3 PM radial and tangential (9 + 9) | split-normal with the published errors |
| HST counts | oMEGACat, F625W < 19, 20 bins over 2–250″, 117 751 stars | Poisson, amplitude profiled |
| Gaia counts | Vasiliev & Baumgardt 2021, G < 17, 12 bins over 300–2400″ | Poisson, amplitude and uniform field profiled |
| Distance prior | — | Gaussian |

The count amplitude converts to a mass per counted star: the best fit has
10.48 M☉ of light-following mass per F625W < 19 star (10.67 in the
2026-09-23 fit). Details of every
product: [Data products](../data/index.md).

## Fit quality

| Term | 2026-10-13 | 2026-09-23 | Gate | Pass |
|---|---|---|---|---|
| \(\chi^2_{\rm kin}\) per bin | 156.4 / 89 = 1.76 | 158.1 / 89 = 1.78 | ≤ 1.3 | no |
| HST PM radial | 2.03 per bin | 2.06 | ≤ 2.0 | no |
| HST PM tangential | 1.89 per bin | 1.95 | ≤ 2.0 | yes |
| MUSE LOS | 1.31 per bin | 1.31 | ≤ 2.0 | yes |
| Gaia PM radial | 1.60 per bin | 1.55 | ≤ 2.0 | yes |
| Gaia PM tangential | 2.39 per bin | 2.43 | ≤ 2.0 | no |
| HST count deviance | 96.9 / 20 = 4.8 per bin | 98.9 | ≤ 2.0 | no |
| Gaia count deviance | 25.2 / 12 = 2.1 per bin | 26.3 | ≤ 2.0 | no |
| Total objective | 280.7 (incl. distance prior 2.2) | 285.6 | — | — |

The two fits use different count likelihoods, so their objectives are not
directly comparable; the table shows that the correction changes little.
The HST count residuals remain structured: the data lie below the model
near 1 pc and in the outermost HST bin (~6 pc), with deviance residuals of
about −4, and above it at 3–4 pc (+4.6).

![Data versus the base model, both fits](../../../plots/base_area_20261013_data.png)

*Counts (HST and Gaia), MUSE line-of-sight and HST/Gaia proper-motion
dispersions against the 2026-10-13 (blue) and 2026-09-23 (orange) fits,
with the earlier Jeans model (dashed). Bottom row: residuals; the count
panel shows signed Poisson deviance residuals. Both fits are evaluated with
the corrected count likelihood.*

## Properties of the fitted model

![Mass components of the base model](../../../plots/base_model_20261013_components.png)

*By component: (a) density; (b) enclosed mass with half-mass radii dotted;
(c) circular velocity; (d) remnant share, local and cumulative; (e)
dynamical-to-stellar mass ratio. Grey: outside the radial range of the
kinematic data (0.1–63 pc).*

![Internal kinematics of the base model](../../../plots/base_model_20261013_kinematics.png)

*(a) Intrinsic radial and tangential dispersions; (b) anisotropy β(r);
(c) stellar surface density with half-mass radii; (d) projected
line-of-sight and proper-motion dispersions from the DF (before the
published rotation is subtracted).*

![Structure of the base-model DF](../../../plots/df_structure_dftwo_counts_area_20261013_observed_no_halo_start3.png)

*(a) Intrinsic moments; (b) anisotropy bias b(q) over the action
distributions of the HST/MUSE and Gaia stars; (c) f(J_r, L).*

Key numbers (current fit): stellar half-mass radius 10.4 pc (projected
7.7 pc), remnant half-mass radius 2.8 pc; remnants 66% of the central
density; \(M_{\rm tot}/M_\star\) from 2.9 in the centre to 1.09 overall;
peak circular velocity 27.0 km/s at 9 pc; central dispersion 21.4 km/s;
β peaks at +0.15 near 8 pc and reaches −0.53 at 100 pc. Interpretation:
[Results](../results.md#base-model-fits).

## What has been tested

- **Convergence**: 4 of 6 starts valid (2 Latin starts infeasible), all four
  reach the same minimum (objectives within 0.01); no parameter at a bound.
- **Numerics**: grid refinement shifts predictions by 0.0005σ; direct
  projection checks agree to ≤ 5×10⁻⁴ in Σ, LOS, PM_R and PM_T.
- **Estimator recovery** (2026-09-23 fit, old count likelihood): 12 mocks drawn from this model with published
  kinematic noise and Poisson counts recover the total mass to ±3% at 20 pc
  and ±5% over 2–40 pc, \(M_\star\) to ±6%, \(M_{\rm rem}\) to ±23%.

Not yet tested: mocks whose truth
breaks the model's own assumptions (mass segregation, rotation, another
remnant shape); calibration of the error bars.

## Built-in assumptions

1. Mass follows light: one stellar DF is both tracer and mass.
2. Isotropic core.
3. One remnant component with a fixed Plummer shape.
4. No central black hole.
5. Spherical; rotation subtracted from published curves.
6. Published error bars taken as they are.

## Reproduce

```bash
source <venv>/bin/activate
export PYTHONPATH=src
python bin/run_compact_df_recovery.py prepare-real --out results/df/<new_label> \
  --two-transition --branches no_halo \
  --counts ocen_counts_hst_f625w19,ocen_counts_gaia_g17 \
  --free-distance 4.8:6.0 --distance-prior 5.43:0.05 \
  --bound stellar.J0=30:800 --bound stellar.J_a=5:500 --bound stellar.b_out=-2:2 \
  --hand-starts 2 --n-starts 4 --probe-workers 5
python bin/run_compact_df_recovery.py run --out results/df/<new_label> --workers 2 --max-calls 600
```

Each run uses 2 jobs × (5 probe workers + 1) = 12 cores and took ~30 min in
September. The batch freezes its source, configuration and inputs; a run
directory is never reused.

Code: `src/ocen_dm/kinematics/regularized_df.py` (DF, remnants, halo),
`df_fit.py` (joint problem), `counts.py` (count likelihood),
`compact_recovery.py` (mocks, residuals), `bin/run_compact_df_recovery.py`
(driver).

## Related results

- Dark-matter variants and the \(\rho_{20}\) and mass-per-star scans:
  [DF observed fit write-up](../../df_observed_fit.pdf).
- Review of next steps: [Codex review, 2026-10-13](../../reviews/codex_mass_model_next_2026-10-13.md).
