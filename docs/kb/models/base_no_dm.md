# Base model: two-transition DF, no dark matter

The reference model for the cluster mass distribution. Every dark-matter
variant is compared with it.

!!! info "Status (2026-10-13)"
    Fitted to the real data on 2026-09-23
    (`results/df/dftwo_counts_20260923`, best job `observed_no_halo_start0`,
    source commit `36bca6f`). The fit is numerically validated but **fails
    its own data-fit gate**. The count likelihood was corrected on 2026-10-13
    (commit `1af8c33`); the model has not yet been refitted with the
    corrected count products.

## Mass components

| Component | Form | Free parameters (range) | Best fit |
|---|---|---|---|
| Stars | Self-consistent spherical action DF (below); both the counted tracer and the mass that follows the light | \(M_\star\) (1–6×10⁶ M☉, log) | 3.00×10⁶ M☉ |
| Remnants | Plummer potential; gravity only, no DF, not a tracer | \(M_{\rm rem}\) (10²–10⁶ M☉, log), \(a_{\rm rem}\) (0.3–20 pc, log) | 2.62×10⁵ M☉, 2.15 pc |
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
| \(J_0\) | envelope scale (pc km/s) | 30–800, log | 54.1 |
| \(\alpha\) | envelope sharpness | 0.6–3 | 0.80 |
| \(b_{\rm out}\) | intermediate anisotropy | −2 to 2 | +0.74 (radial) |
| \(J_a\) | first transition (pc km/s) | 5–500, log | 64.1 |
| \(b_{\rm outer}\) | outer anisotropy | −4 to 2 | −0.33 (tangential) |
| \(J_{\rm outer}\) | second transition (pc km/s) | 2–100 × \(J_a\) | 319 |

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
10.67 M☉ of light-following mass per F625W < 19 star. Details of every
product: [Data products](../data/index.md).

## Fit quality

| Term | Value | Gate | Pass |
|---|---|---|---|
| \(\chi^2_{\rm kin}\) per bin | 158.2 / 89 = 1.78 | ≤ 1.3 | no |
| HST PM radial | 2.07 per bin | ≤ 2.0 | no |
| HST PM tangential | 1.95 per bin | ≤ 2.0 | yes |
| MUSE LOS | 1.31 per bin | ≤ 2.0 | yes |
| Gaia PM radial | 1.55 per bin | ≤ 2.0 | yes |
| Gaia PM tangential | 2.42 per bin | ≤ 2.0 | no |
| HST count deviance | 98.8 / 20 = 4.9 per bin | ≤ 2.0 | no |
| Gaia count deviance | 26.3 / 12 = 2.2 per bin | ≤ 2.0 | no |
| Total objective | 285.6 (incl. distance prior 2.2) | — | — |

The HST count residuals are structured: the model is ~4% low at 2–4 pc,
~4% high at 5–6.5 pc, and low inside 0.1 pc. Replaying the same parameters
with the corrected count products gives HST deviance 98.1 and Gaia 27.2.

## What has been tested

- **Convergence**: 4 of 6 starts valid (2 Latin starts infeasible), all four
  reach the same minimum; no parameter at a bound.
- **Numerics**: grid refinement shifts predictions by 0.0005σ; direct
  projection checks agree to ≤ 5×10⁻⁴ in Σ, LOS, PM_R and PM_T.
- **Estimator recovery**: 12 mocks drawn from this model with published
  kinematic noise and Poisson counts recover the total mass to ±3% at 20 pc
  and ±5% over 2–40 pc, \(M_\star\) to ±6%, \(M_{\rm rem}\) to ±23%.

Not yet tested: a refit with the corrected count products; mocks whose truth
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
