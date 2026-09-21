# Rung 0: analysis of the three isotropic fits

*2026-09-21. Runs `rung0_K1`, `rung0_K2_cored`, `rung0_K2_nfw` in `results/fits/`; comparison
table `results/fits/comparison.md`; figures `plots/fit_posterior_profiles.png`,
`plots/fit_rung0_*_posterior_profiles.png`. Plan: `docs/MODELLING_PLAN.md`; data:
`docs/data_analysis.pdf`; method: `docs/mass_modelling.pdf`.*

## 1. What was run

Rung 0 of the ladder: the simplest defensible model. Isotropic (β ≡ 0), no instrument scale
factors, composite tracer, our own Jeans solver, 89 binned points:

| dataset | points | radial range |
|---|---|---|
| `hst_pm_radial_ours` | 21 | 3.7–342″ (0.1–8.8 pc) |
| `hst_pm_tangential_ours` | 21 | 3.7–342″ |
| `muse_los_dispersion` | 29 | 5.7–290″ |
| `gaia_edr3_ours_radial` | 9 | 357–2149″ (9.2–55 pc) |
| `gaia_edr3_ours_tangential` | 9 | 357–2149″ |

Switches recorded in `run.yaml`: `gaia_errors: raw`, `gaia_rotation: published`,
`isotropic: true`, `no_scales: true`. Three families: K1 (stars + remnants + point mass, 5
parameters), K2-cored and K2-NFW (K1 plus a two-parameter halo, `M_dm_100` and `r_s`).

## 2. Results

| | K1 (no DM) | K2 cored | K2 NFW |
|---|---|---|---|
| ln Z | −190.02 ± 0.48 | **−188.94 ± 0.36** | −190.09 ± 0.35 |
| max ln L | −169.7 | −164.4 | −165.5 |
| χ² at max L (89 points) | 772 | 761 | 763 |
| M★ [M☉] | 2.88 [2.67, 3.03] ×10⁶ | 2.30 [1.68, 2.77] ×10⁶ | 2.34 [1.57, 2.85] ×10⁶ |
| M_rem [M☉] | 3.1 [2.0, 4.8] ×10⁵ | 7.7 [3.9, 13.3] ×10⁵ | 7.1 [3.2, 13.4] ×10⁵ |
| a_rem [pc] | 4.0 [3.2, 4.7] | 5.2 [4.3, 5.7] | 5.1 [4.1, 5.7] |
| M_BH [M☉] | 4.34 [3.90, 4.71] ×10⁴ | 4.65 [4.29, 4.97] ×10⁴ | 4.53 [4.11, 4.88] ×10⁴ |
| M_DM(<100 pc) [M☉] | — | 1.32 [0.19, 2.47] ×10⁶ | 5.3 [0.23, 11.9] ×10⁵ |
| r_s [pc] | — | 261 [46, 663] | 123 [7, 544] |
| D [kpc] | 5.311 [5.291, 5.332] | 5.312 [5.291, 5.333] | 5.313 [5.292, 5.334] |
| 95 % upper limit M_DM(<100 pc) | — | 3.13 ×10⁶ | 1.64 ×10⁶ |
| likelihood calls / wall time | 3.2 ×10⁶ / 4.0 h | 6.6 ×10⁶ / 10.0 h | 14.3 ×10⁶ / 20.6 h |

**Δ ln Z**: cored − K1 = +1.08 ± 0.60; NFW − K1 = −0.07 ± 0.59; cored − NFW = +1.15 ± 0.50.
A Bayes factor of 2.9 is "barely worth mentioning" on the Jeffreys scale. **There is no
evidence for dark matter at rung 0, and the NFW form is indistinguishable from no halo.**

**Comparability verified.** The report flagged the three runs as fitted to different
observations, because the Gaia product's file hash changed on 20 Sep 17:22 between the K1 and
K2 launches. Recomputing K1's χ² at its stored maximum-likelihood point with today's files
reproduces all five per-dataset values to five significant figures, so the change was metadata
only and the evidences are comparable. (`run.yaml` now hashes inputs at launch rather than at
the end of a run, so this cannot recur.)

## 3. Goodness of fit: the χ² is entirely beyond 50″

χ² at the K2-cored maximum-likelihood point, split by projected radius:

| dataset | R < 50″ | 50–150″ | > 150″ |
|---|---|---|---|
| HST radial | 12.9 (11 pts, 1.2/pt) | 49.5 (5) | 90.1 (5) |
| HST tangential | 13.6 (11 pts, 1.2/pt) | 165.8 (5) | 304.4 (5, 61/pt) |
| MUSE LOS | 9.5 (13 pts, 0.7/pt) | 27.5 (10) | 38.9 (6) |

**Inside 50″ the model is an excellent fit to all three datasets.** Every part of the bad χ²
comes from larger radii, and two thirds of the total from the HST tangential component alone.

The misfit is small in physical terms. Fractional residuals (data − model)/data:

| | 108″ | 138″ | 176″ | 224″ | 271″ | 311″ |
|---|---|---|---|---|---|---|
| HST radial | +0.3 % | +1.4 % | +1.2 % | +1.1 % | +0.2 % | −1.2 % |
| HST tangential | −2.0 % | −3.7 % | −5.1 % | −7.5 % | −9.6 % | −10.0 % |
| statistical error | 0.26–0.30 % | 0.24–0.35 % | 0.20–0.49 % | 0.16–0.79 % | 0.21–1.19 % | 0.51–1.61 % |

χ² = 761 corresponds to a model wrong by 1 % in one component and 5–10 % in another. It looks
catastrophic only because the HST errors are 0.16–0.35 % in those bins. In this regime χ² is a
statement about unmodelled systematics and missing degrees of freedom, not about the mass.

## 4. The misfit is in the ratio of the two components, not in the mass

The geometric mean of the two PM components is matched to 2–3 % at 150–300″, so the enclosed
mass is essentially right. What is wrong is σ_T/σ_R. The isotropic model does not predict
exactly 1 for this ratio, because the tangential model has the published rotation subtracted;
the residual excess (data ratio over model ratio) is the part that needs new physics:

| R | 2.2 pc | 3.6 pc | 5.8 pc | 8.0 pc | 9.2 pc | 13.7 pc | 21.4 pc | 26.9 pc | 34.2 pc | 43.6 pc | 55.3 pc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| data σ_T/σ_R | 0.979 | 0.937 | 0.881 | 0.847 | 0.931 | 0.859 | 0.940 | 0.988 | 1.081 | 1.103 | 1.176 |
| model σ_T/σ_R | 0.995 | 0.985 | 0.958 | 0.920 | 0.901 | 0.850 | 0.907 | 0.939 | 0.967 | 0.989 | 0.998 |
| excess | −1.6 % | −4.8 % | −8.0 % | −8.0 % | +3.3 % | +1.1 % | +3.6 % | +5.3 % | **+11.8 %** | **+11.6 %** | **+17.9 %** |

(The first four columns are HST, the rest Gaia EDR3; 9.2 pc onwards is a different instrument,
so the step at 8 → 9 pc should not be over-interpreted.)

**This is the classical anisotropy profile of a tidally limited star cluster**: isotropic in the
core, radially anisotropic in the intermediate region (β ≈ +0.2 to +0.3 at 4–8 pc), turning
tangential in the outskirts (β ≈ −0.2 to −0.4 beyond ~30 pc), because radial orbits are the
ones that reach the tidal boundary and are preferentially stripped. The sign reversal is
visible within the Gaia data alone, where σ_T/σ_R rises monotonically from 0.86 at 11 pc to
1.18 at 55 pc and crosses unity near 27 pc.

The MUSE line-of-sight dispersion is consistent with the same picture: at 290″ (7.5 pc) it is
9.4 % below the model (−4.9 σ), and at large projected radius the line of sight samples the
tangential direction, so radial anisotropy at 4–8 pc must suppress it. That an independent
instrument with independent systematics shows the same deficit at the same radius is the
strongest argument that this is dynamics and not a proper-motion artefact.

### The alternative reading

The tangential component is also where rotation enters, so an error in the rotation curve
produces the same signature. Quantitatively: reconciling the HST tangential data at 224″
requires an extra 5.8 km/s in quadrature on top of the 4.45 km/s the published curve already
supplies there, i.e. a ~64 % larger rotation amplitude at 5.8 pc; and reconciling the outer
Gaia points requires *less* rotation than the model subtracts at 34–55 pc. A single
normalisation error cannot do both — it would need the shape of the rotation curve to be wrong
in opposite directions at 6 and 40 pc. The anisotropy reading needs no such conspiracy, but the
rotation reading is cheap to test (`--gaia-rotation ours`, and the HST equivalent), and should
be, because the two are not mutually exclusive.

## 5. Error inflation cannot rescue the fit

χ² at the same best model with a systematic floor added in quadrature (not refitted):

| floor | 0 | 0.5 % | 1 % | 2 % |
|---|---|---|---|---|
| HST radial | 152 | 38 | 19 | 11 |
| HST tangential | 484 | 272 | 149 | 63 |
| MUSE | 76 | 72 | 62 | 42 |
| Gaia radial | 25 | 24 | 19 | 12 |
| Gaia tangential | 24 | 23 | 21 | 16 |
| **total (89)** | **761** | **428** | **271** | **144** |

A 0.5 % floor fixes the HST radial component almost completely (152 → 38): that misfit really
is error-bar-driven. The tangential component is not: it still contributes 63 of the 144 at a
2 % floor. The shape is wrong, not just the uncertainties. The planned error-model variant
(`--gaia-errors eta`) therefore cannot substitute for the missing anisotropy.

## 6. What the "evidence for dark matter" is actually made of

Δχ² between K1 and K2-cored, by dataset:

| | HST radial | HST tangential | MUSE | Gaia radial | Gaia tangential | total |
|---|---|---|---|---|---|---|
| Δχ² (K2-cored − K1) | −8 | +7 | −1 | **+6** | **−14** | −10 |

The halo buys Δχ² ≈ 10 for two parameters, and the gain comes almost entirely from the outer
Gaia points: the tangential ones improve by 14 and the radial ones worsen by 6. Those are the
same 34–55 pc points that demand β < 0. An isotropic model can only raise both components
together, and raising them helps the tangential data (which sit above the model) more than it
hurts the radial data (which sit below). **In other words, the halo is being recruited as a
proxy for tangential anisotropy in the outskirts.** Until β is free, the two are degenerate and
the ln Z difference of +1.1 cannot be read as a statement about dark matter.

## 7. Other parameters

**Distance.** All three runs give D = 5.311–5.313 ± 0.021 kpc against a prior of
5.43 ± 0.05 — the data pull it down by 2.4 prior σ and tighten it by a factor 2.4. This is a
kinematic distance from matching proper motions in mas/yr to line-of-sight velocities in km/s,
and it is driven mainly by the inner data where the fit is good. It is nevertheless a result
that should be re-checked at rung 1, because anisotropy changes the PM-to-LOS ratio and so can
bias D.

**Central mass.** The total M★ + M_rem + M_BH is 3.23 ×10⁶ (K1) and 3.12 ×10⁶ (K2-cored), i.e.
stable, but the split moves a lot: adding a halo moves 6 ×10⁵ M☉ from the stars to the
remnants. M_BH is 4.3–4.7 ×10⁴ M☉ with an 8 % error in all three runs; the inner fit is good
(χ²/point = 1.2), so this is the most robust number in the table — but the classical degeneracy
between a central point mass and central tangential anisotropy is not tested at rung 0, where
β ≡ 0 everywhere.

**Halo parameters.** r_s is unconstrained in both K2 runs (46–663 pc cored, 7–544 pc NFW),
which is expected: the data end at 55 pc and the halo is a small perturbation there.

## 8. Independent-engine check

Evaluating each best sample with JamPy instead of our Jeans solver:

| run | our solver | JamPy |
|---|---|---|
| K1 | ln L = −169.66 | −170.35 |
| K2-cored | −164.44 | −165.08 |

Agreement to ~0.7 in ln L. The per-dataset shifts are ±18 in χ² between the HST radial and
tangential components and < 1 elsewhere, which is the known difference in how the two engines
project the velocity ellipsoid; it does not affect any conclusion here.

## 9. Sampling cost

K1 took 4 h and 3.2 ×10⁶ likelihood calls; K2-cored 10 h / 6.6 ×10⁶; K2-NFW 20.6 h /
14.3 ×10⁶. The NFW run spent most of its time crawling along the M_DM–r_s degeneracy ridge
(one accepted draw per ~4000 calls for several hours, with 82 % of the evidence still in the
live points after 12 h). For rungs 1 and 2, which add one and three parameters, the NFW family
should be launched with ultranest's step sampler (`step_sampler=True` in `run_nested`) or the
r_s prior narrowed on physical grounds; otherwise the cost will become prohibitive.

## 10. Conclusions

1. Rung 0 is a **good fit inside 50″ and a bad fit outside**, for all three families.
2. The bad fit is a **ratio error between the two PM components**, not a mass error: the
   enclosed mass is matched to 2–3 % where the data are most precise.
3. The pattern — isotropic core, radial at 4–20 pc, tangential beyond ~30 pc — is the
   **textbook anisotropy profile of a tidally limited cluster**, and it is coherent across
   HST, MUSE and Gaia.
4. **A constant β cannot fit it.** Rung 1 will improve matters but cannot reproduce a sign
   reversal; the β(r) family of rung 2, with the turnover parameter that was deferred, is what
   the data are asking for. The trigger for that deferred item has fired.
5. The apparent preference for a cored halo (Δ ln Z = +1.1) is **degenerate with tangential
   anisotropy in the outskirts** and should not be quoted.
6. Error inflation cannot substitute for the missing freedom.

## 11. What to do next, in order

1. **Rung 1** (constant β, 6 and 8 parameters) on the same 89 points — as planned, and now also
   as a diagnostic: the fitted β will be a compromise between +0.25 at 6 pc and −0.3 at 40 pc,
   and the residuals should retain the sign reversal.
2. **Rung 2** (β(r) with a turnover) — promoted from "deferred" to "required" by §4 above.
3. **Rotation variants** (`--gaia-rotation ours` and the HST equivalent) to separate anisotropy
   from an error in the rotation curve. Cheap, and it is the only competing explanation.
4. Only then: the error-model variant, the prior-sensitivity run, the injection test, and the
   one diagnostic run with instrument scales free.
5. Do not quote M_DM limits from rung 0 in any document; the numbers in §2 are conditional on
   isotropy and are superseded by rung 1/2.
