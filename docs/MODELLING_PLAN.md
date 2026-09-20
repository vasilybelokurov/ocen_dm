# Modelling plan: start simple, add only what the residuals demand

*Adopted 2026-09-20.*

## The problem this plan exists to avoid

The data preparation is finished and validated. The modelling machinery has been reviewed
three times. **No fit has been run on the current dataset.** In that state the temptation is
to keep answering reviewers by adding parameters — a ninth anisotropy control, four rotation
nuisances — none of which has been shown to be needed by a single residual. That is premature
optimisation, and it is expensive in exactly the wrong currency: every added parameter makes
the Bayes factor harder to interpret, because the question *is* how much freedom each family
has.

So the plan is a ladder. Climb one rung at a time and stop when the answer stops moving.

## The ladder

The data are **89 points**. Each rung runs the same families on the same observations.

| Rung | Anisotropy | `beta` prior | Scales | K1 | K2 | Points/param (K2) |
|---|---|---|---|---|---|---|
| **0** | isotropic, `beta = 0` fixed | -- | none | 5 | 7 | 12.7 |
| **1** | constant `beta` | `U[-1, 0.5]` | none | 6 | 8 | 11.1 |
| **2** | `beta(r)`, 3 params | `beta_0 ~ U[-1, -0.5]` | none | 8 | 10 | 8.9 |

**Rung 0 is the first experiment.** Isotropic Jeans with stars, remnants and a point mass,
against the same with a cored halo. Five and seven parameters. It is the model everyone
understands, its Bayes factor has the simplest possible meaning, and its residuals tell us
directly whether anisotropy is even needed. Then the NFW halo as a third run, to bracket the
halo shape.

Rung 1 lets the anisotropy float but not vary with radius. Rung 2 lets it vary.

**There is no rung with instrument scales.** A multiplicative scale on an instrument is
degenerate with mass by construction -- one on Gaia EDR3 rescales the outer dispersion, which
is where a halo lives -- and the whole data-preparation effort was spent making the
instruments agree without one (HST/Gaia 1.01 +- 0.12 in the overlap, Pristine/Gaia 1.005 +-
0.012). Adding scales afterwards would say "and if they still disagree, absorb it". If they
still disagree, that is a finding about the data. The MUSE offset in particular is likely
physical (line of sight against proper motion, a different tracer population, equipartition)
and a scale would paper over it rather than model it. Scales are a fudge; removed
2026-09-20 at the user's call.

```
# the first round
ocen-dm fit --family K1       --isotropic --no-scales
ocen-dm fit --family K2-cored --isotropic --no-scales
ocen-dm fit --family K2-nfw   --isotropic --no-scales
# rung 1
ocen-dm fit --family K1       --constant-beta --no-scales
ocen-dm fit --family K2-cored --constant-beta --no-scales
```

Every rung's families build and evaluate the likelihood in 4 ms; a nested run is
$10^5$--$10^6$ evaluations, so expect tens of minutes for rung 0 and a few hours for rung 3.

### Why the anisotropy prior differs between rungs

The An & Evans bound `beta_0 <= -1/2` is a condition **at the centre**. In `beta(r)`,
`beta_0` is genuinely the central value, so the bound applies (rung 2 onwards). A constant
`beta` has no separate centre: imposing `-0.5` on it would force the tangential dispersion to
exceed the radial one by 22 per cent at every radius, which the data never show (projected
`sigma_T/sigma_R` runs 0.84 to 1.18, i.e. `beta` roughly -0.4 to +0.3). So rungs 0 and 1 are
deliberately crude baselines that are **not** required to be DF-realisable at the centre,
and rung 1 gets `U[-1, 0.5]`. This was briefly got wrong on 2026-09-20 and caught in review.

## What decides whether to climb

After each rung, three questions, in order:

1. **Do the residuals show structure the next rung would fix?** If rung 0's residuals are
   flat, anisotropy is not what the data are asking for and rung 1 is unnecessary. If they
   show a radial trend in `sigma_T/sigma_R`, that is the case for climbing.
2. **Has the evidence gap moved?** Compare `Delta ln Z` between rungs. Run-to-run scatter is
   0.5, so movement below ~1.5 is noise. **A gap that is stable across rungs is the finding**:
   it says the answer is not about orbital freedom. A gap that moves a lot is also the
   finding, and a more interesting one.
3. **Is the top rung's model even acceptable?** `chi2` per dataset at the maximum-likelihood
   sample. An evidence comparison between two models that both fit badly is not meaningful.

## One diagnostic, after the first round

Run the best-fitting rung **once** with instrument scales free (`ocen-dm fit ...` without
`--no-scales`) and look at where the scales land. If every `s_inst` sits at 1.00 within a
couple of per cent, that is an independent validation of the data preparation. If one pulls
away, that is a red flag about that instrument and goes back to the data, not into the
model. **Its evidence is never compared with anything.**

Only climb if question 1 says yes.

## Taken now, because they remove freedom rather than adding it

* **`beta_0` prior narrowed to `U[-1, -0.5]` for the radially varying anisotropy (rung 2+).**
  An & Evans give `gamma >= beta + 1/2` for a
  tracer in a point-mass-dominated potential; our tracer is cored and every model carries a
  central point mass, so `-0.5 < beta_0 <= 0` admits models with no non-negative distribution
  function. The counter-argument that the sphere of influence is unresolved **fails**: it
  reaches the innermost datum at 0.11 pc for `M_bh > 7400 Msun`, which is 46 per cent of the
  prior. This makes the model more constrained, not more complicated.
* **The AGAMA positivity claim dropped.** AGAMA clips negative DF values to zero, so a
  successful construction does not certify a positive DF, and its Cuddeford family cannot
  represent the anisotropy we fit anyway. It stays as a limited cross-check with the
  limitation stated.

## Deferred, with the reason written down

Neither is wrong. Both are held back until a baseline shows they are needed.

* **An anisotropy turnover.** One extra control parameter would let `beta(r)` be radial in
  the middle and tangential outside, which is the pattern our projected ratios show (0.93 at
  223 arcsec, 1.11 at 1328). It was motivated by looking at these data, so it would need
  injection tests to show it does not simply absorb a real signal. **Revisit if rung 1 or 2
  leaves an anisotropy-shaped residual.**
* **A shared rotation nuisance.** The rotation-curve uncertainty is coherent across bins and
  currently enters per bin, which understates it. The correct treatment marginalises over the
  published spline's coefficients — but that needs their covariance, and only percentiles are
  published, so it may not be recoverable. **Revisit if the outer tangential bins drive the
  result.**

## Standing rules

* Fits are run only when explicitly asked for, one rung at a time.
* No instrument scale factors in any compared model.
* Evidences are compared only between runs whose recorded provenance matches: dataset keys,
  point count, input hashes, and every switch that changes what the likelihood is shown.
* Every rung's result goes in `JOURNAL.md` with its `Delta ln Z`, its per-dataset `chi2`, and
  a plot of the residuals, before deciding whether to climb.
