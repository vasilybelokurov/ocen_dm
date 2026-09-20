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

The data are **89 points**. Each rung runs the same two families on the same observations.

| Rung | Anisotropy | Scales | K1 params | K2 params | Points per parameter (K2) |
|---|---|---|---|---|---|
| **1** | constant `beta` | none | 6 | 8 | 11.1 |
| **2** | `beta(r)`, 3 params | none | 8 | 10 | 8.9 |
| **3** | `beta(r)` | per instrument | 11 | 13 | 6.8 |

Rung 1 is the most robust pair we can write down and the one where a Bayes factor means
something simple. Rung 3 exists only because instrument scales are known to be degenerate
with halo mass; if a halo survives only when a scale is free, that is a result in itself.

```
ocen-dm fit --family K1       --constant-beta --no-scales
ocen-dm fit --family K2-cored --constant-beta --no-scales
```

## What decides whether to climb

After each rung, three questions, in order:

1. **Do the residuals show structure the next rung would fix?** If rung 1's residuals are
   flat, rung 2 cannot be what the data are asking for.
2. **Has the evidence gap moved?** Compare `Delta ln Z` between rungs. Run-to-run scatter is
   0.5, so movement below ~1.5 is noise. **A gap that is stable across rungs is the finding**:
   it says the answer is not about orbital freedom. A gap that moves a lot is also the
   finding, and a more interesting one.
3. **Is the top rung's model even acceptable?** `chi2` per dataset at the maximum-likelihood
   sample. An evidence comparison between two models that both fit badly is not meaningful.

Only climb if question 1 says yes.

## Taken now, because they remove freedom rather than adding it

* **`beta_0` prior narrowed to `U[-1, -0.5]`.** An & Evans give `gamma >= beta + 1/2` for a
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
* Evidences are compared only between runs whose recorded provenance matches: dataset keys,
  point count, input hashes, and every switch that changes what the likelihood is shown.
* Every rung's result goes in `JOURNAL.md` with its `Delta ln Z`, its per-dataset `chi2`, and
  a plot of the residuals, before deciding whether to climb.
