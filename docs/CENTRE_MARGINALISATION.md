# Marginalising over the cluster centre (plan, tests, results)

Plan agreed 2026-11-03 after the Codex and Fable reviews
(`docs/reviews/{codex,fable}_centre_plan_20261103.md`). The centre is a nuisance
parameter; AvdM10 (oMEGACat pixel 15000, 15000) stays the reporting origin.

## Question

Do the fixed-black-hole DF comparisons (no BH / 8,200 / 4×10⁴ M☉) and the central
count deficit survive marginalisation over the centre?

## Definitions

- **Centre prior.** Isotropic Gaussian about AvdM10, σ_c = 2″ per axis (from the
  F625W 19–22 counts at R < 60″, (+0.5 ± 2.0, −1.4 ± 2.1)″, which are not in the fit).
  Sensitivity: σ_c = 1″ and 3″.
- **Quadrature.** Product Gauss–Hermite grid, n × n nodes, weights w_k (sum 1).
  Default n = 6 (36 centres); convergence check with n = 10.
- **Per-centre data d(c).** Same bin edges and selections as the reference products,
  stars re-assigned about c:
  - HST F625W < 19 counts, all 20 bins (2–250″), coverage recomputed about c;
  - Gaia G < 17 counts, all 12 bins (300–2400″);
  - HST PM radial/tangential dispersions (`hst_pm_*_ours`), the 6 bins with
    r_upper ≤ 15″ (3.24–13.8″), our mixture estimator, nodes recomputed;
  - MUSE LOS dispersion, the 3 bins with r_upper ≤ 15″: published value plus the
    change of our estimator (Gaussian with v₀ + v_rot sin(φ − θ₀), per-star errors,
    `selection_hq_pm_and_los`, which reproduces the published n_stars in all 29 bins
    about AvdM10) between c and AvdM10; published errors and rotation term kept;
  - all other kinematic bins unchanged.
- **Marginal likelihood at fixed model θ.**
  ln L_m(θ) = logsumexp_k [ln w_k + ln L(θ | d(c_k))], with ln L the full
  log-likelihood (split-normal kinematics + Poisson counts with Γ term; count
  amplitudes and the Gaia field profiled per centre), never the deviance.
- **Models.** Saved best fits, not refitted: `dftwo_counts_area_20261013` (no BH),
  `dfbh8200_restart_20261015`, `dfbh40000_restart_20261015` (same data and count likelihood).

## Tests and acceptance

| # | Test | Pass |
|---|---|---|
| T1 | Per-centre builders at zero offset reproduce the frozen run inputs (counts, areas, nodes, weights; HST PM values, errors, nodes; MUSE values) | identical (counts exact, floats ≤ 1e-12 rel.) |
| T2 | ln L(θ \| d(c₀)) equals `DFJointProblem.evaluate(...)["loglike_profiled"]` of the frozen problem | ≤ 1e-9 abs |
| T3 | δ-prior (σ_c → 0) marginal equals T2 | ≤ 1e-9 abs |
| T4 | Quadrature convergence: n = 10 vs n = 6 marginal Δ(−2 ln L) between models | change < 0.5 |
| T5 | Star conservation: total HST counts in 2–250″ change with c only by edge flux (within ±3 √(edge-crossing stars)) | yes |
| T6 | Coverage: per-bin coverage change across the grid | < 0.5 % |
| T7 | Kinematic re-centring beyond 15″: change of HST PM bins 13.8–58.5″ at the ±2σ_c centres | < 0.5 σ per bin (else extend the re-centred range) |
| T8 | Unbinned cross-check (counts only, fixed window 2–58.5″ about AvdM10, λ = a Σ_θ(\|x − c\|), a profiled): marginal Δ(−2 ln L) between models vs the binned counts-only marginal over the same bins | same sign; difference reported |
| T9 | MUSE estimator reproduces the published inner σ_los at AvdM10 | within 0.5 published σ (it only enters as a difference) |

## Decision rule (set before running)

- The 4×10⁴ M☉ count/kinematic penalty is **robust** if the marginal
  −2Δ ln L (BH − no BH) > 5 and its centre-to-centre scatter (weighted s.d.) is
  less than half its value. Then no refit.
- **Ambiguous** otherwise: refit no BH and 4×10⁴ M☉ with the marginalised objective
  (about 30 min per fit on 3–6 cores; launched only after asking).

## Outputs

- `src/ocen_dm/centre_marg.py`, `bin/centre_marginal.py`, `tests/test_centre_marg.py`
- `results/plot_data/centre_marginal_20261103.json`, `plots/centre_marginal_20261103.png`
- Central count deficit: marginalised P(N_inner2 ≤ n_obs).

## Results

Run 2026-11-03 (`bin/centre_marginal.py`, 1 core, 49 min; 209 centres; models fixed at their
saved fits, not refitted). Numbers: `results/plot_data/centre_marginal_20261103.json`. Figures:
`plots/centre_marginal_20261103.png`, `plots/centre_marginal_20261103_data.png`.

**Tests**

| # | Result | Pass |
|---|---|---|
| T1 | counts, areas, nodes, weights and kinematic values identical; HST PM errors to 5e-10 rel. | yes |
| T2, T3 | 7.6e-10 | yes |
| T4 | BH 8200: −0.14; **BH 4e4: +3.9** | **no** (BH 4e4) |
| T5, T6 | pytest | yes |
| T7 | HST PM bins 13.8–58.8″ move by up to **1.85σ** at 4″ offsets | **no**: re-centring must extend beyond 15″ |
| T8 | counts-only window, BH 4e4 − no BH: unbinned +8.6…+9.6, binned +4.8…+8.9 (all σ_c, n) | same sign, yes |
| T9 | our MUSE ML σ is 0.45–0.62σ below published in the 3 inner bins | marginal (only differences used) |

**−2Δ ln L relative to no BH (full likelihood; marginal; s.d. over centres)**

| | AvdM10 only | σ_c = 1″ | 2″ (n = 6) | 2″ (n = 10) | 3″ |
|---|---|---|---|---|---|
| BH 8,200 | −5.5 | −6.0 (0.3) | −5.6 (0.7) | −5.8 (0.7) | −5.7 (0.5) |
| BH 4×10⁴ | −0.5 | −4.4 (2.0) | −5.1 (5.6) | −1.2 (5.4) | +2.1 (3.5) |
| BH 4×10⁴, counts 2–59″ only (unbinned) | — | +9.5 (1.4) | +9.4 (2.8) | +9.6 (2.7) | +8.6 (3.8) |

**Reading**

- **BH 8,200 is stable.** Its gain over no BH (−2Δ ln L ≈ −5.7) does not depend on the centre;
  it comes from the kinematics (counts-only change −0.5 to −0.9).
- **BH 4×10⁴ is ambiguous by the decision rule.** Overall it is neither favoured nor
  disfavoured (−5 to +2), the per-centre scatter (3.5–5.6) exceeds half of any value, and the
  quadrature is not converged (T4). Its count penalty alone is robust (+8.6 to +9.6, unbinned),
  but the inner PM dispersions, which move by 1–3σ when the centre moves 2″, pull the other way:
  off AvdM10 the innermost radial bin rises (AvdM10 0.760 mas/yr; prior-weighted mean 0.81),
  towards the 4×10⁴ model.
- **The central count deficit is not significant once the centre is marginalised.** Two innermost
  bins: P(N ≤ n_obs | no BH) = 0.027 at AvdM10, 0.30–0.37 marginalised (σ_c = 1–3″). For the
  4×10⁴ model 0.001 at AvdM10, 0.05–0.10 marginalised.
- **Not yet valid as a final answer:** T7 shows that the kinematic bins at 15–60″ also resample by
  up to ~2σ, so the re-centred region is too small, and T4 shows the 4×10⁴ marginal is not
  converged. Both must be fixed before a refit decision.

**Centre-marginalised data** (`MDATA=1 python bin/centre_marginal_plot.py` →
`plots/centre_marginal_20261103_mdata.png`). Per bin: prior-weighted mean over the 100 centres
(σ_c = 2″); error² = mean measurement error² + centre scatter² (law of total variance). For
display only: neighbouring bins are correlated through the shared centre, and the PM/MUSE
measurement errors are the AvdM10 ones. Reading: the inner counts become a smooth flat core,
which no BH and 8,200 M☉ follow and 4×10⁴ M☉ over-predicts by ~1σ in every bin from 2″ to 20″;
the inner PM dispersions become a smooth decline from ~0.82 to ~0.775 mas/yr over 4–12″, which
4×10⁴ M☉ follows and no BH / 8,200 M☉ under-predict by ~1σ per bin. The large-BH question is a
counts-versus-kinematics tension of the DF family, not a centre artefact.

## Revision: model-side marginalisation (2026-11-03)

After the Codex and Fable reviews (`docs/reviews/{codex,fable}_recentre_plan_20261103.md`) the data-side
construction above is replaced. Re-binning the stars about every trial centre averages likelihoods of
*different* data summaries; on the first run the centre contributions had an effective sample size of
1.3-3.3 of 36-100 nodes for every model, i.e. each "marginal" was its luckiest re-binning (this is the T4
failure). Now every product stays frozen about AvdM10 and the **model** centre moves
(`src/ocen_dm/centre_model.py`, `DFJointProblem.evaluate(factors=...)`): counts = model surface density
averaged over each bin's covered sky; HST PM = per-star mean of the model variance projected on the data's
(AvdM10) radial/tangential axes; MUSE = per-star mean of the LOS variance; streaming terms fixed (they belong
to the fixed data); Gaia PM bins not moved (checked). c = 0 reproduces the frozen likelihood exactly.

**Mock test** (`bin/centre_model_mock.py`; real star positions and errors, analytic velocity profile steeper
than our DF models, true centre at c, our estimators run in the fixed AvdM10 bins, common random numbers;
`plots/centre_model_mock_20261103_inner.png`, 32 seeds, bins 3.2-17.6"; `..._20261103.png`, 4 seeds, to 36"):

| true centre | HST PM: max |estimated - predicted change| / bin error | MUSE |
|---|---|---|
| (+2, 0)" | 0.06 (s.e. 0.03-0.11) | 0.01 |
| (-1.5, +1.5)" | 0.06 | 0.09 |
| (0, -4)" | 0.48 in the two innermost bins (3.2-5.2"), predicted signal 1.2-1.7; <= 0.17 elsewhere | 0.01 |

Reading: the prediction is right to < 0.1 bin error for offsets up to ~2" (1 sigma_c) in every bin. At 4"
(2 sigma_c) our mixture estimator responds ~30 % less than predicted in the two innermost PM bins (stars at
r_c < 1" with very large model variance are partly absorbed by the field component / down-weighted by the
single-Gaussian fit). At c = 0 the same estimator reads 0.2-0.5 bin errors low in the innermost bins of this
steep mock profile (baseline bias, not centre-specific).
