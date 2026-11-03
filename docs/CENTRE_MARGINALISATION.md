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

(filled in after the run)
