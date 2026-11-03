# Fable review: centre marginalisation plan -- 2026-11-03

Same brief as docs/reviews/codex_centre_plan_20261103.md. Fable ran short read-only checks (raw error-deconvolved dispersions
per centre; not the mixture estimator). Claude verified independently: reference fit uses hst_pm_*_ours; per-star MUSE
columns exist; build_count_profiles.py coverage rule at 30".

## Answers (Fable)
**Q1.** Admissible as a nuisance integral, with corrections: (i) sum (logsumexp), not max; (ii) combine `loglike`
(counts.py:178-186), not `deviance` (df_fit.py:178,183 builds the objective from deviance, whose saturated term differs
between centres); the amplitude is profiled per centre by `compare()` -- state it, or fix it via mass_per_star; (iii) exact
only if the binning information loss is centre-independent. Clean alternative: unbinned Poisson positions in R<60"
(~13k F625W<19 stars), lambda(x) = a Sigma_DF(|x-c|), reusing centre_fit.py:89-117. Coverage: for |dc|<=3" and R<60"
coverage >= 0.99, but rebuild node weights per centre or extend the 30" rule to 60" with a check.

**Q2.** Prior from data not in the fit: F625W 19-22 counts at R<60", (+0.5+-2.0, -1.4+-2.1)" -> isotropic sigma 2" near
AvdM10. Not the joint posterior (aperture-dependent, m=1-degenerate, and its PM part uses the hq stars of the kinematic bins).

**Q3.** Re-centre inner bins: the 0.5% estimate is the smooth-model change; re-centring resamples the noise of 100-400-star
bins. Raw dispersion proxy, hq stars: 3.24-4.13" 0.782+-0.028 -> 0.891+-0.032 mas/yr at the joint centre (0.866 kinematic);
5.25-6.69" 0.836 -> 0.776; 6.69-8.5" 0.777 -> 0.803 (1-3 sigma per bin out to ~8", <1 sigma beyond). MUSE 0-7.9" (55 hq stars)
20.9 -> 19.7 km/s (0.6 sigma). Re-centre HST PM bins with r_upper <~ 15" and MUSE bins < 11.35"; bin MUSE as the full second
moment about v_sys (streaming2 = None).

**Q4.** BH position: not now (influence radii 2.7" for 8,200, 11.5" for 4e4; innermost fitted bins 3.24" HST, 0-7.9" MUSE;
<=1" offset enters at second order). Needed only for a free-BH fit using < 3" or per-star data.

**Q5.** DM at 10-60 pc: 2" is <= 0.7% in radius, second order; irrelevant. BH count penalty: open risk is refit absorption by
the remnant Plummer. Fastest test: evaluate the saved no-BH / 8,200 / 4e4 models on products rebuilt about K prior centres;
if the marginal penalty stays > 5 and the centre scatter < half of it, no refit needed.

**Q6.** (a) Brief point 5 wrong (kinematic HST bins are ours; only MUSE published; per-star MUSE local). (b) The dip: inner
two bins model mu 10.8 and 17.5 vs 6 and 12 observed (12 and 10 at the kinematic centre) -- Poisson resampling of 6-12-star
bins; a p = 0.06 slope should not be called a dip. (c) Deviance comparison at one centre valid; summing across centres not.
(d) "fixed acceptable" holds for the model, not the inner data bins.

## Revised plan (Fable)
1. Prior: 19-22 counts, R<60", counts only; sigma ~2".
2. Products per centre (K = 32, 2-D Gauss-Hermite): HST counts R<60" (coverage per centre), HST PM bins r_upper < 15", MUSE
   < 11.35"; outer shared. Test: delta prior reproduces current products bit-for-bit; coverage change < 0.5% per bin.
3. Objective: per-k loglike, logsumexp. Test: K=1 equals the existing objective to 1e-9; 10x finer K changes the marginal by < 0.1.
4. Fixed-theta evaluation of the three saved BH models; report the marginal penalty and its centre scatter.
5. Only if 4 is ambiguous: refit no-BH and 4e4 with the marginalised objective; mocks at offsets 0-3" with/without BH:
   bias < 0.3 sigma, coverage within 10%.
6. Replace the dip slope by the marginalised two-bin Poisson probability.

## Against (Fable)
The region (R<8", ~700 PM stars, 55 MUSE stars) carries little M_BH information compared with the remnant-Plummer dependence;
the likely outcome is "centre irrelevant, remnant model dominant", obtainable from a +-2" robustness row alone.

## Verified vs assumed (Fable)
Verified: check_count_centre.py, build_count_profiles.py:59-114, counts.py:43-196, df_fit.py:141-185, likelihood.py:274-305,
437-499, hst_profile.py:110-140,168-216, centre_fit.py, catalogue columns (vlos finite for 307,030), saved count-centre JSON,
inner mu, re-centring numbers (raw estimator). Assumed: dftie_restart_20261016 uses the same five datasets; raw proxy tracks
the mixture estimator; m=1 claims from the brief.
