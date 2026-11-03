# Codex review: centre marginalisation plan -- 2026-11-03

Claude verified: reference fit uses hst_pm_*_ours (dftie_restart_20261016 problem files); per-star MUSE vlos/vlos_err/selection_hq_los in catalog_and_selections.dat; build_count_profiles.py full coverage inside 30"; m=1 variant is counts-only (fit_centre.py:50).

## Brief
# omega Cen centre: state of the project and a plan to marginalise over the centre -- critique requested (2026-11-03)

Repository: oCen_dm (read-only for you). Relevant code: src/ocen_dm/centre_fit.py, bin/fit_centre.py, bin/centre_mocks.py,
bin/check_count_centre.py, bin/hst_density_centre.py, bin/build_count_profiles.py, src/ocen_dm/kinematics/hst_profile.py,
src/ocen_dm/kinematics/df_fit.py, docs/kb/data/cluster.md, docs/kb/models/df_bh.md, docs/kb/models/base_no_dm.md.
Earlier reviews of the centre question (by you / the other reviewer): docs/reviews/codex_centre_20261020.md,
docs/reviews/fable_centre_20261020.md. Read them; do not repeat them, build on them.

## Science context
We model omega Cen's mass (stars + Plummer "remnant" component + optional central BH + optional DM halo) with (a) a
self-consistent action DF (AGAMA) fitted to 89 binned kinematic points (published oMEGACat VI HST PM radial/tangential
dispersions and MUSE LOS dispersions, our Gaia PM dispersions outside 300") plus Poisson star counts (HST F625W<19 from
oMEGACat, 20 log bins 2-250", area-weighted; Gaia G<17, 12 bins 300-2400"); (b) spherical Jeans models (JamPy-like) with a
free BH. Centre currently fixed at Anderson & van der Marel 2010 (AvdM10). Units pc, Msun, km/s; 1" = 0.026 pc at 5.36 kpc.

## What has been established (all at fixed model expectations unless stated)
1. Central HST count "dip" at 2-4": for F625W<19, Poisson power-law slope over 2.0-8.5" is +0.33 [+0.16,+0.51] (p 0.06)
   about AvdM10, -0.04..+0.08 (p 0.5-0.9) about the oMEGACat VI kinematic centre. Fainter independent slices (19-23) are
   flat to +-0.07 about either centre. JWST (GO 8322, 4343) sees the same stars; depth/crowding at F625W<19 ruled out.
2. Joint centre fit (centre_fit.py): unbinned counts F625W 19-22 (360k stars) with depth-tied completeness exp(beta_c D),
   plus per-star PM likelihood (374k stars) with an anisotropic velocity ellipse; shared centre; Hessian errors.
   Result (-1.95+-0.38, +1.12+-0.33)" from AvdM10 (east+, north+), 0.7" from the oMEGACat VI kinematic centre.
   Mocks drawn from the model (3 offsets x 15): bias <= 0.13", Hessian errors honest (pull sd 0.6-1.2).
   BUT: the offset grows with fitting aperture (symmetry centres 5-55": (+0.7,-3.0); 60-120": (-1.5,+1.2); 120-200":
   (-2.5,+1.4); counts R<60: (+0.5+-2.0,-1.4+-2.1)), and an m=1 term growing with R absorbs it entirely. The depth term
   changes nothing in the 19-22 slice (beta_c = -0.09). Fable (10-20) found the oMEGACat footprint is depth-asymmetric
   beyond ~100" (W/NW shallow) and argued the W offset is catalogue/exposure-driven. Not resolved.
3. Black hole vs centre: 25 centres within +-2" of AvdM10 plus the kinematic centre (check_count_centre.py, observed counts
   re-binned, model expectations fixed): over the 8 inner count bins (2-13") the DF with a fixed 4e4 Msun BH is worse than no
   BH at every centre by +6 to +19 deviance (median +11.7); 8,200 Msun vs no BH: -0.6 to +1.8. This was a sensitivity scan,
   not a marginalisation.
4. Nothing has been propagated: count products, kinematic bins and all DF/Jeans fits use AvdM10. hst_profile.py fixes
   OCEN_RA/DEC; the DF objective consumes pre-binned radial data.
5. The 89 kinematic bins in the reference fits are the PUBLISHED oMEGACat VI binned profiles (HST PM dispersion and MUSE
   LOS), binned about their adopted centre (oMEGACat adopts AvdM10, per our notes). We also have our own HST PM dispersion
   profile code (hst_profile.py) built from the oMEGACat star catalogue. Whether a public per-star MUSE velocity catalogue
   exists for re-binning: not checked by us.

## Earlier consensus (10-20 reviews)
AvdM10 = reporting origin only. Centre = nuisance, marginalised, for inner measurements: dip tests over a centre prior
(Fable: 1.5"; Codex: an independently estimated distribution, never the centre that makes/removes the dip); IMBH fits with
a free centre (Fable AvdM10 +-1.5" for M <~ 2e4; Codex: also fit BH position separately from the cluster centre); bulk
Jeans/DF: Codex nuisance centre with recomputed radii/projections/selection integrals; Fable: fixed is acceptable for the
bulk (1.5" moves sigma(10") by ~0.5% < 2-3% bin errors) with +-2" robustness rows.

## Proposed plan (Claude) -- critique it
P1. Centre prior. Options: (a) AvdM10 +- 1.5" isotropic Gaussian; (b) the joint-fit posterior (-1.95+-0.38,+1.12+-0.33);
    (c) a broad prior covering both, e.g. centred midway with sigma ~1.5-2". Claude leans to (c) or (a) because (b) may be
    footprint-driven (m=1 degeneracy, point 2).
P2. Marginalise on the DATA side. The models are spherical; a centre shift changes only which stars fall in which
    projected bin. Precompute count products (and our own HST PM dispersion bins) for K centres drawn from the prior
    (K ~ 25-50, or a quadrature grid), only for bins with R < ~60" (outer bins insensitive), outer bins shared. In the
    objective: compute model predictions once per theta, evaluate the likelihood for each centre, L(theta) = sum_k w_k
    L(theta | c_k). Cost: K cheap likelihood evaluations per model evaluation.
P3. Kinematics: switch the inner HST PM dispersion bins from the published profile to our own per-centre bins (same edges);
    keep MUSE published (fixed centre) unless a per-star catalogue is available, and document the inconsistency.
P4. Order of deliverables: (i) the dip test marginalised (cheap); (ii) BH comparison 0 / 8,200 / 4e4 Msun with fixed-BH DF
    refits and a marginalised centre; (iii) Jeans fit with free BH mass and centre as sampled nuisance (2 extra params).
P5. Tests: injection mocks from the fitted DF at a known offset centre (0, 1, 2, 3"), with and without a BH; acceptance:
    BH-mass / Delta-objective bias < 0.3 sigma, coverage of the centre-marginalised interval within +-10% of nominal; check
    that the marginalised result reproduces the fixed-centre result when the prior is a delta function.

## Questions
Q1. Is data-side marginalisation over binned products (P2) statistically correct here, given that the centre is a model
    parameter and alternative binnings share stars? Or should the centre be a sampled/profiled parameter instead (max vs
    sum over k)? Any normalisation subtleties (area-weighted bins near the footprint edge, coverage masks per centre)?
Q2. Which prior (P1), and how to make it independent of the inner data being fitted (avoid double use)?
Q3. Is it worth re-centring the kinematic bins at all, given Fable's 0.5% estimate, or only the inner counts + central
    dispersion bins (R < 10-15")? What about MUSE (P3)?
Q4. BH position vs cluster centre: needed now, or only for a free-BH fit?
Q5. What would change the conclusions (BH count penalty, no-DM result)? Which test would most quickly show the centre is
    irrelevant for our science (DM halo at 10-60 pc, BH at < 0.3 pc)?
Q6. Anything in the current state (points 1-5) you think is wrong or overstated.

Separate what you verified (code/files you opened) from what you assumed. Argue against your own recommendation once.
Keep it under ~900 words, numbered, concrete.

## Answer
1. **Point 5 is wrong:** the frozen reference uses your own HST profiles, `hst_pm_radial_ours` and `hst_pm_tangential_ours`, not the published profiles. P3’s estimator switch has already happened ([snapshot:1](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/results/df/dftwo_counts_area_20261013/observed/data_snapshot.json:1>)). Per-star MUSE inputs also exist locally: `vlos`, `vlos_err`, and `selection_hq_los` ([catalogue:44–75](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/data/raw/omegacat_vi_kinematics/catalog_and_selections.dat:44>)). Their suitability for reproducing the published estimator remains unverified.

2. **Q1: reject P2 as an exact likelihood.** Sharing stars across alternative hypotheses is legitimate; multiplying alternative binnings would double-count them. The deeper problem is comparing likelihoods of different, parameter-dependent summaries. Sampling the centre instead of summing does not fix that.

   Use one fixed sky window and either unbinned positions/Cartesian velocities or fixed sky pixels. For counts,
   \[
   \log L=\sum_i\log\lambda(x_i\mid\theta,c)-\int_W\lambda(x\mid\theta,c)\,dA.
   \]
   Solve the spherical DF once per θ, then interpolate its radial predictions for each centre. Integrate selection over the fixed footprint; include stars crossing the proposed inner exclusion and inner/outer split.

   Marginalisation is `logsumexp(log w + log L)`; profiling is maximisation and answers a different question. Do not exponentiate the existing deviance objective: its discarded constants change with re-binning. The code separately exposes normalized count log-likelihoods and a **profiled**, not nuisance-integrated, likelihood ([counts.py:178–196](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/counts.py:178>), [df_fit.py:176–185](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/df_fit.py:176>)). Area weights must follow each displaced footprint intersection; the builder currently assumes full coverage inside 30″ ([build_count_profiles.py:85–108](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/bin/build_count_profiles.py:85>)).

3. **Q2: use (a) provisionally**, explicitly treating 1.5″ per coordinate as an assumed systematic allowance. Show sensitivity to 1″ and 2″. Avoid (b): the joint fit reuses PMs entering the dynamical analysis ([fit_centre.py:38–53](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/bin/fit_centre.py:38>)). A midpoint chosen from that result is also data-dependent. Obtain an independent constraint from external measurements with assessed covariance, or reserve disjoint control stars and account for shared calibration systematics. Otherwise fit the centre jointly once.

4. **Q3: prioritise counts and PMs inside 15″**, then test extending to 60″. The 0.5% estimate cannot establish insensitivity for BH-induced gradients. Recompute PM directions, uncertainties, selection and radial averaging—not just membership: the loader’s projections depend on the centre, and model averaging uses observed radial nodes ([hst_profile.py:125–140](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/hst_profile.py:125>), [likelihood.py:458–470](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/likelihood.py:458>)). For MUSE, reproduce the original profile first, including its fitted rotation treatment ([likelihood.py:177–186](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/likelihood.py:177>)). If unavailable, predict through its original apertures or omit central MUSE bins as a sensitivity check.

5. **Q4: BH position matters at fixed mass too.** Start with explicitly concentric BH comparisons; assess displacement before interpreting them as general BH constraints. An independently displaced BH breaks spherical symmetry and cannot be implemented merely by reassigning stars. The current potential construction supplies a central point mass without positional parameters ([regularized_df.py:618–655](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/kinematics/regularized_df.py:618>)).

6. **Q5: run a paired robustness experiment first:** refit no-BH/BH and no-halo/halo branches at AvdM10, the joint-fit centre, and several prior-boundary centres; repeat excluding all measurements inside 15″. Compare BH preference, halo limits, and enclosed mass at 10–60 pc. Require shifts below a predeclared scientific tolerance.

   The current count scan holds expectations fixed ([check_count_centre.py:40–56](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/bin/check_count_centre.py:40>)). Refitting, displacement, or relaxing the imposed stellar cusp could change its penalty; that cusp is explicitly a DF-family assumption ([df_bh.md:25–28](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/docs/kb/models/df_bh.md:25>)). Small direct outer-radius effects do not exclude indirect mass–anisotropy changes.

7. **Verified versus assumed / Q6:** I read both earlier reviews, the cited implementation, tests and frozen outputs; I did not execute fits or tests. The numerical joint estimate is recorded, but it is a maximum-likelihood fit with Hessian errors, not a sampled posterior ([centre_fit.py:176–224](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/src/ocen_dm/centre_fit.py:176>)). “m=1 absorbs it entirely” overstates the evidence: that variant fits **counts only**, not the joint likelihood ([fit_centre.py:50–53](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/bin/fit_centre.py:50>)). My prior width and recommended radial cutoffs are assumptions requiring validation.

8. **Strongest argument against my recommendation:** a calibrated approximation using re-binned products might reproduce the scientific answer adequately at much lower implementation cost. A full catalogue likelihood can introduce worse completeness and velocity-model assumptions. Prefer the approximation if comparison against a fixed-data likelihood demonstrates negligible bias.

9. **Checks needed:** delta-prior recovery; centre-grid convergence; fixed-window star conservation; displaced mask integrals; MUSE baseline reproduction; and end-to-end BH/halo injections with asymmetric completeness, PM-error fields and profiles outside the fitted family. Existing centre tests use matched synthetic models, while mocks largely regenerate the fitted family ([tests:23–49](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/tests/test_centre_fit.py:23>), [centre_mocks.py:40–58](</Users/vasilybelokurov/IoA Dropbox/Dr V.A. Belokurov/Code/oCen_dm/bin/centre_mocks.py:40>)). Use hundreds of realizations per regime with binomial coverage intervals; “Δ-objective bias <0.3σ” needs a defined reference distribution. Nothing was modified.