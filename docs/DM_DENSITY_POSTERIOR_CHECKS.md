# Posterior checks of the high-density solutions

Two production runs were launched on 22 September 2026 at 06:52:53 UTC
(07:52:53 BST). Both use the flexible rung-2 cored halo and the original 89
kinematic measurements. The purpose is to check whether the high-density
solutions found by [constrained optimisation](DM_DENSITY_PROFILE_LIMIT.md)
occupy little posterior probability or were insufficiently explored.

| Run | Density at 20 pc | Parameters | Live points | Slice steps | Seed |
|---|---|---:|---:|---:|---:|
| Independent repeat | Free | 12 | 800 | 48 | 1729 |
| Conditional posterior | Fixed at 2 solar masses per cubic parsec | 11 | 800 | 44 | 1730 |

Both use UltraNest 4.5.0, dlogz = 0.25, and one numerical-library thread.
The original run used 400 live points, 24 slice steps, dlogz = 0.5 and seed 42.
The new runs draw their initial live points from the priors; the optimised
solutions are used only for preflight likelihood checks. There is no imposed
production likelihood-call limit. The original stellar-mass floor of one
million solar masses remains in both runs; the lower-floor diagnostic is not
adopted here.

## The conditional prior

Write M = M_DM(<100 pc), s = r_s, and rho = M/A(s). The original prior is
log-uniform in M and s, with independent priors on the other parameters.
Changing variables from M to rho gives

```text
p(rho, s) = p(s) / [rho ln(M_max/M_min)]
```

inside the original mass support. At rho = 2, every allowed scale radius
(3–1,000 pc) remains inside the original halo-mass range (1,000–30 million
solar masses). A conservative integral envelope bounds the derived M100
between 464,803 and 7,266,699 solar masses throughout that radius range.
Consequently p(s | rho=2) is exactly the original log-uniform radius prior,
and all other parameter priors are unchanged. No extra Jacobian factor belongs
in the likelihood. The adapter refuses densities for which full radius support
cannot be established; it never silently changes the conditional prior.

Halo scale radius, stellar and remnant masses, remnant size, black-hole mass,
anisotropy and distance remain free. Halo normalisation is derived at each
sample to maintain the local density constraint. Once both runs finish:

```text
p(rho=2 | data, core) = p_prior(rho=2 | core) Z_conditional / Z_free
p_prior(rho=2 | core) = 1 / [2 ln(30,000,000 / 1,000)]
```

This is a posterior PDF height per unit density, conditional on the cored
model. It is not a probability atom at 2, an integrated upper-tail probability,
or DM/no-DM model odds. The no-DM and cusp runs are not repeated in this batch;
revised model-averaged limits would combine the new result with their existing
evidences and must state that limitation.

## Preserved runs and outputs

The [launch manifest](../results/fits/densitycheck_20260922T065245Z/launch.json)
records commands, process IDs, checksums, environment and job configurations.
The batch directory contains frozen source, the common observation snapshot,
logs and separate worker status files. At launch both workers were confirmed
alive and sampling, with PIDs 77776 and 77778.

Output directories under `results/fits/`:

* `rung2_K2_cored_turnover_freecheck_densitycheck_20260922T065245Z`
* `rung2_K2_cored_turnover_rho2_densitycheck_20260922T065245Z`

Each completed run will contain the standard posterior, summary, model/data
snapshots and radial profiles, plus `density_posterior.ecsv`,
`density_diagnostics.json` and labelled slice-walk statistics. The density table
adds local rho20 and, for the conditional model, the derived halo mass. Once
both finish, the worker writes `density_comparison.json` in the batch directory
with the evidence ratio and implied posterior PDF at density 2.

The launcher refuses existing batch or fit directories. Source and parent-input
checksums are checked by each worker. No earlier fit or active dynamical run was
restarted or overwritten.

## Validation and assessment after completion

All 71 targeted tests passed, including conditional-prior consistency,
halo-density reconstruction, likelihood equivalence, saved-model replay,
sampler configuration and refusal to overwrite existing runs. The production
preflight replayed the archived maximum likelihood, tested 128 prior draws per
model, and recovered the profile-optimisation likelihoods in both models.
A separate conditional smoke run exercised UltraNest, all output products,
saved-run replay and the derived density table. Its 500-call budget was only
an execution test and supplies no scientific inference.

After completion, compare the new free-density posterior and evidence with the
original core run; inspect upper-tail counts, the stellar/remnant mass tradeoff
and slice-walk distances; then compare its density near 2 with the conditional
evidence calculation. Agreement would support the interpretation that the
good high-density fits occupy little integrated probability. Disagreement
would motivate further sampling checks before changing the quoted limit.

The launcher is `bin/run_density_posterior_checks.py`. To regenerate only the
cross-run comparison after both outputs exist:

```sh
PYTHONPATH=src python bin/run_density_posterior_checks.py compare results/fits/densitycheck_20260922T065245Z
```
