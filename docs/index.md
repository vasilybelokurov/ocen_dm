# oCen_dm: dark matter in Omega Centauri

**Question.** Can Omega Centauri's tidal debris distinguish an extended
dark-matter halo from the stars, remnants and central mass that fit its
internal kinematics?

The project has three connected strands:

1. **Cluster mass model.** Infer the present mass distribution from the
   internal kinematics and star counts with action-based distribution
   functions (DFs). The [base no-DM model](kb/models/base_no_dm.md) is the
   reference for every dark-matter variant.
2. **Dark-matter survival and tails.** Estimate how much dark matter could
   survive the progenitor's disruption and what the tidal tails of no-DM and
   DM models look like (N-body).
3. **Stream 54 (Fimbulthul).** Model the stream in a barred Milky Way with the
   LMC, as the external test of the cluster's mass and orbit.

A joint kinematic-and-tail inference is the long-term goal and is not built.

## Status (2026-10-17)

| Strand | Where it stands | Read |
|---|---|---|
| Cluster mass model | A two-transition DF with Poisson star counts fits the 89 kinematic bins without dark matter (χ²_kin = 156), but fails its own data-fit gate. The data fix the total mass profile (±3% at 20 pc) but not the split between stars and dark mass. Fits with a fixed central black hole (8,200 and 4×10⁴ M☉) also fail the gate and do not measure its mass. A refit with completeness-corrected counts and a tied HST/Gaia amplitude ratio changes the predictions by a few per cent (χ²_kin = 159). With those counts a free cored halo goes to zero; Jeans fits with NFW or cored halos show no preference for dark matter with either tracer. The ρ₂₀ profile, mass-per-star and mock scans still use the old (schema-1) counts. | [Base model](kb/models/base_no_dm.md), [DF + black hole](kb/models/df_bh.md), [DM fits](kb/models/dm_fits.md), [write-up (PDF)](df_observed_fit.pdf) |
| DM survival and tails | N-body tails of a no-DM and a DM model on Omega Cen's orbit agree to 10% in kinematics; the DM model loses 1.3–1.6× fewer stars. Neither difference is testable with published tail data. | [Results](kb/results.md) |
| Stream 54 | Barred models beat the best axisymmetric model by about 3,300–4,100 in ln L, but none reproduces the steep proper motions of the turn. The LMC shifts the distal proper motions by a few tenths of a mas/yr but is absorbed by the cluster's PM. A paper draft compares the barred models with the axisymmetric null. Fits to the binned tracks disagree with the star-by-star likelihood. That likelihood (legacy score) did not converge with particle number; it now uses fixed χ bins with a jackknife bias correction, and a 9-D emulator training set is running. | [Stream status](kb/stream_status.md), [likelihood versions](kb/models/stream_likelihood_versions.md) |

## How to navigate

- **[Getting started](kb/getting_started.md)**: environment, CLI, tests, conventions.
- **[Data products](kb/data/index.md)**: every dataset, its source, how it is
  processed and where it is used; [registered sources](kb/data/sources.md);
  [known gaps](kb/data/gaps.md).
- **[Models](kb/models/index.md)**: the base no-DM DF model, the dark-matter
  variants and the earlier Jeans ladder.
- **[Results](kb/results.md)**: current findings by strand, with the reports
  they rest on.
- **[Literature](kb/literature.md)**: every paper the project cites, by topic,
  mirrored in the public [SciX library](https://scixplorer.org/public-libraries/RSXq3-NtTTCApvfPUbdwrg).
- **[Reviews](reviews/index.md)**: external (Codex) reviews of code and plans.

Code: <https://github.com/vasilybelokurov/ocen_dm>
