# Dynamical experiments

The `ocen nbody` commands prepare and evolve two independent experiment families:

* **Remnant survival:** a compact DM distribution around a rigid or live stellar nucleus.
* **Progenitor disruption:** an extended DM halo and stellar body around that nucleus.

Both use a spherical equilibrium in the combined satellite potential. Live runs use
gyrfalcON self-gravity for the particles. The default nucleus is an analytic Plummer
potential; `nucleus.live: true` replaces it with sampled, self-gravitating Plummer stars.
The nucleus contributes gravity once in either representation. Frozen controls
use AGAMA orbit integration in the same initial satellite potential, held fixed about
the same nucleus track. Particle masses in the frozen runs are diagnostic weights.

The first 24 completed pilot evolutions are assessed in the
[21 September batch analysis](DYNAMICAL_BATCH_ANALYSIS.md), including the matched
isolation references, repaired softening control and progenitor resolution checks.

The [long-term evolution plan](LIFETIME_EVOLUTION.md) extends these pilots to
multi-Gyr histories and makes the surviving stellar cluster an observational
constraint. The [active lifetime batch](LIFETIME_BATCH_REPORT.md) adds a responsive
nucleus, measured stellar centres, physical heating estimates and numerical gates.

Rigid-nucleus experiments prescribe its orbit. A live nucleus responds to particle
self-gravity and the imposed host, and its centre can depart from the reference orbit.
The code does not couple orbital decay to a live host wake or physical dynamical
friction. A self-consistent accretion history needs those additional dynamics.

## Commands and products

The staged lifetime driver is `bin/run_lifetime_batch.py`. Its `init` action creates
a new batch with archived source and observational inputs; `execute` runs up to four
workers by default and advances only after numerical checks. The current batch is
`results/dynamics/lifetime_20260921T162000Z`; **do not initialise or execute it again**.
Inspect `batch.json`, `queue.log`, the per-case logs and gate JSON files. Refresh
its report with `PYTHONPATH=src python bin/report_lifetime_batch.py
results/dynamics/lifetime_20260921T162000Z` on one command line.

Configure a live nucleus with, for example,
`nucleus: {live: true, mass_msun: 3550000, scale_pc: 5.4, n_particles: 100000,
softening_pc: 0.1, beta: -0.1}`. The signed-DF validation applies to its stars as
well as the DM. `orbit.end_at_present: true` rewinds and replays the nominal
present phase-space point for class 1 or class 3; class 3 must end at the stored
bar history's exact endpoint. This option checks the reference orbit, while
the live stellar endpoint displacement is measured separately.

From the project root, using the existing environment:

```bash
source ~/Work/venvs/.venv/bin/activate
export PYTHONPATH=src
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
python -m ocen_dm.cli nbody --help

# A short integration check, explicitly configured at low particle count.
python -m ocen_dm.cli nbody prepare configs/dynamics/smoke.yaml --label smoke_01
python -m ocen_dm.cli nbody run results/dynamics/smoke_01 --engine frozen
python -m ocen_dm.cli nbody run results/dynamics/smoke_01 --engine live
python -m ocen_dm.cli nbody plot \
  results/dynamics/smoke_01/evolution_frozen \
  results/dynamics/smoke_01/evolution_live --output plots/nbody_smoke_01.png
```

`prepare` performs the DF checks, samples particles and saves the full orbit, without
launching an evolution. `run` evolves those saved inputs and computes diagnostics.
The same prepared model can have both live and frozen evolutions. Existing preparation
directories and evolution labels are refused. A fresh `--label` on `run` preserves an
earlier attempt; there is no automatic resume.

The preparation contains the resolved `config.yaml`, `initial_local.npz`, `orbit.npz`,
potential files, DF and density checks, `equilibrium.json`, a source archive and `run.yaml`.
The manifest hashes every prepared input. Evolution refuses missing or altered inputs.
Host assets are copied into the run, including the class-3 bar histories; subsequent
integration does not read the original external potential files.

Each evolution has its own manifest and source archive. Live runs also record executable
paths and hashes, exact arguments, compilation and integration logs, `energy.log`, and
the `snapshots.nemo` stream. Frozen snapshots are stored in `snapshots.h5`. Both produce
`diagnostics.ecsv`. Interrupted or failed runs retain their inputs and failure status.
If integration finishes but analysis fails, `nbody analyse <evolution-directory>` can
retry the analysis after the problem is corrected.

## Initial conditions and matched comparisons

The density of each live component is

\[
\rho(r)=\rho_0 (r/a)^{-\gamma}(1+r/a)^{\gamma-b}
         \exp[-(r/r_{\rm cut})^2].
\]

Here `gamma` is the central slope and `outer_slope` is b. A `gamma: 0` core has finite
central density; its radial form differs from a Plummer core. The cutoff is a smooth
**density** taper. It is not an instantaneous Jacobi boundary or a DF energy cut.
The templates specify finite total masses or aperture masses, not cosmological M200.

The rigid nucleus has mass 3.55e6 Msun and Plummer scale 5.4 pc. Each component's DF
is inverted in the potential of that nucleus plus **all** particle components.
The nucleus is counted once. Constant `beta` belongs to that particle component;
the stellar anisotropy in the rung-1 fits is a separate parameter.

Before sampling, preparation computes a signed constant-beta inversion independently
of AGAMA, which clips negative DF values internally. A negative contribution exceeding
0.5% of the absolute inversion integral fails preparation. The signed energy grid and
its radial domain are saved: passing this finite grid is a numerical check, not a
proof over all phase space. Preparation also requires density recovery within 5% and
beta recovery within 0.03 on a separate radial grid. These are rejection tolerances;
the measured errors are recorded for each model.

| Templates | Matching convention |
|---|---|
| `remnant_cusp.yaml`, `remnant_core.yaml` | M_DM(<20 pc)=1e5 Msun, a=20 pc, cutoff=100 pc and beta_DM=0; total DM mass differs |
| `progenitor_cusp.yaml`, `progenitor_core.yaml` | Same total DM mass 1e9 Msun, stellar mass 1e8 Msun, scales, cutoffs and beta_DM=-0.3; central DM mass differs |

The broad **isotropic** progenitor core fails the signed DF check in this particular
nucleus-plus-dwarf potential. The beta_DM=-0.3 version passes, so the paired progenitor
templates both use that value. Their stellar body has beta=0. This does not imply that
every cored halo requires beta=-0.3.

All four templates start in isolation. They are pilot specifications, not established
resolution choices. Inspect `equilibrium.json` for mass ratios, particle mass ranges and
effective particle counts in the central apertures before choosing production resolution.

## Particle refinement and numerical controls

Without `refinement_radius_pc`, a component uses equal masses. With refinement, AGAMA
samples the DF multiplied by

\[
p(r_{\rm peri})=\max[p_{\min},\min(1,(r_{\rm protect}/r_{\rm peri})^2)].
\]

Particle weights are divided by p and normalized to the specified component mass.
Pericentres are computed in the initial **combined spherical potential**. This protects
orbits that pass through the centre even when sampled near apocentre. It avoids the
phase-dependent bias of selecting particles by their instantaneous radius alone.
Particle masses remain fixed during evolution. Softening scales as m^(1/3), with
`softening_pc` the minimum sampled softening; the kernel is falcON P1. The analytic
Plummer nucleus is an external potential, not a softened falcON particle.

Tides can subsequently change pericentres. Diagnostics therefore report the mass fraction
inside the protection radius contributed by particles heavier than ten times the lightest
particle of that component. Repeat with a larger protection radius, lower probability
floor, more particles and smaller softening to assess contamination and convergence.

The continuum DFs use unsoftened gravity. Isolation runs must measure the response to
particle sampling and softened self-gravity before interpreting tidal evolution. Vary
particle count, softening and time step separately, and compare to the matched frozen
control. The low-count smoke runs verify integration and I/O, not long-term equilibrium
or converged remnant densities.

## Host, orbit and clock

User-facing units are pc for satellite scales, kpc for Cartesian positions, km/s for
velocities, Msun for mass and Myr for time. The integrators use kpc/(km/s), explicitly
converted with 1 unit = 977.7922216807892 Myr. The same G is supplied to gyrfalcON and
the AGAMA external-potential plugin.

* `isolation`: stationary nucleus at the origin, no host.
* `class1`: McMillan17. By default, start at the nearest past apocentre of the nominal
  present-day Omega Cen orbit and integrate forward. `initial_xv` can supply a different
  six-element Cartesian state in the same frame.
* `class3`: Hunter24 with the existing growing, slowing bar. Rewind the nominal present
  orbit to `start_time_myr`, then evolve forward from that state. The absolute bar clock
  is preserved; its present epoch is 8 natural time units. Preparation requires the
  existing `oCen_bar` potential files, or the location specified by `OCEN_BAR_DIR`.
* `plummer_host`: a portable analytic host for integration tests.

For a class-1 pilot, copy a template and change its orbit block, for example:

```yaml
orbit:
  kind: class1
  duration_myr: 100.0
  sample_step_myr: 0.02
```

Inspect the saved track to select an apocentre-to-apocentre interval when a single-passage
comparison is required. Increasing the duration provides repeated passages. Class 3 uses
`start_time_myr` in its original bar clock, not lookback time. A model starting at 2000 Myr
therefore encounters the bar state at 2000 Myr, not its initial growth state.

Trajectories store all six Cartesian coordinates in float64. Moving potentials use Hermite
interpolation with the velocities as derivatives. Preparation compares the reconstructed
acceleration with the host force and rejects errors above 1%. It also checks that saved
host assets reproduce the original forces. The old radius-only orbit summaries are not
accepted as simulation inputs. Reordering a backward track would preserve velocity signs.

The live maximum step is rounded downward to a binary step; duration is rounded upward
to a whole number of these steps. Actual values are saved. Frozen controls use adaptive
AGAMA integration over the same time interval. Both run in the host's inertial frame.

## Diagnostics and interpretation

For DM and stars separately, diagnostics include total and bound-mass proxies, M(<r),
shell density, effective particle numbers, radial dispersion and beta at 3, 10, 20, 35,
50 and 70 pc. The radii are configurable. Density uses a shell from r exp(-0.1) to
r exp(0.1); the estimate and its effective particle count must be read together.
Missing velocity moments are NaN when fewer than three effective particles are present.

The bound-mass proxy iteratively removes positive-energy particles from a self-excluded
spherical estimate of the particle potential, adding the analytic nucleus only when
it is rigid. It measures energy relative to the prescribed rigid nucleus or the
shrinking-sphere centre and inner mean velocity of the live stars. It is approximate for
distorted systems and is not a conserved Jacobi energy. Total-aperture and bound-aperture
measurements are both retained. No particles are removed from the simulation.

Live-nucleus diagnostics include bound stellar mass, 3D and three projected half-mass
radii, velocity dispersions and offsets from the reference orbit. Projected half-mass
radius is a half-light radius only under constant M/L. The batch also saves initial
and final projected kinematic comparisons in `stellar_observables.ecsv`, with missing
values in bins below 30 effective particles. Complete annular selections and Cartesian
views make these diagnostic comparisons; no observational survival cut is imposed.

The instantaneous tidal scale uses the largest eigenvalue of the host force gradient
plus the centrifugal matrix for omega = x cross v / |x|^2. It excludes the Euler and
velocity-dependent Coriolis terms and is an explicitly defined diagnostic, not an escape
surface. The implementation recovers the circular spherical coefficient 3-dlnM/dlnR.

Energy-change columns use the **fixed initial satellite potential**, with mass-weighted
medians grouped by component and initial radius. They are suitable for the frozen shock comparison and a reference diagnostic
for live runs; they do not equal the energy change in the evolving live potential.

## NEMO setup and validation

The existing `bin/nemo_env.sh` configures this machine's NEMO installation. An alternative
script can be selected with `--nemo-env` or `OCEN_NEMO_ENV`. A run-local no-op manipulator
requests particle IDs during input; `give=k` alone does not retain them in this build.
It changes no forces or particle properties. It is built with `OCEN_NEMO_CXX` (default
`g++-15`) and `OCEN_NEMO_PRECISION` (default `SINGLE`, matching the installed build).
These must match another installation's compiler ABI and precision.

New live evolutions also build a run-local falcON library using
`bin/build_nbody_kernel.py`. The installed single-precision P1 Taylor expansion can
overflow an unnecessary higher derivative at small separations and softenings. This
caused the 0.1-pc-softening cusp control to abort with a NaN position. The repair
evaluates the same coefficients as
`D_n = T_n * [1 + (2n+1) * (epsilon^2/2)/(r^2+epsilon^2)]`, avoiding the overflowing
intermediate. It preserves the P1 force law, masses, softening and timestep.
Coefficients that genuinely exceed the native precision produce an explicit error.

Only `kernel.o` is replaced in the run-local library; the installed NEMO files are
untouched. The builder checks the upstream P1 source block and links the installed
objects, so those build objects must be available and ABI-compatible. Its compiler
defaults to `clang++` on macOS and `g++` elsewhere; `OCEN_NEMO_LIBRARY_CXX` can select
the compiler used for that installation. An unfamiliar kernel implementation fails
the source check and needs review. This workflow has been validated on this Mac.

`runtime/kernel_build.json` and the evolution manifest record the patch, compiler,
source/object hashes and library checksum. The library and patched source stay with
the evolution; subsequent NEMO calls verify the library checksum. The source archive
includes the builder and coefficient helper. A compact 64-particle regression case
crashes with the original library and passes with the repair, with forces agreeing
with direct P1 summation to better than 0.01%.

Snapshot parsing checks unique IDs, particle masses, monotonic times and the final epoch.
Times are read from full-precision table columns because NEMO's standalone time header
prints only six significant digits. Frozen controls also use float64 trajectory output.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OCEN_RUN_NEMO_TESTS=1 \
  python -m pytest tests/test_dynamics.py tests/test_nbody_kernel.py -q
```

The tests include the analytic Plummer DF, circular tidal-radius limits, rejection of an
inadmissible isotropic core, weighted sampling, density and anisotropy recovery, portable
host assets and absolute clocks, frozen energy conservation, and live NEMO units and IDs.
The NEMO integration test is opt-in for environments without the local compiled tools.

The implementation follows the [AGAMA N-body example](https://github.com/GalacticDynamics-Oxford/Agama/blob/master/py/example_nbody_simulation.py)
and [gyrfalcON interface](https://teuben.github.io/nemo/man_html/gyrfalcON.1.html), with local
checks of the installed APIs. Hard DF energy truncation with potential iteration,
mass-dependent friction, responsive/live nuclei, and fitting the present-day baryonic
potential into the IC family remain separate extensions.
