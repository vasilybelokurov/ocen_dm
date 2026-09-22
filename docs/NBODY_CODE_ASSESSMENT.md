# N-body code for the host-disruption simulations: assessment (2026-09-20)

> **Status, 2026-09-21:** this is a design/tool assessment, not a production simulation.
> The current first experiments are [N0–N6](SECTION11_TESTS.md), beginning with a consistent
> tidal radius and truncated equilibrium. Hardware timings below are historical measurements.
> See the [audit](CODE_ANALYSIS_AUDIT.md) for limitations that remain unresolved.

> **Implementation update:** use the in-repository
> [dynamical experiment pipeline](DYNAMICAL_EXPERIMENTS.md), rather than the old
> `satellite_experiment` scaffold. It includes joint equilibrium sampling, live and
> frozen drivers, portable host inputs and identity-preserving diagnostics.


## What exists locally

| item | state |
|---|---|
| `~/Work/Code/satellite_experiment` | scaffold from 2026-07-03: AGAMA host (analytic Miyamoto-Nagai + Dehnen + NFW placeholders), single-component Plummer satellite ICs, gyrfalcON driver via falcON's `Combined` analytic acceleration, snapshot diagnostics, animation. Only pipeline-test values; no two-component dwarf, no nucleus, no bar, no friction. |
| NEMO + gyrfalcON | built at `~/Work/src/nemo_satellite_experiment` (falcON 3.7, 2025-10-07, gcc). Runs once the library path is fixed (`bin/nemo_env.sh`). |
| AGAMA 1.0.152 | in the venv; provides potentials, self-consistent multi-component ICs, and the NEMO acceleration plug-in (`agama.so`). |
| pyfalcon, Gadget-4, Arepo, Bonsai, PKDGRAV | not installed. |
| hardware | Apple M4 Pro, 10 P + 4 E cores, 48 GB. gyrfalcON is single-threaded. |

## Measured

* gyrfalcON self-gravity: 10^5 particles, 64 steps: 5.3 s (0.08 s/step); 10^6 particles, 16 steps: 11 s (0.7 s/step).
* AGAMA plug-in with the **rotating Hunter+2024 bar** (`MWPotentialHunter24_rotating.ini`), 10^4-particle cluster, 0.1 time units, kmax=8: runs, 0.3 s. Time-dependent potentials (decelerating bar from `bar_migration.py`) are supported by the same route (INI with `scale`/`rotation` modifiers referencing text files).
* Three environment defects fixed without rebuilding anything: relative dylib install names (set `DYLD_LIBRARY_PATH` to `$FALCON/lib` and `$FALCON/utils/lib`); SIP strips `DYLD_*` when launching via `/usr/bin/time`; LLVM libomp leaves a stale `dlerror()` that NEMO's `loadobj` misreads (stub dylib preloaded). All in `bin/nemo_env.sh`.

## Requirements of our problem

1. Host: static axisymmetric (classes 1, 2) **and** time-dependent barred (class 3, slow-bar variants) potentials; semi-analytic dynamical friction for class 2 (falcON has no friction; add as a `UniformAcceleration`-type time series on the satellite's own orbit, with an explicit treatment of the corresponding bulk acceleration; a translating frame alone does not supply physical friction).
2. Satellite: nucleus (3.6e6 Msun, r_h ~ 7 pc) + stars (1e8-1e9) + DM halo (1e10-1e11): three components, four decades in mass, three in size. Needs multi-mass particles with individual softening (falcON supports `eps<0`) and block time steps (`Nlev`, `fea`); the light nucleus particles will suffer heating from heavy DM particles unless the mass ratio is kept modest (~10-30) inside the nucleus's region.
3. 8-10 Gyr, pericentres to 0.4 kpc; nucleus crossing time ~0.3 Myr -> smallest step ~0.02 Myr. With block steps only the ~2e4 nucleus particles take it.
4. Output: bound mass of the nucleus + DM within its Jacobi radius vs time; needs snapshots every ~50 Myr.

## Options

| option | pros | cons | verdict |
|---|---|---|---|
| **gyrfalcON + AGAMA plug-in** (works now) | fast multipole, individual softening, block steps, any AGAMA potential incl. the decelerating bar, NEMO tool chain, satellite_experiment scaffold reusable for orbit/IC/diagnostics | serial; no friction (work-around above); NEMO snapshot handling | **recommended for the first round** |
| pyfalcon (Python force calculator) | same tree code from Python; trivial to add friction and the time-dependent AGAMA potential; full control | must be built from GitHub; Python leapfrog overhead per step (fine for 1e5-1e6 particles, ~0.1-1 s/step); no block steps unless written | second choice; good for the friction runs |
| Gadget-4 + AGAMA patch | parallel (MPI, 10 cores), adaptive steps, mature | needs building (open-mpi, gsl present; hdf5/fftw not), patching; heavier tooling | only if serial falcON is too slow |
| AGAMA "restricted N-body" (`example_tidal_stream.py`) | minutes; nucleus potential updated on the fly | not self-gravitating at large; unsuitable for "how much DM stays bound" | quick look-ahead only |

## Cost estimate for gyrfalcON

Particle budget: nucleus 2e4 x 180 Msun, stars 1e5 x 1e3-1e4, DM 1e6 x 1e4-1e5 (total ~1.1e6). Largest step 1 Myr (1e-3 units) -> 1e4 full-tree steps x 0.7 s = 2 h; nucleus block at 0.02 Myr: 5e5 partial steps at ~0.02 s = 3 h. **~5-8 h per 10-Gyr run**, single core, so 10 runs in parallel on the 10 P-cores. A 1e5-particle pilot runs in under an hour.

## Recommendation

Use the new `ocen_dm.dynamics` driver with gyrfalcON and the AGAMA plug-in. Its
McMillan17 and Hunter24 experiments use prescribed nucleus orbits without friction.
Changing to a translating frame does not supply physical dynamical friction.
Class-2 historical reconstructions still need a specified drag model coupled to
the evolving bound mass, including a treatment of particles that escape. A Python
force-driver interface remains an option when that per-particle control is needed.
