# Compact DF implementation and numerical tests

The compact stellar model is implemented in
[`regularized_df.py`](../src/ocen_dm/kinematics/regularized_df.py).
The independent generalized lowered-isothermal control is implemented in
[`lowered_isothermal.py`](../src/ocen_dm/kinematics/lowered_isothermal.py).
Both predict the same surface-density and three projected second-moment
quantities used by the existing binned likelihood. Existing DoublePowerLaw
models and saved configurations keep their original interpretation.

## Regularized stellar family

The stellar DF is

\[
f=A\exp[-(L_c/J_0)^\alpha],\qquad
b(q)=b_{\rm out}\frac{q^2}{q^2+J_a^2},\quad q=J_r+L.
\]

The default model has five stellar coordinates: mass, action scale, envelope
index, anisotropy amplitude and transition action. An optional second
transition adds an outer amplitude and a larger transition action. Neither
amplitude is identified with the velocity-moment anisotropy beta(r).

The code follows constant-DF contours in log(q) and c=L/q. It tabulates their
circular endpoints and uses a cubic interpolation whose radial boundary
preserves f_Jr/f_L=2. Action-space quadrature includes the signed-Jphi volume
factor and checks its normalization by doubling the integration order.
DF evaluation is positive by construction; no negative-value clipping is
used. Unbound phase-space points return zero.

The exact spherical orbital-frequency reference is the default. AGAMA maps
actions to an orbit, then its direct action/frequency integrator supplies
Omega_r/Omega_t. Direct integration avoids a regression case in which the
interpolated Hamiltonian derivative becomes negative on a very eccentric,
weakly bound orbit. The analytic radial limit is used at L=0, and the circular
epicycle limit at Jr/L<1e-5 avoids singular period integrals when the turning
points coincide. The latter has an error of order Jr/L. This local numerical
limit is distinct from the optional global `frequency_mode="epicycle"`
approximation. That approximation is retained for comparison, not adopted as
an exactly isotropic model.

Every self-gravity iteration rebuilds the contour map and its normalization
in the updated total potential. The final density is checked between radial
grid points with twice the velocity quadrature. Convergence includes stellar
density as well as total force, so a dominant halo cannot conceal an
unconverged tracer. Relative density closure is measured where density exceeds
1e-12 of its peak or stellar mass per log-radius exceeds 1e-8 of Mstar.
Only the negligible numerical tail below 1e-30 of the peak is continued for
logarithmic density/moment interpolation. Its contribution is included in the
independent stellar mass-closure check; the DF itself is not floored.

The implementation requires alpha>1/2. In a harmonic core at r=0,
Lc is proportional to v^2, so f(v)-f(0) is proportional to |v|^(2 alpha).
This bound makes its first velocity derivative tend to zero at the origin.
Finite normalization alone would require only alpha>0. Central black holes
are not supported by this finite-central-potential prescription.

## Prescribed dark densities

The halo uses the exact coordinate rho20=rho_DM(20 pc), a scale radius, and
either gamma=0 or gamma=1 as separate families. Its outer taper is
exp[-(r/r_t)^2], normalized at 20 pc along with the inner profile. The taper
radius is fixed in the supplied examples. This is deliberately a new density
implementation: the older Jeans/DF templates use a different taper, and are
unchanged. No-DM is a separate model with rho20=0. The optional remnant density
is Plummer with its own mass and radius.

The dark components enter Poisson gravity without dark orbital parameters.
This constructs a positive equilibrium **stellar** DF. It does not establish
phase-space admissibility of the prescribed halo/remnant densities, dynamical
stability or survival in the Galactic tidal field.

## Independent lowered-isothermal control

The control implements

\[
f(E,L)=A\exp[-L^2/(2r_a^2s^2)]
 E_\gamma(g,(E_t-E)/s^2),\qquad E<E_t,
\]

and zero otherwise. This is the generalized Michie/LIMEPY family, implemented
with direct velocity quadrature and an embedded spherical Poisson solve.
It does not call the external LIMEPY package, and it is not the Osipkov-Merritt
energy-shift mock in `df_mock.py`.

The five stellar controls are central dimensionless depth W0, velocity scale,
central stellar density, anisotropy radius and truncation index. Stellar mass,
DF amplitude and energy cutoff follow from the solution. Dark gravity enters
the Poisson integration itself. Tests compare the isotropic King limit to the
independent existing mock solver and check Jeans balance in embedded cored
and cusped halos. The supplied examples have finite anisotropy radii; a large
fixed radius supplies an explicitly isotropic control.

## Running and reproducing the checks

All commands below use the project virtual environment and one numerical
thread. They create fresh result directories; they do not launch posterior
sampling.

```sh
PYTHONPATH=src OMP_NUM_THREADS=1 MPLCONFIGDIR=/private/tmp/ocen-df-mpl \
  /Users/vasilybelokurov/Work/venvs/.venv/bin/python -m pytest \
  tests/test_regularized_df.py tests/test_lowered_isothermal.py \
  tests/test_positive_df.py tests/test_df_capacity.py tests/test_df_mass_recovery.py -q

/Users/vasilybelokurov/Work/venvs/.venv/bin/python \
  bin/validate_regularized_df.py --out results/diagnostics/CHOOSE_A_FRESH_NAME
```

The validation driver records source copies and hashes, analytical checks,
constant-DF contour accuracy, self-consistency, independent projections and
Jeans balance, different starting potentials, doubled numerical grids and
the numerical shift in the actual saved observational bins. Its JSON report
and iteration log live in the chosen result directory. PNGs go to `plots/`;
their arrays and provenance go to `results/plot_data/`.

Six unfitted starting specifications are in `configs/df/`: the
`regularized_*` and `lowered_isothermal_*` files, each with no-DM, cored and
cusped variants. Their fit bounds are explicit **optimizer bounds**, not
adopted probability priors. The existing `bin/run_df_pilot.py` accepts either
family through its usual `--config` argument. Light normalization remains a
profiled photometric zero point in this diagnostic likelihood. The driver
preserves data snapshots and stores the new source files with each run.

## Scope of the present validation

These checks establish numerical implementation and representative
equilibria. They do not establish that the compact family fits Omega Cen or
can recover dark density without bias. Matched and independent mock recovery,
prior and outer-taper sensitivity, tracer-selection/equipartition checks,
and subsequent rotation/flattening tests remain separate research stages.

The next stage is the [matched-family compact DF recovery experiment](COMPACT_DF_RECOVERY.md),
with noiseless no-DM and injected-halo mocks, two starting points, and separate
checks of numerical accuracy, observable residuals, and physical mass recovery.
