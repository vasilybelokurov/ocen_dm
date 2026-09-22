# Frequency-mapper failure found during compact DF recovery

The two refined mock truths passed their resolution checks, but a two-evaluation
fitting smoke test failed when rebuilding its selected model at finer resolution.
The original error was `integrateGL: order is too high (not implemented)` in
AGAMA's direct action/frequency integration.

The frequency integrator received nonfinite velocities from the preceding
action mapper. The requested spherical actions were Jr=9.911050890959185e-6
and L=2.4039457821260765e-5 pc km/s. They were represented as Jz=L, Jphi=0,
with all three angles set to 1 radian. The mapper returned x=y=0,
z=0.001495402583200513 pc, and velocities (-Infinity, -Infinity, NaN).

The installed AGAMA source `actions_spherical.cpp`, in
`mapPointFromActionAngles`, computes a polar-plane velocity component using
division by cylindrical radius R. The pole R=0 is a coordinate singularity.
The source and captured nonfinite coordinates identify this as the immediate
cause of the downstream integration failure.

Representing the same Jr and L as Jz=0, Jphi=L puts the orbit in the equatorial
plane of the spherical potential. A replay then returned finite coordinates,
recovered the input actions, and gave Omega_r/Omega_t approximately 1.999999.
Changing the orbital plane preserves the spherical invariants and frequencies.
The frequency calculation now uses that equatorial representative and rejects
nonfinite mapped coordinates explicitly before calling the frequency integrator.

The angular-momentum regression failed before the fix with relative error
0.00042647 and passes afterwards at a tolerance of 1e-6. All 64 focused tests
pass, and the complete refined-model reproduction that originally raised the
exception now exits successfully. The original v1 mocks and failed smoke test
are retained; a fresh v2 batch will use the corrected source.

## Preserved evidence

- [Original probe](../bin/diagnostics/compact_frequency_probe_20260923.py):
  the full-model failure-isolation script, initially placed in `/tmp` and now
  retained in version control. Its input is now the tracked
  [smoke-test checkpoint](../tests/fixtures/compact_frequency/recovery_checkpoint.json).
- [Standalone replay](../bin/diagnose_compact_frequency.py): compares orbital
  planes using the tracked [potential](../tests/fixtures/compact_frequency/potential.ini)
  and [action/angle record](../tests/fixtures/compact_frequency/orbit.json).
- Original experiment outputs remain in
  `results/diagnostics/compact_frequency_failure_20260923/` and
  `results/diagnostics/compact_recovery_smoke_20260923_v1/`.

The exported potential has finite decimal precision. Exporting and reloading
changes the old polar output to a finite point close to the axis; the standalone
replay does not reproduce the original nonfinite values exactly. The action
record preserves this distinction. The full-model probe supplies the original
reproduction path with the original source revision `29b4694`.
