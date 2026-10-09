# JH2 Validation Log

## Baseline environment

- Peridigm source commit: `f25c4928e850b2aa0012bc5d3b8b1ec5c78c8903`
- Runtime image digest:
  `sha256:58b5acc40a4f2bc97dbffa783ed650953e6e65356cb0bf6706a3338eccc1a9b5`
- Runtime-reported Peridigm version: `1.5.0 (Dev)`
- Runtime compiler: GCC 12.2.0
- Runtime MPI: Open MPI 4.1.4
- Host: Ubuntu 26.04 LTS under WSL2

## Unmodified runtime baseline

The `Compression_2x1x1` verification input completed five explicit time steps
with exit code zero using the pinned runtime image.  The stable time step was
`5.78881e-4 s`, and the requested time step was `1.0e-5 s`.

The full tensile example was not selected as the iterative baseline because it
assembles a global tangent with approximately 30 million nonzeros.  It will be
used only for milestone validation when appropriate.

## First JH2 implementation milestone (2026-09-25)

### Build

- Development image: `localhost/peridigm-jh2-dev:f25c492`.
- Base runtime is pinned by digest; the build stage uses Debian Bookworm,
  GCC 12.2.0, CMake 3.25.1, Open MPI 4.1.4, and the runtime's Trilinos 13.1.
- A compatibility link for `libdl.so` is required because Trilinos records the
  pre-glibc-2.34 linker path while Bookworm folds libdl into libc.
- The complete Peridigm source and all test executables build successfully with
  C++17.
- The standalone JH2 kernel and material-point tests also compile with
  `-Wall -Wextra -Wpedantic -Werror`.

### Material-point and wrapper tests

`utJH2MaterialPoint` verifies thirteen paths:

1. analytical elastic response;
2. zero-increment state preservation;
3. hydrostatic compression;
4. plastic return, damage, energy loss, and bulking pressure against
   hard-coded reference values;
5. pure-shear radial return;
6. strain-rate strengthening;
7. monotonic multi-step history variables;
8. persistent hydrostatic-tension failure;
9. non-positive-J rejection;
10. invalid-history rejection;
11. full-tensor shear-rate contraction;
12. parameter validation;
13. independent single-increment comparison of the Eq. (27)--(28) energy
    difference and trapezoidal plastic work.

`utPeridigm_JH2CorrespondenceMaterial` verifies the DataManager field wiring,
analytical elastic update, preloaded rigid rotation, and four-point affine
compression.  Both CTest targets pass.  The pre-existing
`utPeridigm_ElasticMaterial` and `Compression_2x1x1_np1` tests also pass.

### End-to-end prescribed compression

The 384-point `jh2_uniaxial_compression` example completed ten explicit steps
with `dt = 1.0e-7 s`; the reported stable step was `7.29403e-5 s`.  The Exodus
database contains 11 states, 9 nodal variables, and 39 element variables.  At
the final state the imposed deformation has `J = 0.99`; all inspected JH2
history and diagnostic values are finite, and the hydrostatic-tension flag is
zero.

The same example completed with two MPI processes after partitioning the mesh.
The two files were merged with SEACAS EPU and compared with the single-process
database using Exodiff.  With a `1.0e-11` relative tolerance and a `1.0e-4`
floor for near-zero roundoff, Exodiff reported that the files are the same.

### Scientific limitations still open

- The production path uses Zhu and Zhao's literal energy-difference form for
  Eqs. (27)--(28).  The alternative trapezoidal plastic-work definition has not
  yet been implemented as a selectable option or compared over multi-step
  paths.
- `JH2_Damage = 1` degrades strength to the fractured surface but does not
  delete Peridigm bonds.  Hydrostatic tensile failure zeros point stress and is
  persistent, but its coupling to bond failure remains a separate milestone.
- Plastic-regime objectivity under combined stretch and rotation, formal
  time-step convergence, published benchmark curves, experimental calibration,
  and stabilization sensitivity are not yet complete.  Therefore this
  milestone establishes software integration and deterministic behavior, not
  predictive validation for blasting simulations.

## Objectivity and energy-definition diagnostics (2026-10-09)

### Preloaded rigid rotation

The Peridigm wrapper test now initializes a one-point state with the diagonal
unrotated Cauchy stress `(10, -5, 2) MPa`, nonzero damage, plastic strain, and
bulking pressure.  A zero rate-of-deformation and a 30-degree rotation about
the z axis are then applied.  The JH2 update preserves the unrotated stress and
all history variables.  Peridigm's correspondence rotation routine produces
`R sigma R^T` to a relative tolerance of `1.0e-13`.  This verifies the intended
single-rotation call chain independently of the multi-step test below.

### Four-point affine compression

A tetrahedral four-point neighborhood, with three independent bonds per point,
is subjected to the homogeneous fields
`F = diag(0.9995, 1, 1)` and `L = diag(-500, 0, 0) 1/s`.  Coordinates and
velocities are assigned as `y = F X` and `v = L y`; the test then calls the
complete `CorrespondenceMaterial::computeForce()` path rather than supplying
the deformation gradient or rate-of-deformation directly.

For all four points, the shape-tensor calculation recovers `F`, the velocity
gradient recovers the prescribed unrotated rate, and JH2 returns the analytical
elastic stress.  The deformation Jacobian is `0.9995`, damage remains zero,
the assembled internal force is self-equilibrated, and the hourglass force
vanishes within `1.0e-12` of the characteristic force scale.  The absolute
tolerances for theoretical zero components are `1.0e-14` for `F` and rotation,
`1.0e-10 1/s` for rate, and `1.0e-6 Pa` for stress.  This establishes affine
reproduction for one homogeneous increment.

### Multi-step stretch and rotation

The same four-point neighborhood is advanced through the homogeneous motion

`F(t) = Rz(omega*t) * diag(exp(-100*t), 1, 1)`,

with `omega = 1.0e4 rad/s`.  At every increment, coordinates use `y = F X` and
velocities use `v = (Fdot F^-1) y`.  The test performs 100 state updates with
`dt = 1.0e-7 s`, giving a final rotation of `0.1 rad` and an axial compression
of approximately `0.1%`.

At every step and all four points, the computed deformation gradient and
Jacobian reproduce the prescribed motion, damage remains zero, internal force
is self-equilibrated, and affine hourglass force remains within `1.0e-12` of
the characteristic force scale.  The computed rotation, unrotated rate,
unrotated elastic stress, and spatially rotated stress agree with the continuum
reference within relative tolerances of `2.0e-4`, `3.0e-4`, `3.0e-4`, and
`5.0e-4`, respectively.

An initial 10-step trial used `dt = 1.0e-6 s` while preserving the same total
motion.  Its final errors were approximately `0.1%` in rotation, `0.11%` in
unrotated stress, and `0.31%` in spatial stress, exceeding the intended test
tolerances.  Refinement to 100 steps reduced the errors sufficiently for all
checks to pass.  This behavior is consistent with time-discretization error in
Peridigm's Flanagan--Taylor rotation update; it is not evidence of a second
stress rotation in JH2.  A formal convergence-rate study and a plastic-regime
counterpart remain open.

### Dissipated-energy definitions

For the established one-increment plastic reference point, the production
Eq. (27)--(28) energy-difference definition gives
`1.9444855600160017e5 J/m^3`.  The independently evaluated trapezoidal plastic
work,

`0.5 * (trial equivalent stress + yield strength) * plastic strain increment`,

gives `8.7692853506050073e5 J/m^3`, or 4.509822819 times the production value.
The test fixes both reference values and confirms that they are not silently
interchanged.  No production parameter or constitutive path was changed;
selecting the alternative remains contingent on source justification and a
multi-step sensitivity study.

Both JH2 CTest targets pass after these additions.  The wrapper target now
contains the analytical elastic update, preloaded rigid rotation, and
four-point single- and multi-increment affine-motion tests.
