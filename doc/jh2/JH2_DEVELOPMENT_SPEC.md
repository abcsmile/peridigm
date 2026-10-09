# JH2 Correspondence Material Development Specification

## Scope

The first milestone implements a JH2 constitutive material for Peridigm's
correspondence formulation.  It includes pressure-dependent strength,
strain-rate dependence, radial return, continuous damage, and bulking
pressure.  Bond failure, explosive loading, contact, and new stabilization
schemes are outside this milestone.

## References

1. Johnson and Holmquist (1994), JH-2 ceramic damage model.
2. Zhu and Zhao (2021), peridynamic modelling of blasting-induced rock
   fractures, Eqs. (19)--(36).

## Sign and kinematic conventions

- Cauchy stress is positive in tension.
- JH2 pressure is positive in compression:
  `p = -tr(sigma)/3`.
- `J = det(F)` and mass conservation gives `rho/rho0 = 1/J`.
- The compression measure is `mu = 1/J - 1`.
- The equivalent strain rate is
  `sqrt(2/3 * d_dev:d_dev)`.
- The Peridigm `CorrespondenceMaterial` base class supplies the unrotated
  rate-of-deformation and rotates the updated unrotated Cauchy stress back to
  the spatial frame.  The JH2 material must not rotate stress a second time.

## State variables

| Variable | Unit | Time storage | Meaning |
|---|---:|---|---|
| Unrotated Cauchy stress | Pa | TWO_STEP | Objective stress state managed by the base class |
| JH2 damage | 1 | TWO_STEP | Continuous strength-degradation variable in [0, 1] |
| Equivalent plastic strain | 1 | TWO_STEP | Accumulated effective plastic strain |
| Bulking pressure increment | Pa | TWO_STEP | Irreversible pressure increment from damaged plastic work |
| Hydrostatic tensile failure | 1 | TWO_STEP | Persistent failure flag set when `p* <= -T*` |

Diagnostic fields include equivalent stress, yield strength, EOS pressure,
equivalent strain rate, and deformation Jacobian.

## Peridigm material parameters

The input deck must provide every parameter explicitly.  No JH2 material
constant is assigned a silent default.

```yaml
Material Model: "JH2 Correspondence"
Density: 2660.0
Bulk Modulus: 21.9e9
Shear Modulus: 25.7e9
Hourglass Coefficient: 0.01
JH2 Intact Strength Constant: 1.9
JH2 Intact Strength Exponent: 0.62
JH2 Fractured Strength Constant: 0.25
JH2 Fractured Strength Exponent: 0.62
JH2 Strain Rate Coefficient: 0.005
JH2 Damage Constant D1: 0.005
JH2 Damage Exponent D2: 0.7
JH2 Energy Conversion Factor: 0.5
Pressure at HEL: 3.7e9
Strength at HEL: 1.2e9
JH2 Tensile Pressure: 54.0e6
JH2 Reference Strain Rate: 1.0
JH2 EOS K1: 21.9e9
JH2 EOS K2: -4500.0e9
JH2 EOS K3: 300000.0e9
```

`JH2 EOS K1` is required to equal `Bulk Modulus`, following Eq. (25).
The wrapper rejects non-finite parameters, invalid signs, and energy
conversion factors outside `[0, 1]`.

The hydrostatic JH2 tensile condition `p* <= -T*` sets the material-point
stress to zero and persistently marks `JH2_Tensile_Failure = 1`. It does not alter
Peridigm bond damage in this milestone.  Coupling this point failure flag to a
bond-failure model is a separate validation task.

## Update sequence

1. Validate `dt`, material parameters, state values, and `J`.
2. Compute strain increment from the unrotated rate-of-deformation.
3. Form the elastic trial stress, Eq. (19).
4. Compute trial deviatoric stress and equivalent stress, Eq. (20).
5. Compute `mu`, EOS pressure, and normalized pressure, Eqs. (25)--(26).
6. Evaluate intact and fractured normalized strengths, Eqs. (22)--(23).
7. Interpolate the current yield strength using damage, Eq. (21).
8. If elastic, accept the trial stress and retain history variables.
9. If plastic, perform radial return, Eq. (31).
10. Compute the effective plastic strain increment, Eqs. (32)--(34).
11. Update damage and fracture strain, Eqs. (29)--(30).
12. Update bulking pressure and reconstruct the unrotated stress,
    Eqs. (25)--(28), (35).

## Open modelling decision

Zhu and Zhao's Eqs. (27)--(28) define the energy loss used by the bulking
pressure update, while an alternative implementation uses trapezoidal plastic
work.  These definitions are not numerically identical.  They are compared as
independent material-point calculations: for the current single-increment
plastic reference, they give `1.94448556e5 J/m^3` and `8.76928535e5 J/m^3`,
respectively.  The production path retains the literal energy-difference
definition; making the alternative selectable requires a separate scientific
decision and multi-step sensitivity study.

## Acceptance tests

Completed in the first implementation milestone:

- zero increment and elastic analytical response;
- hydrostatic compression and tension;
- pure shear yield and radial return;
- strain-rate comparison;
- monotonic damage, plastic strain, and bulking pressure;
- invalid-parameter, invalid-history, and non-positive-J rejection;
- comparison with hard-coded independent material-point reference values;
- Peridigm DataManager wrapper update;
- preloaded rigid-rotation objectivity through the Peridigm stress-rotation
  call chain;
- single-increment comparison of the Eq. (27)--(28) energy loss and the
  trapezoidal plastic-work diagnostic;
- one- and two-process consistency for the prescribed compression example.

Required before using the model for calibrated blast predictions:

- an independent implementation of the complete multi-step loading path;
- an affine-deformation test at the Peridigm level;
- multi-step sensitivity comparison of the Eq. (27)--(28) energy loss and
  trapezoidal plastic-work alternatives;
- coupling of material-point failure to a documented bond-failure rule;
- comparison with published JH2 material tests and experimental rock data.
