# JH2 fully prescribed compression

This small verification example applies one percent uniaxial compressive
strain over `1.0e-6 s` to the existing fully prescribed correspondence mesh.
All displacement degrees of freedom are prescribed, so the example isolates
the constitutive update from global equilibrium and rigid-body modes.

The imposed condition is uniaxial **strain**, not uniaxial stress. It is a
high-rate material verification case and is not intended to represent a
laboratory compression test.

Run the input from this directory with a Peridigm executable containing the
`JH2 Correspondence` material:

```bash
Peridigm jh2_uniaxial_compression.yaml
```

The mesh is reused from
`test/verification/ElasticCorrespondenceFullyPrescribedTension` to avoid
duplicating a binary Exodus file.
