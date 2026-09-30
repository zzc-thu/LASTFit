# Restart Data

[Case guide](../README.md)

The `SF_Results*.flowsfg` and `SF_Results*.shksfg` files form matching
flow/shock pairs for the 4 x 2 MPI decomposition and 101 x 61 x 20 grid.

Copy only these numerical files, not this README, into the run directory's
`RESU/`. The supplied configuration enables continuation.

Restart records are compiler-dependent Fortran unformatted data. Retain the
matching grid and MPI partition; cross-toolchain portability is not established.
