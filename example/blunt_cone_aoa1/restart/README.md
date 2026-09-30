# Restart Data

[Case guide](../README.md)

The `SF_Results*.flowsfg` and `SF_Results*.shksfg` files form matching
flow/shock pairs for the 4 x 2 MPI decomposition and 120 x 151 x 40 grid.

Copy only these numerical files, not this README, into the run directory's
`RESU/`. The default configuration starts a new calculation;
set `IF_Continue_Calculate=1` in a staged copy to read this checkpoint.

Restart records are compiler-dependent Fortran unformatted data. Retain the
matching grid and MPI partition; cross-toolchain portability is not established.
