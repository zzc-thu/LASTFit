# Base-Flow Restart

[Case guide](../README.md)

The `SF_Results*.flowsfg` and `SF_Results*.shksfg` files form matching
flow/shock pairs for the 1 x 1 MPI decomposition and 81 x 51 x 8 grid.

Copy only these numerical files, not this README, into the run directory's
`RESU_STEADY/`. The LNS solver requires this nonlinear base flow even when
perturbation continuation is disabled. Perturbation restarts are separate
files under `LNSResults/`.

Restart records are compiler-dependent Fortran unformatted data. Retain the
matching grid and MPI partition; cross-toolchain portability is not established.
