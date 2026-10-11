# Restart Data

[Case guide](../README.md)

`unsteady/` contains four flow/shock file pairs for the 90 x 81 x 160
compact acoustic setup. Copy `*.flowsfg` and `*.shksfg` into the run's
`RESU/` directory and enable `IF_Continue_Calculate=1` for continuation.
Use the `2 x 2` MPI decomposition. These files match an archived steady
snapshot byte for byte; the historical `unsteady/` path is retained for
compatibility. The default acoustic run also stages them into `RESU_STEADY/`.
No convergence guarantee is inferred from that archive provenance.

The default configuration starts a new calculation. The separate 192-rank
production restart is not distributed here; the four-rank checkpoint cannot
replace it. Restart files use compiler-dependent Fortran unformatted records.
