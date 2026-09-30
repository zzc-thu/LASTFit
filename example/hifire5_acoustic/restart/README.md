# Restart Data

[Case guide](../README.md)

`unsteady/` contains four flow/shock file pairs for the 90 x 81 x 160
compact acoustic setup. Copy `*.flowsfg` and `*.shksfg` into the run's
`RESU/` directory and enable `IF_Continue_Calculate=1` for continuation.

The default configuration starts a new calculation. The separate 192-rank
production restart is not distributed here; the four-rank checkpoint cannot
replace it. Restart files use compiler-dependent Fortran unformatted records.
