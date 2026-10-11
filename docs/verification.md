# Verification Scope

[Overview](../README.md) | [Usage guide](usage.md) | [Examples](../example/README.md)

## Inputs and Build

Six configurations and 44 flow/shock input files pass sequential-record,
record-shape, interior finite-value and density/pressure positivity checks.
The 79 input-file checksums and reference-grid PVTS links have been verified.
Some supplied nonlinear checkpoints contain non-finite physical exterior halos;
these are distinguished from owned cells by the input checker.

All 37 compilation units were built from an empty directory using Intel
Fortran Classic 2021.11 and Intel MPI 2021.11 on Windows. The Linux Makefile
has not been executed in this test environment.

## Runtime Tests

Tests use the supplied grids, MPI partitions and physical settings with
shortened iteration limits and output intervals.

| Configuration | MPI ranks | Completed steps | Result |
|---|---:|---:|---|
| Circular cylinder | 8 | 10 | Normal exit, finite fields, valid VTK references |
| Parabolic nonlinear | 1 | 10 | Normal exit, finite fields, valid VTK references |
| Parabolic LNS | 1 | 100 | Normal exit, finite fields and all four history levels |
| Blunt cone | 8 | 10 | Normal exit, finite fields, valid VTK references |
| HIFiRE acoustic | 4 | 10 | Normal exit, finite fields, valid VTK references |
| HIFiRE compact steady | 1 | 10 | Normal exit, finite fields, valid VTK references |

Nonlinear density, pressure and temperature remain positive in the checked
output fields. Unsteady checkpoints contain matching physical-time metadata
on every rank. These tests establish short-run execution, not convergence.

## LNS Consistency Tests

Tests use the supplied 81 x 51 x 8 parabolic base flow.

| Check | Observed result |
|---|---|
| Stored history | All four flow/shock history levels are finite |
| Zero forcing | Exactly zero interior flow and shock perturbations |
| Amplitude scaling | Maximum absolute difference between the full response and twice the half-amplitude response: 2.43e-15 |
| Restart | Ten steps continued to twenty match the uninterrupted twenty-step history exactly |
| MPI partition | Two-rank partitions in either direction match the single-rank owned histories exactly |
| Forcing modes | Fast acoustic, slow acoustic, entropy and vortical modes complete ten-step LNS and nonlinear runs with finite fields |
| Diagnostic build | Owned-cell LNS passes bounds checks and floating-point traps enabled after base-flow loading |
| Restart guards | Missing metadata, incompatible time step/grid and non-finite histories are rejected |

The optional result checker requires Python and NumPy:

```bash
python tools/check_lns_results.py run/parabolic_lns/LNSResults
python tools/check_lns_results.py run/resumed/LNSResults --reference run/continuous/LNSResults
```

Six synthetic tests cover the input checker; five cover LNS record checks
and comparisons. These test the utilities, not the physical accuracy of the solver.

## Validation Limits

Grid/time convergence, complete AD Taylor consistency, long-time stability
and production-grid manuscript reproduction remain unverified. The compact
parabolic grid differs from the 181 x 161 x 8 manuscript comparison grid.

The MPI consistency tests use spanwise-homogeneous base flow and forcing;
they do not validate nonzero azimuthal-wavenumber three-dimensional LNS physics.
Diagnostic floating-point traps apply to the LNS phase, not all initialization
paths. Short-run consistency is not an AD Taylor test. HIFiRE demonstrates
nonlinear three-dimensional acoustic response, not controlled three-dimensional
LNS validation.

No Linux runtime, cross-compiler checkpoint portability, converged benchmark
accuracy or MPI scaling guarantee is inferred from these tests.
Record the exact source commit, compiler, MPI version and run configuration
when reporting numerical results.
