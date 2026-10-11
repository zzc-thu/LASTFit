# Verification Status

[Overview](../README.md) | [Build and execution](usage.md) | [Examples](../example/README.md)

## Inputs and Build

Six active configurations and 44 flow/shock files pass sequential-record,
record-shape, interior finite-value and density/pressure positivity checks.
Example SHA-256 manifests and supplied PVTS references have been checked.
The input checker has six synthetic regression tests. These checks establish
input integrity, not numerical accuracy.

The recovered parabolic nonlinear checkpoint includes all four required
records. The HIFiRE acoustic checkpoint uses its matching `2 x 2` partition.
The incompatible historical E4 configuration is isolated under `docs/archive/`.
Some original nonlinear checkpoints have non-finite physical exterior halos;
these are distinguished from owned cells. Original numerical data have not
been replaced with manufactured values.

All 37 compilation units were built from an empty directory using Intel
Fortran Classic 2021.11 and Intel MPI 2021.11 on Windows. A separate harness
reads the Makefile prerequisites; missing Fortran module dependencies have
been corrected. The Linux Makefile itself has not been executed here.

Unused implicit flux matrices are no longer allocated for RK4 or LNS;
the Krylov basis is allocated only for GMRES. The CFL diagonal remains
allocated for every solver. Ordered blocking sends replace buffered planar
output, removing a fixed 800 MB allocation per MPI rank. This is an allocation
change, not a measured parallel-performance result.

## Runtime Checks

The final test binary was used with the distributed settings except for
shortened iteration limits and output intervals:

| Configuration | MPI ranks | Completed steps | Output check |
|---|---:|---:|---|
| Circular cylinder | 8 | 10 | Finite fields and valid VTK blocks/references |
| Parabolic nonlinear | 1 | 10 | Finite fields and valid VTK blocks/references |
| Parabolic LNS | 1 | 100 | Finite fields and all four stored history levels |
| Blunt cone, original grid | 8 | 10 | Finite fields and valid VTK blocks/references |
| HIFiRE acoustic, original grid | 4 | 10 | Finite fields and valid VTK blocks/references |
| HIFiRE compact steady | 1 | 10 | Finite fields and valid VTK blocks/references |

Nonlinear output density, pressure and temperature were checked for positive
values. A reduced 60 x 21 x 40 blunt-cone test completed ten steps. It is a
separate startup test, not a grid-convergence study.

All six active configurations now pass the short-run checks above, including
normal process exit, completed solver steps and finite final fields. Nonlinear
unsteady checkpoints include matching physical-time metadata on every rank;
LNS checks include all four saved history levels. These are runtime checks,
not converged benchmark comparisons.

Earlier 180-second tests were incomplete under severe memory pressure.
Repeating the original-grid blunt-cone and HIFiRE acoustic tests sequentially
with a 30-minute limit completed ten steps in 662.36 s and 1127.41 s,
respectively. The grid, MPI partition and physical parameters were unchanged.
These elapsed times describe a heavily loaded 16 GB Windows workstation;
they are not performance benchmarks or production-runtime guarantees.

## LNS Repairs and Regression

An initial failure was traced to uninitialized transformed viscous mass
fluxes. Further diagnostics found incomplete workspace initialization,
boundary formulas evaluated on undefined exterior metric storage, forcing
evaluated on the perturbed rather than fixed base-flow shock, and forcing
included in both primal and tangent inputs. These paths have been corrected.
Boundary formulas evaluate owned cells before halo exchange; stencil halo
communication is retained.

All four forcing modes now have consistent state and time-derivative inputs.
The interpolation search has a defined upper-end fallback. Nonlinear RK4
uses the same mode definitions and the appropriate singular-grid update.

Checks on the supplied 81 x 51 x 8 parabolic base flow include:

| Check | Observed result |
|---|---|
| LNS stored history | All four levels, including exterior storage, finite |
| Zero forcing | Exactly zero interior flow and shock perturbations |
| Amplitude scaling | Full response minus twice half-amplitude response: maximum absolute difference 2.43e-15 |
| Restart | Ten steps continued to twenty match uninterrupted twenty-step history exactly |
| MPI partition | `2 x 1` and `1 x 2` match single-rank owned histories exactly |
| Forcing modes | Fast acoustic, slow acoustic, entropy and vortical modes complete ten-step LNS and nonlinear runs with finite fields |
| Diagnostic build | Owned-cell LNS passes bounds checks and floating-point traps enabled after base-flow loading |
| Restart guards | Missing metadata, wrong time step, wrong grid and non-finite history are rejected |

The previous ghost-storage and strict-restart failures are superseded by
these repairs, not passing results for unmodified earlier commits.
Record the source commit with every run.

`tools/check_lns_results.py` requires NumPy, audits every stored history record
and compares restart, partition and amplitude-scaling results. Five synthetic
regression tests cover its record checks and comparisons. For example:

```bash
python tools/check_lns_results.py run/parabolic_lns/LNSResults
python tools/check_lns_results.py run/resumed/LNSResults --reference run/continuous/LNSResults
python tools/test_check_lns_results.py
```

## Remaining Validation

The preparation tool retains `--allow-unverified` for LNS because grid/time
convergence, complete AD Taylor consistency, long-time stability and
production-grid manuscript reproduction remain unverified. The compact
parabolic example is not the 181 x 161 x 8 production comparison dataset.

The MPI checks use spanwise-homogeneous base flow and forcing; they do not
validate nonzero azimuthal-wavenumber three-dimensional LNS physics.
Diagnostic floating-point traps apply to the LNS phase, not all initialization
paths. Short-run consistency is not an AD Taylor test. The manuscript HIFiRE
result demonstrates nonlinear three-dimensional acoustic response, not a
controlled three-dimensional LNS validation.

No Linux runtime, cross-compiler checkpoint portability, converged benchmark
accuracy or MPI scaling guarantee is inferred.

## Release Status

Version 0.1.0 denotes the source version in the citation metadata. No formal
tag, GitHub Release or archival DOI is published by this update. The existing
BSD-3-Clause license is retained. Record the exact commit when using the code.
