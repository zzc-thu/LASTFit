# Changelog

[Overview](../README.md) | [Verification status](verification.md)

## Unreleased: 0.1.0 Preparation

### Software and Examples

- Initialize LNS workspaces and evaluate boundary operators on owned cells
  before halo exchange; remove undefined exterior-geometry arithmetic.
- Keep linear forcing on the fixed base-flow shock with unforced primal
  inputs, and supply consistent derivatives for all four forcing modes.
- Correct the LNS interpolation fallback and the singular-grid RK4 updates.
- Allocate implicit matrices and Krylov storage only for the solvers that
  use them; remove the fixed per-rank buffered-output allocation.
- Correct planar MPI gathering and missing Fortran module prerequisites.
- Add finite-history, restart, partition and amplitude-scaling checks with
  synthetic regression tests; distinguish these from accuracy validation.
- Recovered the complete compact nonlinear checkpoint without altering its
  original flow or shock values; supplied the missing coordinate records.
- Corrected the compact HIFiRE MPI partition and base-flow staging commands.
- Added nonlinear physical-time and LNS time-step/grid/partition metadata.
- Reject non-finite interior LNS checkpoints rather than silently saving them.
- Corrected VTK coordinate block lengths and ordinary serial/single-partition
  output extents; removed duplicate finalization deallocations.
- Added example-input checks and isolated short-run staging tools.
- Moved the incompatible E4 configuration into historical documentation.

- MPI-parallel, perfect-gas shock-fitting source with steady, nonlinear unsteady,
  and AD-generated LNS modes.
- Five examples containing configurations, reference grids, restart data, and
  representative outputs.
- GNU Make configuration with selectable compiler and external BLAS/LAPACK.
- LASTFit startup banner, citation metadata, and existing license file.

### Documentation

- Recorded passing original-grid short runs for all six active configurations,
  including ten-step blunt-cone and HIFiRE acoustic runs after extending the
  bounded test under severe memory pressure; retained accuracy-validation limits.
- Concise repository overview with a linked example index and source map.
- Complete Bash staging commands for the primary configurations and the
  HIFiRE compact steady setup.
- Explicit distinction between reference grids, runtime paths, and
  nonlinear/LNS restart requirements.
- Documented production-grid exclusions, historical configuration compatibility,
  and verification limits.

### Scope

This entry records preparation of a research snapshot, not a tagged release.
Part II chemistry and multigrid extensions are outside the current scope.
Full manuscript time histories and numerical/performance verification records
are not included merely by documenting their requirements.
