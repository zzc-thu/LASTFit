# Changelog

[Overview](../README.md) | [Verification status](verification.md)

## Unreleased: 0.1.0 Preparation

### Software and Examples

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
