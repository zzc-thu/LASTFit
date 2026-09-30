# Changelog

All notable changes to LASTFit Part I are recorded in this file.

## [0.1.0] - 2026-09-30

### Added

- MPI-parallel perfect-gas shock-fitting solver source.
- Nonlinear steady and unsteady Navier-Stokes modes.
- AD-generated linearized Navier-Stokes mode.
- Five documented benchmark and demonstration cases with representative data.
- Portable GNU Make configuration with external BLAS/LAPACK linkage.
- Release metadata, checksums, citation information, and BSD-3-Clause license.

### Known limitations

- The repository contains compact release examples, not every production-grid
  time history used in the manuscript.
- Reacting five-species and eleven-species models and multigrid extensions are
  planned for Part II and are not included in this release.
