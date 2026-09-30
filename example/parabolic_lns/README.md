# Parabolic leading edge: AD-LNS

This case exercises the automatic-differentiation-generated linearized
Navier-Stokes solver on the same moving-shock formulation used by the nonlinear
calculation.

## Included configuration

- Grid: 81 x 51 x 8
- Mach number: 15
- Reynolds number: 6026.6
- LNS time step: 6e-6
- Acoustic amplitude used for comparison: 5e-4

## Release status

The configuration, one-piece grid, and one-rank nonlinear base-flow restart are
included. A complete local LNS record containing 71 parallel VTK time records
has been identified, but the full time sequence is not committed here. The
manuscript currently states a 181 x 161 x 8 grid, so the release data and
manuscript description must be reconciled.

`output/` contains the final available LNS VTK record. A complete figure
reproduction package must also include the matching nonlinear time histories
and a relative-path script that regenerates the LNS/nonlinear amplitude and
phase comparisons.
