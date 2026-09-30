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

The configuration is associated with a complete local LNS record containing 71
parallel VTK time records and flow/shock restart data. Those binary files are
not committed here yet. The manuscript currently states a 181 x 161 x 8 grid,
so the release data and manuscript description must be reconciled.

A complete reproduction package must also include the matching nonlinear
time histories and a relative-path script that regenerates the LNS/nonlinear
amplitude and phase comparisons.
