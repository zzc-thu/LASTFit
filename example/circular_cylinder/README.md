# Circular cylinder

This case evaluates the steady shock-fitting solution for Mach 5.73 viscous
flow over a circular cylinder. The three-dimensional solver is run with a
spanwise-uniform, periodic solution to exercise the three-dimensional storage,
metric, boundary-condition, and output paths.

## Included configuration

- Grid: 101 x 61 x 20
- Reynolds number: 2050
- Wall condition: isothermal
- Analysis: steady nonlinear Navier-Stokes

The grid dimensions match the case documented in the manuscript.

## Release status

`Config.cfg` has been checked against the available calculation. The distributed
grid, converged flow/shock restart, and one final ParaView snapshot still need
to be copied from that same run and assigned checksums. Do not use perturbation
or LNS files found beside the legacy calculation: representative files in
those directories belong to a parabolic-leading-edge grid.

The reduced release output should contain the wall-pressure samples and the
temperature field used for the validation figure.
