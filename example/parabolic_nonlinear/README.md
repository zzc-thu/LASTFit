# Parabolic leading edge: nonlinear acoustic forcing

This case computes the nonlinear time-domain response of the Mach 15
parabolic-leading-edge flow to a planar fast-acoustic disturbance.

## Included configuration

- Grid declared by `Config.cfg`: 81 x 51 x 8
- Reynolds number: 6026.6
- Acoustic amplitude: 5e-4
- Streamwise wave number: 15
- Analysis: nonlinear unsteady Navier-Stokes

## Release status

The configuration, one-piece parabolic grid, and matching one-rank flow/shock
restart are included. This is not yet a manuscript reproduction case: the
manuscript currently reports a 181 x 161 x 8 grid, while available legacy
nonlinear perturbation output uses another grid. Do not combine those files
until their provenance has been established.

`output/` contains a representative final unsteady field with the same
81 x 51 x 8 extent as the included configuration. The completed manuscript
reproduction case must additionally include the wall-normal time histories,
Fourier extraction script, and sampling-window metadata.
