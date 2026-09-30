# HIFiRE-5-type elliptic cone

This example is the non-axisymmetric three-dimensional application of the
steady and nonlinear acoustic-forcing workflow.

## Included configurations

- `Config.cfg`: unsteady acoustic setup, 90 x 81 x 160
- `Config_steady.cfg`: steady setup, 90 x 41 x 40
- Mach number: 6
- Unit Reynolds number: 10.2e6
- Cone half-angle: 7 degrees

## Release status

These are provisional setup records rather than the final manuscript
reproduction package. The manuscript reports a 160 x 151 x 160 grid. The raw
E-2, E-3, and E-4 acoustic histories and their harmonic-analysis products must
be restored from the original data disk before this example can be marked
complete.

The final package should contain a matching steady restart, selected wall
histories, first-harmonic amplitude and phase fields, spectra, and the script
that defines the Fourier normalization and sampling window.
