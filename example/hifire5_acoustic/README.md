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

The 90 x 81 x 160 unsteady grid and its four-rank base-flow restart are
included. `Config.cfg` uses four circumferential MPI partitions to match those
files. The one-piece 90 x 41 x 40 steady grid is also included. These remain
provisional setup records rather than the final manuscript reproduction
package, which reports a 160 x 151 x 160 grid. The raw E-2, E-3, and E-4
acoustic histories and their harmonic-analysis products must be restored from
the original data disk before this example can be marked complete.

`output/steady/` contains the available steady VTK field. No matching nonlinear
acoustic VTK field is present on the currently mounted data disks. The final
package still requires selected wall histories, first-harmonic amplitude and
phase fields, spectra, and the script that defines the Fourier normalization
and sampling window.
