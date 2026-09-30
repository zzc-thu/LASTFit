# HIFiRE-5-type elliptic cone

This example is the non-axisymmetric three-dimensional application of the
steady and nonlinear acoustic-forcing workflow.

## Included configurations

- `Config.cfg`: unsteady acoustic setup, 90 x 81 x 160
- `Config_steady.cfg`: steady setup, 90 x 41 x 40
- `Config_E4_fine.cfg`: archived E-4 manuscript run, 160 x 151 x 160
- Mach number: 6
- Unit Reynolds number: 10.2e6
- Cone half-angle: 7 degrees

## Release status

The 90 x 81 x 160 unsteady grid and its four-rank base-flow restart are
included. `Config.cfg` uses four circumferential MPI partitions to match those
files. The one-piece 90 x 41 x 40 steady grid is also included. These are
legacy setup records and must not be presented as the 160 x 151 x 160
manuscript calculation.

`Config_E4_fine.cfg` is the configuration recovered with the manuscript E-4
dataset. It records `epsilon=5e-4`, `k_infty=3700`, `dt=1e-7`, an output
interval of 2000 steps, and a 12 x 16 MPI layout. It uses the model numbering
of the historical executable (`ModelType=2`), which differs from the current
release configuration; retain it as provenance unless compatibility with a
new executable has been tested.

`output/manuscript_fine/` now contains the complete wall first-harmonic field,
time-domain diagnostics, late-periodic amplitude and phase products, and the
E-2/E-3/E-4 amplitude and spectral comparisons used to check the frequency
normalization. The corresponding scripts are in `postprocess/` and no longer
contain machine-specific absolute paths.

The full fine-grid base-flow field, one complete E-4 perturbation snapshot,
and full three-dimensional harmonic field were verified in the raw archive,
but are approximately 281 MB, 360 MB, and 320 MB, respectively. They should be
distributed through the versioned data archive/DOI rather than duplicated in
the ordinary Git tree. The fine-grid initial `INIT` dataset has not yet been
located and remains the main missing input artifact.
