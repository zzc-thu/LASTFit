# Reduced output

`steady/` contains `Result_00000000.pvts` and its referenced VTS piece. Its
stored VTK extent is 0:90, 1:41, 1:41, which differs from the initial-grid index
range in `Config_steady.cfg` and is preserved here exactly as generated.

`manuscript_fine/wall_harmonic/` contains a complete 192-piece PVTS wall field
at `f0 = 666.666666667`, including pressure-harmonic amplitude, phase, real,
and imaginary components. The smoothed single-wall VTS file is included for
publication rendering. All 192 files referenced by
`Harmonic_f0_666p67_Wall_3D.pvts` are present.

The other `manuscript_fine/` directories contain reduced CSV data, summaries,
and figures generated from the E-2 (`epsilon=5e-2`), E-3 (`epsilon=5e-3`), and
E-4 (`epsilon=5e-4`) nonlinear acoustic histories:

- `analysis_unsteady/`: E-4 time histories and wall RMS diagnostics;
- `analysis_late_periodic/`: late-window fundamental amplitude and phase;
- `amplitude_comparison/`: response normalized by the input amplitude and
  second-harmonic scaling;
- `modal_comparison/`: Fourier, POD, and DMD comparison products.

The archived sampling interval is `0.0002`, corresponding to a sampling
frequency of `5000`. Spectral amplitudes in the comparison products are
normalized by the imposed disturbance amplitude, not independently by the
maximum at each station.

The full three-dimensional base-flow, perturbation, and harmonic PVTS bundles
remain external archive products because each is several hundred megabytes.
