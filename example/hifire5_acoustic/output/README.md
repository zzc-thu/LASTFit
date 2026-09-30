# Reduced output

`steady/` contains `Result_00000000.pvts` and its referenced VTS piece. Its
stored VTK extent is 0:90, 1:41, 1:41, which differs from the initial-grid index
range in `Config_steady.cfg` and is preserved here exactly as generated.

`manuscript_fine/wall_harmonic/` contains a complete 192-piece PVTS wall field
at `f0 = 666.666666667`, including pressure-harmonic amplitude, phase, real,
and imaginary components. All 192 files referenced by
`Harmonic_f0_666p67_Wall_3D.pvts` are present. Only this compact field required
by the manuscript is retained; derived figures, CSV tables, smoothed display
files, POD/DMD products, and post-processing scripts are excluded.

The full three-dimensional base-flow, perturbation, and harmonic PVTS bundles
remain external archive products because each is several hundred megabytes.
