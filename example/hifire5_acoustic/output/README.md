# Reference Output

[Case guide](../README.md)

| Dataset | PVTS header | Contents |
|---|---|---|
| Compact steady field | [Result_00000000.pvts](steady/Result_00000000.pvts) | One VTS piece |
| Production wall harmonic | [Harmonic_f0_666p67_Wall_3D.pvts](manuscript_fine/wall_harmonic/Harmonic_f0_666p67_Wall_3D.pvts) | All 192 referenced VTS pieces |

The wall-harmonic dataset stores pressure-harmonic amplitude, phase, real part,
and imaginary part at `f0=666.666666667`. Open each PVTS header with its
referenced pieces at their relative paths.

Full production-grid base flows, perturbation histories, and three-dimensional
harmonic volumes are not distributed in Git. The wall field supports inspection
of amplitude and phase maps, not reconstruction of the full time series.
