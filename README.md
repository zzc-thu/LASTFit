# LASTFit Part I

An MPI-parallel, three-dimensional shock-fitting solver for perfect-gas
hypersonic base flows and receptivity.

LASTFit combines nonlinear Navier-Stokes and automatic-differentiation-generated
linearized Navier-Stokes (AD-LNS) calculations within a shared moving-shock
framework, using high-order finite differences on structured grids.

## Build and Run

Requires GNU Make, an MPI Fortran compiler, and BLAS/LAPACK.

```bash
make -C src
```

The executable is `src/SFSolver`. See the [build and run guide](docs/usage.md#run)
for configuration, restart staging, and MPI process counts.

## Examples

| Case | Analysis |
|---|---|
| [Circular cylinder](example/circular_cylinder/) | Steady base flow |
| [Parabolic leading edge](example/parabolic_nonlinear/) | Nonlinear acoustic response |
| [Parabolic leading edge, LNS](example/parabolic_lns/) | AD-LNS response |
| [Blunt cone at incidence](example/blunt_cone_aoa1/) | Three-dimensional base flow |
| [HIFiRE-5 elliptic cone](example/hifire5_acoustic/) | Three-dimensional acoustic response |

Each case contains configuration, reference grids, restart files, and representative
outputs. See the [case index](example/README.md) for grid sizes and data scope.
These compact examples do not include all manuscript production histories.

## Documentation

- [Build, execution, and numerical scope](docs/usage.md)
- [Source modules](src/README.md)
- [Changelog](docs/CHANGELOG.md)
- [Citation metadata](CITATION.cff)

Part I covers a calorically perfect gas. Five-species and eleven-species models
and multigrid extensions are planned for Part II.

Distributed under the [BSD-3-Clause license](LICENSE).
Report reproducible problems through [GitHub Issues](https://github.com/zzc-thu/LASTFit/issues).
