# LASTFit Part I

```text
 _        _    ____ _____ _____ _ _
| |      / \  / ___|_   _|  ___(_) |_
| |     / _ \ \___ \ | | | |_  | | __|
| |___ / ___ \ ___) || | |  _| | | |_
|_____/_/   \_\____/ |_| |_|   |_|\__|

Version 0.1.0
LAST Group, Tsinghua University
Developed by Zhichao Zhu & Youcheng Xi
```

An MPI-parallel, three-dimensional shock-fitting solver for perfect-gas
hypersonic base flows and receptivity.

LASTFit combines nonlinear Navier-Stokes calculations and
automatic-differentiation-generated linearized Navier-Stokes (AD-LNS)
evolution in a shared moving-shock framework with high-order finite differences.

## Build and Run

Requires a Linux toolchain with GNU Make, MPI Fortran, and BLAS/LAPACK.

```bash
make -C src
```

The executable is `src/SFSolver`. The [usage guide](docs/usage.md) covers
compiler options, isolated run directories, and restart staging.

## Examples

| Case | Analysis | MPI ranks |
|---|---|---:|
| [Circular cylinder](example/circular_cylinder/) | Steady base flow | 8 |
| [Parabolic leading edge](example/parabolic_nonlinear/) | Nonlinear acoustic response | 1 |
| [Parabolic leading edge, LNS](example/parabolic_lns/) | AD-LNS response | 1 |
| [Blunt cone at incidence](example/blunt_cone_aoa1/) | Three-dimensional base flow | 8 |
| [HIFiRE-5-type elliptic cone](example/hifire5_acoustic/) | Steady / nonlinear acoustic response | 1 / 4 |

Each case contains a configuration, reference grid, restart data, and
representative outputs. See the [case index](example/README.md) for data scope.
The examples do not include all manuscript production histories.

## Documentation

- [Build and execution](docs/usage.md)
- [Source modules](src/README.md)
- [Verification status](docs/verification.md)
- [Changelog](docs/CHANGELOG.md)
- [Citation metadata](CITATION.cff)

Version 0.1.0 identifies the current research snapshot; no formal release tag
is published. The repository currently contains a [BSD-3-Clause license](LICENSE).
Part I covers a calorically perfect gas. Five-species and eleven-species models
and multigrid extensions are planned for Part II, not included here.

Report reproducible problems through [GitHub Issues](https://github.com/zzc-thu/LASTFit/issues).
