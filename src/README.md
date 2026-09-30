# Source Code

`src/` contains the LASTFit Part I Fortran source and the default GNU Make
build. Run `make` in this directory, or `make -C src` from the repository root,
to create `SFSolver`.

## Module Map

| Area | Principal files |
|---|---|
| Program and input | `main.f90`, `Read_Para.f90`, `Module_Constant.f90` |
| State and parallel decomposition | `Module_CFD.f90`, `Module_Parallel.f90` |
| Grid and metrics | `Initial_Grid.f90`, `Singular_Initial_Grid.f90`, `Calculate_Jaco.f90` |
| Flow initialization | `Initial_Fields.f90`, `Singular_Initial_Fields.f90` |
| Spatial discretization | `Calculate_Flux.f90`, `Module_FD_1st.f90`, `Module_FD_5th.f90` |
| Shock fitting | `Module_SFitting.f90`, `Calculate_ShockAcDeri.f90` |
| Time advancement | `Module_RK4.f90`, `Module_DPLUR.f90`, `Module_Krylov.f90` |
| Linearized solver | `Module_LNS.f90` and the generated tangent routines |
| Output and diagnostics | `Module_OutputParaView.f90`, `OutputResults.f90`, `Module_Monitor.f90` |

The Makefile is the authoritative list of compilation units for the default
release. Other source variants retained in this directory are not compiled
unless they are added explicitly to that dependency graph.

## Build Variables

| Variable | Default | Purpose |
|---|---|---|
| `FC` | `mpif90` | MPI Fortran compiler wrapper |
| `FFLAGS` | GNU release flags | Compilation and preprocessing flags |
| `LDFLAGS` | empty | Linker-only flags |
| `LIBS` | `-llapack -lblas` | Numerical libraries |

Command-line assignments override the defaults:

```bash
make FC=mpiifx FFLAGS="-O3 -cpp -heap-arrays" LIBS="-qmkl"
```

No BLAS/LAPACK implementation is vendored. See `LAPACK/README.md`.

## Runtime Paths

The executable reads `Config.cfg` from its current working directory and uses
fixed relative output paths: `INIT/`, `CheckFiles/`, `JACO/`, `RESU/`,
`RESU_STEADY/`, `LNSResults/`, and `Pert/`. Create these directories before
launching outside `src/`. The repository root README gives a complete staging
example.

## Automatic Differentiation

The generated tangent routines required by the LNS solver are versioned with
the release. TAPENADE is needed only to regenerate them after changing the
corresponding nonlinear operations. See `AutoDiff/README.md`.
