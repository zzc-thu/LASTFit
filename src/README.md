# Source Code

[Overview](../README.md) | [Build and execution](../docs/usage.md)

The Makefile builds `SFSolver` from the default Fortran compilation units.
Use `make -C src` from the repository root with an initialized MPI toolchain.

## Module Map

| Area | Principal files |
|---|---|
| Program and input | [main.f90](main.f90), [Read_Para.f90](Read_Para.f90), [Module_Constant.f90](Module_Constant.f90) |
| Flow state and MPI | [Module_CFD.f90](Module_CFD.f90), [Module_Parallel.f90](Module_Parallel.f90) |
| Grid and metrics | [Initial_Grid.f90](Initial_Grid.f90), [Singular_Initial_Grid.f90](Singular_Initial_Grid.f90), [Calculate_Jaco.f90](Calculate_Jaco.f90) |
| Initialization | [Initial_Fields.f90](Initial_Fields.f90), [Singular_Initial_Fields.f90](Singular_Initial_Fields.f90) |
| Spatial discretization | [Calculate_Flux.f90](Calculate_Flux.f90), [Module_FD_1st.f90](Module_FD_1st.f90), [Module_FD_5th.f90](Module_FD_5th.f90) |
| Shock fitting | [Module_SFitting.f90](Module_SFitting.f90), [Calculate_ShockAcDeri.f90](Calculate_ShockAcDeri.f90) |
| Time advancement | [Module_RK4.f90](Module_RK4.f90), [Module_DPLUR.f90](Module_DPLUR.f90), [Module_Krylov.f90](Module_Krylov.f90) |
| Linearized evolution | [Module_LNS.f90](Module_LNS.f90) and generated tangent sources |
| Output and diagnostics | [Module_OutputParaView.f90](Module_OutputParaView.f90), [OutputResults.f90](OutputResults.f90), [Module_Monitor.f90](Module_Monitor.f90) |

The [Makefile](Makefile) is the authoritative list of default compilation units.
Alternative source variants are not necessarily enabled or validated. In
particular, source files with multigrid or harmonic-analysis names do not
establish those capabilities as supported Part I modes.

## Build Variables

| Variable | Default | Purpose |
|---|---|---|
| `FC` | `mpif90` | MPI Fortran wrapper |
| `FFLAGS` | `-O3 -cpp -ffree-line-length-0 -fcray-pointer` | Compile flags |
| `FLFLAGS` | `$(FFLAGS)` | Compile-rule alias |
| `LDFLAGS` | Empty | Linker flags |
| `LIBS` | `-llapack -lblas` | External libraries |

Command-line assignments override defaults. `FFLAGS` is used for compilation,
not automatically passed to the final link; supply link options through
`LDFLAGS` or `LIBS`. No BLAS/LAPACK implementation is vendored; see the
[dependency note](LAPACK/README.md).

## Execution

The executable reads `Config.cfg` from its working directory. An isolated run
needs `INIT/`, `CheckFiles/`, `JACO/`, `RESU/`, `RESU_STEADY/`,
`LNSResults/`, and `Pert/`. The [usage guide](../docs/usage.md#run) and
[case READMEs](../example/README.md) provide complete commands.

## Automatic Differentiation

The default LNS build includes generated tangent sources. TAPENADE is required
only for regeneration after modifying differentiated operations. See the
[AD support note](AutoDiff/README.md).
