# Build and Execution

[Overview](../README.md) | [Examples](../example/README.md) | [Verification](verification.md)

## Requirements

The documented target is Linux with GNU Make, an MPI Fortran compiler wrapper,
an MPI launcher, and external BLAS/LAPACK. ParaView is optional for inspecting
VTK output. Generated tangent routines are included, so TAPENADE is not needed
to build or run the LNS solver.

The commands below use Bash and start from the repository root. Compiler
commands are configuration examples, not a claim of successful testing on
every toolchain; see the [verification status](verification.md).

## Build

The default Makefile uses `mpif90` and links `-llapack -lblas`:

```bash
make -C src
```

For another GNU MPI wrapper:

```bash
make -C src clean
make -C src FC=mpifort
```

For Intel oneAPI, initialize the compiler/MPI environment first, then use:

```bash
make -C src clean
make -C src FC=mpiifx FFLAGS="-O3 -cpp -heap-arrays" LIBS="-qmkl"
```

The executable is `src/SFSolver`. Use a serial Make invocation; the documented
build procedure does not assume that all Fortran module dependencies support
`make -j`. See the [source map](../src/README.md) for compilation units.

`make -C src clean` removes compiler products, not numerical output.
`make -C src clean-output` deletes files from the runtime directories under
`src/`; it is a separate, destructive cleanup operation.

## Run

`SFSolver` reads `Config.cfg` from its working directory. Run each case in a
separate directory and retain the repository data as reference files.

For the circular-cylinder continuation:

```bash
set -e
make -C src
case_dir="$PWD/example/circular_cylinder"
run_dir="$PWD/run/circular_cylinder"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
cp "$case_dir"/restart/*.flowsfg "$case_dir"/restart/*.shksfg "$run_dir/RESU/"
(
  cd "$run_dir"
  mpirun -np 8 ./SFSolver
)
```

Each [case README](../example/README.md#case-index) gives a complete staging
command for its own configuration. The MPI process count is `npx0 * npz0`.
On a cluster, replace `mpirun` with the site-supported launcher.

### Restart Rules

| Analysis | Configuration control | Data directory in the run |
|---|---|---|
| Nonlinear fresh start | `IF_Continue_Calculate=0` | No nonlinear restart required |
| Nonlinear continuation | `IF_Continue_Calculate=1` | All rank-matched flow/shock files in `RESU/` |
| LNS with new perturbations | `AnalysisType=3`, `IF_Continue_LNS=0` | Nonlinear base flow in `RESU_STEADY/` |
| LNS continuation | `IF_Continue_LNS=1` | Base flow plus compatible perturbation restart in `LNSResults/` |

LNS base-flow data are required even when perturbation continuation is disabled.
A separate `IF_Continue_Calculate=1` also activates the nonlinear initialization
restart path. The included LNS configuration sets that control to zero.

Do not change grid dimensions or MPI partitioning when using an existing
restart. Restart files use Fortran unformatted records; portability between
compiler/runtime environments has not been established.

The supplied iteration limits are research-run controls, not short test
budgets. A shortened run requires changing the time-step limit and output
cadence in a staged copy of `Config.cfg`; it cannot establish convergence or
reproduce a manuscript spectrum. No few-minute runtime is asserted.

## Data Layout

```text
case/
  Config.cfg
  grid/
  restart/
  output/
  README.md
  SHA256SUMS
```

The solver generates its grid from `Config.cfg`. The PVTS files in `grid/`
are reference artifacts, not runtime mesh inputs. The case `output/` directory
stores representative fields and is not a live output destination.

| Runtime directory | Contents |
|---|---|
| `INIT/` | Generated grids |
| `CheckFiles/` | Diagnostics and monitoring data |
| `JACO/` | Grid metrics and Jacobian diagnostics |
| `RESU/` | Nonlinear fields and flow/shock restarts |
| `RESU_STEADY/` | Input base flow for LNS |
| `LNSResults/` | Linearized fields and perturbation restarts |
| `Pert/` | Nonlinear perturbation fields |

## Inspect and Verify

Open a `.pvts` header in ParaView with all referenced `.vts` pieces present
at their relative paths. Individual pieces contain only one MPI subdomain.

Verify the included configuration and numerical files from a case directory:

```bash
(cd example/circular_cylinder && sha256sum -c SHA256SUMS)
```

Documentation is excluded from these manifests. Matching checksums confirm
file integrity, not numerical correctness.

## Scope and Citation

Part I uses a calorically perfect gas, Sutherland-law viscosity, structured
shock-fitted grids, no-slip isothermal or adiabatic walls, and CPU-based MPI.
Reacting chemistry, thermal nonequilibrium, general multi-block coupling,
unstructured meshes, and GPU execution are outside this snapshot.

Use the [citation metadata](../CITATION.cff) and record the exact commit used.
No archival DOI or formal release tag is assigned. Bug reports should include
the commit, compiler/MPI versions, configuration, rank count, and relevant logs.
