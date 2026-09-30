# LASTFit Part I

Three-dimensional MPI shock fitting for perfect-gas hypersonic base flows and
receptivity calculations.

LASTFit advances the flow field together with a fitted bow-shock boundary. The
post-shock domain therefore remains smooth and can be discretized with
high-order finite differences without representing the shock as a smeared
capturing layer. Part I integrates steady and nonlinear unsteady Navier-Stokes
calculations with an automatic-differentiation-generated linearized
Navier-Stokes (LNS) solver in the same moving-shock framework.

The current research release is **v0.1.0**. It contains the perfect-gas solver,
source-generated grids, MPI decomposition, representative validation data, and
compact examples. Reacting five-species and eleven-species models and
multigrid extensions are planned for Part II and are not included here.

## Capabilities

- Body-fitted shock fitting for two-and-a-half-dimensional and fully
  three-dimensional structured grids.
- Steady and time-accurate nonlinear Navier-Stokes calculations.
- Fast- and slow-acoustic forcing for receptivity studies.
- AD-generated LNS evolution of the coupled flow, grid, and shock-motion
  system.
- Fifth-order upwind inviscid discretization and sixth-order central viscous
  discretization.
- Explicit Runge-Kutta and implicit solver paths.
- MPI decomposition in the streamwise and circumferential directions.
- Distributed VTK output for ParaView.

## Repository Layout

```text
LASTFit/
  src/                 Fortran source and Makefile
  src/AutoDiff/        AD support files for shock acceleration
  example/             release cases and representative data
  CITATION.cff         software citation metadata
  CHANGELOG.md         release history
  LICENSE              BSD-3-Clause license
  VERSION              release version
```

Each case under `example/` uses the following public layout:

```text
case_name/
  Config.cfg
  grid/
  restart/
  output/
  README.md
  SHA256SUMS
```

`grid/` and `output/` contain reference artifacts. The solver generates its
grid from `Config.cfg`; it does not read the PVTS file in `grid/` as an input.
`restart/` contains files that can be staged into a run directory when a case
continues from a saved base flow.

## Requirements

- A POSIX-like Linux or macOS environment.
- GNU Make.
- An MPI Fortran compiler wrapper such as `mpif90`, `mpifort`, `mpiifort`, or
  `mpiifx`.
- A BLAS/LAPACK implementation providing the standard double-precision LAPACK
  interface.
- An MPI launcher such as `mpirun` or the scheduler-specific equivalent.
- ParaView or another VTK reader for optional visualization.

TAPENADE is not required to build the release because the generated tangent
routines are included in `src/`. It is needed only when regenerating those
routines from modified nonlinear source.

## Build

The default build uses an MPI-enabled GNU Fortran wrapper and system
BLAS/LAPACK:

```bash
make -C src
```

The executable is written to `src/SFSolver`. Compiler and library settings can
be overridden without editing the Makefile. For example, an Intel oneAPI build
can be requested with:

```bash
make -C src \
  FC=mpiifx \
  FFLAGS="-O3 -cpp -heap-arrays" \
  LIBS="-qmkl"
```

To remove compiler products while retaining generated solver data:

```bash
make -C src clean
```

`make -C src clean-output` also empties the runtime output directories under
`src/` and should be used only when those files are no longer needed.

## Run

`SFSolver` reads `Config.cfg` from its current working directory. It also uses
fixed relative directory names, so each calculation should run in an isolated
directory containing:

```text
run-directory/
  Config.cfg
  SFSolver
  INIT/
  CheckFiles/
  JACO/
  RESU/
  RESU_STEADY/
  LNSResults/
  Pert/
```

The MPI process count must equal `npx0 * npz0` in `Config.cfg`. The following
example stages the circular-cylinder case, whose configuration continues from
an eight-rank restart:

```bash
make -C src

case_dir="$PWD/example/circular_cylinder"
run_dir="$PWD/run/circular_cylinder"

mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
cp "$case_dir"/restart/* "$run_dir/RESU/"

cd "$run_dir"
mpirun -np 8 ./SFSolver
```

For an LNS case, copy its base-flow restart into `RESU_STEADY/` instead of
`RESU/`. A new calculation with `IF_Continue_Calculate=0` does not require a
flow restart. Consult the case README before staging a run.

## Runtime Outputs

Depending on the selected analysis mode, a calculation writes to:

- `INIT/`: generated grid in distributed VTK format;
- `JACO/`: grid metrics and Jacobian diagnostics;
- `RESU/`: nonlinear fields and flow/shock restart files;
- `Pert/`: nonlinear perturbation fields;
- `LNSResults/`: linearized disturbance fields and restarts;
- `CheckFiles/`: diagnostic slices and monitoring data.

The `output/` directory inside each release example is a read-only reference
snapshot, not the directory used by a live calculation.

## Examples

| Case | Mode | Release grid | MPI ranks | Included reference output |
|---|---|---:|---:|---|
| Circular cylinder | Steady nonlinear | 101 x 61 x 20 | 8 | Final distributed field |
| Parabolic leading edge | Nonlinear acoustic | 81 x 51 x 8 | 1 | Representative unsteady field |
| Parabolic leading edge | AD-LNS | 81 x 51 x 8 | 1 | Final representative LNS field |
| Blunt cone, 1 degree angle of attack | Steady nonlinear | 120 x 151 x 40 | 8 | Three-dimensional field |
| HIFiRE-5-type elliptic cone | Steady and nonlinear acoustic | Multiple, documented in case README | 1-192 | Steady field and compact wall-harmonic field |

These are release assets, not a claim that every production-grid time history
from the manuscript is stored in Git. Large histories and full three-dimensional
HIFiRE fields are intended for a versioned data archive.

## Data Integrity

Each example has a `SHA256SUMS` file covering its configuration and numerical
data. Verify a case from its directory with:

```bash
sha256sum -c SHA256SUMS
```

The PVTS headers in `output/` reference only pieces committed with the same
case. Case READMEs document any production-grid data that are intentionally
outside this repository.

## Numerical Scope

Part I is limited to:

- a calorically perfect gas;
- Sutherland-law viscosity;
- structured, shock-fitted, single-domain grids;
- no-slip isothermal or adiabatic walls;
- CPU-based MPI execution.

Reacting chemistry, thermal nonequilibrium, unstructured meshes, general
multi-block coupling, and GPU acceleration are outside this release.

## Citation

Citation metadata are provided in `CITATION.cff`. Until the accompanying CPC
article receives final bibliographic information, cite the software release as:

```bibtex
@software{lastfit_part1_0_1_0,
  author  = {Zhu, Zhichao and Xi, Youcheng and Fu, Song},
  title   = {LASTFit Part I: A Three-Dimensional Shock-Fitting Solver for
             Perfect-Gas Hypersonic Base Flows and Receptivity},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/zzc-thu/LASTFit}
}
```

## License

LASTFit Part I is distributed under the BSD-3-Clause license. See `LICENSE`.

## Support

Use the GitHub issue tracker for reproducible bug reports and release-data
questions. Include the commit, compiler and MPI versions, `Config.cfg`, MPI
rank count, and the shortest input that reproduces the problem.
