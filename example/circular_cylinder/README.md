# Circular Cylinder

[Case index](../README.md) | [Usage guide](../../docs/usage.md)

Steady viscous shock fitting for Mach 5.73 flow over a circular cylinder.
The solution is spanwise uniform and periodic; three-dimensional grid and MPI
operations remain active.

## Configuration

| Parameter | Value |
|---|---|
| Analysis | Steady nonlinear Navier-Stokes |
| Grid | 101 x 61 x 20 |
| MPI decomposition | 4 x 2 (8 ranks) |
| Mach number | 5.73 |
| Reynolds parameter (`Re_Ref`) | 2050 |
| Wall | Isothermal |
| Nonlinear continuation | Enabled |

## Included Data

- [Configuration](Config.cfg): runtime namelist.
- [Grid](grid/): eight-piece reference grid.
- [Restart](restart/): eight flow/shock file pairs.
- [Output](output/): distributed steady field.
- [Checksums](SHA256SUMS): configuration and numerical-file integrity.

## Run

Run the following Bash commands from the repository root after initializing
the MPI/compiler environment. The source configuration remains unchanged.

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

The solver writes generated grids to `INIT/` and nonlinear fields and restart
data to `RESU/`.

## Data Scope

The grid dimensions match the manuscript cylinder benchmark. Digitized
literature curves and figure-generation files are not distributed in this case.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md) for checks performed
and remaining numerical validation.
