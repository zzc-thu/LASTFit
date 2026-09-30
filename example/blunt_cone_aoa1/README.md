# Blunt Cone at 1 Degree Angle of Attack

[Case index](../README.md) | [Usage guide](../../docs/usage.md)

Asymmetric three-dimensional steady shock fitting and polar-axis treatment for
Mach 6 flow over a spherical-nosed cone.

## Configuration

| Parameter | Value |
|---|---|
| Analysis | Steady nonlinear Navier-Stokes |
| Grid | 120 x 151 x 40 |
| MPI decomposition | 4 x 2 (8 ranks) |
| Mach number | 6 |
| Reynolds parameter (`Re_Ref`) | 10000 |
| Cone half-angle | 5 degrees |
| Angle of attack | 1 degree |
| Nonlinear continuation | Disabled |

## Included Data

- [Configuration](Config.cfg): runtime namelist.
- [Grid](grid/): eight-piece reference grid.
- [Restart](restart/): eight flow/shock file pairs for optional continuation.
- [Output](output/): distributed three-dimensional steady field.
- [Checksums](SHA256SUMS): configuration and numerical-file integrity.

## Run

Run the following Bash commands from the repository root after initializing
the MPI/compiler environment. The source configuration remains unchanged.

```bash
set -e
make -C src
case_dir="$PWD/example/blunt_cone_aoa1"
run_dir="$PWD/run/blunt_cone_aoa1"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
(
  cd "$run_dir"
  mpirun -np 8 ./SFSolver
)
```

The default configuration starts a new calculation. For continuation, copy
`restart/*.flowsfg` and `restart/*.shksfg` into the staged `RESU/` directory,
then set `IF_Continue_Calculate=1` in the staged configuration.

## Data Scope

The example grid is 120 x 151 x 40. It differs from the 30 x 101 x 40 grid
used for the manuscript wall-pressure comparison. Literature curves and
publication plotting data are not distributed in this case.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md) for checks performed
and remaining numerical validation.
