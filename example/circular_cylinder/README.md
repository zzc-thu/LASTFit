# Circular Cylinder

This case exercises steady viscous shock fitting for Mach 5.73 flow over a
circular cylinder. The solution is spanwise uniform and periodic, while the
three-dimensional storage, metric, MPI, boundary-condition, and output paths
remain active.

## Configuration

| Parameter | Value |
|---|---:|
| Analysis mode | Steady nonlinear Navier-Stokes |
| Grid | 101 x 61 x 20 |
| MPI decomposition | 4 x 2 (8 ranks) |
| Mach number | 5.73 |
| Reynolds number | 2050 |
| Wall condition | Isothermal |
| Continuation | Enabled |

## Included Data

- `Config.cfg`: runtime namelist.
- `grid/`: complete eight-piece reference grid.
- `restart/`: eight flow files and eight fitted-shock files.
- `output/`: complete eight-piece final VTK field.
- `SHA256SUMS`: checksums for all configuration and numerical files.

## Run

Stage the restart under `RESU/`, then launch eight MPI ranks:

```bash
run_dir=/path/to/prepared/run
cp restart/* "$run_dir/RESU/"
cd "$run_dir"
mpirun -np 8 ./SFSolver
```

The executable regenerates the grid in `INIT/` and writes nonlinear fields and
updated restart files to `RESU/`.

## Data Scope

The grid dimensions match the circular-cylinder case described in the
manuscript. The repository includes the solver field used for release
inspection; digitized literature curves and publication plotting files are not
part of this solver example.
