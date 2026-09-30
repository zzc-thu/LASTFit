# Blunt Cone at 1 Degree Angle of Attack

This case exercises asymmetric three-dimensional steady shock fitting and the
polar-axis treatment for Mach 6 flow over a spherical-nosed cone.

## Configuration

| Parameter | Value |
|---|---:|
| Analysis mode | Steady nonlinear Navier-Stokes |
| Grid | 120 x 151 x 40 |
| MPI decomposition | 4 x 2 (8 ranks) |
| Mach number | 6 |
| Reynolds number | 10000 |
| Cone half-angle | 5 degrees |
| Angle of attack | 1 degree |
| Continuation | Disabled |

## Included Data

- `Config.cfg`: runtime namelist for a new steady calculation.
- `grid/`: complete eight-piece reference grid.
- `restart/`: matching eight-rank checkpoint retained as a reference.
- `output/`: complete eight-piece three-dimensional VTK field.
- `SHA256SUMS`: checksums for all configuration and numerical files.

## Run

The supplied configuration starts a new calculation and therefore does not
read the files under `restart/`:

```bash
run_dir=/path/to/prepared/run
cd "$run_dir"
mpirun -np 8 ./SFSolver
```

To continue from the supplied checkpoint, copy `restart/*` to `RESU/` and set
`IF_Continue_Calculate=1` in the staged `Config.cfg`.

## Data Scope

The release assets form a consistent 120 x 151 x 40 case. They are separate
from the lower-resolution grid quoted for the manuscript wall-pressure
comparison; reference curves and plotting data for that comparison are not
included here.
