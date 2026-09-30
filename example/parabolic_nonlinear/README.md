# Parabolic Leading Edge: Nonlinear Acoustic Forcing

This case advances the nonlinear time-domain response of a Mach 15
parabolic-leading-edge flow to a planar fast-acoustic disturbance.

## Configuration

| Parameter | Value |
|---|---:|
| Analysis mode | Unsteady nonlinear Navier-Stokes |
| Grid | 81 x 51 x 8 |
| MPI decomposition | 1 x 1 (1 rank) |
| Mach number | 15 |
| Reynolds number | 6026.6 |
| Acoustic amplitude | 5e-4 |
| Streamwise wave number | 15 |
| Fixed time step | 6e-5 |
| Continuation | Enabled |

## Included Data

- `Config.cfg`: runtime namelist.
- `grid/`: one-piece reference grid.
- `restart/`: matching nonlinear flow and fitted-shock restart.
- `output/`: one representative unsteady VTK field.
- `SHA256SUMS`: checksums for all configuration and numerical files.

## Run

Stage the restart under `RESU/` and launch one MPI rank:

```bash
run_dir=/path/to/prepared/run
cp restart/* "$run_dir/RESU/"
cd "$run_dir"
mpirun -np 1 ./SFSolver
```

Nonlinear fields are written to `RESU/`; perturbation fields sampled at the
configured output interval are written to `Pert/`.

## Data Scope

This compact 81 x 51 x 8 release case is internally consistent. It is not the
181 x 161 x 8 production grid cited for the manuscript comparison, and the
repository does not include the full production time history or Fourier input
series.
