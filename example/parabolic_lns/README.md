# Parabolic Leading Edge: AD-LNS

This case exercises the automatic-differentiation-generated linearized
Navier-Stokes solver about a fitted-shock nonlinear base flow.

## Configuration

| Parameter | Value |
|---|---:|
| Analysis mode | Linearized Navier-Stokes |
| Grid | 81 x 51 x 8 |
| MPI decomposition | 1 x 1 (1 rank) |
| Mach number | 15 |
| Reynolds number | 6026.6 |
| Acoustic amplitude | 5e-4 |
| LNS time step | 6e-6 |
| LNS continuation | Disabled |

`Config.cfg` uses the parameter names and `AnalysisType=3` expected by the
current `Read_Para.f90` namelist.

## Included Data

- `Config.cfg`: runtime namelist.
- `grid/`: one-piece reference grid.
- `restart/`: converged nonlinear base-flow and fitted-shock restart.
- `output/`: final representative record from the available LNS sequence.
- `SHA256SUMS`: checksums for all configuration and numerical files.

## Run

The LNS solver reads its nonlinear base flow from `RESU_STEADY/`:

```bash
run_dir=/path/to/prepared/run
cp restart/* "$run_dir/RESU_STEADY/"
cd "$run_dir"
mpirun -np 1 ./SFSolver
```

Linearized fields and disturbance restarts are written to `LNSResults/`.

## Data Scope

This compact 81 x 51 x 8 release case is intended for software inspection and
regression. It is not the 181 x 161 x 8 production-grid comparison described in
the manuscript, and the complete LNS/nonlinear time-series comparison is not
stored in this Git repository.
