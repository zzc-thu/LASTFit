# Parabolic Leading Edge: Nonlinear Acoustic Response

[Case index](../README.md) | [Usage guide](../../docs/usage.md)

Nonlinear time-domain response of a Mach 15 parabolic-leading-edge flow to
planar fast-acoustic forcing.

## Configuration

| Parameter | Value |
|---|---|
| Analysis | Unsteady nonlinear Navier-Stokes |
| Grid | 81 x 51 x 8 |
| MPI decomposition | 1 x 1 (1 rank) |
| Mach number | 15 |
| Reynolds parameter (`Re_Ref`) | 6026.6 |
| Acoustic amplitude | 5e-4 |
| Streamwise wave number | 15 |
| Fixed time step | 6e-5 |
| Nonlinear continuation | Enabled |

## Included Data

- [Configuration](Config.cfg): runtime namelist.
- [Grid](grid/): one-piece reference grid.
- [Restart](restart/): one nonlinear flow/shock file pair.
- [Checksums](SHA256SUMS): configuration and numerical-file integrity.

## Run

The restored checkpoint contains conservative variables and all three
coordinate records. Its conservative record and shock file are identical
to the previously distributed checkpoint. The matching compact base-flow
snapshot is shared with `parabolic_lns`; the grid, wall coordinates and
nominal physical parameters match. This does not establish convergence.

Run the following Bash commands from the repository root after initializing
the MPI/compiler environment. The source configuration remains unchanged.

```bash
set -e
make -C src
case_dir="$PWD/example/parabolic_nonlinear"
base_dir="$PWD/example/parabolic_lns/restart"
run_dir="$PWD/run/parabolic_nonlinear"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
cp "$case_dir"/restart/SF_Results*.flowsfg "$run_dir/RESU/"
cp "$case_dir"/restart/SF_Results*.shksfg "$run_dir/RESU/"
cp "$base_dir"/SF_Results*.flowsfg "$run_dir/RESU_STEADY/"
cp "$base_dir"/SF_Results*.shksfg "$run_dir/RESU_STEADY/"
# For time-continuous continuation, also copy the matching SF_Results*.time files.
(
  cd "$run_dir"
  mpirun -np 1 ./SFSolver
)
```

Nonlinear fields and restarts are written to `RESU/`; perturbation fields are
written to `Pert/` at the configured cadence.
Copy checkpoint `.time` files into `RESU/` when continuing an unsteady run.
Without time metadata, legacy files are treated as a warm start at zero time
and the solver prints a warning.

## Data Scope

The compact grid is 81 x 51 x 8, whereas the manuscript comparison uses
181 x 161 x 8. Full production time histories and Fourier input series are not
distributed in this case.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md) for checks performed
and remaining numerical validation.
