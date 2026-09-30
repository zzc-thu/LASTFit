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
- [Output](output/): representative unsteady field.
- [Checksums](SHA256SUMS): configuration and numerical-file integrity.

## Run

Run the following Bash commands from the repository root after initializing
the MPI/compiler environment. The source configuration remains unchanged.

```bash
set -e
make -C src
case_dir="$PWD/example/parabolic_nonlinear"
run_dir="$PWD/run/parabolic_nonlinear"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
cp "$case_dir"/restart/*.flowsfg "$case_dir"/restart/*.shksfg "$run_dir/RESU/"
(
  cd "$run_dir"
  mpirun -np 1 ./SFSolver
)
```

Nonlinear fields and restarts are written to `RESU/`; perturbation fields are
written to `Pert/` at the configured cadence.

## Data Scope

The compact grid is 81 x 51 x 8, whereas the manuscript comparison uses
181 x 161 x 8. Full production time histories and Fourier input series are not
distributed in this case.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md) for checks performed
and remaining numerical validation.
