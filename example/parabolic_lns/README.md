# Parabolic Leading Edge: AD-LNS Response

[Case index](../README.md) | [Usage guide](../../docs/usage.md)

Small-amplitude disturbance evolution about a nonlinear shock-fitted base flow
using the automatic-differentiation-generated linearized Navier-Stokes solver.

## Configuration

| Parameter | Value |
|---|---|
| Analysis | AD-LNS (`AnalysisType=3`) |
| Grid | 81 x 51 x 8 |
| MPI decomposition | 1 x 1 (1 rank) |
| Mach number | 15 |
| Reynolds parameter (`Re_Ref`) | 6026.6 |
| Acoustic amplitude | 5e-4 |
| LNS time step | 6e-6 |
| Nonlinear initialization continuation | Disabled |
| Perturbation continuation (`IF_Continue_LNS`) | Disabled |

## Included Data

- [Configuration](Config.cfg): runtime namelist.
- [Grid](grid/): one-piece reference grid.
- [Restart](restart/): one nonlinear base-flow/shock file pair.
- [Output](output/README.md): representative linearized field and runtime paths.
- [Checksums](SHA256SUMS): configuration and numerical-file integrity.

## Run

Run the following Bash commands from the repository root after initializing
the MPI/compiler environment. The source configuration remains unchanged.

```bash
set -e
make -C src
case_dir="$PWD/example/parabolic_lns"
run_dir="$PWD/run/parabolic_lns"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
cp "$case_dir"/restart/*.flowsfg "$case_dir"/restart/*.shksfg "$run_dir/RESU_STEADY/"
(
  cd "$run_dir"
  mpirun -np 1 ./SFSolver
)
```

The base flow is read from `RESU_STEADY/` even though perturbation continuation
is disabled. Linearized fields and perturbation restarts are written to
`LNSResults/`.

For this analysis, `LNS_dt=6e-6` controls integration; `FixedDeltaTime` and
the nonlinear CFL ramp do not set the LNS time step. Perturbation continuation
requires all four saved time levels, `IterationNumber.dat`, and
`RestartMetadata.dat`, in addition to the steady base flow. The metadata
record the next step, `LNS_dt`, grid and MPI partition. Incompatible or
metadata-free histories are rejected. Metadata must match the calculation
that produced the checkpoint. Keep the base
flow and forcing parameters unchanged as well.

## Data Scope

The supplied base flow supports short-run LNS checks at `LNS_dt=6e-6`.
The [verification guide](../../docs/verification.md) reports runtime,
restart, MPI partition and amplitude-scaling tests and their limitations.
Convergence and complete AD consistency validation remain incomplete.
Record the source commit and do not mix checkpoints from different operators.

The compact grid is 81 x 51 x 8, whereas the manuscript LNS/nonlinear comparison
uses 181 x 161 x 8. The complete production-grid comparison time series are not
distributed in this case.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md) for checks performed
and remaining numerical validation.
