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
- [Output](output/): representative LNS field.
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
`RestartMetadata.dat`, in addition to the steady base flow. New metadata
record the next step, `LNS_dt`, grid and MPI partition. Incompatible or
metadata-free legacy histories are rejected; never fabricate metadata for
an archive without first verifying its original settings. Keep the base
flow and forcing parameters unchanged as well.

## Data Scope

The initial short-run audit failed at step 4 because the transformed viscous
mass flux and its tangent were read without initialization. These components
are now explicitly zero, consistent with the perfect-gas mass equation.
With the supplied base flow and `LNS_dt=6e-6`, one-rank Windows runs reached
10 and 100 steps with finite interior histories and valid, finite VTK output.
A zero-forcing run preserved exactly zero interior flow and shock perturbations.

This is a limited runtime check, not a completed numerical validation.
Stored physical-boundary ghost cells still contain non-finite values.
A 10-step run restarted to step 20 differed from an uninterrupted run by
approximately 1.12e-12 in the interior histories, slightly exceeding the
diagnostic absolute tolerance of 1e-12. This strict restart check remains open.
The preparation tool retains `--allow-unverified` for diagnostic use.

The compact grid is 81 x 51 x 8, whereas the manuscript LNS/nonlinear comparison
uses 181 x 161 x 8. The complete production-grid comparison time series are not
distributed in this case.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md) for checks performed
and remaining numerical validation.
