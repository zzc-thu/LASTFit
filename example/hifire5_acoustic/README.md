# HIFiRE-5-Type Elliptic Cone

[Case index](../README.md) | [Usage guide](../../docs/usage.md)

Non-axisymmetric three-dimensional steady flow and nonlinear fast-acoustic
response on an elliptic cone.

## Configurations

| Configuration | Analysis | Grid | MPI ranks | Scope |
|---|---|---:|---:|---|
| [Config.cfg](Config.cfg) | Nonlinear acoustic | 90 x 81 x 160 | 4 | Compact example |
| [Config_steady.cfg](Config_steady.cfg) | Steady nonlinear | 90 x 41 x 40 | 1 | Compact example |

The compact configurations use Mach 6, `Re_Ref=10.2e6`, and a 7-degree cone half-angle.
The compact acoustic example specifies `epsilon=5e-4` and `k_infty=15`.
The compact acoustic MPI decomposition is `2 x 2`, matching the supplied
four-rank checkpoint and reference grid.
The production record instead specifies `k_infty=3700`, fixed time step
`1e-7`, and output every 2000 steps. The configurations are not interchangeable.

The incompatible [E4 production configuration](../../docs/archive/README.md)
is stored separately as a historical record, not a runnable example.

## Included Data

- [Grids](grid/): separate compact steady and unsteady reference grids.
- [Restart](restart/): four-rank flow/shock checkpoint for the compact setup.
- [Checksums](SHA256SUMS): configuration and numerical-file integrity.

## Run

The primary configuration starts a new four-rank acoustic calculation.
The files retained under the historical path `restart/unsteady/` match an
archived `RESU_STEADY` snapshot byte for byte. They provide the compact
`90 x 81 x 160`, `2 x 2` base-flow input; the commands below also copy them
to `RESU_STEADY/`. Their provenance does not establish convergence.
The distributed `Config_steady.cfg` uses a different grid and is not that
base flow. Run these Bash commands from the repository root.

```bash
set -e
make -C src
case_dir="$PWD/example/hifire5_acoustic"
run_dir="$PWD/run/hifire5_acoustic"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
cp "$case_dir"/restart/unsteady/SF_Results*.flowsfg "$run_dir/RESU_STEADY/"
cp "$case_dir"/restart/unsteady/SF_Results*.shksfg "$run_dir/RESU_STEADY/"
(
  cd "$run_dir"
  mpirun -np 4 ./SFSolver
)
```

For compact continuation, copy `restart/unsteady/*.flowsfg` and
`restart/unsteady/*.shksfg` into the staged `RESU/` directory, then set
`IF_Continue_Calculate=1` in the staged configuration. Keep `npx0=2` and
`npz0=2`. The supplied historical checkpoint has no time metadata and is
only a warm start, not a time-continuous acoustic restart. New checkpoints
also write per-rank `.time` files; copy those files with the flow/shock files.

For the one-rank steady calculation, stage `Config_steady.cfg` as
`Config.cfg` in a separate run directory and launch one MPI rank:

```bash
set -e
make -C src
case_dir="$PWD/example/hifire5_acoustic"
run_dir="$PWD/run/hifire5_steady"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config_steady.cfg" "$run_dir/Config.cfg"
(
  cd "$run_dir"
  mpirun -np 1 ./SFSolver
)
```

## Data Scope

Only the compact reference grids and matching restart/base-flow inputs are
distributed. Generated flow fields and postprocessed spectra or harmonic maps
are excluded. The compact example does not reproduce the production-grid spectra.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md).
