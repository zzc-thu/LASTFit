# HIFiRE-5-Type Elliptic Cone

[Case index](../README.md) | [Usage guide](../../docs/usage.md)

Non-axisymmetric three-dimensional steady flow and nonlinear fast-acoustic
response on an elliptic cone.

## Configurations

| Configuration | Analysis | Grid | MPI ranks | Scope |
|---|---|---:|---:|---|
| [Config.cfg](Config.cfg) | Nonlinear acoustic | 90 x 81 x 160 | 4 | Compact example |
| [Config_steady.cfg](Config_steady.cfg) | Steady nonlinear | 90 x 41 x 40 | 1 | Compact example |
| [Config_E4_fine.cfg](Config_E4_fine.cfg) | Nonlinear acoustic | 160 x 151 x 160 | 192 | Production-run record only |

All configurations use Mach 6, `Re_Ref=10.2e6`, and a 7-degree cone half-angle.
The compact acoustic example specifies `epsilon=5e-4` and `k_infty=15`.
The production record instead specifies `k_infty=3700`, fixed time step
`1e-7`, and output every 2000 steps. The configurations are not interchangeable.

`Config_E4_fine.cfg` uses model numbering from a different source revision and
requires a 192-rank restart that is not distributed here. It is retained for
provenance, not as a runnable configuration for the current source.

## Included Data

- [Grids](grid/): separate compact steady and unsteady reference grids.
- [Restart](restart/): four-rank flow/shock checkpoint for the compact setup.
- [Output](output/): compact steady field and production-grid wall-harmonic field.
- [Checksums](SHA256SUMS): configuration and numerical-file integrity.

## Run

Run these Bash commands from the repository root. The primary configuration
starts a new four-rank acoustic calculation without reading the checkpoint.

```bash
set -e
make -C src
case_dir="$PWD/example/hifire5_acoustic"
run_dir="$PWD/run/hifire5_acoustic"
mkdir -p "$run_dir"/{INIT,CheckFiles,JACO,RESU,RESU_STEADY,LNSResults,Pert}
cp src/SFSolver "$run_dir/"
cp "$case_dir/Config.cfg" "$run_dir/"
(
  cd "$run_dir"
  mpirun -np 4 ./SFSolver
)
```

For compact continuation, copy `restart/unsteady/*.flowsfg` and
`restart/unsteady/*.shksfg` into the staged `RESU/` directory, then set
`IF_Continue_Calculate=1` in the staged configuration.

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

The wall-harmonic field supports inspection of the manuscript amplitude and
phase maps. Full production-grid base-flow fields, unsteady histories, and
three-dimensional harmonic volumes are not distributed in Git. The compact
example does not reproduce the production-grid spectra.

Iteration limits are research-run settings, not a short runtime guarantee.
See the [verification status](../../docs/verification.md).
