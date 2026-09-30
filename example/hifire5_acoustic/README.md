# HIFiRE-5-Type Elliptic Cone

This directory records the non-axisymmetric three-dimensional steady and
nonlinear acoustic-forcing workflow for a HIFiRE-5-type elliptic cone.

## Configurations

| File | Purpose | Grid | MPI ranks | Status |
|---|---|---:|---:|---|
| `Config.cfg` | Compact nonlinear acoustic setup | 90 x 81 x 160 | 4 | Current release setup |
| `Config_steady.cfg` | Compact steady setup | 90 x 41 x 40 | 1 | Current release setup |
| `Config_E4_fine.cfg` | E-4 manuscript-run record | 160 x 151 x 160 | 192 | Historical provenance |

All three configurations use Mach 6, unit Reynolds number 10.2e6, and a
7-degree cone half-angle. The E-4 record specifies `epsilon=5e-4`,
`k_infty=3700`, `FixedDeltaTime=1e-7`, and an output interval of 2000 steps.

`Config_E4_fine.cfg` was recovered with a historical executable whose model
numbering differs from the current source. It also requires the external
192-rank base-flow restart. Treat it as provenance; do not substitute the
included four-rank grid or restart.

## Included Data

- `grid/steady/`: one-piece compact steady reference grid.
- `grid/unsteady/`: four-piece compact unsteady reference grid.
- `restart/unsteady/`: matching four-rank flow and fitted-shock restart.
- `output/steady/`: compact steady VTK field.
- `output/manuscript_fine/wall_harmonic/`: complete 192-piece wall harmonic
  used by the manuscript discussion.
- `SHA256SUMS`: checksums for all configurations and numerical files.

## Run the Compact Unsteady Setup

The primary `Config.cfg` starts a new four-rank calculation:

```bash
run_dir=/path/to/prepared/run
cd "$run_dir"
mpirun -np 4 ./SFSolver
```

To use the included checkpoint, copy `restart/unsteady/*` to `RESU/` and set
`IF_Continue_Calculate=1` in the staged configuration.

## Data Scope

The full 160 x 151 x 160 base-flow field, one E-4 perturbation snapshot, and
the full three-dimensional harmonic field are approximately 281 MB, 360 MB,
and 320 MB, respectively, and remain outside the Git tree. The repository
retains the compact wall-harmonic field required to inspect the manuscript
amplitude and phase result.
