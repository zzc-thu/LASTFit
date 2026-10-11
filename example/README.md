# Example Cases

[Overview](../README.md) | [Build and execution](../docs/usage.md) | [Verification](../docs/verification.md)

The five examples cover steady base flows, nonlinear acoustic response, and
AD-LNS evolution. Configurations, reference grids, restart or base-flow inputs,
and representative solver output are supplied. Run calculations in isolated directories.

## Case Index

| Case | Analysis | Grid in primary configuration | MPI ranks | Input files for default run |
|---|---|---:|---:|---|
| [Circular cylinder](circular_cylinder/) | Steady nonlinear | 101 x 61 x 20 | 8 | `restart/` to `RESU/` |
| [Parabolic leading edge](parabolic_nonlinear/) | Nonlinear acoustic | 81 x 51 x 8 | 1 | `restart/` to `RESU/` |
| [Parabolic leading edge, LNS](parabolic_lns/) | AD-LNS | 81 x 51 x 8 | 1 | `restart/` to `RESU_STEADY/` |
| [Blunt cone at incidence](blunt_cone_aoa1/) | Steady nonlinear | 120 x 151 x 40 | 8 | None; checkpoint provided |
| [HIFiRE-5-type elliptic cone](hifire5_acoustic/) | Nonlinear acoustic | 90 x 81 x 160 | 4 | `restart/unsteady/` to `RESU_STEADY/` |

HIFiRE also includes a one-rank steady configuration with a separate grid.

## Organization

```text
case/
  Config.cfg       Runtime namelist
  grid/            Reference VTK grid
  restart/         Flow and shock restart data
  output/          Representative solver fields and output documentation
  README.md        Parameters, run commands, and data scope
  SHA256SUMS       Configuration and numerical-file checksums
```

The grid is generated from the configuration. Runtime directory names differ
from this reference-data layout; use the complete commands in the case READMEs.
The `output/` directory contains one representative PVTS dataset with all
referenced VTS pieces. It is not the live output destination for a calculation.

## Integrity and Scope

Each manifest covers the supplied configuration and numerical data, excluding
documentation. Run `sha256sum -c SHA256SUMS` from the selected case directory.
Open grid PVTS headers, with every referenced piece available, to inspect meshes.

The MPI decomposition must match the restart. The supplied inputs are not
complete manuscript reproduction datasets. Preserve them and modify only
staged copies during a run. Postprocessed spectra, harmonic maps and
plotting datasets are excluded from this repository.
