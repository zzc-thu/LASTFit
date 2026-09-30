# Example Cases

[Overview](../README.md) | [Build and execution](../docs/usage.md) | [Verification](../docs/verification.md)

The five examples cover steady base flows, nonlinear acoustic response, and
AD-LNS evolution. Case data are representative solver records, not complete
reproductions of every manuscript figure.

## Case Index

| Case | Analysis | Grid in primary configuration | MPI ranks | Restart for default run |
|---|---|---:|---:|---|
| [Circular cylinder](circular_cylinder/) | Steady nonlinear | 101 x 61 x 20 | 8 | `restart/` to `RESU/` |
| [Parabolic leading edge](parabolic_nonlinear/) | Nonlinear acoustic | 81 x 51 x 8 | 1 | `restart/` to `RESU/` |
| [Parabolic leading edge, LNS](parabolic_lns/) | AD-LNS | 81 x 51 x 8 | 1 | `restart/` to `RESU_STEADY/` |
| [Blunt cone at incidence](blunt_cone_aoa1/) | Steady nonlinear | 120 x 151 x 40 | 8 | None; checkpoint provided |
| [HIFiRE-5-type elliptic cone](hifire5_acoustic/) | Nonlinear acoustic | 90 x 81 x 160 | 4 | None; checkpoint provided |

HIFiRE also includes a one-rank steady configuration and a separate 192-rank
production-run configuration record. The latter is not a self-contained
runnable example; see its [case README](hifire5_acoustic/README.md).

## Organization

```text
case/
  Config.cfg       Runtime namelist
  grid/            Reference VTK grid
  restart/         Flow and shock restart data
  output/          Representative reference fields
  README.md        Parameters, run commands, and data scope
  SHA256SUMS       Configuration and numerical-file checksums
```

The grid is generated from the configuration. Runtime directory names differ
from this reference-data layout; use the complete commands in the case READMEs.

## Integrity and Scope

Each manifest covers the supplied configuration and numerical data, excluding
documentation. Run `sha256sum -c SHA256SUMS` from the selected case directory.
Open PVTS headers, with every referenced piece available, to inspect fields.

The MPI decomposition must match the restart. Compact grids and representative
outputs cannot replace production-grid time histories for convergence or
spectral analysis. All examples retain their original numerical data; only
the staged configuration and runtime outputs should be modified during a run.
