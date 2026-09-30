# LASTFit examples

Each example has a stable, descriptive directory name and follows the same
layout:

```text
case_name/
  Config.cfg
  grid/
  restart/
  output/
  README.md
```

`Config.cfg` is the solver input. `grid/` stores the distributed initial grid,
`restart/` stores matching flow and fitted-shock restart files, and `output/`
stores only the reduced data needed to check or reproduce the documented
result. Large time histories should be archived in the release dataset rather
than committed repeatedly to the source repository.

## Case index

| Directory | Capability | Configuration in this repository | Release status |
|---|---|---|---|
| `circular_cylinder` | Steady shock fitting and three-dimensional infrastructure | 101 x 61 x 20 | Configuration checked; binary data pending |
| `parabolic_nonlinear` | Nonlinear fast-acoustic forcing | 81 x 51 x 8 | Provisional configuration; manuscript-grid data pending |
| `parabolic_lns` | AD-generated LNS disturbance evolution | 81 x 51 x 8 | Configuration checked; reduced LNS data pending |
| `blunt_cone_aoa1` | Asymmetric three-dimensional base flow and polar treatment | 120 x 151 x 40 | Provisional configuration; manuscript-grid provenance pending |
| `hifire5_acoustic` | Non-axisymmetric steady and nonlinear acoustic response | 90 x 81 x 160 unsteady; 90 x 41 x 40 steady | Provisional configurations; final acoustic dataset pending |

The directories are deliberately explicit about incomplete data. A case should
be marked reproducible only after its configuration, grid, restart, executable
revision, reduced output, plotting command, and checksums are tied to the same
run.

## Data policy

- Keep source and small text inputs in Git.
- Keep one reduced, inspectable result in `output/` when practical.
- Store large time histories in a versioned release archive with a DOI.
- Do not commit files copied from a different geometry or grid merely to fill a
  directory.
- Do not document machine-specific absolute paths.

Each case README records the remaining work and the expected contents of its
data directories.
