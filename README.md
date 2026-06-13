# LASTFit

```
 _        _    ____ _____ _____ _ _
| |      / \  / ___|_   _|  ___(_) |_
| |     / _ \ \___ \ | | | |_  | | __|
| |___ / ___ \ ___) || | |  _| | | |_
|_____/_/   \_\____/ |_| |_|   |_|\__|

Three-dimensional high-order shock-fitting solver for hypersonic flows
```

LASTFit is a Fortran 90 software package for smooth hypersonic base-flow
computation and acoustic-forcing unsteady-field analysis using a high-order
shock-fitting finite-difference method.

The code treats the bow shock as a sharp moving computational boundary governed
by the Rankine-Hugoniot relations. This keeps the region between the wall and
the shock smooth, allowing high-order finite differences to be used for
boundary-layer stability and receptivity calculations without shock smearing.


## Main Features

- Three-dimensional body-fitted structured-grid shock-fitting formulation.
- Nonlinear steady and unsteady Navier-Stokes solver for hypersonic flows.
- Explicit and implicit time-advancement options for base-flow computation.
- Acoustic-forcing unsteady simulations for receptivity studies.
- Linearized Navier-Stokes shock-fitting solver obtained by perturbing the full
  coupled flow, grid, and shock-motion system.
- TAPENADE-generated tangent Fortran routines for LNS Jacobian operations.
- High-order finite differences: fifth-order upwind discretization for
  inviscid terms and sixth-order central discretization for viscous terms.
- MPI domain decomposition for CPU-based parallel calculations.
- ParaView-compatible VTK output for flow, shock, wall-history, and
  perturbation-field visualization.

## Validation and Demonstration Cases

The accompanying manuscript and figures document the following benchmark and
demonstration cases:

- Hypersonic viscous flow over a circular cylinder.
- Steady and unsteady flow over a parabolic leading edge.
- Linearized shock-fitting calculation over a parabolic leading edge.
- Three-dimensional blunt cone at 1 degree angle of attack.
- HIFiRE-5-type elliptic cone with steady base flow and fast-acoustic
  receptivity response.

The reported outputs include wall pressure, temperature contours, vorticity,
pressure derivatives, perturbation amplitudes, wall-pressure histories,
Fourier amplitudes and phases, and spectra.

## Requirements

LASTFit is intended for Linux-based high-performance computing environments.
The recommended release configuration is:

- Fortran 90 compiler with MPI support, for example Intel oneAPI Fortran with
  an MPI Fortran compiler wrapper.
- MPI runtime.
- BLAS/LAPACK or Intel MKL.
- VTK-compatible output workflow for post-processing.
- TAPENADE, only when regenerating the tangent routines. The generated tangent
  Fortran files should be included in the release source tree and compiled by
  the supplied Makefile.

GPU acceleration is not used.

## Repository Layout

A release source tree is expected to use a layout similar to:

```text
LASTFit/
  README.md
  LICENSE
  Makefile
  src/                 Fortran source files
  src_tapenade/        TAPENADE-generated tangent routines
  cases/               benchmark namelist files and grids
  scripts/             plotting and post-processing scripts
  Figure/              manuscript and README figures
  doc/                 manuscript or additional documentation
```

This working directory currently contains the CPC manuscript, validation
figures, and post-processing scripts. Before public release, add the solver
source tree, example input files, build files, and license information.

## Compiling

Load the compiler and MPI environment first. For Intel oneAPI, this is usually:

```bash
source /opt/intel/oneapi/setvars.sh
```

Then build the executable:

```bash
make clean
make
```

The released `Makefile` should document the selected compiler wrapper,
optimization flags, MPI settings, and BLAS/LAPACK or MKL linkage. If TAPENADE
tangent files are regenerated, rebuild the corresponding generated Fortran
objects before linking the LNS executable.

## Running a Case

Runtime options are read from a Fortran namelist input file. A typical MPI run
has the form:

```bash
mpirun -np <nproc> ./LASTFit < case.nml
```

or, if the executable reads the case file name from the command line:

```bash
mpirun -np <nproc> ./LASTFit case.nml
```

Replace `LASTFit`, `<nproc>`, and `case.nml` with the executable name, processor
count, and case-control file used in the released version.

## Input File

Each case should define the numerical configuration through namelist entries,
including:

- Geometry or model type.
- Analysis type: steady base flow, nonlinear unsteady shock-fitting, or LNS
  disturbance calculation.
- Grid dimensions and MPI partition.
- Freestream Mach number, Reynolds number, Prandtl number, wall temperature,
  and gas parameters.
- Wall condition: isothermal or adiabatic.
- Time-step, final time, restart, and residual-control settings.
- Spatial scheme and temporal advancement method.
- Acoustic-disturbance amplitude, frequency, incidence direction, and phase.
- Output interval and selected field/history variables.

Benchmark input files should be placed under `cases/` in the public release so
that the manuscript figures can be reproduced.

## Output Files

LASTFit calculations may write:

- Steady base-flow fields.
- Shock height, shock velocity, and shock-acceleration histories.
- Perturbation fields from nonlinear or linearized unsteady calculations.
- Wall-pressure and wall-temperature histories.
- Fourier amplitudes and phase distributions.
- Residual histories.
- ParaView-compatible VTK/PVTK files.

Post-processing scripts in this repository generate the validation figures used
in the manuscript, including wall-pressure maps, wall-normal profiles,
spectra, and harmonic amplitude/phase maps.

## Current Restrictions

The documented cases use:

- Perfect-gas thermodynamics.
- Sutherland-law viscosity.
- Structured shock-fitted grids.
- Single-domain calculations.
- No-slip isothermal or adiabatic wall boundary conditions.

Chemical nonequilibrium, thermal nonequilibrium, real-gas effects,
unstructured meshes, and general multi-block coupling are outside the current
documented scope.

## How to Cite

If you use LASTFit, please cite the software manuscript:

```bibtex
@article{Zhu_LASTFit,
  title   = {LASTFit: A Three-Dimensional High-Order Shock-Fitting Software Package for Smooth Hypersonic Base Flows and Acoustic-Forcing Unsteady Fields},
  author  = {Zhu, Zhichao and Xi, Youcheng and Fu, Song},
  journal = {Computer Physics Communications},
  year    = {to appear}
}
```

Update the bibliographic information after CPC submission or publication.

## License

The open-source license should be specified before public release. Add a
`LICENSE` file and update this section accordingly.

## Contact

For questions about LASTFit, please contact the authors listed in the software
manuscript or open an issue in the public repository after release.
