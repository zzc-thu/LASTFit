# Automatic-Differentiation Support

This directory contains source fragments used in the automatic-differentiation
workflow for shock acceleration:

- `CalculateShockAc.f90`: nonlinear shock-acceleration operation;
- `CalculateShockU0.f90`: associated shock-state operation.

The tangent routines required by the default LNS build are already generated
and stored in the parent `src/` directory, including `Calculate_AcDeri.f90`,
`Calculate_UshkDeri_CV*.f90`, `Complete_Flux_un_d.f90`,
`CompleteFlux_D_mpi.f90`, and `LNS_shockBC_d.f90`. TAPENADE is therefore not a
build-time dependency for the published release.

Regenerating tangent code is a developer workflow. Generated files should be
reviewed and tested before replacing the versioned routines in `src/`.
