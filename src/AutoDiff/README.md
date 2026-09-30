# Automatic-Differentiation Support

[Source map](../README.md) | [Usage guide](../../docs/usage.md)

## Input Sources

- [CalculateShockAc.f90](CalculateShockAc.f90): shock-acceleration operation.
- [CalculateShockU0.f90](CalculateShockU0.f90): associated shock-state operation.

## Generated Sources

The default build compiles tangent routines in the parent source directory:
`Calculate_AcDeri.f90`, `Calculate_UshkDeri_CV*.f90`,
`Complete_Flux_un_d.f90`, `CompleteFlux_D_mpi.f90`, and
`LNS_shockBC_d.f90`. Their compilation is defined in the [Makefile](../Makefile).

TAPENADE is not a build-time dependency. Regeneration is a separate developer
operation; a complete, validated regeneration command is not supplied here.
Changes to generated routines require operator-level consistency checks and
LNS/nonlinear response comparisons before replacing the distributed sources.
