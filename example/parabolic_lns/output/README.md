# Calculation Output

[Case guide](../README.md)

Open [LNSResults_00071.pvts](LNSResults_00071.pvts) with its VTS piece
to inspect the representative linearized field. It is not a full time history,
an operator-consistency test or a regression acceptance dataset.
No postprocessed data are supplied.

In the isolated run directory,
the solver writes linearized fields and perturbation checkpoints to
`LNSResults/`, generated grids to `INIT/`, and monitoring data to `CheckFiles/`.
The nonlinear base-flow files in `RESU_STEADY/` are inputs, not LNS results.
