# Calculation Output

[Case guide](../README.md)

Open [Result_00000000.pvts](steady/Result_00000000.pvts) with its VTS piece
to inspect the compact steady field on the 90 x 41 x 40 grid. It is not an
acoustic time history or a convergence acceptance dataset. Postprocessed
spectra, harmonic maps and production-grid results are not supplied.

In each isolated run directory,
the solver writes nonlinear fields and checkpoints to `RESU/`, generated grids
to `INIT/`, and monitoring data to `CheckFiles/`. The acoustic calculation
also writes disturbance fields to `Pert/`; its `RESU_STEADY/` files are
base-flow inputs. Use separate run directories for steady and acoustic cases.
