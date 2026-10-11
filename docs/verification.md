# Verification Status

[Overview](../README.md) | [Build and execution](usage.md) | [Examples](../example/README.md)

## Repository Checks

The following checks were performed for this documentation revision:

| Check | Result |
|---|---|
| Example SHA-256 manifests | Included configuration and numerical files match their checksums |
| Distributed VTK references | Every PVTS piece reference resolves to a supplied file |
| Configuration/documentation review | Grid sizes, primary MPI counts, continuation modes, and restart paths checked against configuration and source |
| Documentation navigation | Local Markdown links checked for existing targets |
| Documented Bash commands | Syntax checked with Bash; commands not executed as solver runs |

File-integrity checks do not establish the accuracy of the solver. A supplied
reference field is not a validated regression tolerance.

## Numerical and Runtime Scope

The manuscript reports benchmark comparisons and a small-amplitude
LNS/nonlinear response comparison. It does not report observed grid or
time-step convergence order, complete AD Taylor tests, or MPI performance
measurements. The HIFiRE example demonstrates nonlinear three-dimensional
acoustic response, not a controlled three-dimensional LNS verification.

The Linux Makefile and run instructions have been reviewed but have not been
compiled or executed in this revision environment. No runtime, memory,
restart-portability, or supported-compiler guarantee is inferred from the
presence of configuration and reference data.

## Local Restart Repairs (2026-10-11)

The source was compiled with Intel Fortran Classic 2021.11 and Intel MPI
2021.11 on Windows using a separate test harness, not the Linux Makefile.
A six-step nonlinear diagnostic and a three-step run continued to step six
agreed within 2.4e-22 in the interior conservative variables and 1.1e-19 in
shock velocity; shock height and saved coordinates were identical. Physical
time metadata were identical. Exterior halos were excluded from that
comparison because some contain non-finite values.

The diagnostic used the complete parabolic LNS base-flow checkpoint as a
nonlinear initial condition. It did not use the distributed nonlinear
checkpoint and is not manuscript validation. See the later data recovery.

The HIFiRE compact configuration now uses the checkpoint-compatible `2 x 2`
partition. The supplied checkpoint was traced to an archived steady-base-flow
snapshot with identical files and can be staged into `RESU_STEADY/`.
The parabolic LNS base flow passes file-shape and interior positivity checks,
but the short LNS run produces non-finite perturbations. LNS numerical
correctness and time-continuous LNS restart have not been established.

These checks do not establish grid/time convergence, AD consistency, Linux
support, or production-case performance.

After the repairs, a ten-step cylinder run on eight ranks exited normally;
16 VTK files had valid block lengths and finite values. A ten-step LNS run
on one rank also exited normally and saved next iteration 11, but its VTK
field contained 456,512 non-finite values. All 310,590 values in its raw t0
perturbation conservative record were non-finite. That is a failed numerical
test despite successful process cleanup and checkpoint metadata output.

## External Data Recovery (2026-10-11)

The complete parabolic nonlinear checkpoint was recovered from an archive.
Its conservative-variable record and shock file exactly match the previously
distributed single-record checkpoint; all three coordinate records are now
included. The compact parabolic LNS base flow has identical wall coordinates
and nominal physical settings and is shared for perturbation output.
With the recovered nonlinear checkpoint, a ten-step compact run exited
normally and both output VTK files passed block-length and finite-value
checks. This is a startup/output test, not convergence or spectral validation.

A 192-rank HIFiRE checkpoint was located externally and passed record-size
and interior positivity checks. It is not added to Git and does not resolve
the historical geometry-numbering mismatch in the E4 configuration.
Archived LNS histories contain finite interior perturbations, but their VTK
time series use different grid sizes. They are not a verified continuation
package for the current source or a substitute for the manuscript dataset.

## Current Repair Checks (2026-10-11)

Six active configurations and 44 flow/shock files pass sequential-record,
interior finite-value and density/pressure positivity checks. Stored exterior
halos in cylinder and HIFiRE data still include non-finite values; original
numerical data have not been replaced with invented values. The incompatible
E4 configuration is now isolated under `docs/archive/`.

The input-check tool passes six synthetic regression tests covering valid
records, missing coordinates, missing shock files, non-finite interior flow,
non-positive pressure and corrupt record trailers. These are input checks,
not accuracy or convergence tests.

The current compact nonlinear ten-step run exits normally; two VTK files
have valid block lengths and finite values. Rebuilding loaded LNS base-flow
metrics does not resolve its numerical failure: the first detected bad
interior cell is (63,51,1) on rank 0 at step 4, on the shock boundary.
New guards stop the run without saving the invalid perturbation checkpoint.
This LNS example remains a failing diagnostic, not a validated regression;
the staging tool requires `--allow-unverified` for it.

New LNS restart metadata verify the time step, grid, MPI partition and next
iteration. Legacy histories without verified metadata are rejected. This
safety check is not proof of a passing time-continuous LNS restart.
Four runtime guard regressions pass: missing metadata, wrong time step,
wrong grid and non-finite saved history are all rejected with nonzero exit.

## Release Status

Version 0.1.0 denotes the research snapshot in the citation metadata.
No formal version tag, GitHub Release, or archival DOI is published by this
documentation revision. Record the exact commit when using the source.

The existing BSD-3-Clause license file is retained. A formal release remains
separate from this documentation update and requires author approval.
