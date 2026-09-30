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

## Versioning

Version 0.1.0 denotes the research snapshot in the citation metadata.
No formal version tag, GitHub Release, or archival DOI is published by this
documentation revision. Record the exact commit when using the source.

The existing BSD-3-Clause license file is retained. A formal release remains
separate from this documentation update and requires author approval.
