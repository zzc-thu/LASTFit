# Optional Utilities

[Overview](../README.md) | [Usage guide](../docs/usage.md)

These Python utilities are optional. They are not used to compile or execute
the Fortran solver and do not contain reference result datasets.

| Script | Purpose | Requirements |
|---|---|---|
| `check_examples.py` | Check configurations and flow/shock input records | Python 3.9+, NumPy |
| `prepare_case.py` | Stage an isolated short run with matching input files | Python 3.9+ |
| `check_lns_results.py` | Check generated LNS histories and compare runs | Python 3.9+, NumPy |
| `test_check_examples.py` | Regression tests for the input checker | Python 3.9+, NumPy |
| `test_check_lns_results.py` | Regression tests for the LNS checker | Python 3.9+, NumPy |

Run commands from the repository root; see the usage guide for staging and
the [verification guide](../docs/verification.md) for LNS checks.
The solver can also be run directly using the commands in each case README.
