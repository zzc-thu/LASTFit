# BLAS/LAPACK Dependency

[Source map](../README.md) | [Build guide](../../docs/usage.md#build)

No BLAS/LAPACK source or binaries are distributed in this directory. The
default Makefile links an external implementation using `-llapack -lblas`.

Set `LIBS` for the implementation installed on the target system. For Intel
oneAPI with MKL, the documented configuration is:

```bash
make -C src FC=mpiifx FFLAGS="-O3 -cpp -heap-arrays" LIBS="-qmkl"
```

Run this command from the repository root after initializing the oneAPI
environment. Compiler examples are not a tested-platform guarantee.
