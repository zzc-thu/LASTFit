# BLAS/LAPACK Dependency

LASTFit does not vendor BLAS or LAPACK source or binaries. The default Makefile
links the system libraries with `-llapack -lblas`.

Override `LIBS` when using another implementation, for example Intel MKL:

```bash
make FC=mpiifx FFLAGS="-O3 -cpp -heap-arrays" LIBS="-qmkl"
```

Keeping the numerical library external avoids machine-specific library paths
and allows each HPC environment to select its optimized implementation.
