# Reference Grid

[Case guide](../README.md)

Open [Initial_grid.pvts](Initial_grid.pvts) with its one referenced VTS
piece to inspect the 81 x 51 x 8 grid. The MPI decomposition is 1 x 1.

These files are reference output, not runtime mesh inputs. The solver generates
its grid from the configuration and writes live grids to `INIT/`.
