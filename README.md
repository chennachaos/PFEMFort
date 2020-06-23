# ParallelFEM-Fortran
Parallel programming for Finite Element Analysis using FORTRAN and PETSc

This code includes several parallel implementations for the Lapalace equation and linear elasticity in 2D and 3D. This code has been tested on HPC machines at Swansea University.

Please report if you find any issues/bugs.

## Compilation steps
* Enter `bin` directory using `cd bin`
* Modify the makefile as per the location of the PETSc library files in your system
* Execute `make <exename>`. This will also compile all the dependencies.
* For example, to build `tetrapoissonparallelimpl1` use the command `make tetrapoissonparallelimpl1`


