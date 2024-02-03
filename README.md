[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20450.svg)](http://dx.doi.org/10.5281/zenodo.20450)

# HubHam

Hubbard Hamiltonian 
===================

(under GNU GENERAL PUBLIC LICENSE v2)

This program can perform Exact diagonalization calculations of various types of
model Hamiltonians. It is especially optimized for the Hubbard Hamiltonian
type model Hamiltonians. The core feature which the program is specialized for
is the adressing of determinant in an efficient manner to quickly construct the
Hamiltonian non-zero matrix-elements. Once the Hamiltonian is constructed in 
its sparse format, it is stored in distributed memory for all linear algebra
operations.

The main work of diagonalizing the Hamiltonian is performed using PETSc and
SLEPc helper functions. These functions return the eigenvectors which are 
not stored to disk by default due to their large size. 

This project also contains subroutines which analyze the wavefunction in 
its distributed memory form and calculates the various observables. The
output of the program are the energies and the various observables such as 
the total Spin, various Spin-Spin correlation functions, and one-and two-body
density matrices.

Note: Caluculation of S2 can be deactivated to save memory and perform 
larger matrix diagonalizations.

_Dependencies_
---------------

  1. [PETSc](https://www.mcs.anl.gov/petsc/documentation/installation.html) 

  2. [SLEPc](http://slepc.upv.es/documentation/instal.htm)

  3. [IGraph](http://igraph.org/c/)

  4. [slaterlib](http://github.com/scemama/slaterlib)

_Compiling_
------------

  1. Configure with explicit paths of installation directories of PETSc and SLEPc.

```shell
./autogen.sh
./configure CFLAGS="-I{PATH-TO-INCLUDE-SLEPC-PETSC-FILES} " LDFLAGS="-L{PATH-FOR-LIBRARY-FILES} -lpetsc -lslepc -lpetsc -lslater -ligraph -lm"  CC=mpicc
```


  2. Make the executable

```shell
make
make install
```


_Using HubHam_
---------------

  1. The HubHam program requires an input file which 
   has is in the `graphml` format. A few sample files
   have been provided in the `data/` directory.


  2. running HubHam

```shell
mpiexec -n [nprocs] ./bin/ex1 -eps_nev 10 -f data/c8h8_mma.graphml 
```

_Publications using this code_
-------------------------------

  1. High-Spin Chains and Crowns from Double-Exchange Mechanism [doi:10.3390/cryst6040039](http://www.dx.doi.org/10.3390/cryst6040039)
