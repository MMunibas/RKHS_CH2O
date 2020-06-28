# RKHS_CH2O
Potential energy surfaces for formaldehyde using ab initio energies and gradients

**Requirements**

(1) LAPACK library

**Calculating the kernel coefficients using the training energies and gradients**

A program file (training.f90) is given and it can be compiled as

`gfortran training.f90 -llapack

Run as ./a.out and provide training.dat file with training geometries, energies and gradients. It will write the coefficients in coeff.dat file.

**Calculating the energies and gradinets using the coefficients**

A program file (training.f90) is given and it can be compiled as

`gfortran rkhs_pes.f90

Run as ./a.out and provide the geometry of CH2O for which the energy and gradiets to be computed in inp.xyz file. Keep the coeff.dat file in the same directory.
