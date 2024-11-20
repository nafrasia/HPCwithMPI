# HPCwithMPI

## Description
This is an introduction to High Performance Computing in C/C++. Massage Passing Interface (MPI), Basic Linear Algebra Subprograms (BLAS), and Linear Algebra Package (LAPACK) are utilized to solve mathematical problems.
Problems such as
  - Computing moments of a distribution
  - Matrix-Vector Multiplication
  - Matrix-Matrix Multiplication
  - 1D and 2D Poisson Problem using Jacobi method, and Conjugate Gradient

## Setup
All of the programs are written in C/C++ and need to be compiled before running. This requires C, C++, and MPI compilers. You can install these compilers using the following commands

### Linux
#### C/C++
`sudo apt install gcc g++`

#### OpenMPI
`sudo apt install openmpi-bin libopenmpi-dev` 

#### MPICH
`sudo apt install mpich libmpich-dev`

## Acknowledgement
This project was completed as part of the SC9505R by Dr. Colin Denniston at University of Western Ontario. I would like to thank Dr. Denniston for their guidance and support.

## License
This project is licensed under the GNU GPLv3 license. Please see the LICENSE file for details
