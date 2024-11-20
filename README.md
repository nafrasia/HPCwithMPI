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

### Windows
You can use Windows Subsystem for Linux (WSL) to set up a linux environment and then use the instructions above to set up the C, C++, and MPI compilers. See https://learn.microsoft.com/en-us/windows/wsl/install for information on how to install WSL.

An alternative would be to use C/C++ compiler software like Visual Studio Code. For VS code, you can download from

https://code.visualstudio.com/download

Make sure to install the C/C++ extension. You also need to install msmpi to run MPI code.

## Boost
Boost is a set of libraries for the C++ programming language that provides support for tasks and structures such as linear algebra, pseudorandom number generation, multithreading, image processing, regular expressions, and unit testing. It contains 164 individual libraries (as of version 1.76)[Wikipedia]. For some of the programs, we use boost arrays. You can install the boost libraries using

`sudo apt install libboost-all-dev`

## BLAS AND LAPACK
Install LAPACK and BLAS using:

`sudo apt install liblapack-dev`

These packages need to be explicitly linked during compilation using `-lblas` and `-llapack`

## Acknowledgement
This project was completed as part of the SC9505R by Dr. Colin Denniston at University of Western Ontario. I would like to thank Dr. Denniston for their guidance and support.

## License
This project is licensed under the GNU GPLv3 license. Please see the LICENSE file for details
