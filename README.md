# Computational-Linear-Algebra

This Git repository contains three main folders:
  1. eigensolvers
  2. lu_factorisation
  3. sparse_cholesky

## eigensolvers

This folder contains three main general eigenproblem solvers each using a different technique. The input for all the solvers are tridiagonal matrices corresponding to a mass and stiffness matrix. All solvers output the eigenvalues and eigenvectors of the system

The solvers are as follows

### qr_iteration

- Performs Cholesky decomposition on M
- Converts the problem into Hessenberg form
- Uses the QR iteration technique with Gram-Schmidt to solve the eigenproblem

### rayleigh_quotient_iteration

- Performs Cholesky decomposition on M and converts to standard eigenvalue form
- Uses the Rayleigh quotient technique with Householder deflation to solve the eigenproblem

### subspace_iteration_ritz

- Performs a tridiagonal LU decomposition
- Uses the subspace iteration method to find the eigenvalues and eigenvectors on a projected RITZ subspace

### tridiagonal_eigensensitivty

This contains scripts of an exmaple use of an eigenvalue sensitivity analysis for a tridiagonal matrix

## lu_factorisation

This main folder contains functions which perform lu decomposition using complete and partial pivoting. An example use case is provided in main

## sparse_cholesky

This contains functions which creates a random tridiagonal matrix and performs a sparse cholesky decomposition

