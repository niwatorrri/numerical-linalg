##############################################################
###   Numerical Linear Algebra Codes                       ###
###   Chapter 4: Classical Iterations for Linear Systems   ###
###   Author: Hongyi Zhang                                 ###
###   Date: 2017/11/19                                     ###
##############################################################

A total of 2 functions and 1 script are included.

Functions:

solvelower.m:
- Solving Ly=b where L lower triangular
- Arguments: L, b
- Return: y

itersolve.m:
- Solving Ax=b with classical iterations
- Arguments: A, b, method, omega
  - method in {‘Jacobi’, ‘G-S’, ‘SOR’}
  - omega: relaxation factor for SOR, 1 by default
- Return: x, rho, k
  - rho: the spectral radius of the iteration matrix
  - k: the number of iterations till convergence, maxiter otherwise

Script:

itersolve_test.m:
- Test script for the assignment problem
- Solving a discretized difference equation with iterations
