# Chapter 1: Direct Methods for Linear Systems

### Chapter overview

This chapter basically contains:
- Solving triangular linear systems
- Solving general linear systems: LU decomposition (with pivoting)
- Symmetric positive definite systems: Cholesky and LDL^T decompositions

### Code description

A total of 3 scripts and 6 functions are included.

Functions:

solvelower.m:
- Solve a lower triangular system Ly=b
- Arguments: Matrix L and vector b
- Return value: Solution y

solveupper.m:
- Solve an upper triangular system Ux=y
- Arguments: Matrix U and vector y
- Return value: Solution b

solve.m:
- Solve a linear system Ax=b with LU decomposition
- Arguments: Matrix A, vector b and pivoting option pivot
  - pivot = ‘none’: no pivoting
  - pivot = ‘col’: partial pivoting
- Return value: Solution x
- Internal call: solvelower.m, solveupper.m

cholesky.m:
- Return Cholesky factor of a symmetric positive definite matrix
- Argument: Matrix A
- Return value: Factor L

ldlt.m:
- Return LDL^T decomposition of a symmetric positive definite matrix
- Argument: Matrix A
- Return value: Matrix L and vector D

solvespd.m:
- Solve a symmetric positive definite system Ax=b
- Arguments: Matrix A, vector b and mode option
  - mode = ‘chol’: solving with Cholesky decomposition
  - mode = ‘ldlt’: solving with LDL^T decomposition
- Return value: Solution x
- Internal call: solvelower.m, solveupper.m, cholesky.m, ldlt.m

Scripts:

solve_test.m:
- Test script for assignment problem 1

solvespd_test1.m:
- Test script for assignment problem 2(1)

solvespd_test2.m:
- Test script for assignment problem 2(2)
