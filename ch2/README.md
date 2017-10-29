# Chapter 2: Error Analysis for Linear Systems

### Chapter overview
This chapter basically includes:
- Sensitivity analysis: condition number
- Numerical stability analysis: round-off error
- Accuracy estimation of a numerical solution

### Code description
A total of 3 functions and 2 scripts are included.

Functions:

gauss.m:
- Solving Ax=b using Gaussian elimination
- Arguments: A, b, pivot = {‘none’, ‘col’}
- Return: x

invnorm.m:
- Estimate L-infty norm of the inverse of A
- Argument: A
- Return: L-infty norm of A^{-1}

acc.m:
- Estimate the numerical solution accuracy for Ax=b
- Arguments: A, b, true x(optional)
- Return: Numerical error and accurate error (-1 if x not given)

Scripts:

infcond_test.m:
- Test script for assignment problem 1
- Estimate condition number for Hilbert matrices

acc_test.m:
- Test script for assignment problem 2
- Estimate solution accuracy and compare with the actual one
