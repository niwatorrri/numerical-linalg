########################################################
###   Numerical Linear Algebra Codes                 ###
###   Chapter 3: The Least Squares Problem           ###
###   Author: Hongyi Zhang                           ###
###   Date: 2017/11/06                               ###
########################################################

A total of 3 functions and 2 scripts are included.

Functions:

house.m:
- Return a Householder reflection Hx = alpha * e1
- Argument: x
- Return: [beta, v], H = I - beta * v * v’

qrfact.m:
- Return the QR factorization of a given matrix A
- Argument: A
- Return: [V, R] where V stores the Householder vectors

qrls.m:
- Solving the least squares problem ||Ax-b||_2 with QR
- Arguments: A, b
- Return: x

Scripts:

qrls_test.m:
- Test script for assignment problem 1
- Solving least squares with QR

qrsolve_test.m:
- Test script for assignment problem 2
- Solving linear systems with QR
