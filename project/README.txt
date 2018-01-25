########################################################
###   Numerical Linear Algebra Course Project        ###
###   Solving 2D Poisson Equation with MGCG Method   ###
###   Author: Hongyi Zhang                           ###
###   Date: 2018/01/25                               ###
########################################################

A total of 4 main functions, 1 auxiliary function and 5 test scripts are included.

Main functions:

gauss.m:
% Gaussian elimination for symmetric banded matrix
% Input: Equation Ax=b and the bandwidth of A
% Output: x

symgs.m:
% Symmetric Gauss-Seidel iteration
% Input: Equation Ax=f, current iterate u, mode (line or point)
% Output: u updated by one iterate of symmetric g-s

multigrid.m:
% Multigrid scheme
% Input: Equation Au=f, chosen smoother (line-gs / point-gs)
% Output: u

mgcg.m:
% Multigrid preconditioned conjugate gradient
% Input: Equation Au=b, initial guess x, mode for multigrid
% Output: solution u and number of iterations cnt

Auxiliary function:

vec2mat.m:
% Convert a vector to matrix
% Input: vector u of dimension n^2
% Output: matrix U of dimension n

Test scripts:

test1_gauss.m:
% Test for Gaussian elimination

test2_mgcg_linegs.m:
% Test for MGCG with line-gs as the smoother

test2_mgcg_pointgs.m:
% Test for MGCG with point-gs as the smoother

test3_dst.m:
% Test for DST method

plot_l2error.m:
% Plot l2-error of MGCG method
