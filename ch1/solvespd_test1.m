n = 100;
A = diag(ones(n-1, 1), -1) + eye(n) * 10 + diag(ones(n-1, 1), 1);
rand('seed', 5);
b = rand(n, 1) * 200 - 100;

format long;
tic; x1 = solvespd(A, b, 'chol'); toc
tic; x2 = solvespd(A, b, 'ldlt'); toc
sol = A \ b;
norm(x1 - sol, inf)
norm(x2 - sol, inf)