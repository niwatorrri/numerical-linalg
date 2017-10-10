n = 40;
A = hilb(n);
b = sum(A, 2);
sol = ones(n, 1);

format long;
tic; x1 = solve(A, b, 'none'); toc
tic; x2 = solve(A, b, 'col'); toc
tic; x3 = solvespd(A, b, 'chol'); toc
tic; x4 = solvespd(A, b, 'ldlt'); toc
tic; x5 = A \ b; toc
norm(x1 - sol, inf)
norm(x2 - sol, inf)
norm(x3 - sol, inf)
norm(x4 - sol, inf)
norm(x5 - sol, inf)
