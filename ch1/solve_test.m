n = 84;
A = diag(ones(n-1, 1) * 8, -1) + eye(n) * 6 + diag(ones(n-1, 1), 1);
b = [7; ones(n-2, 1) * 15; 14];
format long;
sol = ones(84, 1);
tic; x1 = solve(A, b, 'none'); toc
tic; x2 = solve(A, b, 'col'); toc
norm(x1 - sol, inf)
norm(x2 - sol, inf)