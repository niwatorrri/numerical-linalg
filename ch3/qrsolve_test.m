n = 84;
A = diag(ones(n-1, 1) * 8, -1) + eye(n) * 6 + diag(ones(n-1, 1), 1);
b = [7; ones(n-2, 1) * 15; 14];
sol = ones(n, 1);
tic; x = qrls(A, b); toc
norm(x - sol, inf)

n = 100;
A = diag(ones(n-1, 1), -1) + eye(n) * 10 + diag(ones(n-1, 1), 1);
rand('seed', 5);
b = rand(n, 1) * 200 - 100;
sol = A \ b;
tic; x = qrls(A, b); toc
norm(x - sol, inf)

n = 40;
A = hilb(n);
b = sum(A, 2);
sol = ones(n, 1);
tic; x = qrls(A, b); toc
norm(x - sol, inf)