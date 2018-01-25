% Test for Gaussian elimination
N = 64;

for e = [1, 1e-1, 1e-3, 1e-5]
    T = gallery('tridiag', N-1, -1, 2, -1) * N^2;
    A = kron(T, speye(N-1)) + e * kron(speye(N-1), T);

    x = repmat(1:N-1, 1, N-1)' / N;
    y = repelem(1:N-1, N-1)' / N;
    u = sin(pi*x) .* sin(pi*y);
    f = pi^2 * (1+e) * u;

    tic
    sol = gauss(A, f, N-1);
    toc
    norm(u - sol) / N
    norm(A \ f - sol)
end