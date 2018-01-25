% Test for MGCG with point-gs as the smoother
for N = [128, 256]
    for e = [1e-1, 1e-3, 1e-5, 1e-7]
        T = gallery('tridiag', N-1, -1, 2, -1) * N^2;
        A = kron(T, speye(N-1)) + e * kron(speye(N-1), T);

        x = repmat(1:N-1, 1, N-1)' / N;
        y = repelem(1:N-1, N-1)' / N;
        u = sin(pi*x) .* sin(pi*y);
        f = pi^2 * (1+e) * u;
        u0 = ones((N-1)^2, 1);

        tic;
        [sol, cnt] = mgcg(A, f, u0, 'point');
        time = toc;
        err = norm(u - sol) / N;
        fprintf('N = %d, eps = %e, time = %f, n_iter = %d, err = %e\n', ...
                N, e, time, cnt, err);
    end
end