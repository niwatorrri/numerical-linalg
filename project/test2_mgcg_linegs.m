% Test for MGCG with line-gs as the smoother
fd = fopen('result_mgcg.txt', 'w');

for N = [32, 64, 128, 256, 512, 1024]
    for e = [1, 1e-1, 1e-3, 1e-5, 1e-7]
        T = gallery('tridiag', N-1, -1, 2, -1) * N^2;
        A = kron(T, speye(N-1)) + e * kron(speye(N-1), T); % row-major order
        % A = e * kron(T, speye(N-1)) + kron(speye(N-1), T); % column-major order

        x = repmat(1:N-1, 1, N-1)' / N;
        y = repelem(1:N-1, N-1)' / N;
        u = sin(pi*x) .* sin(pi*y);
        f = pi^2 * (1+e) * u;
        u0 = ones((N-1)^2, 1);

        tic;
        [sol, cnt] = mgcg(A, f, u0, 'line');
        time = toc;
        err = norm(u - sol) / N;
        fprintf(fd, ['N = %d, eps = %e, time = %f, n_iter = %d, ' ...
                     'l2-error = %e\n'], N, e, time, cnt, err);
        mesh((1:N-1)/N, (1:N-1)/N, vec2mat(sol - u));
    end
end

fclose(fd);