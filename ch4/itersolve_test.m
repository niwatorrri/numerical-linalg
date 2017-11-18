N = 100;
n = N - 1;
h = 1 / N;
a = 0.5;
eps = [1, 0.1, 0.01, 0.0001];

for j = 1:4
    % Accurate solution %
    x = h * (1:n);
    f = @(x)(a*x + (1-a) * (1-exp(-x/eps(j))) ./ (1-exp(-1/eps(j))));
    y = f(x');

    % Numerical solution %
    A = diag(ones(n-1, 1) * eps(j), -1) - ...
        eye(n) * (2 * eps(j) + h) + ...
        diag(ones(n-1, 1) * (eps(j) + h), 1);
    b = ones(n, 1) * a * h^2;
    b(n) = b(n) - eps(j) - h;

    [sol, rho, k] = itersolve(A, b, 'Jacobi');
    fprintf('Jacobi & %.4f & %d & %.4f \\\\\n', rho, k, norm(sol - y));
    [sol, rho, k] = itersolve(A, b, 'G-S');
    fprintf('G-S & %.4f & %d & %.4f \\\\\n', rho, k, norm(sol - y));
    rhoB = max(abs(eig( diag(diag(A).^(-1)) * A - eye(n) )));
    omega = 2 / (1 + sqrt(1 - rhoB ^ 2));
    [sol, rho, k] = itersolve(A, b, 'SOR', omega);
    fprintf('SOR & %.4f & %d & %.4f \\\\\n', rho, k, norm(sol - y));
end
