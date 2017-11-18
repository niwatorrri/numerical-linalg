function [sol, rho, k] = itersolve(A, b, method, omega)
    % Solve Ax=b with classical iterations
    % Return: sol - numerical solution
    %           q - L2-norm of the iteration matrix
    %           k - the iteration number
    
    % Default relaxation factor %
    if nargin < 4
        omega = 1;
    end

    % Decompose the matrix A=L+D+U %
    L = tril(A, -1);
    D = diag(diag(A));
    U = triu(A, 1);
    wb = omega * b;
    
    % Set the iteration matrix %
    if strcmp(method, 'Jacobi')
        N = D;
        M = D - A;
    else  % Gauss-Seidel or SOR %
        N = D + omega * L;
        M = D - omega * (A - L);
    end
    rho = max(abs(eig(N \ M)));
    q = norm(N \ M);
    c = q / (1 - q);
    if c < 0; c = 5e3; end
    
    % Start iteration %
    x = zeros(size(A, 1), 1);
    maxiter = 100000;
    for k = 1:maxiter
        x1 = solvelower(N, M * x + wb);
        if c * norm(x1 - x) < 5e-5
            break
        end
        x = x1;
    end
    
    sol = x;
    if k == maxiter
        fprintf('Not convergent after max iterations\n');
    end
end