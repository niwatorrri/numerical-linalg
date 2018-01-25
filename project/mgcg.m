function [sol, cnt] = mgcg(A, b, x, mode)
% Multigrid preconditioned conjugate gradient
% Input: Equation Au=b, initial guess x, mode for multigrid
% Output: solution u and number of iterations cnt

    N = size(b, 1);
    r = b - A * x;
    z = multigrid(A, r, mode);
    r0 = eye(N, 1);
    z0 = eye(N, 1);
    p = zeros(N, 1);
    res = norm(r);
    cnt = 0;
    
    while norm(r) / res > 1e-6
        cnt = cnt + 1;
        
        % Update direction p
        s = r' * z;
        beta = s / (r0' * z0);
        p = z + beta * p;
        
        % Update solution x
        Ap = A * p;
        alpha = s / (p' * Ap);
        x = x + alpha * p;
        
        % Update residual r (and z)
        r0 = r; z0 = z;
        r = r - alpha * Ap;
        z = multigrid(A, r, mode);
    end    
    sol = x;
end