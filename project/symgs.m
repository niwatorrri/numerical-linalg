function sol = symgs(A, f, u, mode)
% Symmetric Gauss-Seidel iteration
% Input: Equation Ax=f, current iterate u, mode (line or point)
% Output: u updated by one iterate of symmetric g-s

    % Line GS
    if strcmp(mode, 'line')
        n = sqrt(size(f, 1));        
        i = 1:n;
        
        for j = 1:1:n
            u(i) = A(i,i) \ (f(i) - A(i,:) * u + A(i,i) * u(i));
            i = i + n;
        end
        for j = n:-1:1
            i = i - n;
            u(i) = A(i,i) \ (f(i) - A(i,:) * u + A(i,i) * u(i));
        end
        
    % Point GS
    elseif strcmp(mode, 'point')
        n = size(f, 1);
        
        % Solve (D+L)x=r
        r = f - A * u;
        L = tril(A);
        for k = 1:n-1
            m = min(k+1, n);
            r(k) = r(k) / L(k, k);
            r(k+1:m) = r(k+1:m) - r(k) * L(k+1:m, k);
        end
        r(n) = r(n) / L(n, n);
        u = u + r;
        
        % Solve (D+U)x=r
        r = f - A * u;
        U = triu(A);
        for k = n:-1:2
            m = max(k-1, 1);
            r(k) = r(k) / U(k, k);
            r(m:k-1) = r(m:k-1) - r(k) * U(m:k-1, k);
        end
        r(1) = r(1) / U(1, 1);
        u = u + r;
    end
    sol = u;
end