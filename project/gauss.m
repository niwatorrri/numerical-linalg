function sol = gauss(A, b, band)
% Gaussian elimination for symmetric banded matrix
% Input: Equation Ax=b and the bandwidth of A
% Output: x

    n = size(b, 1);
    
    % LU decomposition
    for k = 1:n-1
        m = min(k+band, n);
        A(k+1:m, k) = A(k+1:m, k) / A(k, k);
        A(k+1:m, k+1:m) = A(k+1:m, k+1:m) - A(k+1:m, k) * A(k, k+1:m);
    end
    
    % Solve Ly=b
    L = speye(n) + tril(A, -1);
    for k = 1:n-1
        m = min(k+band, n);
        b(k) = b(k) / L(k, k);
        b(k+1:m) = b(k+1:m) - b(k) * L(k+1:m, k);
    end
    b(n) = b(n) / L(n, n);
    
    % Solve Ux=y
    U = triu(A);
    for k = n:-1:2
        m = max(k-band, 1);
        b(k) = b(k) / U(k, k);
        b(m:k-1) = b(m:k-1) - b(k) * U(m:k-1, k);
    end
    b(1) = b(1) / U(1, 1);
    sol = b;
end