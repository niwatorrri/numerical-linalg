function x = qrls(A, b)
    % Solving least squares problem ||Ax-b||_2 with QR factorization
    % Also applicable to solving linear system Ax=b
    
    % Factorize A = QR = [Q1,Q2]*[R1,0]' and compute Q'b
    [m, n] = size(A);
    [V, R] = qrfact(A);
    for j = 1:min(n, m-1)
        v = V(j:m, j);
        beta = 2 / (v' * v);
        b(j:m) = b(j:m) - beta * (v' * b(j:m)) * v;
    end
    
    % Solve R1 * x = Q1' * b
    b = b(1:n);
    for j = n:-1:2
        b(j) = b(j) / R(j, j);
        b(1:j-1) = b(1:j-1) - b(j) * R(1:j-1, j);
    end
    b(1) = b(1) / R(1, 1);
    x = b;
end