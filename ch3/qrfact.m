function [V, R] = qrfact(A)
    % Return QR factorization of a given matrix A
    % V stores v's for Householder vectors, R stores upper triangular part
    
    [m, n] = size(A);
    V = zeros(m, n);
    for j = 1:min(n, m-1)
        [beta, v] = house(A(j:m, j));
        A(j:m, j:n) = A(j:m, j:n) - (beta * v) * (v' * A(j:m, j:n));
        V(j:m, j) = v(1:m-j+1);
    end
    R = triu(A);
end