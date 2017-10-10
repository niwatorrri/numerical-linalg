function [L, D] = ldlt(A)
    %% Return LDL^T factorization of a symmetric positive definite matrix %%
    n = size(A, 1);
    v = zeros(n, 1);
    for j = 1:n
        v = A(j,:)' .* diag(A);
        A(j, j) = A(j, j) - A(j, 1:j-1) * v(1:j-1);
        A(j+1:n, j) = (A(j+1:n, j) - A(j+1:n, 1:j-1) * v(1:j-1)) / A(j, j);
    end
    L = eye(n) + tril(A, -1);
    D = diag(A);
end