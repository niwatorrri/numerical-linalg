function sol = gauss(A, b, pivot)
    %% Solve linear system Ax=b with Gaussian elimination %%
    n = size(A, 1);
    if nargin < 3
        pivot = 'col';
    end
    
    %% LU decomposition %%
    for k = 1:n-1
        if strcmp(pivot, 'col')
            [~, idx] = max(abs(A(k:n, k)));
            idx = idx + k - 1;
            A([k,idx], :) = A([idx,k], :);
            b([k,idx]) = b([idx,k]);
        end
        A(k+1:n, k) = A(k+1:n, k) / A(k, k);
        A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - A(k+1:n, k) * A(k, k+1:n);
    end
    
    %% Solve Ly=b %%
    L = eye(n) + tril(A, -1);
    for j = 1:n-1
        b(j) = b(j) / L(j, j);
        b(j+1:n) = b(j+1:n) - b(j) * L(j+1:n, j);
    end
    b(n) = b(n) / L(n, n);
    
    %% Solve Ux=y %%
    U = triu(A);
    for j = n:-1:2
        b(j) = b(j) / U(j, j);
        b(1:j-1) = b(1:j-1) - b(j) * U(1:j-1, j);
    end
    b(1) = b(1) / U(1, 1);
    
    %% Output the solution %%
    sol = b;
end