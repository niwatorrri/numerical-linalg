function sol = solve(A, b, pivot)
    %% Solve linear system Ax=b with Gaussian elimination %%
    n = size(A, 1);
    if nargin < 3
        pivot = 'none';
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
    b = solvelower(eye(n) + tril(A, -1), b);
    
    %% Solve Ux=y %%
    b = solveupper(triu(A), b);
    
    %% Output the solution %%
    sol = b;
end