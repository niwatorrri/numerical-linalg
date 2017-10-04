function sol = solvespd(A, b, mode)
    %% Solve linear system Ax=b with A symmetric positive definite %%
    if nargin < 3
        mode = 'chol';
    end
    
    %% Solve with Cholesky decomposition %%
    if strcmp(mode, 'chol')
        L = cholesky(A);
        b = solvelower(L, b);
        b = solveupper(L', b);
    end
    
    %% Solve with LDL^T decomposition %%
    if strcmp(mode, 'ldlt')
        [L, D] = ldlt(A);
        b = solvelower(L, b);
        b = solveupper(L', b ./ D);
    end
    sol = b;
end