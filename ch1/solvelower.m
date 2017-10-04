function sol = solvelower(L, b)
    %% Solve Ly=b where L lower triangular %%
    n = size(L, 1);
    for j = 1:n-1
        b(j) = b(j) / L(j, j);
        b(j+1:n) = b(j+1:n) - b(j) * L(j+1:n, j);
    end
    b(n) = b(n) / L(n, n);
    sol = b;
end