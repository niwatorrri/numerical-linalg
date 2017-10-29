function [numerr, accerr] = acc(A, b, x)
    %% Estimated relative error %%
    tildex = gauss(A, b);
    r = A * tildex - b;
    numerr = cond(A, inf) * norm(r, inf) / norm(b, inf);
    
    %% Accurate relative error if given real x %%
    if nargin == 3
        accerr = norm(x - tildex, inf) / norm(x, inf);
    else
        accerr = -1;
    end
end