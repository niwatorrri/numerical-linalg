function u = multigrid(A, f, mode)
% Multigrid scheme
% Input: Equation Au=f, chosen smoother (line-gs / point-gs)
% Output: u

    v1 = 1;
    v2 = 1;
    n2 = sqrt(size(f, 1));
    n1 = (n2-1)/2;
    u = zeros(size(f, 1), 1);
    
    if n2 < 8
        % Solve directly
        u = gauss(A, f, 2 * n2);
        return
    else
        % Pre-smoothing
        for k = 1:v1
            u = symgs(A, f, u, mode);
        end
        
        % Restriction
        I = sparse(n1, n2);
        for j = 1:n1
            I(j, 2*j-1:2*j+1) = [1/4, 1/2, 1/4];
        end
        I = kron(I, I);
        
        newr = I * (f - A * u);
        newA = I * A * I' * 4;
        e = multigrid(newA, newr, mode);
        
        % Prolongation
        newe = I' * (e * 4);
        
        % Correction
        u = u + newe;
        
        % Post-smoothing
        for k = 1:v2
            u = symgs(A, f, u, mode);
        end
    end
end