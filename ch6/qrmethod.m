function schur = qrmethod(A)
    % Compute eigenvalues for unsymmetric matrix A
    % Assumption: A is upper Hessenberg
    % Output: diagonal identical to A's Schur normal form
    
    n = size(A, 1);
    H = A;
    tol = 1e-12;
    
    % QR iteration
    % Reduction to an upper quasi-triangular matrix
    p = 1; q = n + 1;
    while q > 3
        while q > 3
            if abs(H(q-2, q-3)) < tol; q = q - 2; H(q, q-1) = 0;
            elseif abs(H(q-1, q-2)) < tol; q = q - 1; H(q, q-1) = 0;
            else break
            end
        end
        if q > 3
            p = q - 1;
            while p > 1 && abs(H(p, p-1)) > tol
                p = p - 1;
            end
            H(p:q-1, p:q-1) = francis(H(p:q-1, p:q-1));
        end
    end
    
    % Check if 2x2 blocks have real eigenvalues
    q = 2;
    while q <= n
        while q <= n && abs(H(q, q-1)) < tol
            q = q + 1;
        end
        if q <= n
            p = q - 1;
            s = H(p,p) + H(q,q);
            t = H(p,p) * H(q,q) - H(p,q) * H(q,p);
            delta = s^2 - 4*t;
            if delta >= 0
                H(p,p) = (s + sqrt(delta)) / 2;
                H(q,q) = (s - sqrt(delta)) / 2;
                H(p,q) = 0; H(q,p) = 0;
            end
        end
        q = q + 1;
    end
    schur = H;
end

function H1 = francis(H)
    % Francis QR step
    
    n = size(H, 1);
    m = n - 1;
    
    s = H(m,m) + H(n,n);
    t = H(m,m) * H(n,n) - H(m,n) * H(n,m);    
    x = H(1,1)^2 + H(1,2) * H(2,1) - s * H(1,1) + t;
    y = H(2,1) * (H(1,1) + H(2,2) - s);
    z = H(2,1) * H(3,2);
    
    for k = 0:n-3
        [beta, v] = house([x, y, z]');
        q = max(1, k);
        H(k+1:k+3, q:n) = H(k+1:k+3, q:n) - (beta * v) * (v' * H(k+1:k+3, q:n));
        r = min(k+4, n);
        H(1:r, k+1:k+3) = H(1:r, k+1:k+3) - (H(1:r, k+1:k+3) * (beta * v)) * v';
        x = H(k+2, k+1);
        y = H(k+3, k+1);
        if k < n-3
            z = H(k+4, k+1);
        end
    end
    
    [beta, v] = house([x, y]');
    H(n-1:n, n-2:n) = H(n-1:n, n-2:n) - (beta * v) * (v' * H(n-1:n, n-2:n));
    H(1:n, n-1:n) = H(1:n, n-1:n) - (H(1:n, n-1:n) * (beta * v)) * v';
    H1 = H;
end

function [beta, v] = house(x)
    % Return a Householder transformation H=I-bvv'
    % Keep the first entry positive and vanish others
    
    m = size(x, 1);
    x = x / norm(x); % L2-norm normalization
    sigma = x(2:m)' * x(2:m);
    v = [1; x(2:m)];
    
    if sqrt(sigma) < eps(1)
        beta = 0;
    else
        alpha = sqrt(x(1)^2 + sigma);
        if x(1) <= 0
            v(1) = x(1) - alpha;
        else
            v(1) = -sigma / (x(1) + alpha);
        end
        beta = 2 / (sigma + v(1)^2);
    end
end