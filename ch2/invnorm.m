function ans = invnorm(A)
    n = size(A, 1);
    x = ones(n, 1) / n;
    while 1
        w = gauss(A', x);
        z = gauss(A, sign(w));
        if abs(norm(z, inf) - z' * x) < 1e-10
            ans = norm(w, 1);
            return;
        else
            [~, j] = max(abs(z));
            x = zeros(n, 1); x(j) = 1;
        end
    end
end