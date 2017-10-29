for n = 5:30
    A = diag(ones(n, 1)) - tril(ones(n, n), -1);
    A(:, n) = 1;
    x = rand(n, 1) * 200 - 100;   % Uniform in [-100,100]
    [numerr, accerr] = acc(A, A * x, x);
    fprintf('&%d&%e&%e&\\\\\n', n, numerr, accerr);
end