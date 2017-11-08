function [beta, v] = house(x)
    % Return a Householder transformation H=I-bvv'
    % Keep the first entry positive and vanish others
    
    m = size(x, 1);
    x = x / norm(x, inf);
    sigma = x(2:m)' * x(2:m);
    v = [1; x(2:m)];
    if abs(sigma) < eps(1)
        beta = 2;
    else
        alpha = sqrt(x(1)^2 + sigma);
        if x(1) <= 0
            v(1) = x(1) - alpha;
        else
            v(1) = -sigma / (x(1) + alpha);
        end
        beta = 2 / (sigma + v(1) ^ 2);
    end
end