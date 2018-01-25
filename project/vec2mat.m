function U = vec2mat(u)
% Convert a vector to matrix
% Input: vector u of dimension n^2
% Output: matrix U of dimension n

    n = sqrt(size(u, 1));
    U = zeros(n, n);
    
    for i = 1:n
        U(i,:) = u((i-1)*n+1:i*n);
    end
end

