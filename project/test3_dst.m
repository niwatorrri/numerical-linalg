% Test for DST method
N = 128;
e = 1e-3;

T = gallery('tridiag', N-1, -1, 2, -1) * N^2;
l = ones(N-1, 1);
x = (1:N-1)' * l' / N;
y = l * (1:N-1) / N;
u = sin(pi*x) .* sin(pi*y);
F = pi^2 * (1+e) * u / N^2;

% DST method routine
tic
Fp = idst(dst(F)');
d = 4 * sin((1:N-1) * pi / (2 * N)).^2;
Up = Fp ./ (d' * l' + e * l * d);
U = dst(idst(Up')');
toc
mesh(x, y, U - u);