load('data1.mat')
n = size(y, 1);
A = [ones(n,1), t, t.^2];
x = qrls(A, y)
norm(A * x - y)

load('data2.mat')
n = size(price, 1);
A = [ones(n,1), data];
x = qrls(A, price)
norm(A * x - price)
rho = corrcoef(A * x, price);
R2 = rho(1,2)^2