function [c,ceq] = constraint(z)

alfa = 0.5;
beta = 20;
N = 40;
nx = 6;
nu = 2;
lambda_t = 2*pi/3;

c1 = @(lambda,e) alfa*exp(-beta.*(lambda - lambda_t).^2) - e;

c = c1(z(1:nx:nx*N),z(5:nx:N*nx));

ceq = [];