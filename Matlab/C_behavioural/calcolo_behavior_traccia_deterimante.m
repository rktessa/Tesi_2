clc;
clear all;
close all;

% Equilibrium
syms x y z k_1 k_2 lambda_1 lambda_2 
eq_H(x,y) = -k_1*x*y - k_2*(1-x-y)*x + lambda_1*y + lambda_2 * (1-x-y);
eq_C(x,y) = +k_1*x*y - lambda_1*y;
% Equilibrium point
X = 1; Y = 0; %All Heedles in case 1
% The Jacobian of this equation in this case is:
Jac(x, y) = jacobian([eq_H, eq_C], [x, y])
% Computing the Trace
Trac(x,y) = trace(Jac(x,y))
% Determinant
d(x,y) = det(Jac)