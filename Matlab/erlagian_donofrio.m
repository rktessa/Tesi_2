% clc;
% close all;
% clear all;
%% Erlagian function


% vector x
x = 0:0.1:20;
a = 0.5; % delay 1/2 in days? 
n = 1; % exponentially fading memory
% n=2; % humped Erlagian memory
pdf_erl =   erl(x, a, 1);
pdf_erl_2 = erl(x, a, 2);

figure(1)
hold on
plot(x,pdf_erl,'LineWidth',2)
plot(x,pdf_erl_2,'LineWidth',2)



function [y] = erl(x,a, n)
% mean delay T= n/a
% standard deviation sigma = sqrt(n)/a
y = (a^n /factorial((n-1))) .* x.^(n-1) .* exp(-a.*x);

end