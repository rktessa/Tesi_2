clc;
clear all;
close all;

%% sirs nullcline
% call the class functions container
obj = functionsContainer;

fig = 0;
s = linspace(0,1,100);
beta = 0.5;
delta = 0.3;
gamma = 0.22;
R0 = beta/gamma
time =1500;


fig = fig +1;
figure(fig)
hold on 
plot(s, first_null(s,beta,delta))
xline((gamma/beta), LineWidth= 1.5)
legend('s nullcline', 'i nullcline')
hold off

%% simulation
[taxis,xaxis,yaxis, zaxis] = obj.SIRS(beta,gamma,delta,time);

fig = fig +1;
figure(fig)
hold on
plot(taxis,xaxis, 'linewidth',1.5 )
plot(taxis,yaxis, 'linewidth',1.5 )
plot(taxis,zaxis, 'linewidth',1.5 )
title("SIRS model")
legend('S','I','R')

%% function 
function y = first_null(s,beta, delta)
    y = (- delta* (s-1))./(beta.*s + delta);

end