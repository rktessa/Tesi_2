clc; 
close all; 
clear;

x = [0, 0, 0];
b = 0.5; % probability of contracting the disease

y = ode45(@epidemia, [0 60], [495; 5; 0]);  

N = x(1) + x(2) + x(3); % a fixed population for 500 individuals
g = 0.15; % recovery rate
xdot(1) = - b*x(2)*x(1)/N; % Susceptible differential equation
xdot(2) = b*x(2)*x(1)/N - g*x(2); % Infectious differential equation
xdot(3) = g*x(2); % Recovered differential equation
x0 = [495; 5; 0]; % Initial 495 Susceptible and 5 Infectious
t = linspace(0, 60, 600); % 60-week duration (600 data points)
plot(t, y);
legend ('Susceptible', 'Infectious', 'Recovered');

function xdot = epidemia(x, t)
xdot = zeros(3, 1);
end