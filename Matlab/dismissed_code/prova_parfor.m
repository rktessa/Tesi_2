clc; 
clear all;
close all;

%user parameters
N = 45742000; %total population
I0 = 10; %initial infected population
r0 = 12.9 %reproductive value
R0 = 0;%initial recovered population
i_period = 9; %duration of infectious period
i_period = 5; %duration of infectious period
tspan = [1, 70]; %length of time to run simulation over
%rate constant calculations
mu = 1 / i_period %recovery rate
N = 10000;
beta=0.5/N;
%beta = r0 * mu / N %infection rate
S0 = N - I0 - R0 %initial susceptible population
N0 = I0 + R0 + S0; %total population
%---------------------------------------------------
%feeding parameters into function
pars = [beta, mu, N, r0];
y0 = [S0 I0 R0];

%Running SIR model function 
%using the ode45 function to perform intergration 
[t,y] = ode45(@sir_rhs, tspan, y0, [], pars);
figure()
hold on
plot(t,y(:,1), 'b');
plot(t,y(:,2), 'r');
plot(t,y(:,3), 'g');
xlim(tspan);
ylabel('Population (n)');
xlabel('Time (days)');
legend('Susceptible','Infected','Recovered');


function f = sir_rhs(t,y,pars)
f = zeros(3,1);
f(1) = -pars(1)*y(1)*y(2);
f(2) = pars(1)*y(1)*y(2) - pars(2)*y(2);
f(3) = pars(2) * y(2);
end
