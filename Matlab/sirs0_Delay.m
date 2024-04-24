close all;
clc;
%% Readme!
% To run this code Optimization toolbox is necessary
% Test of solve a SIRS model with a Runge Kutta second order and ODE45
% matlab solve (that is a built in six-stage, fifth-order, Runge-Kutta
% method). Comparison of the result and analysis of SIR behavior with
% different beta and gamma values. 
%% Load covid data
M = readmatrix('covid_italia.csv'); % read a csv file
total_infected = M(:,7);
% consider firs 15 days
t0t_15d = total_infected(1:15)./60e6; %normalized


%%  SIRS simulation model with a delay
% try replicating the experiment of the article
% call the class functions container
obj = functionsContainer;
% Simulation parameters
% Uso gli stessi dell'articolo

gamma=0.119;
beta= 0.392; %Così dovrei avere R0 = 3.2
delta = 1/205;  % 7 mesi, nell'arco temporale di 15 giorni non ha un impatto
Ro = beta/gamma
time =24;
D = 8.6; %delay between case reported and number of infected
I0 = 19/60e6;

%add delay to vector of tot
tot2 = [ zeros(round(D),1); t0t_15d ];

% For plots
i = 0;


[taxis,xaxis,yaxis, zaxis] = obj.SIRS(beta,gamma,delta,time);


% i = i+1;
% figure(i)
% hold on
% plot(taxis,xaxis, 'r', 'linewidth',2.0 )
% plot(taxis,yaxis, 'b', 'linewidth',2.0 )
% plot(taxis,zaxis, 'm', 'linewidth',2.0 )
% title("SIR MODEL Runge Kutta")
% legend('S_{RK}','I_{RK}','R_{RK}')
% txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)],['delta: ' num2str(delta)]};
% text(20,0.75, txt)

i = i+1;
figure(i)
hold on

plot(taxis,yaxis,  'linewidth',2.0 )
plot(taxis,tot2, 'o' )
title("SIR MODEL Runge Kutta wit delay of 8 days")
legend('I_{RK}')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)],['delta: ' num2str(delta)]};
text(20,tot2(18), txt)

