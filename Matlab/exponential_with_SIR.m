close all;
clear all;
clc;
%% SIR and exponenial fitting
% try to calculate the value of x that  exprime the exponential coefficient
% of a R curve in a simple SIR model

% Simulation parameters
N = 1;
time =1000;
dt=0.001; %un millesimo di secondo di dt

% Create object to call function from the class
Obj = functionsContainer;

res1 = Obj.func1(2);
% For plots
i = 0;

beta_f= 0.5/N;
gamma_f=0.37;
[i] = Obj.SIR_RK_figure(N,beta_f,gamma_f,time,dt,i);

