%% Try to read Covid italia .csv file
clc;
clear all;
close all; 
%%
M = readmatrix('covid_italia.csv'); % read a csv file
number_infection = M(:,7);

i = i+1;
figure(i)
plot(M(:,7)/60e6, LineWidth=2)
legend(' Totali positivi normalizzati covid italia')

% Valori dei coefficienti per alfa
infected = 0:1:100;
p = 1;
q = 1;
beta = 1/20000;
beta_1_1 = 1.1*beta;
beta_0_9 = 0.9*beta;
%beta = 1/10;
%beta = 1,
alfa0 = 1;

for l = 1:length(number_infection)
    alfa_2(l) = alpha_function(p,q,alfa0,beta,number_infection(l));
    alfa_21(l) = alpha_function(p,q,alfa0,beta_1_1,number_infection(l));
    alfa_29(l) = alpha_function(p,q,alfa0,beta_0_9,number_infection(l));
end
% Normalize the value of alfa
max_alfa2 = max(alfa_2);
alfa_2 = alfa_2./max_alfa2;

max_alfa21 = max(alfa_21);
alfa_21 = alfa_21./max_alfa21;

max_alfa29 = max(alfa_29);
alfa_29 = alfa_29./max_alfa29;

i = i+1;
figure(i)
hold on
plot(alfa_2,LineWidth=2)
plot(alfa_21,LineWidth=2)
plot(alfa_29,LineWidth=2)
title("\alpha during covid in taly")
legend("\alpha", "\alpha +10 \%", "\alpha -10\%")
xlabel("iteration")

for l = 1:length(infected)
    alfa(l) = alpha_function(1,1,alfa0,beta,infected(l));  
    alfa3(l) = alpha_function(1,2,alfa0,beta,infected(l)); 
    alfa4(l) = alpha_function(2,2,alfa0,beta,infected(l)); 
end


i = i+1;
figure(i)
hold on
plot(infected,alfa,LineWidth=2)
plot(infected,alfa3,LineWidth=2)
plot(infected,alfa4,LineWidth=2)
legend("Andamento di \alpha in funzione di I")
xlabel("valore di I")
%alpha(I) risk perception che satura
function [alpha] = alpha_function(p,q, alpha_0,beta, i_rate)
alpha = (alpha_0*i_rate.^p)/(1+beta*i_rate.^q);
end