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

% Simulate a SIR
i = 0;
beta_f= 0.5/N;
gamma_f=0.37;
[i,t_vec, S1, I1,R1 ] = Obj.SIR_RK_figure(N,beta_f,gamma_f,time,dt,i);

% Find the max value of Infected
[max_I1,index_I1] = max(I1);
% Gradient of I
gradient_I1 = gradient(I1); %derivata prima
% Ci sono due flessi
[max_grad_I1,inx_flesso_1] = max(gradient_I1); %primo
[min_grad_I1,inx_flesso_2] = min(gradient_I1); %primo

% Plot del gradiente di I
i = i+1;
figure(i)
hold on
plot(gradient_I1, 'b', 'linewidth',1.0 )
title("gradient of I function")

i = i+1;
figure(i)
hold on
plot(t_vec,I1, 'b', 'linewidth',1.0 )
plot(t_vec(inx_flesso_1),I1(inx_flesso_1), 'r*' )
plot(t_vec(inx_flesso_2),I1(inx_flesso_2), 'r*' )
title("I function")


% Plot in Log scale

i = i+1;
figure(i)
hold on
plot(log(S1), 'r', 'linewidth',1.0 )
plot(log(I1), 'b', 'linewidth',1.0 )
plot(log(R1), 'g', 'linewidth',1.0 )
title("SIR MODEL")
legend('S','I','R')
txt = {['beta: ' num2str(beta_f)],['gamma: ' num2str(gamma_f)]};
text(9000, -20,txt)

%% Using LS to interpolate
% switch from max time to inflection time to entire sample
l = inx_flesso_1;
% l = index_I1;
% l = length(I1);


I1_log = log(I1(2:l)); %calulate the interpolation to the first inflection point
% I1_log = I1(2:l);

t_pol = t_vec(2:l); %time vector used for calculations
% The model I want to fit is y = exp(ax+b), y = exp(cx^2+ax+b)


[a_I, b_I] = LS_linear( t_pol, I1_log)
f_inter = t_pol*a_I + b_I;
%% Trying with matlab polinomial fitting 
p = polyfit(t_pol, I1_log,2);

%% Plot the Recovered on increasing phase

i = i+1;
figure(i)
hold on
%plot(t_vec(1:1000),R1(1:1000), 'g', 'linewidth',1.0 )
plot(t_pol, I1_log,  'g*')
plot(t_pol,polyval(p,t_pol), 'r', 'linewidth',1.0)
plot(t_pol, f_inter, 'b', 'linewidth',1.0 )
title("I MODEL")
legend('I', 'polyfit','I fitted')

%% Function 
function [a, b] = LS_linear(x, y)
    %with x and y data points 
    m_x = mean(x); % mean
    m_y = mean(y);
    v_x = var(x);  % variance
    v_y = var(y);
    c_xy = cov(x,y); % covariance
    rms_x_2 = rms(x)^2; % rms square 
    rms_y_2 = rms(y)^2;
    m_xy = mean(x.*y);

    M = [ 1, m_x;
          m_x, rms_x_2];
    v = [m_y; m_xy];

    sol = inv(M)*v;

    b = sol(1);
    a = sol(2);
end