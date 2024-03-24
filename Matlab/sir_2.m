close all;
clear all;
clc;
%% Readme!
% To run this code Optimization toolbox is necessary
% Test of solve a SIR model with a Runge Kutta second order and ODE45
% matlab solve (that is a built in six-stage, fifth-order, Runge-Kutta
% method). Comparison of the result and analysis of SIR behavior with
% different beta and gamma values. 

%%  SIR simulation model


% Simulation parameters
N = 1;
gamma=1/8.9;
beta=gamma*2.01/N;
R0 = beta/gamma
time =1000;
dt=0.01; %un millesimo di secondo di dt

% For plots
i = 0;

%% Runge Kutta solution
%Heun's method
a1 = 1/2; a2=1/2; p1 =1; q11 =1;
[taxis,xaxis,yaxis,zaxis,i] = SIR(a1,a2,p1,q11,N,beta,gamma,time,dt,i);



%% ODE solution
[taxisODE,xaxisODE,yaxisODE,zaxisODE,i] = SIR_ODE(N,beta,gamma,time,dt,i);

%% Comparison between Runge Kutta and ODE
i = i+1;
figure(i)
hold on
plot(taxis,xaxis, 'r', 'linewidth',2.0 )
plot(taxis,yaxis, 'g', 'linewidth',2.0 )
plot(taxis,zaxis, 'b', 'linewidth',2.0 )
plot(taxisODE,xaxisODE, 'color', [0.6350 0.0780 0.1840], 'linewidth',2.0 )
plot(taxisODE,yaxisODE, 'color', [0.4660 0.6740 0.1880], 'linewidth',2.0 )
plot(taxisODE,zaxisODE, 'color', [0 0.4470 0.7410], 'linewidth',2.0 )
title("SIR MODEL Runge Kutta vs ODE15")
legend('S_{RK}','I_{RK}','R_{RK}', 'S_{ode}','I_{ode}','R_{ode}')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)]};
    text(600,0.5,txt)


%Also in this case the solutions are very close
% Calculate the percentage difference
perc_diff_S = zeros(length(xaxis),1);
perc_diff_I = zeros(length(xaxis),1);
perc_diff_R = zeros(length(xaxis),1);

for h =1 :length(xaxis)
    perc_diff_S(h) = percentage_difference(xaxis(h),xaxisODE(h));
    perc_diff_I(h) = percentage_difference(yaxis(h),yaxisODE(h));
    perc_diff_R(h) = percentage_difference(zaxis(h),zaxisODE(h));
    
end

%Plot of the difference of the 2 solutions
i = i+1;
figure(i)
hold on
plot(taxis,perc_diff_S, 'color', [0.6350 0.0780 0.1840], 'linewidth',2.0 )
plot(taxis,perc_diff_I, 'color', [0.4660 0.6740 0.1880], 'linewidth',2.0 )
plot(taxis,perc_diff_R, 'color', [0 0.4470 0.7410], 'linewidth',2.0 )
title("Percentage Difference of SIR model between Runge Kutta 2^{nd} order and ODE45")
legend('S_{dif}','I_{dif}','R_{dif}')   


%% Runge Kutta second order solution of SIR      
function [taxis,xaxis,yaxis,zaxis,i] = SIR(a1,a2,p1,q11,N,beta,gamma,time,dt,i)
    
    x = N-N/100000; % susceptible
    y = N/100000; % infected
    z = 0; % recovered
    t = 0;
    cnt=0;
    %Array creation and inititialization
    taxis=[]; taxis(1) =0; 
    xaxis=[]; xaxis(1) = x;
    yaxis=[]; yaxis(1) = y;
    zaxis=[]; zaxis(1) = z;
    while t < time
        if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
            taxis = cat(2,taxis,t);
            xaxis = cat(2,xaxis,x);
            yaxis = cat(2,yaxis,y);
            zaxis = cat(2,zaxis,z);
        end
        % step 1
        kx1 = - beta*x*y;
        ky1 = beta*x*y - gamma*y;
        % step 2
        t2 = t+p1*dt;
        x2 = x + q11*kx1*dt;
        y2 = y + q11*ky1*dt;
        kx2 = - beta*x2*y2;
        ky2 = beta*x2*y2 - gamma*y2;
        % update
        x = x + (a1*kx1+a2*kx2)*dt;
        y = y + (a1*ky1+a2*ky2)*dt;
        z = N - x - y;
        t = t + dt;
        cnt = cnt + 1;


    end
    
    
%     i = i+1;
%     figure(i)
%     hold on
%     plot(taxis,xaxis, 'r', 'linewidth',1.0 )
%     plot(taxis,yaxis, 'b', 'linewidth',1.0 )
%     plot(taxis,zaxis, 'g', 'linewidth',1.0 )
%     title("SIR MODEL")
%     legend('S','I','R')

end

%% Solve the problem with Matlab ODE functions
function [taxis,xaxis,yaxis,zaxis,i] = SIR_ODE(N,beta,gamma,time,dt,i)
    
    s0= N-N/100000; % susceptible
    i0 = N/100000; % infected
    r0 = 0; % recovered
    tspan = 0:dt:time;
    y0 = [s0,i0,r0];
    pars = [beta, gamma];
    
    %ode45 solver input: function to solve, t_span, initial condition 

    [t,y] = ode15s(@sir_rhs, tspan, y0, [], pars);
    
    cnt=0;
    p  =0;
    %Array creation and inititialization
    taxis=[]; taxis(1) =t(1); 
    xaxis=[]; xaxis(1) = y(1,1);
    yaxis=[]; yaxis(1) = y(1,2);
    zaxis=[]; zaxis(1) = y(1,3);
    
    
    while p < length(t) %Scorro tutta la soluzione
        if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
            
            taxis = cat(2,taxis,t(cnt));
            xaxis = cat(2,xaxis,y(cnt,1));
            yaxis = cat(2,yaxis,y(cnt,2));
            zaxis = cat(2,zaxis,y(cnt,3));
        end
    cnt = cnt + 1;  
    p = p+1;
    end
   

end

%function used to the ODE45 solver to compute correctly the system
%evolution
function f = sir_rhs(t,y,pars)
f = zeros(3,1);
f(1) = -pars(1)*y(1)*y(2);
f(2) = pars(1)*y(1)*y(2) - pars(2)*y(2);
f(3) = pars(2) * y(2);
end




% Function that calculate the percentage difference between two values
function [perc_diff] = percentage_difference(v1,v2)
perc_diff = abs(v1-v2)/((v1+v2)/2) *100;
end


