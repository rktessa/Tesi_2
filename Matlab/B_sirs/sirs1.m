close all;
clear; 
clc;
%% Readme!
% To run this code Optimization toolbox is necessary
% Test of solve a SIRS model with a Runge Kutta second order and ODE45
% matlab solve (that is a built in six-stage, fifth-order, Runge-Kutta
% method). Comparison of the result and analysis of SIR behavior with
% different beta and gamma values. 
%%  SIRS simulation model 
% with beta dependent on I, simulating the effect of 
% government policies

% Simulation parameters
N = 1;

beta = 0.40;
gamma = 1/9;
delta = 1/90;
Ro = beta/gamma
time =1500;
dt=0.01; %un millesimo di secondo di dt

% For plots
i = 0;

%%  Trying different Runge Kutta coefficients
%Heun's method
a1 = 1/2; a2=1/2; p1 =1; q11 =1;
[taxis,xaxis,yaxis, zaxis,beta_vec,i] = SIRS(a1,a2,p1,q11,N,beta,gamma, delta,time,dt,i);


%% ODE solution
[taxisODE,xaxisODE,yaxisODE, zaxisODE,i] = SIRS_ODE(N,beta,gamma, delta,time,dt,i);

%% Comparison between Runge Kutta and ODE
i = i+1;
figure(i)
hold on
plot(taxis,xaxis, 'r', 'linewidth',2.0 )
plot(taxis,yaxis, 'b', 'linewidth',2.0 )
plot(taxis,zaxis, 'm', 'linewidth',2.0 )

plot(taxisODE,xaxisODE, 'color', [0.6350 0.0780 0.1840], 'linewidth',1.0 )
plot(taxisODE,yaxisODE, 'color', [0 0.4470 0.7410], 'linewidth',1.0 )
plot(taxisODE,zaxisODE, 'color', [0.4660 0.6740 0.1880], 'linewidth',1.0 )
title("SIR MODEL Runge Kutta vs ODE15s")
legend('S_{RK}','I_{RK}','R_{RK}', 'S_{ode}','I_{ode}','R_{ode}')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)],['delta: ' num2str(delta)]};
    text(800,0.75, txt)

%%
i = i+1;
figure(i)
hold on
plot(taxis,yaxis, 'b', 'linewidth',2.0 )
plot(taxisODE, yaxisODE, 'r', LineWidth=2.0)
title("SIRS MODEL Runge Kutta")
legend('I_{RK}', 'I_{ode}')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)],['delta: ' num2str(delta)]};
    text(600,0.08,txt)    
%%
i = i+1;
figure(i)
hold on
plot(taxis,xaxis, 'r', 'linewidth',2.0 )
plot(taxis,yaxis, 'b', 'linewidth',2.0 )
plot(taxis,zaxis, 'm', 'linewidth',2.0 )
title("SIR MODEL Runge Kutta with beta varying" )
legend('S_{RK}','I_{RK}','R_{RK}')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)],['delta: ' num2str(delta)]};
    text(800,0.75, txt)

%% Evoluzione di beta nel corso della pandemia
i = i+1;
figure(i)
plot(taxis, beta_vec./gamma, LineWidth=1.5)
title("Evolution of Ro during the course of the epidemic")


%% Calculate the percentage difference
% perc_diff_S = zeros(length(xaxis),1);
% perc_diff_I = zeros(length(xaxis),1);
% 
% 
% for h =1 :length(xaxis)
%     perc_diff_S(h) = percentage_difference(xaxis(h),xaxisODE(h));
%     perc_diff_I(h) = percentage_difference(yaxis(h),yaxisODE(h));
%     perc_diff_R(h) = percentage_difference(zaxis(h),zaxisODE(h));
% 
% end
% 
% %Plot of the difference of the 2 solutions
% i = i+1;
% figure(i)
% hold on
% plot(taxis,perc_diff_S, 'color', [0.6350 0.0780 0.1840], 'linewidth',1.5 )
% plot(taxis,perc_diff_I, 'color', [0.4660 0.6740 0.1880], 'linewidth',1.5 )
% plot(taxis,perc_diff_R, 'color', [0 0.4470 0.7410], 'linewidth',1.5 )
% title("Percentage Difference of SIR model between Runge Kutta 2^{nd} order and ODE15s")
% legend('S_{dif}','I_{dif}','R_{dif}')   


%% Runge Kutta second order solution of SIRS  with beta state varying     
function [taxis,xaxis,yaxis, zaxis,beta_vec, i] = SIRS(a1,a2,p1,q11,N,beta,gamma, delta,time,dt,i)
    
    x = N-100/60e6; % susceptible
    y = 100/60e6;  % infected
    z = 0; % recovered

    t = 0;
    cnt=0;
    cnt2=0; %contatore per cambio policy
    %Array creation and inititialization
    taxis=[]; taxis(1) =0; 
    xaxis=[]; xaxis(1) = x;
    yaxis=[]; yaxis(1) = y;
    zaxis=[]; zaxis(1) = z;
    beta_vec = []; beta_vec(1) = beta;
    beta_0 = gamma;
    while t < time
        % Control the I value
        if mod(cnt2,30/dt) == 0 % ogni 30 giorni posso cambiare policy 
            if y <= 0.0005 
                beta= beta_0*2.38*N;
            elseif y > 0.0005 && y < 0.0025
                beta= beta_0*1.3*N;
            else 
                beta= beta_0*0.6*N;
            end
        end
        cnt2 = cnt2 +1;
        if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
            taxis = cat(2,taxis,t);
            xaxis = cat(2,xaxis,x);
            yaxis = cat(2,yaxis,y);
            zaxis = cat(2,zaxis,z);
            beta_vec = cat(2,beta_vec,beta);
            
        end
        % step 1
        kx1 = - beta*x*y + delta*z;
        ky1 = beta*x*y - gamma*y;
        kz1 = + gamma*y - delta*z;
        % step 2
        t2 = t+p1*dt;
        x2 = x + q11*kx1*dt;
        y2 = y + q11*ky1*dt;
        z2 = z + q11*kz1*dt;

        kx2 = - beta*x2*y2 + delta*z2;
        ky2 = beta*x2*y2 - gamma*y2;
        kz2 = gamma*y2 - delta*z2;
        % update
        x = x + (a1*kx1+a2*kx2)*dt;
        y = y + (a1*ky1+a2*ky2)*dt;
        z = z + (a1*kz1+a2*kz2)*dt;
       
        t = t + dt;
        cnt = cnt + 1;


    end
    
end

%% Solve the problem with Matlab ODE functions
function [taxis,xaxis,yaxis, zaxis,i] = SIRS_ODE(N,beta,gamma,delta, time,dt,i)
    
    s0 = N-200/60e6; % susceptible
    i0 = 200/60e6;  % infected
    r0 = 0; 
    tspan = 0:dt:time;
    y0 = [s0,i0, r0];
    pars = [beta, gamma, delta];
    
    %ode45 solver input: function to solve, t_span, initial condition 

    [t,y] = ode45(@sirs_rhs, tspan, y0, [], pars);
    
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
function f = sirs_rhs(t,y,pars)
f = zeros(3,1);
f(1) = -pars(1)*y(1)*y(2)+ pars(3)*y(3);
f(2) = pars(1)*y(1)*y(2) - pars(2)*y(2);
f(3) = pars(2)*y(2) - pars(3)*y(3);
end


% function f = sirs_rhs(t,y,pars)
% f = zeros(3,1);
% f(1) = -pars(1)*y(1)*y(2);
% f(2) = pars(1)*y(1)*y(2) - pars(2)*y(2);
% f(3) = pars(2)*y(2);
% end
% Function that calculate the percentage difference between two values
function [perc_diff] = percentage_difference(v1,v2)
perc_diff = abs(v1-v2)/((v1+v2)/2) *100;
end

