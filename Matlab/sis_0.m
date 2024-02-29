close all;
clear all;
clc;
%% Readme!
% To run this code Optimization toolbox is necessary
% Test of solve a SIR model with a Runge Kutta second order and ODE45
% matlab solve (that is a built in six-stage, fifth-order, Runge-Kutta
% method). Comparison of the result and analysis of SIR behavior with
% different beta and gamma values. 

%%  SIS simulation model

% Simulation parameters
N = 1;

gamma=1/3;
beta= gamma*2*N; %Così dovrei avere R0 = 2
Ro = beta/gamma
time =1000;
dt=1; %un millesimo di secondo di dt

% For plots
i = 0;

%%  Trying different Runge Kutta coefficients
%Heun's method
a1 = 1/2; a2=1/2; p1 =1; q11 =1;
[taxis,xaxis,yaxis,i] = SIS(a1,a2,p1,q11,N,beta,gamma,time,dt,i);


%% ODE solution
[taxisODE,xaxisODE,yaxisODE,i] = SIS_ODE(N,beta,gamma,time,dt,i);

%% Comparison between Runge Kutta and ODE
i = i+1;
figure(i)
hold on
plot(taxis,xaxis, ':r', 'linewidth',1.0 )
plot(taxis,yaxis, ':g', 'linewidth',1.0 )
plot(taxisODE,xaxisODE, 'color', [0.6350 0.0780 0.1840], 'linewidth',1.0 )
plot(taxisODE,yaxisODE, 'color', [0.4660 0.6740 0.1880], 'linewidth',1.0 )
title("SIR MODEL Runge Kutta vs ODE45")
legend('S','I', 'S_{ode}','I_{ode}','R_{ode}')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)]};
    text(50,6000,txt)


%Also in this case the solutions are very close
% Calculate the percentage difference
perc_diff_S = zeros(length(xaxis),1);
perc_diff_I = zeros(length(xaxis),1);


for h =1 :length(xaxis)
    perc_diff_S(h) = percentage_difference(xaxis(h),xaxisODE(h));
    perc_diff_I(h) = percentage_difference(yaxis(h),yaxisODE(h));
    
    
end

%Plot of the difference of the 2 solutions
% i = i+1;
% figure(i)
% hold on
% plot(taxis,perc_diff_S, 'color', [0.6350 0.0780 0.1840], 'linewidth',1.0 )
% plot(taxis,perc_diff_I, 'color', [0.4660 0.6740 0.1880], 'linewidth',1.0 )
% plot(taxis,perc_diff_R, 'color', [0 0.4470 0.7410], 'linewidth',1.0 )
% title("Percentage Difference of SIR model between Runge Kutta 2^{nd} order and ODE45")
% legend('S_{dif}','I_{dif}','R_{dif}')   

%% Analisys of epidemic spread varying beta and gamma

% The Initial condition are the same in all the case, but beta and gamma
% coefficients vary in an interval [0,1] 

% d = 10;
% step = 0.01;
% [I_max, T_imax,R0] = multipleSIR(d,step);    
% 
% %% Plot the results
% beta_pl = (0.01:step:(d*step))/N;
% gamma_pl = 0.01:step:(d*step);
% [xx, yy] = meshgrid(gamma_pl,beta_pl);
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,I_max)
% colorbar;
% colormap default;
% xlabel('gamma') 
% ylabel('beta') 
% title('Max value of I for beta and gamma')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,T_imax)
% colorbar;
% colormap default;
% xlabel('gamma') 
% ylabel('beta') 
% title('Time of I peak for beta and gamma')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,R0)
% colorbar;
% colormap default;
% xlabel('gamma') 
% ylabel('beta') 
% title('R_0 for beta and gamma')
% 

%% Runge Kutta second order solution of SIS      
function [taxis,xaxis,yaxis,i] = SIS(a1,a2,p1,q11,N,beta,gamma,time,dt,i)
    
    x = N-N/100000; % susceptible
    y = N/100000;  % infected
    
    t = 0;
    cnt=0;
    %Array creation and inititialization
    taxis=[]; taxis(1) =0; 
    xaxis=[]; xaxis(1) = x;
    yaxis=[]; yaxis(1) = y;
    
    while t < time
        if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
            taxis = cat(2,taxis,t);
            xaxis = cat(2,xaxis,x);
            yaxis = cat(2,yaxis,y);
            
        end
        % step 1
        kx1 = - beta*x*y + gamma*y;
        ky1 = beta*x*y - gamma*y;
        % step 2
        t2 = t+p1*dt;
        x2 = x + q11*kx1*dt;
        y2 = y + q11*ky1*dt;
        kx2 = - beta*x2*y2 + gamma*y2;
        ky2 = beta*x2*y2 - gamma*y2;
        % update
        x = x + (a1*kx1+a2*kx2)*dt;
        y = y + (a1*ky1+a2*ky2)*dt;
       
        t = t + dt;
        cnt = cnt + 1;


    end
    
    
%     i = i+1;
%     figure(i)
%     hold on
%     plot(taxis,xaxis, 'r', 'linewidth',1.0 )
%     plot(taxis,yaxis, 'b', 'linewidth',1.0 )
%     title("SIR MODEL")
%     legend('S','I','S')

end

%% Solve the problem with Matlab ODE functions
function [taxis,xaxis,yaxis,i] = SIS_ODE(N,beta,gamma,time,dt,i)
    
    s0 = N-N/100000; % susceptible
    i0 = N/100000;  % infected
    tspan = 0:dt:time;
    y0 = [s0,i0];
    pars = [beta, gamma];
    
    %ode45 solver input: function to solve, t_span, initial condition 

    [t,y] = ode45(@sis_rhs, tspan, y0, [], pars);
    
    cnt=0;
    p  =0;
    %Array creation and inititialization
    taxis=[]; taxis(1) =t(1); 
    xaxis=[]; xaxis(1) = y0(1,1);
    yaxis=[]; yaxis(1) = y0(1,2);
    
    
    
    while p < length(t) %Scorro tutta la soluzione
        if mod(cnt,10) == 0 && cnt ~=0 %every 100 millisecond I save the result
            
            taxis = cat(2,taxis,t(cnt));
            xaxis = cat(2,xaxis,y(cnt,1));
            yaxis = cat(2,yaxis,y(cnt,2));
           
        end
    cnt = cnt + 1;  
    p = p+1;
    end
    
%      i = i+1;
%     figure(i)
%     hold on
%     plot(taxis,xaxis, 'r', 'linewidth',1.0 )
%     plot(taxis,yaxis, 'b', 'linewidth',1.0 )
%     title("SIS MODEL")
%     legend('S','I')

end

%function used to the ODE45 solver to compute correctly the system
%evolution
function f = sis_rhs(t,y,pars)
f = zeros(2,1);
f(1) = -pars(1)*y(1)*y(2)+ pars(2)*y(2);
f(2) = pars(1)*y(1)*y(2) - pars(2)*y(2);
end


function [I_max, T_imax,R0] = multipleSIR(d,step)
    I_max = zeros(d,d);
    T_imax = zeros(d,d);
    R0 = zeros(d,d);
    N = 10000;
    time =70;
    dt=0.001; %un millesimo di secondo di dt
    a1 = 1/2; a2=1/2; p1 =1; q11 =1;
    
    for i = 1:d
        beta = i*step/N;
        for j = 1:d
            gamma = j*step;
            
            %Esecuzione della trasformazione, variando beta e gamma
            [taxis,xaxis,yaxis,zaxis,i] = SIR(a1,a2,p1,q11,N,beta,gamma,time,dt,i);
            
            %Compute max value of infected
            [max_i, index] = max(yaxis);
            time_i_max = taxis(index);
            I_max(i,j) = max_i;
            T_imax(i,j) = time_i_max;
            R0(i,j) = (beta/gamma)*N;
            
        end
    end
end

% Function that calculate the percentage difference between two values
function [perc_diff] = percentage_difference(v1,v2)
perc_diff = abs(v1-v2)/((v1+v2)/2) *100;
end


function [i] = SIR_RK_figure(N,beta,gamma,time,dt,i)
    a1 = 1/2; a2=1/2; p1 =1; q11 =1;
    x = N-10; % susceptible
    y = 10;  % infected
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
    
    
    i = i+1;
    figure(i)
    hold on
    plot(taxis,xaxis, 'r', 'linewidth',1.0 )
    plot(taxis,yaxis, 'b', 'linewidth',1.0 )
    plot(taxis,zaxis, 'g', 'linewidth',1.0 )
    title("SIR MODEL")
    legend('S','I','R')
    txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)]};
    text(50,6000,txt)
end

