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
beta=0.45;
R0 = beta/gamma
time =300;
dt=0.01; %un millesimo di secondo di dt
infected_zero = 1/10e7;
% max number of infected
i_max = 1 - gamma/beta - log(R0*(1-infected_zero))*gamma/beta
s_maxx =  gamma/beta

% For plots
i = 0;

%% Runge Kutta solution
%Heun's method
a1 = 1/2; a2=1/2; p1 =1; q11 =1;
[taxis,xaxis,yaxis,zaxis,infected_zero] = SIR(a1,a2,p1,q11,N,beta,gamma,time,dt,infected_zero);
% 2
[taxis2,xaxis2,yaxis2,zaxis2,infected_zero2] = SIR(a1,a2,p1,q11,N,0.35,gamma,time,dt,infected_zero);
% 3
[taxis3,xaxis3,yaxis3,zaxis3,infected_zero3] = SIR(a1,a2,p1,q11,N,0.3,gamma,time,dt,infected_zero);
% 4
[taxis4,xaxis4,yaxis4,zaxis4,infected_zero4] = SIR(a1,a2,p1,q11,N,0.25,gamma,time,dt,infected_zero);
% 5 
[taxis5,xaxis5,yaxis5,zaxis5,infected_zero5] = SIR(a1,a2,p1,q11,N,0.19,gamma,time,dt,infected_zero);

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

%% Figure tesi
i = i+1;
figure(i)
box on
hold on
% yyaxis left
plot(taxis,xaxis3, 'linewidth',1.1 )
plot(taxis,yaxis3, 'linewidth',1.1)
plot(taxis,zaxis3, 'linewidth',1.1)
xlabel("t[days]");
ylabel("S[t], I[t], R[t]");
yyaxis right
plot(taxis,(beta/gamma).*xaxis3./N, 'linewidth',1.1)
yline(1, "--")
ylim([0 7])
ylabel("R_0[t]");
title("SIR MODEL")
legend('S','I','R', 'R_0[t]', Orientation='horizontal', Location='southoutside')
txt = {['\beta = ' num2str(beta)],['\gamma = ' num2str(gamma)]};
    text(220,2.9,txt)
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', ['sir_con_rt.pdf'])


i = i+1;
figure(i)
box on
hold on
plot(taxis,xaxis,'.','Color', [0 0.4470 0.7410],'HandleVisibility','off' )
plot(taxis,yaxis, 'linewidth',1.1,'Color', [0 0.4470 0.7410],'DisplayName','\beta = 0.45' )
plot(taxis,zaxis, 'LineStyle', "--",'Color', [0 0.4470 0.7410],'HandleVisibility','off')
%2
plot(taxis,xaxis2,'.','Color', [0.8500 0.3250 0.0980],'HandleVisibility','off')
plot(taxis,yaxis2, 'linewidth',1.1,'Color', [0.8500 0.3250 0.0980],'DisplayName','\beta = 0.35')
plot(taxis,zaxis2, 'LineStyle', "--",'Color', [0.8500 0.3250 0.0980],'HandleVisibility','off')
%3
plot(taxis,xaxis3,'.','Color', [0.4940 0.1840 0.5560],'HandleVisibility','off')
plot(taxis,yaxis3, 'linewidth',1.1,'Color', [0.4940 0.1840 0.5560],'DisplayName','\beta = 0.3')
plot(taxis,zaxis3, 'LineStyle', "--",'Color', [0.4940 0.1840 0.5560],'HandleVisibility','off')
%4
plot(taxis,xaxis4,'.','Color', [0.4660 0.6740 0.1880],'HandleVisibility','off')
plot(taxis,yaxis4, 'linewidth',1.1,'Color', [0.4660 0.6740 0.1880],'DisplayName','\beta = 0.25')
plot(taxis,zaxis4, 'LineStyle', "--",'Color',[0.4660 0.6740 0.1880],'HandleVisibility','off')
%5
plot(taxis,xaxis5,'.','Color', [0.6350 0.0780 0.1840],'HandleVisibility','off')
plot(taxis,yaxis5, 'linewidth',1.1,'Color', [0.6350 0.0780 0.1840],'DisplayName','\beta = 0.19')
plot(taxis,zaxis5, 'LineStyle', "--",'Color',[0.6350 0.0780 0.1840],'HandleVisibility','off')
title("SIR MODEL")
legend("AutoUpdate","off",  Orientation='horizontal', Location='southoutside')
txt = {['\gamma = ' num2str(gamma)],['1/\gamma = ' num2str(1/gamma)]};
text(230,0.25,txt)
xlabel("t[days]");
ylabel("S[t], I[t], R[t]");
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
 print(gcf,'-dpdf', ['sir_multipli_beta.pdf'])

%% Runge Kutta second order solution of SIR      
function [taxis,xaxis,yaxis,zaxis,i] = SIR(a1,a2,p1,q11,N,beta,gamma,time,dt,i)
    
    x = N-i; % susceptible
    y = i; % infected
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


