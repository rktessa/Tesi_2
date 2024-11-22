close all;
clear; 
clc;
%% Readme!
% To run this code Optimization toolbox is necessary

%%  SIRS simulation model 
% with beta dependent on I, simulating the effect of 
% government policies

% Simulation parameters
N = 1; % Population normalised
beta = 0.40;
gamma = 1/9;
delta = 1/90;
Ro = beta/gamma
time =400;
dt=0.01; %un millesimo di secondo di dt

% For plots
i = 0;

%%  Trying different Runge Kutta coefficients
%Heun's method
a1 = 1/2; a2=1/2; p1 =1; q11 =1;
[taxis,xaxis,yaxis, zaxis,beta_vec,i] = SIRS(a1,a2,p1,q11,N,beta,gamma,delta,time,dt,i);

%% Comparison between Runge Kutta and ODE
i = i+1;
figure(i)
hold on
plot(taxis,xaxis, 'linewidth',3.0, 'color', [0.6350 0.0780 0.1840] )
plot(taxis,yaxis, 'linewidth',3.0, 'color', [0 0.4470 0.7410] )
plot(taxis,zaxis, 'linewidth',3.0, 'color', [0.4660 0.6740 0.1880] )
title("SIRS MODEL")
xlabel("Time [t]")
ylabel("S[t], I[t], R[t]")
legend('S', 'I', 'R', Location='southoutside', Orientation='horizontal')
txt = {['R_0= ' num2str(round(beta/(gamma+delta),3))],['\beta= ' num2str(round(beta,3))],['\gamma= ' num2str(round(gamma,3))],['\delta= ' num2str(round(delta,3))],...
    ['S_\infty= ' num2str(round(xaxis(end),3))], ['I_\infty= ' num2str(round(yaxis(end),3))], ['R_\infty= ' num2str(round(zaxis(end),3))]};
dim = [.93 .85 .1 .1];
annotation('textbox',dim, ...
    'String',txt,'EdgeColor','none')
hold off
fontsize(24,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 21 15]);
set(gcf, 'PaperSize', [25 15]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', "sirs_figure")
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
   
    while t < time
    
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

