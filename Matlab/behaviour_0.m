close all;
clear;
clc;
%%  Opinion simulation model

% Simulation parameters
                                                                                            
k1 =  1/5; % from Ca to Co
k2 = 1/2; % from Ca to Ag

lambda_1 = 1/50; %fatigue to mantain Co behaviour
lambda_2 = 1/10; %fatigue to mantain Ag behaviour
time = 2000;
R1 = k1/lambda_1
R2 = k2/lambda_2
% For plots
i = 0;

%% Runge Kutta solution
%using Heun's method 
a1 = 1/2; a2=1/2; p1 =1; q11 =1;
[taxisRK,xaxisRK,yaxisRK,zaxisRK,awareness] = Behaviour_RK(k1,k2,lambda_1,lambda_2,time);

% [rCo,lt,ut,ll,ul] = risetime(xaxisRK,taxisRK,PercentReferenceLevels=[10 90]);

    i = i+1;
    figure(i)
    hold on
    plot(taxisRK,xaxisRK,  'linewidth',1.0 )
    plot(taxisRK,yaxisRK,  'linewidth',1.0 )
    plot(taxisRK,zaxisRK, 'linewidth',1.0 )
    % plot(lt,ll, 'o', 'color', [0.4940 0.1840 0.5560])
    % plot(ut,ul, 'o', 'color', [0.4940 0.1840 0.5560])
    % yline([ll ul],'--',{'low level','upper level'})
    title("Behaviour model Runge Kutta")
    legend('Careless','Compliant','Against')
    txt = {['k1: ' num2str(k1)],['k2: ' num2str(k2)],['lambda_1: ' num2str(lambda_1)],['lambda_2: ' num2str(lambda_2)]};
    text(2000,0.75,txt)

%% Varying alpha and alpha_prime parameters
% The Initial condition are the same in all the case, but alpha and
% alpha_prime coefficients vary in an interval [0.01,1] 

% d and step used for computations
% d = 3;
% step = 0.005;
% [Co_max, A_max,max_gr_Co, rCo_m, rA_m] = multiple_Opinion(d,step);    

%% Plot the results

% alpha_pl = (0.005:step:(step*d))/N;
% alpha_prime_pl = (0.005:step:(step*d))/N;
% [xx, yy] = meshgrid(alpha_pl,alpha_prime_pl);
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,Co_max)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Max value of Compliants for alpha and alpha_{prime}')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,A_max)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Max value of Anti for alpha and alpha_{prime}')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,rCo_m)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Rise time of Compliant for alpha and alpha_{prime}')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,rA_m)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Rise time of Anti for alpha and alpha_{prime}')

%% Runge Kutta second order solution of Opinion     
function [taxis,xaxis,yaxis,zaxis,awareness] = Behaviour_RK(k1,k2,lambda_1,lambda_2,time)
    
    a1 = 1/2; a2=1/2; p1 =1; q11 =1;
    x = 1-400/60e6; % Careless
    y = 200/60e6;  % Compliant
    z = 200/60e6; % Against
    dt = 0.01;
    t = 0;
    cnt=0;
    awareness = 0;
    %Array creation and inititialization
    taxis=[]; taxis(1) = 0; 
    xaxis=[]; xaxis(1) = x;
    yaxis=[]; yaxis(1) = y;
    zaxis=[]; zaxis(1) = z;
    while t < time
        if mod(cnt,100) == 0 && cnt ~=0 %every 100 iterations I save the result
            taxis = cat(2,taxis,t);
            xaxis = cat(2,xaxis,x);
            yaxis = cat(2,yaxis,y);
            zaxis = cat(2,zaxis,z);
        end
        % step 1
        kx1 =  -k1*x*y - k2*x*z + lambda_1*y + lambda_2*z;
        ky1 =   k1*x*y - lambda_1*y;
        kz1 =   k2*x*z - lambda_2*z;
        % step 2
        t2 = t+p1*dt;
        x2 = x + q11*kx1*dt;
        y2 = y + q11*ky1*dt;
        z2 = z + q11*kz1*dt;
        kx2 =  -k1*x2*y2 - k2*x2*z2 + lambda_1*y2 + lambda_2*z2;
        ky2 =   k1*x2*y2 - lambda_1*y2;
        kz2 =   k2*x2*z2 - lambda_2*z2;
        % update
        x = x + (a1*kx1+a2*kx2)*dt;
        y = y + (a1*ky1+a2*ky2)*dt;
        z = z + (a1*kz1+a2*kz2)*dt;
        t = t + dt;
        cnt = cnt + 1;
    end
end