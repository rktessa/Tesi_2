close all;
clear all;
clc;

%%  Opinion simulation model

% Simulation parameters
N = 10000;
alpha=0.27/N;
alpha_prime=0.28/N;

alpha=5/N;
alpha_prime=3/N;
lambda = 2000/N;

kappa = 7.9/N;
time =2500;
dt=0.001; %un centesimo/millesimo di secondo di dt

% For plots
i = 0;

%% Runge Kutta solution
%using Heun's method 
a1 = 1/2; a2=1/2; p1 =1; q11 =1;
[taxisRK,xaxisRK,yaxisRK,zaxisRK,i] = Opinion_RK(a1,a2,p1,q11,N,alpha,alpha_prime,time,dt,i);

[rCo,lt,ut,ll,ul] = risetime(xaxisRK,taxisRK,PercentReferenceLevels=[10 90]);

    i = i+1;
    figure(i)
    hold on
    plot(taxisRK,xaxisRK, 'r', 'linewidth',1.0 )
    plot(taxisRK,yaxisRK, 'b', 'linewidth',1.0 )
    plot(taxisRK,zaxisRK, 'g', 'linewidth',1.0 )
    %plot(taxisODE, 10*exp(taxisODE/rCo),'color', [0.6350 0.0780 0.1840], 'linewidth',1.0 )
    plot(lt,ll, 'o', 'color', [0.4940 0.1840 0.5560])
    plot(ut,ul, 'o', 'color', [0.4940 0.1840 0.5560])
    yline([ll ul],'--',{'low level','upper level'})
    title("OPINION MODEL Runge Kutta")
    legend('Co','A','Ca')
    txt = {['alpha: ' num2str(alpha)],['alpha_{prime}: ' num2str(alpha_prime)],['rise time Co: ' num2str(rCo)]};
    text(70,6000,txt)


%% ODE solution
[taxisODE,xaxisODE,yaxisODE,zaxisODE,i] = Opinion_ODE(N,alpha,alpha_prime,lambda,time,dt,i);

[rCo,lt,ut,ll,ul] = risetime(xaxisODE,taxisODE,PercentReferenceLevels=[10 90]);

    i = i+1;
    figure(i)
    hold on
    plot(taxisODE,xaxisODE, 'r', 'linewidth',1.0 )
    plot(taxisODE,yaxisODE, ':b', 'linewidth',1.0 )
    plot(taxisODE,zaxisODE, 'g', 'linewidth',1.0 )
    %plot(taxisODE, 10*exp(taxisODE/rCo),'color', [0.6350 0.0780 0.1840], 'linewidth',1.0 )
    plot(lt,ll, 'o', 'color', [0.4940 0.1840 0.5560])
    plot(ut,ul, 'o', 'color', [0.4940 0.1840 0.5560])
    yline([ll ul],'--',{'low level','upper level'})
    title("OPINION MODEL Ode45")
    legend('Co','A','Ca')
    txt = {['alpha: ' num2str(alpha)],['alpha_{prime}: ' num2str(alpha_prime)],['lambda: ' num2str(lambda)],['rise time Co: ' num2str(rCo)]};
    text(70,6000,txt)
    
    
%% Varying alpha and alpha_prime parameters
% The Initial condition are the same in all the case, but alpha and
% alpha_prime coefficients vary in an interval [0.01,1] 

% d and step used for computations
d = 3;
step = 0.005;
[Co_max, A_max,max_gr_Co, rCo_m, rA_m] = multiple_Opinion(d,step);    

%% Plot the results

alpha_pl = (0.005:step:(step*d))/N;
alpha_prime_pl = (0.005:step:(step*d))/N;
[xx, yy] = meshgrid(alpha_pl,alpha_prime_pl);

i =i+1; 
figure(i)
surf(xx,yy,Co_max)
colorbar;
colormap default;
xlabel('alpha_{prime}') 
ylabel('alpha') 
title('Max value of Compliants for alpha and alpha_{prime}')

i =i+1; 
figure(i)
surf(xx,yy,A_max)
colorbar;
colormap default;
xlabel('alpha_{prime}') 
ylabel('alpha') 
title('Max value of Anti for alpha and alpha_{prime}')

i =i+1; 
figure(i)
surf(xx,yy,rCo_m)
colorbar;
colormap default;
xlabel('alpha_{prime}') 
ylabel('alpha') 
title('Rise time of Compliant for alpha and alpha_{prime}')

i =i+1; 
figure(i)
surf(xx,yy,rA_m)
colorbar;
colormap default;
xlabel('alpha_{prime}') 
ylabel('alpha') 
title('Rise time of Anti for alpha and alpha_{prime}')

%
%% Runge Kutta second order solution of Opinion     
function [taxis,xaxis,yaxis,zaxis,i] = Opinion_RK(a1,a2,p1,q11,N,alpha,alpha_prime,time,dt,i)
    
    
    x = 10;  % Compliant
    y = 10;  % Anti
    z = N-20;% Careless
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
        kx1 =  alpha*x*z;
        ky1 = alpha_prime*y*z;
        kz1 = - alpha*x*z - alpha_prime*y*z;
        % step 2
        t2 = t+p1*dt;
        x2 = x + q11*kx1*dt;
        y2 = y + q11*ky1*dt;
        z2 = z + q11*kz1*dt;
        kx2 = alpha*x2*z2;
        ky2 = alpha_prime*y2*z2;
        kz2 = - alpha*x2*z2 - alpha_prime*y2*z2;
        % update
        x = x + (a1*kx1+a2*kx2)*dt;
        y = y + (a1*ky1+a2*ky2)*dt;
        z = z + (a1*kz1+a2*kz2)*dt;
        t = t + dt;
        cnt = cnt + 1;


    end
end



%% Solve the problem with Matlab ODE functions
function [taxis,xaxis,yaxis,zaxis,i] = Opinion_ODE(N,alpha,alpha_prime,lambda,time,dt,i)
    
    Ca = N-20; % careless
    Co = 10;  % compliant
    A = 10; % anti social rules
    tspan = [0:dt:time];
    y0 = [Co,A,Ca];
    pars = [alpha, alpha_prime,lambda,N];
    
    %ode45 solver input: function to solve, t_span, initial condition 

    [t,y] = ode45(@opi_eqs, tspan, y0, [], pars);
    
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
    
%     i = i+1;
%     figure(i)
%     hold on
%     plot(taxis,xaxis, 'r', 'linewidth',1.0 )
%     plot(taxis,yaxis, 'b', 'linewidth',1.0 )
%     plot(taxis,zaxis, 'g', 'linewidth',1.0 )
%     title("OPINION MODEL")
%     legend('Co','A','Ca')
%     txt = {['alpha: ' num2str(alpha)],['alpha_{prime}: ' num2str(alpha_prime)]};
%     text(70,6000,txt)

end

% Function used to the ODE45 solver to compute correctly the system
%evolution
function f = opi_eqs(t,y,pars)
% 1- alpha - alpha' < 1 Condition
f = zeros(3,1);
f(1) = pars(1)*y(1)*y(3) - pars(3)*y(1);                       % Co' = alpha *Co*Ca -lambda_1 Co
f(2) = pars(2)*y(2)*y(3) - pars(3)*y(2) ;                       % A'  = alpha'*A*Ca - lambda_1 A
f(3) = - pars(1)*y(1)*y(3) - pars(2)*y(2)*y(3)+pars(3)*y(1)+pars(3)*y(2); % Ca' = -alpha*Co*Ca -alpha'*A*Ca 
end


function [Co_max, A_max,max_gr_Co, rCo_m, rA_m] = multiple_Opinion(d,step)
    
    Co_max = zeros(d,d);
    A_max = zeros(d,d);
    max_gr_Co = zeros(d,d);
    rCo_m = zeros(d,d);
    rA_m = zeros(d,d);
    N = 10000;
    time =2500;
    dt=0.001; %un millesimo di secondo di dt
    a1 = 1/2; a2=1/2; p1 =1; q11 =1;
    
    for i = 1:d
        alpha = i*step/N;
        for j = 1:d
            alpha_prime = j*step/N;
            
            %Esecuzione della trasformazione, variando beta e gamma
            
            % Usando ODE
            %[taxis,xaxis,yaxis,zaxis,i] = Opinion_ODE(N,alpha,alpha_prime,time,dt,i);
            
            % Usando Runge Kutta
            [taxis,xaxis,yaxis,zaxis,i] = Opinion_RK(a1,a2,p1,q11,N,alpha,alpha_prime,time,dt,i);
            %Compute max value of infected
            max_Co = max(xaxis);
            max_A = max(yaxis);
            Co_max(i,j) = max_Co
            A_max(i,j) = max_A;
            gr_Co = gradient(xaxis); %Calcolo il valore massimo del gradiente di ogni curva
            max_gr_Co(i,j) = max(gr_Co);
            rCo_m(i,j) = risetime(xaxis,taxis,PercentReferenceLevels=[10 90]);
            rA_m(i,j) = risetime(yaxis,taxis,PercentReferenceLevels=[10 90]);
            
        end
    end
end


