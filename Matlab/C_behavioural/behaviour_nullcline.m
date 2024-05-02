clc;
clear all;
close all;
addpath('../Matlab')
%% NULLCLINE PLOT OF BEHAVIOURAL MODEL
x = linspace(0,1,100);
R1 = 1.2;
R2 = 1.2;
obj = functionsContainer;
% Simulation parameters

k1 =  1/4; %from Ca to Co
k2 = 1/2; %from Ca to Ag

lambda_1 = k1/R1; %fatigue to mantain Co behaviour
lambda_2 = k2/R2; %fatigue to mantain Ag behaviour


final_Ca = lambda_1/k1
final_Ca = lambda_2/k2
if R1 > R2
    figure(1)
    hold on
    plot(x,x_nullcline(x, R1))
    plot(x,y_nullcline(x, R2))
    xline((lambda_1/k1))
    legend('x^ = 0', 'y^ = 0', 'z^ = 0')
    hold off
else
    figure(1)
    hold on
    plot(x,x_nullcline(x, R1))
    plot(x,y_nullcline(x, R2))
    xline((lambda_2/k2))
    yline((1-1/R2), 'r')
    yline((1-1/R1), 'b')
    legend('x nullcline', 'y nullcline', 'z_nullcline', 'orizzontal R2','orizzontal R1' )
    hold off
end

time = 500;
Co_zero = 0.001*60e6;
Ag_zero = 0.001*60e6;
[taxisRK,xaxisRK,yaxisRK,zaxisRK,awareness] = obj.Behaviour_RK_ic(k1,k2,lambda_1,lambda_2,time,Co_zero, Ag_zero);
    
figure(2)
hold on
plot(taxisRK,xaxisRK, 'linewidth',1.3 )
plot(taxisRK,yaxisRK, 'linewidth',1.3 )
plot(taxisRK,zaxisRK, 'linewidth',1.3 )
legend('Careless', 'Compliant', 'Against')


%% Caso interessante per i coefficienti da recuperare
% k1 =  1/5; % from Ca to Co
% k2 = 1/2; % from Ca to Ag
% 
% lambda_1 = 1/50; %fatigue to mantain Co behaviour
% lambda_2 = 1/10; %fatigue to mantain Ag behaviour
% R1 = k1/lambda_1
% R2 = k2/lambda_2
% 
% 
% figure(3)
% hold on
% plot(x,x_nullcline(x, R1))
% plot(x,y_nullcline(x, R2))
% xline((lambda_1/k1))
% yline((1-1/R2), 'r')
% yline((1-1/R1), 'b')
% legend('x nullcline', 'y nullcline', 'z_nullcline', 'orizzontal R2','orizzontal R1' )
% hold off


% figure(3)
% y = 0.667;
% plot(x,x_prime(x,y, R1))
% legend('segno x dot funzione di y e R1')
% hold off
% 
% 
% figure(4)
% xx = 0.5;
% y = x; % vario lungo asse y con x fisso
% plot(x,y_prime(xx,y, R2))
% legend('segno x dot funzione di y e R2')
% hold off
%% FUNCTIONS
function y = x_nullcline(x, R1)
    y = -x - 1/R1 +1;
end

function y = y_nullcline(x, R2)
    y = -x - 1/R2 +1;
end


function yy = x_prime(x,y, R1)
    yy = -x * (R1+1) + R1*(1-y);
end


function yy = y_prime(x,y, R2)
    yy = - y* (R2+1) + R2*(1-x);
end