clc;
clear all;
close all;
%% NULLCLINE PLOT OF BEHAVIOURAL MODEL
obj = functionsContainer;
x = linspace(0,1,1000);
% Simulation parameters
R1 = 1.2;
R2 = 1.2;

k1 =  1/4; %from Ca to Co
k2 = 1/2; %from Ca to Ag
lambda_1 = k1/R1; %fatigue to mantain Co behaviour
lambda_2 = k2/R2; %fatigue to mantain Ag behaviour

final_Ca = lambda_1/k1

figure(1)
hold on
plot(x,x_nullcline(x,k1,k2,lambda_1,lambda_2))
xline((lambda_1/k1))
legend('x nullcline', 'y nullcline')
hold off



% R1 = k1/lambda_1;
% R2 = k2/lambda_2;
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
function y = x_nullcline(x, k1,k2,lambda_1,lambda_2)
    y = (x.*(k2-k2*x+lambda_2)-lambda_2)./(x.*(k2-k1)+(lambda_1-lambda_2));
end


function yy = x_prime(x,y, R1)
    yy = -x * (R1+1) + R1*(1-y);
end


function yy = y_prime(x,y, R2)
    yy = - y* (R2+1) + R2*(1-x);
end