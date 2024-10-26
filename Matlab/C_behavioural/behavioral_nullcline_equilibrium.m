clc;
clear all;
close all;
%% Behavior model simulation for thesis
% 24 Ottobre 2024


addpath(genpath('..\'))
obj = functionsContainer;
x = linspace(0,1,1000);
y = linspace(0,1,20);
fig = 0;
time = 800;

%% I CASE: B1 < 1 and B2 < 1
% Simulation parameters
B1 = 0.89;
B2 = 0.65;
lambda_1 = 1/31; %fatigue to mantain C behaviour
lambda_2 = 1/40; %fatigue to mantain A behaviour
k1 = B1*lambda_1 ; %from H to C
k2 = B2*lambda_2; %from H to A

den = (lambda_2-lambda_1)/(k2-k1)

% k1 = 1/7; %from Ca to Co
% k2 = 1/2; %from C to Ag
% lambda_1 = k1/B1; %fatigue to mantain Co behaviour
% lambda_2 = k2/B2; %fatigue to mantain Ag behaviour
%
C_zero = 20e6;
A_zero = 20e6;
[taxisRK,HaxisRK,CaxisRK,AaxisRK,awareness] = obj.Behaviour_RK_ic(k1,k2,lambda_1,lambda_2,time,C_zero, A_zero);

% Equilibrium
syms H C A 
eq_H(H,C) = -k1*H*C - k2*(1-H-C)*H + lambda_1*C + lambda_2 * (1-H-C);
eq_C(H,C) = +k1*H*C - lambda_1*C;
% Plot of the nullcline
% Calculation of system vector field 
u = zeros(length(y), length(y));
v = zeros(length(y), length(y));
for i = 1:numel(y)
    for j = 1:numel(y)
        u(i,j) = eq_H(y(i),y(j));
        v(i,j) = eq_C(y(i),y(j));
    end
end
u2 = u; u2 = transpose(u2);
v2 = v; v2 = transpose(v2);

% figure
fig = fig+1;
figure(fig)
hold on
plot(x,x_nullcline(x,k1,k2,lambda_1,lambda_2), LineWidth=1.3)
xline((lambda_1/k1), LineWidth= 1.3)
quiver(y,y,u2,v2,'r')
xlabel("H")
ylabel("C")
legend('x nullcline', 'y nullcline', 'vector field', Orientation='horizontal', Location='southoutside')
hold off
xlim([0 1])
ylim([0 1])
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 23 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', ['nullcline_B1_B2_less_1.pdf'])



fig = fig+1;
figure(fig)
hold on
plot(taxisRK,HaxisRK, 'linewidth',1.4 )
plot(taxisRK,CaxisRK, 'linewidth',1.4 )
plot(taxisRK,AaxisRK, 'linewidth',1.4 )
xlabel("t[days]");
ylabel("H[t], C[t], A[t]");
ylim([0,1])
legend('Heedless', 'Compliant', 'Against', Orientation='horizontal', Location='southoutside')
txt = {['B_1=' num2str(B1)], ['B_2=' num2str(B2)], ['k_1=' num2str(round(k1,2))],['k_2=' num2str(round(k2,2))],['\lambda_1=' num2str(round(lambda_1,2))],['\lambda_2=' num2str(round(lambda_2,2))]};    
dim = [.91 .8 .1 .1];
annotation('textbox',dim, ...
    'String',txt,'EdgeColor','none')
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 23 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
% print(gcf,'-dpdf', ['behavior_B1_B2_less_1.pdf'])

%% II CASE: B1 = B2 and greather than 1
% Simulation parameters
B1 = 10;
B2 = 10;
lambda_1 = 1/30; %fatigue to mantain C behaviour
lambda_2 = 1/20; %fatigue to mantain A behaviour
k1 = B1*lambda_1 ; %from H to C
k2 = B2*lambda_2; %from H to A
%
time = 100;
C_zero = 100;
A_zero = 100;
[taxisRK,HaxisRK,CaxisRK,AaxisRK,awareness] = obj.Behaviour_RK_ic(k1,k2,lambda_1,lambda_2,time,C_zero, A_zero);



%% III CASE: B1 > 1 and B1 > B2
% Simulation parameters
B1 = 7;
B2 = 3;
lambda_1 = 1/30; %fatigue to mantain C behaviour
lambda_2 = 1/30; %fatigue to mantain A behaviour
k1 = B1*lambda_1 ; %from H to C
k2 = B2*lambda_2; %from H to A
%
time = 200;
C_zero = 100;
A_zero = 100;
[taxisRK,HaxisRK,CaxisRK,AaxisRK,awareness] = obj.Behaviour_RK_ic(k1,k2,lambda_1,lambda_2,time,C_zero, A_zero);


%% IV CASE: B1 > 1 and B1 > B2
% Simulation parameters
B1 = 7;
B2 = 3;
lambda_1 = 1/30; %fatigue to mantain C behaviour
lambda_2 = 1/5; %fatigue to mantain A behaviour
k1 = B1*lambda_1 ; %from H to C
k2 = B2*lambda_2; %from H to A
%
time = 300;
C_zero = 100;
A_zero = 100;
[taxisRK,HaxisRK,CaxisRK,AaxisRK,awareness] = obj.Behaviour_RK_ic(k1,k2,lambda_1,lambda_2,time,C_zero, A_zero);

%% Function section

function y = x_nullcline(x, k1,k2,lambda_1, lambda_2)
y = (x.*(k2-k2*x+lambda_2)-lambda_2)./(x.*(k2-k1)+(lambda_1-lambda_2));
end

function y = y_nullcline(Ca, k1,k2,lambda_1, lambda_2)
y = (Ca.*(k1-k1*Ca+lambda_1)-lambda_1)./(Ca.*(k1-k2)-lambda_1+lambda_2);
end

function y = x_nullcline2(x, k1,k2,lambda_1, lambda_2)
    R1 = k1/lambda_1; 
    y = -x - 1/R1 +1;
end

function y = y_nullcline2(x, k1,k2,lambda_1, lambda_2)
    R2 = k2/lambda_2;
    y = -x - 1/R2 +1;
end





