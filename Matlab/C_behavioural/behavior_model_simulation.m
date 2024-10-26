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
time = 500;

%% I CASE: B1 < 1 and B2 < 1
% Simulation parameters
B1 = 0.89;
B2 = 0.45;
lambda_1 = 1/30; %fatigue to mantain C behaviour
lambda_2 = 1/40; %fatigue to mantain A behaviour
k1 = B1*lambda_1 ; %from H to C
k2 = B2*lambda_2; %from H to A
%
C_zero = 20e6;
A_zero = 20e6;
[taxisRK,HaxisRK,CaxisRK,AaxisRK,awareness] = obj.Behaviour_RK_ic(k1,k2,lambda_1,lambda_2,time,C_zero, A_zero);

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
print(gcf,'-dpdf', ['behavior_B1_B2_less_1.pdf'])

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
print(gcf,'-dpdf', ['behavior_B1_equal_B2.pdf'])


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
print(gcf,'-dpdf', ['behavior_B1_mag_B2_k1_mag_k2.pdf'])

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
print(gcf,'-dpdf', ['behavior_B1_mag_B2_k2_mag_k1.pdf'])





