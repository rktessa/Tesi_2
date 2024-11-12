clc;
%clear all;
close all;
%% Behavior model simulation for thesis
% save('nullcline_data')
%load('nullcline_data')
% 24 Ottobre 2024
addpath(genpath('..\'))
obj = functionsContainer;
x = linspace(0,1,30);
y = linspace(0,1,20);
fig = 0;
%% I CASE: B1 < 1 and B2 < 1
%% II CASE: B1 = B2 and greather than 1
%% III CASE: B1 > 1 and B1 > B2
%% IV CASE: B1 > 1 and B1 > B2, lambda2 > lambda1
%% V CASE: B1 < 1 B2 > 1
% Simulation parameters

caso = 5;

[lambda_1,lambda_2,k1,k2, B1, B2, C_zero,A_zero, time,title, X,Y] = scenario(caso);
[taxisRK,HaxisRK,CaxisRK,AaxisRK,awareness] = obj.Behaviour_RK_ic(k1,k2,lambda_1,lambda_2,time,C_zero, A_zero);

syms  xa ya
eq_y(xa,ya) = -k2 *xa^2+(k2*ya -k1*ya -k2 + lambda_2)*xa+ (lambda_1*ya-lambda_2*ya + lambda_2);

eq_y(lambda_1/k1, ya)

% Equilibrium
syms H C A 
% Reduced system
eq_H(H,C) = -k1*H*C - k2*(1-H-C)*H + lambda_1*C + lambda_2 * (1-H-C);
eq_C(H,C) = +k1*H*C - lambda_1*C;
% FULL system of equation
eq_H2(H,C,A) = -k1*H*C - k2*A*H + lambda_1*C + lambda_2 * A;
eq_C2(H,C,A) = +k1*H*C - lambda_1*C;
eq_A2(H,C,A) = +k2*H*A - lambda_2*A;

% The Jacobian of this equation in this case is:
Jac(H, C) = jacobian([eq_H, eq_C], [H, C])
% Computing the Trace
Trac(H,C) = trace(Jac(H,C))
% Determinant
d(H,C) = det(Jac)
% X = lambda_1/k1; Y = 1-lambda_1;
value1 = X;
value2 = Y
trace_J = Trac(value1, value2)
det_J =  d(value1, value2)

% y = x_nullcline(, k1,k2,lambda_1, lambda_2)
%% 

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

% figure nullcline
fig = fig+1;
figure(fig)
hold on
plot(x,x_nullcline(x,k1,k2,lambda_1,lambda_2), LineWidth=3)
xline((lambda_1/k1), 'Color', [0.4660 0.6740 0.1880],LineWidth= 3)
% xline((lambda_2/k2), 'Color', [0.4660 0.6740 0.1880],LineWidth= 3)
yline(0, 'Color', [0.4660 0.6740 0.1880],LineWidth= 3)
quiver(y,y,u2,v2,'r', LineWidth=2)
plot(1,0, 'o', 'Color',[0.4940 0.1840 0.5560], 'LineWidth',5)
plot(lambda_2/k2,0,"d", 'Color',[0.4940 0.1840 0.5560], 'LineWidth',5)
plot(lambda_1/k1,1-lambda_1/k1,"o", 'Color',[0.4940 0.1840 0.5560], 'LineWidth',5)
% xline((lambda_2-lambda_1)/(k2-k1),"--",{'vertical', 'asymptote'},'LineWidth', 3,Color='magenta')
xlabel("x = Heedless")
ylabel("y = Compliant")
legend('x nullcline','', 'y nullcline', 'vector field','equilibrium pts', Orientation='horizontal', Location='southoutside')
txt = {['1/B_1=' num2str(round(lambda_1/k1,3))], ['1/B_2=' num2str(round(lambda_2/k2,3))]};    
dim = [.91 .8 .1 .1];
annotation('textbox',dim, ...
    'String',txt,'EdgeColor','none')
if lambda_1/k1 <1
    txt2 = {'1/B_1'};
    txtA = {'A'};
    pos = [lambda_1/k1+0.1  0.6 0.1 0.1];
    annotation('textbox',pos, ...
        'String',txt2,'Color', [0.4660 0.6740 0.1880], 'EdgeColor', 'none', 'FontSize',15, 'Rotation',90)
    posA = [lambda_1/k1+0.05  1-lambda_1/k1-0.1 0.1 0.1];
    annotation('textbox',posA, ...
        'String',txtA,'Color', [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FontSize',15 )
end
if lambda_2/k2 <1
    txtC = {'C'};
    posC = [lambda_2/k2+0.1  0.2 0.1 0.1];
    annotation('textbox',posC, ...
        'String',txtC,'Color', [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FontSize',15 ) 
end
txtB = {'B'};
posB = [0.9  0.2 0.1 0.1];
annotation('textbox',posB, ...
    'String',txtB,'Color', [0.4940 0.1840 0.5560], 'EdgeColor', 'none', 'FontSize',15 ) 
hold off
xlim([0 1])
ylim([0 1])
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 22 22]);
set(gcf, 'PaperSize', [24 24]); % dimension on x axis and y axis resp.
 print(gcf,title,'-dpdf')

%FIGURA SIMULAZIONE
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




%% Function section

function y = x_nullcline(x, k1,k2,lambda_1, lambda_2)
y = (x.*(k2-k2*x+lambda_2)-lambda_2)./(x.*(k2-k1)+(lambda_1-lambda_2));
end

function y = y_nullcline(Ca, k1,k2,lambda_1, lambda_2)
y = (Ca.*(k1-k1*Ca+lambda_1)-lambda_1)./(Ca.*(k1-k2)-lambda_1+lambda_2);
end
% Da usare nel caso di B1 == B2
function y = y_nullcline2(z, k1,k2,lambda_1, lambda_2)
    R1 = k1/lambda_1; 
    y = 1 - z - 1/R1;
end

function y = z_nullcline2(z, k1,k2,lambda_1, lambda_2)
    R2 = k2/lambda_2;
    y = 1-z - 1/R2;
end

function [lambda_1,lambda_2,k1,k2, B1, B2, C_zero,A_zero, time,title,X,Y] = scenario(caso)
    if caso == 1
    B1 = 0.89;
    B2 = 0.45;
    lambda_1 = 1/30; %fatigue to mantain C behaviour
    lambda_2 = 1/40; %fatigue to mantain A behaviour
    k1 = B1*lambda_1 ; %from H to C
    k2 = B2*lambda_2; %from H to A
    % Population initial condition
    C_zero = 20e6;
    A_zero = 20e6;
    time = 500;
    title = 'Pr_nullcline_B1_B2_less_1.pdf';
    % Equilibrium point
    X = 1; Y = 0; %All Heedles in case 1
    end
    if caso == 2
            B1 = 8.5;
            B2 = 8.5;
            lambda_1 = 1/25; %fatigue to mantain C behaviour
            lambda_2 = 1/30; %fatigue to mantain A behaviour
            k1 = B1*lambda_1 ; %from H to C
            k2 = B2*lambda_2; %from H to A
            time = 2000;
            % Population initial condition
            C_zero = 100;
            A_zero = 100;     
            title = 'Pr_nullcline_B1_B2_equal.pdf';
            % Equilibrium point
            X = 1/B1; Y = 0; %All Heedles in case 1
    end
    if caso == 3
        B1 = 7;
        B2 = 3;
        lambda_1 = 1/30; %fatigue to mantain C behaviour
        lambda_2 = 1/30; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        %  Population initial condition
        time = 200;
        C_zero = 100;
        A_zero = 100;     
        title = 'Pr_nullcline_B1_mag_B2.pdf';
        X = 1/B1; Y = 1-X;
    end
    if caso == 4
        B1 = 7;
        B2 = 3;
        lambda_1 = 1/30; %fatigue to mantain C behaviour
        lambda_2 = 1/5; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        %  Population initial condition
        time = 300;
        C_zero = 100;
        A_zero = 100;
        title = 'Pr_nullcline_B1_mag_B2_lambda2_mag.pdf';
        X = 1/B1; Y = 1-X;
    end
     if caso == 5
        B1 = 0.6;
        B2 = 4;
        lambda_1 = 1/20; %fatigue to mantain C behaviour
        lambda_2 = 1/30; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        %  Population initial condition
        time = 2000;
        C_zero = 20e6;
        A_zero = 20e6;
        title = 'Pr_nullcline_B1_less_B2.pdf';
        X = 1/B1; Y = 1-X;
    end

end