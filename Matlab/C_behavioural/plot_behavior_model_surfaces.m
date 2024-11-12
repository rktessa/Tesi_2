clc;
close all;
%%  Plot of the surface resulting from the 3 equations
x = linspace(0,1,20);
fig = 0;
caso =2; 

[lambda_1,lambda_2,k1,k2, B1, B2, C_zero,A_zero, time,title,X,Y] = scenario(caso);

syms H C A
% Reduced system
eq_H(H,C) = -k1*H*C - k2*(1-H-C)*H + lambda_1*C + lambda_2 * (1-H-C);
eq_C(H,C) = +k1*H*C - lambda_1*C;
% FULL system of equation
eq_H2(H,C,A) = -k1*H*C - k2*A*H + lambda_1*C + lambda_2 * A;
eq_C2(H,C,A) = +k1*H*C - lambda_1*C;
eq_A2(H,C,A) = +k2*H*A - lambda_2*A;

% [X,Y,Z] = ndgrid(x,x,x) ;
% F =eq_H2(X,Y,Z) ;
% F2 =eq_C2(X,Y,Z) ;
% F3 =eq_A2(X,Y,Z) ;


[X2,Y2] = ndgrid(x,x) ;
F_second =eq_H(X2,Y2) ;
F2_second =eq_C(X2,Y2) ;

fig = fig+1;
figure(fig)
hold on

surf(X2,Y2, double(F_second));
colorbar;
surf(X2,Y2, double(F2_second));
% shading interp
view(45,45)
legend('H surface', 'C surface')
colormap("jet")
xlabel('x = Heedless')
ylabel('y = Compliant')
zlabel('Gradient module')
hold off
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 26 23]);
set(gcf, 'PaperSize', [28 24]); % dimension on x axis and y axis resp.
 print(gcf,title,'-dpdf')

% 
% fig = fig+1;
% figure(fig)
% text(X2(ind),Y2(ind),'Z1 = 0','fontweight','bold','color','k')
% ind2 = F2_second == 0;
% plot(X2(ind2),Y2(ind2))
% text(X2(ind2),Y2(ind2),'Z2 = 0','fontweight','bold','color','k')
% contour(X2,Y2,F_second,[0 0],'magenta','linewidth',4)
% % fig = fig+1;
% figure(fig)
% hold on
% for i = 1:length(x)
%     surf(X(:,:,i),Y(:,:,i),Z(:,:,i), double(F(:,:,i))) ;
% end
% colorbar;
% colormap default;
% view(45,45)
% hold off
% 
% fig = fig+1;
% figure(fig)
% hold on
% for i = 1:length(x)
%     surf(X(:,:,i),Y(:,:,i),Z(:,:,i), double(F2(:,:,i))) ;
% end
% colorbar;
% colormap default;
% view(45,45)
% hold off
% 
% fig = fig+1;
% figure(fig)
% hold on
% for i = 1:length(x)
%     surf(X(:,:,i),Y(:,:,i),Z(:,:,i), double(F3(:,:,i))) ;
% end
% colorbar;
% colormap default;
% view(45,45)
% hold off
%% FUNCTIONS
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
           title = 'Surface_nullcline_B1_equal_B2.pdf'
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
        B2 = 1;
        lambda_1 = 1/30; %fatigue to mantain C behaviour
        lambda_2 = 1/30; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        %  Population initial condition
        time = 2000;
        C_zero = 20e6;
        A_zero = 20e6;
        title = 'Pr_nullcline_B1_less_B2_equal_1.pdf';
        X = 1/B1; Y = 1-X;
    end

end