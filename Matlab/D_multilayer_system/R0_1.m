%% R_0 estimation I version
% In this code is find and calculated the value of the R_0 index. 
% The method used to estimate is the one describe in Arino_2007 work. 
clc;
clear;
close all;
%% Initialization
syms psi SC SA SH IC IA rho epsilon k3 k4 lambda3 lambda4 beta gamma A C SC0 SA0 SH0 
fig = 0;
% To calculate the R0 the following matrices are defined:
    Pi = [ rho, 1 0; 0, 0, 1];
    D =  [ 1, 0, 0; 0, 1, 0; 0, 0, 1];
    y =  [ SC; SH; SA];
    b =  [ epsilon, 1];
    x =  [ IC; IA];  
    V = [  k4*A + lambda3 + gamma, -psi*k3*C - lambda4; 
          -k4*A - lambda3, psi * k3*C + gamma + lambda4 ];
    y0= [SC0; SH0; SA0];


%% Definizione del R0
R_0(psi, rho, epsilon, k3, k4, lambda3, lambda4, beta, gamma, A, C, SC0, SA0, SH0) = beta .* b * inv(V) * Pi * D * y0;
R_0_simplified(psi, rho, epsilon, k3, k4, lambda3, lambda4, beta, gamma, A, C, SC0, SA0, SH0) = beta/gamma *R_0;
% Diamo dei valori ai parametri
rho1 = 0.65; % protezione da infezione
epsilon1 = 0.15; % gli IC che vanno a infettare in giro
k41 = 0.243; k31 = 0.48;
lambda41 = 0.143; lambda31 = 0.143;
gamma1 = 0.35; beta1 = 0.40;
% Population division
    SC01 = 50/60e6; SA01 = 50/60e6; SH01 = 1-SC01 - SA01;
    IC01 = 10/60e6; IA01 = 10/60e6;
    RC01 = 0;  RA01 = 0;
    C1 = SC01 + IC01+RC01; A1 = SA01 +IA01 +RA01;
 
%% Varying initial conditions

% Multiple plots varying initial number of Compliant and Against, respectively the SC0 and SA0 variables. This change affect also the A and C groups. Several heatmaps are plotted using different awareness value.
SC0i = linspace(50/60e6,30e6/60e6,10); %interval from 50 to 10 milion 
SA0i = linspace(50/60e6,30e6/60e6,10); 

psi = 1; %not considered the effect in this first analysis
R0_initial_SC_and_SA = zeros(length(SC01),length(SA01), length(psi));
for i = 1:length(SC0i)
    for j = 1:length(SA0i)
        for k = 1: length(psi)
            SH01 = 1-SC0i(i)-SA0i(j);
            A1 = SA0i(j)+ 50/60e6 + 0;
            C1 = SC0i(i) + 50/60e6 + 0;

            % x1 = lambda31 + A1*k41;
            % x2 = lambda41 + awk(k) *C1*k31;
            % ratio1 = (gamma1+x1+epsilon1*x2)/(x1+x2+gamma1);
            % ratio2 = (x1+epsilon1*(x2+gamma1))/(x1+x2+gamma1);
            % somma = SA0j(j)*ratio1 + (SH01 + rho1*SC0i(i))*ratio2;
            R0_initial_SC_and_SA(i,j,k) = R_0(psi(k), rho1, epsilon1, k31, k41, lambda31, lambda41, beta1, gamma1, A1, C1, SC0i(i), SA0i(j), SH01);
        end
     end 
end
% filename = "R0_variandoSC_SA,.mat";
% save(filename, "R0_initial_SC_and_SA", '-v7.3');
% load("R0_variandoSC_SA,.mat")
% Prossimi casi, plot funzione di awareness SC e C, SC,SA C e A  più egli altri coeeficienit, anche se solo roh e epsilon hannno un valore compartamentale. 
% Varying SC0 and SA0

% SC0i = linspace(50/60e6,10e6/60e6); %interval from 50 to 10 milion 
% SA0j = linspace(50/60e6,10e6/60e6); 
% awk = linspace(0.1,100,10);
%% Equilibrium plot of R1 and R2

%1,1
d1 = 1;
fig = fig+1;
figure(fig)
box on
[xx, yy] = meshgrid(SC0i,SA0i);
aga = reshape(R0_initial_SC_and_SA(:,:,d1),10,10,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "R_0 function of initial SC and SA";
title(txt3)
xlabel('initial susceptible compliant') 
ylabel('initial susceptible against') 
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
print(gcf,'-dsvg', ['r0_epi_behav_SC_SA.svg'])



%% Varying SC0 e rho
% Multiple plots varying initial number of Compliant and Against, respectively the SC0 and SA0 variables. This change affect also the A and C groups. Several heatmaps are plotted using different awareness value.
SC0i = linspace(50/60e6,30e6/60e6,10); %interval from 50 to 10 milion 
rhoj = linspace(0,1,10); 
SA01 = 50/60e6;
psi = 1; %not considered in this initial case
R0_initial_SC_and_SA = zeros(length(SC01),length(SA01), length(psi));
for i = 1:length(SC0i)
    for j = 1:length(rhoj)
        for k = 1: length(psi)
            SH01 = 1-SC0i(i)-SA01;
            A1 = SA01+ 50/60e6 + 0;
            C1 = SC0i(i) + 50/60e6 + 0;

            % x1 = lambda31 + A1*k41;
            % x2 = lambda41 + awk(k) *C1*k31;
            % ratio1 = (gamma1+x1+epsilon1*x2)/(x1+x2+gamma1);
            % ratio2 = (x1+epsilon1*(x2+gamma1))/(x1+x2+gamma1);
            % somma = SA0j(j)*ratio1 + (SH01 + rho1*SC0i(i))*ratio2;
            R0_initial_SC_and_SA(i,j,k) = R_0(psi(k), rhoj(j), epsilon1, k31, k41, lambda31, lambda41, beta1, gamma1, A1, C1, SC0i(i), SA01, SH01);
        end
     end 
end



%1,1
d1 = 1;
fig = fig+1;
figure(fig)
box on
[xx, yy] = meshgrid(SC0i,rhoj);
aga = reshape(R0_initial_SC_and_SA(:,:,d1),10,10,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "R_0 function of initial SC and efficacy of NPIs";
title(txt3)
xlabel('initial susceptible compliant') 
ylabel('\rho ') 
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', ['r0_epi_behav_SC_rho.pdf'])