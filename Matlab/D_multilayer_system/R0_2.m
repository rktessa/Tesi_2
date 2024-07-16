%% R_0 estimation II version
% In this code is find and calculated the value of the R_0 index. 
% The method used to estimate is the one describe in Arino_2007 work. 
clc;
clear;
close all;
%% Initialization
fig = 0;

%% Esempio 1          
% Diamo dei valori ai parametri
rho1 = 0.65; % protezione da infezione
epsilon1 = 0.15; % gli IC che vanno a infettare in giro
k41 = 0.243; k31 = 0.48;
lambda41 = 0.143; lambda31 = 0.143;
gamma1 = 0.35;  
beta1 = 0.40;
% Population division
    SC01 = 50/60e6; SA01 = 50/60e6; 
    IC01 = 10/60e6; IA01 = 10/60e6;
    RC01 = 0;  RA01 = 0;
    phi1 = 0.9091;
    % Calculation of the resulting R_0_epi_behav value
    R_0_alt = Reproductive_rate(SC01, SA01,IC01, IA01, RC01, RA01, gamma1, lambda31, lambda41,k31, k41, phi1, epsilon1, rho1 )
    % Total R_0_behav_epi
    R_0_alt_tot = R_0_alt*beta1/gamma1

    
%% Varying initial conditions plot with different SC0 and SA0
%% Population change rapidly behavior, small lambdas
% Caso E3 > E4 >1 and population change rapidly behavior
E3 = 3; E4 = 1.4;
lambda3 = 1/7; lambda4 = 1/7;
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_SC_and_SA, fig] = varying_SC0_SA0(fig,k3,lambda3,k4,lambda4);
% Caso E3 < E4 >1 and population change rapidly behavior
E3 = 2; E4 = 20;
lambda3 = 1/7; lambda4 = 1/7;
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_SC_and_SA, fig] = varying_SC0_SA0(fig,k3,lambda3,k4,lambda4);
% Caso E3 = E4 >1 and population change rapidly behavior
E3 = 2.4; E4 = 2.4;
lambda3 = 1/7; lambda4 = 1/7;
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_SC_and_SA, fig] = varying_SC0_SA0(fig,k3,lambda3,k4,lambda4);

%% Population stubborn
% Caso E3 > E4 >1 and population change rapidly behavior in Compliant but
% is stubborn in Against group
E3 = 3; E4 = 1.4;
lambda3 = 1/7; lambda4 = 1/40;
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_SC_and_SA, fig] = varying_SC0_SA0(fig,k3,lambda3,k4,lambda4);
% Caso E3 < E4 >1 and population change rapidly behavior
E3 = 2; E4 = 20;
lambda3 = 1/7; lambda4 = 1/40;
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_SC_and_SA, fig] = varying_SC0_SA0(fig,k3,lambda3,k4,lambda4);
% Caso E3 = E4 >1 and population change rapidly behavior
E3 = 2.4; E4 = 2.4;
lambda3 = 1/7; lambda4 = 1/40
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_SC_and_SA, fig] = varying_SC0_SA0(fig,k3,lambda3,k4,lambda4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Varying initial conditions plot with different rho and epsilon
%% Population change rapidly behavior, small lambdas
% Caso E3 > E4 >1 and population change rapidly behavior
E3 = 3; E4 = 100.4;
lambda3 = 1/7; lambda4 = 1/7;
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_rho_and_epsilon, fig] = varying_epsilon_rho(fig,k3,lambda3,k4,lambda4);

% Caso E3 >> E4 >1 and population change rapidly behavior
E3 = 30; E4 = 1.4;
lambda3 = 1/7; lambda4 = 1/7;
k3 = E3*lambda3; k4 = E4*lambda4;
[R0_initial_rho_and_epsilon, fig] = varying_epsilon_rho(fig,k3,lambda3,k4,lambda4);

%% Funzioni
%% Definizione del R0 epi_behav senza R0_epi = beta/gamma
% Forma alternativa del R_0 raccogliendo e portando fuori un beta/gamma per
% poter fare analisi di sensitività escludendo quel valore epidemiologico.
% Così ottengo una analisi più generica, meno basata sulla malattia.
% Comunque non riesco a togliere un gamma che rimane a combnarsi con
% qualche parametro. 

function R_0_alt = Reproductive_rate(SC0, SA0,IC0, IA0, RC0, RA0, gamma, lambda3, lambda4,k3, k4, phi, epsilon, rho )
    
    SH0 = 1- SC0 - SA0 - IA0 - IC0 - RC0 - RA0; %semplifico il modello eliminando una variabile
    C = SC0+IC0+RC0;
    A = SA0+IA0+RA0;
    x1 = lambda3+ A*k4;
    x2 = lambda4+phi*C*k3;
   
    R_0_alt = SA0*((gamma+x1+epsilon*x2)/(x1+x2+gamma)) + (SH0+ rho*SC0)*((x1+epsilon*(x2+gamma))/(x1+x2+gamma));

end

function [R0_initial_SC_and_SA, fig] = varying_SC0_SA0(fig,k3,lambda3,k4,lambda4)
    % Multiple plots varying initial number of Compliant and Against,
    % respectively the SC0 and SA0 variables. 
    % This change affect also the A and C groups. Several heatmaps are plotted using different awareness value.
    SC0i = linspace(50/60e6,29.9e6/60e6,100); %interval from 50 to 30 milion 
    SA0j = linspace(50/60e6,29.9e6/60e6,100); 
    IC_0 = 10/60e6; IA_0 = 10/60e6;
    RC_0 = 0;  RA_0 = 0;
    gamma_0 = 1/8.5; 
    R0_initial_SC_and_SA = zeros(length(SC0i),length(SA0j));
    epsilon_0 = 0.3; % rate of I_c NOT staying at home
    rho_0 = 0.35; % rate of protection to infection
    phi_1 = 1; %in this sensitivity The Awareness is not considered initially in the study
    for i = 1:length(SC0i)
        for j = 1:length(SA0j)
            R0_initial_SC_and_SA(i,j) = Reproductive_rate(SC0i(i), SA0j(j),IC_0, IA_0, RC_0, RA_0, gamma_0, lambda3, lambda4,k3, k4, phi_1, epsilon_0, rho_0 );
        end
    end 
    %Plot
    d1 = 1;
    fig = fig+1;
    figure(fig)
    box on
    [xx, yy] = meshgrid(SC0i,SA0j);
    %aga = reshape(R0_initial_SC_and_SA(:,:,d1),10,10,[]);
    zz = R0_initial_SC_and_SA; zz = transpose(zz);
    surface(xx,yy,zz, 'edgecolor','none')
    colorbar;
    colormap default;
    txt3 = "Epi-behav R_0 with: "+ 'E_3 = '+ num2str(k3/lambda3) +', E_4 = '+ num2str(k4/lambda4)+', \lambda_3 = '+ num2str(lambda3)+', \lambda_4 = '+ num2str(lambda4);
    title(txt3)
    xlabel('S_C_0') 
    ylabel('S_A_0') 
    fontsize(20,"points")
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 24 15]);
    set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
    text1 = "ro_epi-behav_SC_SA_num_"+ num2str(fig); 
    print(gcf,'-dpdf', [text1])
end


function [R0_initial_rho_and_epsilon, fig] = varying_epsilon_rho(fig,k3,lambda3,k4,lambda4)
    % Multiple plots varying initial number of Compliant and Against,
    % respectively the SC0 and SA0 variables. 
    % This change affect also the A and C groups. Several heatmaps are plotted using different awareness value.
    rho_i = linspace(0,1,100);  % rate of protection to infection
    epsilon_j = linspace(0,1,100); % rate of I_c NOT staying at home
    SC0 = 50/60e6;
    SA0 = 50/60e6;
    IC_0 = 10/60e6; IA_0 = 10/60e6;
    RC_0 = 0;  RA_0 = 0;
    gamma_0 = 1/8.5; 
    R0_initial_rho_and_epsilon = zeros(length(rho_i),length(epsilon_j));
    
    phi_1 = 1; %in this sensitivity The Awareness is not considered initially in the study
    for i = 1:length(rho_i)
        for j = 1:length(epsilon_j)
            R0_initial_rho_and_epsilon(i,j) = Reproductive_rate(SC0, SA0,IC_0, IA_0, RC_0, RA_0, gamma_0, lambda3, lambda4,k3, k4, phi_1, epsilon_j(j), rho_i(i) );
        end
    end 
    %Plot
    d1 = 1;
    fig = fig+1;
    figure(fig)
    box on
    [xx, yy] = meshgrid(rho_i,epsilon_j);
    %aga = reshape(R0_initial_SC_and_SA(:,:,d1),10,10,[]);
    zz = R0_initial_rho_and_epsilon; zz = transpose(zz);
    surface(xx,yy,zz, 'edgecolor','none')
    colorbar;
    colormap default;
    txt3 = "Epi-behav R_0 with: "+ 'E_3 = '+ num2str(k3/lambda3) +', E_4 = '+ num2str(k4/lambda4)+', \lambda_3 = '+ num2str(lambda3)+', \lambda_4 = '+ num2str(lambda4);
    title(txt3)
    xlabel('\rho') 
    ylabel('\epsilon') 
    fontsize(20,"points")
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 24 15]);
    set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
    text1 = "ro_epi-behav_rho_epsilon_num_"+ num2str(fig); 
    print(gcf,'-dpdf', [text1])
end