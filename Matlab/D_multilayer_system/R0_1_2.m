%% R_0 estimation 1_second version
% In this code is find and calculated the value of the R_0 index. 
% The method used to estimate is the one describe in Arino_2007 work. 
% RESULT: A prima vista i vari casi non differiscono molto. 
%% Initialization
 fig = 0;
% Diamo dei valori ai parametri
rho1 = 0.65; % protezione da infezione
epsilon1 = 0.15; % gli IC che vanno a infettare in giro
k41 = 0.243; k31 = 0.48;
lambda41 = 0.143; lambda31 = 0.143;
gamma1 = 1/9; beta1 = 0.40;
% Population division
    SC01 = 50/60e6; SA01 = 50/60e6; SH01 = 1-SC01 - SA01;
    IC01 = 10/60e6; IA01 = 10/60e6;
    RC01 = 0;  RA01 = 0;
    C1 = SC01 + IC01+RC01; A1 = SA01 +IA01 +RA01;
 
%% Varying initial conditions
%% I CASE: B1 < 1 and B2 < 1
%% II CASE: B1 = B2 and greather than 1
%% III CASE: B1 > 1 and B1 > B2
%% IV CASE: B1 > 1 and B1 > B2, lambda2 > lambda1
caso = 1;
[lambda_3,lambda_4,k_3,k_4, B3, B4,title_1,title_2] = scenario(caso);
fig = plot_cases_SC0_SA0(fig,epsilon1, k_3, k_4, lambda_3, lambda_4, beta1, gamma1,rho1, title_1);
fig = plot_cases_rho_SC0(fig,epsilon1, k_3, k_4, lambda_3, lambda_4, beta1, gamma1, title_2);

%% function section

function [lambda_3,lambda_4,k3,k4, B3, B4,title_1,title_2] = scenario(caso)
    if caso == 1
        B3 = 0.89;
        B4 = 0.45;
        lambda_3 = 1/30; %fatigue to mantain C behaviour
        lambda_4 = 1/40; %fatigue to mantain A behaviour
        k3 = B3*lambda_3 ; %from H to C
        k4 = B4*lambda_4; %from H to A
        title_1 = 'HMap_SC0_SA0_B1_B2_less_1.pdf';
        title_2 = 'HMap_SC0_rho_B1_B2_less_1.pdf';
    end
    if caso == 2
        B3 = 8.5;
        B4 = 8.5;
        lambda_3 = 1/25; %fatigue to mantain C behaviour
        lambda_4 = 1/30; %fatigue to mantain A behaviour
        k3 = B3*lambda_3 ; %from H to C
        k4 = B4*lambda_4; %from H to A
        title_1 = 'HMap_SC0_SA0_B1_B2_equal.pdf';   
        title_2 = 'HMap_SC0_rho_B1_B2_equal.pdf';  
    end
    if caso == 3
        B3 = 7;
        B4 = 3;
        lambda_3 = 1/30; %fatigue to mantain C behaviour
        lambda_4 = 1/30; %fatigue to mantain A behaviour
        k3 = B3*lambda_3 ; %from H to C
        k4 = B4*lambda_4; %from H to A
        title_1 = 'HMap_SC0_SA0_B1_mag_B2.pdf';
        title_2 = 'HMap_SC0_rho_B1_mag_B2.pdf'; 
    end
    if caso == 4
        B3 = 7;
        B4 = 3;
        lambda_3 = 1/30; %fatigue to mantain C behaviour
        lambda_4 = 1/5; %fatigue to mantain A behaviour
        k3 = B3*lambda_3 ; %from H to C
        k4 = B4*lambda_4; %from H to A
        %  Population initial condition
        title_1 = 'HMap_SC0_SA0_B1_mag_B2_lambda2_mag.pdf';
        title_2 = 'HMap_SC0_rho_B1_mag_B2_lambda2_mag.pdf';
    end
end


%% Varying SC0 and SA0
function fig = plot_cases_SC0_SA0(fig,epsilon1, k31, k41, lambda31, lambda41, beta1, gamma1,rho1, titolo)
    syms psi2 SC SA SH IC IA rho epsilon k3 k_4 lambda3 lambda4 beta2 gamma2 A C SC0 SA0 SH0 
    % To calculate the R0 the following matrices are defined:
    Pi = [ rho, 1 0; 0, 0, 1];
    D =  [ 1, 0, 0; 0, 1, 0; 0, 0, 1];
    y =  [ SC; SH; SA];
    b =  [ epsilon, 1];
    x =  [ IC; IA];  
    V = [  k_4*A + lambda3 + gamma2, -psi2*k3*C - lambda4; 
          -k_4*A - lambda3, psi2 * k3*C + gamma2 + lambda4 ];
    y0= [SC0; SH0; SA0];
    %% Definizione del R0
    R_0(psi2, rho, epsilon, k3, k_4, lambda3, lambda4, beta2, gamma2, A, C, SC0, SA0, SH0) = beta2 .* b * inv(V) * Pi * D * y0;
    R_0_simplified(psi2, rho, epsilon, k3, k_4, lambda3, lambda4, beta2, gamma2, A, C, SC0, SA0, SH0) = gamma2/beta2 *R_0;

    % Multiple plots varying initial number of Compliant and Against, respectively the SC0 and SA0 variables. This change affect also the A and C groups. Several heatmaps are plotted using different awareness value.
    SC0i = linspace(50/60e6,30e6/60e6,10); %interval from 50 to 10 milion 
    SA0i = linspace(50/60e6,30e6/60e6,10); 
    psi = 1; %not considered the effect in this first analysis
    IC01 =50/60e6; IA01 = 50/60e6;
    R0_initial_SC_and_SA = zeros(length(SC0i),length(SA0i), length(psi));
    for i = 1:length(SC0i)
        for j = 1:length(SA0i)
            for k = 1: length(psi)
                SH01 = 1-SC0i(i)-SA0i(j);
                A1 = SA0i(j)+ IA01 + 0;
                C1 = SC0i(i) + IC01 + 0;
                % x1 = lambda31 + A1*k41;
                % x2 = lambda41 + awk(k) *C1*k31;
                % ratio1 = (gamma1+x1+epsilon1*x2)/(x1+x2+gamma1);
                % ratio2 = (x1+epsilon1*(x2+gamma1))/(x1+x2+gamma1);
                % somma = SA0j(j)*ratio1 + (SH01 + rho1*SC0i(i))*ratio2;
                R0_initial_SC_and_SA(i,j,k) = R_0_simplified(psi(k), rho1, epsilon1, k31, k41, lambda31, lambda41, beta1, gamma1, A1, C1, SC0i(i), SA0i(j), SH01);
            end
         end 
    end 
    %Plot
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
    txt3 = "E_0^{s} function of initial SC and SA";
    title(txt3)
    xlabel('SC_0') 
    ylabel('SA_0') 
    fontsize(20,"points")
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 24 15]);
    set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
    print(gcf,'-dpdf', titolo)
end

%% Varying SC0 e rho
function fig = plot_cases_rho_SC0(fig, epsilon1, k31, k41, lambda31, lambda41, beta1, gamma1, titolo)
     syms psi2 SC SA SH IC IA rho epsilon k3 k_4 lambda3 lambda4 beta2 gamma2 A C SC0 SA0 SH0 
    % To calculate the R0 the following matrices are defined:
    Pi = [ rho, 1 0; 0, 0, 1];
    D =  [ 1, 0, 0; 0, 1, 0; 0, 0, 1];
    y =  [ SC; SH; SA];
    b =  [ epsilon, 1];
    x =  [ IC; IA];  
    V = [  k_4*A + lambda3 + gamma2, -psi2*k3*C - lambda4; 
          -k_4*A - lambda3, psi2 * k3*C + gamma2 + lambda4 ];
    y0= [SC0; SH0; SA0];
    %% Definizione del R0
    R_0(psi2, rho, epsilon, k3, k_4, lambda3, lambda4, beta2, gamma2, A, C, SC0, SA0, SH0) = beta2 .* b * inv(V) * Pi * D * y0;
    R_0_simplified(psi2, rho, epsilon, k3, k_4, lambda3, lambda4, beta2, gamma2, A, C, SC0, SA0, SH0) = beta2/gamma2 *R_0;
    % Multiple plots varying initial number of Compliant and Against, respectively the SC0 and SA0 variables. This change affect also the A and C groups. Several heatmaps are plotted using different awareness value.
    SC0i = linspace(50/60e6,30e6/60e6,10); %interval from 50 to 10 milion 
    rhoj = linspace(0,1,10); 
    SA01 = 50/60e6; IC01 =50/60e6; IA01 = 50/60e6;
    psi = 1; %not considered in this initial case
    R0_initial_SC_and_SA = zeros(length(SC0i),length(rhoj), length(psi));
    for i = 1:length(SC0i)
        for j = 1:length(rhoj)
            for k = 1: length(psi)
                SH01 = 1-SC0i(i)-SA01;
                A1 = SA01+ IA01 + 0;
                C1 = SC0i(i) + IC01 + 0;
                R0_initial_SC_and_SA(i,j,k) = R_0_simplified(psi(k), rhoj(j), epsilon1, k31, k41, lambda31, lambda41, beta1, gamma1, A1, C1, SC0i(i), SA01, SH01);
            end
         end 
    end
    %Figure
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
    txt3 = "E_0^{s} function of initial SC and efficacy of NPIs";
    title(txt3)
    xlabel('initial susceptible compliant') 
    ylabel('\rho ') 
    fontsize(20,"points")
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 24 15]);
    set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
    print(gcf,'-dpdf', titolo)
end