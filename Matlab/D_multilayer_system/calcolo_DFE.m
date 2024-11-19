clc;
clear;
close all;
%% Set the scenario; first number scenario, second what F/V matrix compute
caso =[2,2];
%% Run the code to se R_0, B1, B2 and the three resulting E_0 values
[lambda_1,lambda_2,k_1,k_2,time,title_fig, beta,gamma,delta_v,rho_v,epsilon_v,psi_v, SH0_val,SC0_val,SA0_val,IC0_val,IA0_val,RA0_val,RC0_val]=scenario(caso(1));
phi_n = 0.5;
DFE = equilibri_FDE(psi_v,rho_v,epsilon_v, k_1, k_2, lambda_1, lambda_2, beta, gamma, delta_v,phi_n)

E_0 = zeros(1,size(DFE,1));
for i = 1:size(DFE,1)
    SC_0 = DFE(i,1); SA_0 = DFE(i,3); SH_0 = DFE(i,2); RC_0 = DFE(i,4); RA_0 = DFE(i,5); IC_0 =0; IA_0 = 0; RC_0 = 0; RA_0 = 0;
    [E_0(i), R_0, B_1, B_2,FV] = calcolo_E_0_Van(caso(2),beta, gamma,psi_v, rho_v, epsilon_v, k_1, k_2,lambda_1, lambda_2, SH_0, SC_0,SA_0, IC_0,IA_0,delta_v,RA_0,RC_0);
end

% B_1
% B_2
% R_0
E_0
%% Function Disease Free Equilibrium
function FDE = equilibri_FDE(psi_v,rho_v,epsilon_v, k_1, k_2, lambda_1, lambda_2, beta, gamma, delta_v,phi_n_val)
    syms psi2 SC SA RA RC SH IC IA rho epsilon k1 k2 k3 k4 k5 k6 lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 beta2 gamma2 SC0 SH0 SA0 RC0 RA0 delta phi_n
        
        %  F = [epsilon*beta2*(rho*SC+SH)+ psi2*k3*IA,   beta2*(rho*SC+SH)+psi2*k3*(SC+IC+RC)+lambda4;
        %      beta2*epsilon*SA+ lambda3+k4*(SA+IA+RA), beta2*SA+ k4*IC];
        % 
        % V = [lambda3+k4*(SA+IA+RA)+gamma2, k4*IC;
        %     psi2*k3*IA, psi2*k3*(SC+IC+RC)+lambda4+gamma2];
        
        psi2 =psi_v;  rho=rho_v; epsilon=epsilon_v; 
        k1 = k_1; k2 = k_2; k3=k_1; k4=k_2; k5 = k_1; k6 = k_2;
        lambda1 = lambda_1; lambda2= lambda_2; lambda3=lambda_1; lambda4=lambda_2; lambda5 = lambda_1; lambda6 =lambda_2;
        beta2=beta; gamma2=gamma; delta = delta_v; phi_n = phi_n_val;
        % DFE
        IA = 0; IC = 0;
       
        % FV = F/V;
        % 
        eq1 = - beta2*SH*(epsilon*IC+IA) - k1*psi2*SH*(SC+IC+RC) -k2*SH*(SA+IA+RA) +lambda1*SC + lambda2*SA + delta*(1-phi_n)*RC == 0;
        eq2 = k1*psi2*SH*(SC+IC+RC) - lambda1*SC -beta2*rho*SC*(epsilon*IC+IA) + delta*phi_n*RC==0;
        eq3 = k2*SH*(SA+IA+RA) - lambda2*SA - beta2*SA*(epsilon*IC+IA) + delta*RA==0;
        eq4 = beta2*rho*SC*(epsilon*IC+IA) + beta2*SH*(epsilon*IC+IA) - k4 *IC*(SA+IA+RA) + k3*psi2*IA*(SC+IC+RC) - lambda3*IC- gamma2*IC + lambda4*IA==0;
        eq5 = +beta2*SA*(epsilon*IC+IA) + k4*IC*(SA+IA+RA)+lambda3*IC - k3*psi2*IA*(SC+IC+RC)-gamma2*IA - lambda4*IA==0;
        eq6 = gamma2 * IC - k6*RC*(SA+IA+RA) + psi2*k5*RA*(SC+IC+RC) - delta*RC + lambda6*RA - lambda5 * RC==0;
        eq7 = + gamma2*IA - psi2*k5*RA*(SC+IC+RC)+ k6*RC*(SA+IA+RA) -delta*RA - lambda6*RA + lambda5 * RC==0;
        eq8 = SC + SA + SH + IA + IC + RC + RA ==1;
    
        eq1 = subs(eq1); eq2 = subs(eq2); eq3 = subs(eq3); 
        eq4 = subs(eq4); eq5 = subs(eq5); eq6 = subs(eq6); 
        eq7 = subs(eq7); eq8= subs(eq8); 
        sol = solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7, eq8],[SC,SA,SH,RC,RA]);
        
        
        FDE = [double(sol.SC), double(sol.SH), double(sol.SA), double(sol.RC),  double(sol.RA)];
        % soluzioni(4,:)
end
  
%% Function "Scenari"
function [lambda_1,lambda_2,k1,k2,time,title_fig, beta,gamma,delta,rho,epsilon,psi, SH0,SC0,SA0,IC0,IA0,RA0,RC0] = scenario(caso)
    % Parameters equal in all situations
    beta = 0.40;
    gamma = 1/9;
    delta = 1/90;
    % rho = 0.65;
    % epsilon = 0.15;
    psi = 1;
    rho = 0.15;
    epsilon = 0.15;
    if caso == 1
        B1 = 0.89;
        B2 = 0.45;
        lambda_1 = 1/30; %fatigue to mantain C behaviour
        lambda_2 = 1/40; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        % Population initial condition
        SC0 = 50/60e6;
        SA0 = 50/60e6;
        IC0 = 10/60e6;
        IA0 = 10/60e6;
        RC0 = 0;
        RA0 = 0;
        SH0 = 1- SC0-SA0-IC0-IA0;
        time = 3000;
        title_fig = 'epi_behav_sim_B1_B2_less_1.pdf';
    end
    if caso == 2
            B1 = 8.5;
            B2 = 8.5;
            lambda_1 = 1/40; %fatigue to mantain C behaviour
            lambda_2 = 1/400; %fatigue to mantain A behaviour
            k1 = B1*lambda_1 ; %from H to C
            k2 = B2*lambda_2; %from H to A
            time = 1000;
            % Population initial condition
            SC0 = 50/60e6;
            SA0 = 50/60e6;
            IC0 = 10/60e6;
            IA0 = 10/60e6;
            RC0 = 0;
            RA0 = 0;
            SH0 = 1- SC0-SA0-IC0-IA0;
            title_fig = 'epi_behav_sim_B1_B2_equal.pdf';
    end
    if caso == 3
        B1 = 7;
        B2 = 3;
        lambda_1 = 1/30; %fatigue to mantain C behaviour
        lambda_2 = 1/30; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        %  Population initial condition
        time = 1000;
        SC0 = 50/60e6;
        SA0 = 50/60e6;
        IC0 = 0;
        IA0 = 0;
        RC0 = 0;
        RA0 = 0;
        SH0 = 1- SC0-SA0-IC0-IA0;
        title_fig = 'epi_behav_sim_B1_mag_B2.pdf';
    end
    if caso == 4
        B1 = 7;
        B2 = 3;
        lambda_1 = 1/30; %fatigue to mantain C behaviour
        lambda_2 = 1/30; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        %  Population initial condition
        time = 1000;
        SC0 = 0;
        SA0 = lambda_2/k2;
        IC0 = 0;
        IA0 = 0;
        RC0 = 0;
        RA0 = 0;
        SH0 = 1- SC0-SA0-IC0-IA0;
        title_fig = 'epi_behav_sim_B1_mag_B2_lambda2_mag.pdf';
    end
    if caso == 5
        B1 = 3;
        B2 = 7;
        lambda_1 = 1/30; %fatigue to mantain C behaviour
        lambda_2 = 1/2; %fatigue to mantain A behaviour
        k1 = B1*lambda_1 ; %from H to C
        k2 = B2*lambda_2; %from H to A
        %  Population initial condition
        time = 1000;
        SC0 = 500/60e6;
        SA0 = 500/60e6;
        IC0 = 10/60e6;
        IA0 = 10/60e6;
        RC0 = 0;
        RA0 = 0;
        SH0 = 1- SC0-SA0-IC0-IA0;
        title_fig = 'epi_behav_sim_B2_mag_B1.pdf';
    end
end
%% Function TNG
function [E_0, R_0, B_1, B_2,FV] = calcolo_E_0_Van(caso, beta_v, gamma_v,psi_v, rho_v, epsilon_v, k3_v, k4_v,lam3_v, lam4_v, SH0, SC0,SA0, IC0,IA0,delta_v,RA0,RC0 )
    syms psi2 SC SA RA RC SH IC IA rho epsilon k1 k2 k3 k4 k5 k6 lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 beta2 gamma2 
    
    if caso == 1
        % Matrix F
        F = [epsilon*beta2*(rho*SC+SH)+ psi2*k3*IA,   beta2*(rho*SC+SH)+psi2*k3*(SC+IC+RC)+lambda4;
         beta2*epsilon*SA+ lambda3+k4*(SA+IA+RA), beta2*SA+ k4*IC];
        % Matrix V
        V = [lambda3+k4*(SA+IA+RA)+gamma2, k4*IC;
        psi2*k3*IA, psi2*k3*(SC+IC+RC)+lambda4+gamma2];
    elseif caso == 2
        % Matrix F
        F = [epsilon*beta2*(rho*SC+SH)+ psi2*k3*IA,   beta2*(rho*SC+SH)+psi2*k3*(SC+IC+RC);
         beta2*epsilon*SA+k4*(SA+IA+RA), beta2*SA+ k4*IC];
        % Matrix V
        V = [lambda3+k4*(SA+IA+RA)+gamma2, k4*IC-lambda4;
        psi2*k3*IA-lambda3, psi2*k3*(SC+IC+RC)+lambda4+gamma2];
    elseif caso == 3
        % Matrix F
        F = [epsilon*beta2*(rho*SC+SH),   beta2*(rho*SC+SH);
             beta2*epsilon*SA ,            beta2*SA];
        % Matrix V
       V = [lambda3+k4*(SA+IA+RA)+gamma2-psi2*k3*IA, k4*IC-lambda4-psi2*k3*(SC+IC+RC);
           +psi2*k3*IA-lambda3- k4*(SA+IA+RA),  psi2*k3*(SC+IC+RC)-k4*IC+lambda4+gamma2];
       
        % V = [lambda3+k4*(SA+IA+RA)+gamma2-psi2*k3*IA,  k4*IC-lambda4-psi2*k3*(SC+IC+RC);
        %     +psi2*k3*IA-lambda3- k4*(SA+IA+RA),           psi2*k3*(SC+IC+RC)+lambda4+gamma2 -k4*IC]
      
    end
    % Substitute symbols with variables
    psi2 =psi_v;  rho=rho_v; epsilon=epsilon_v; 
    k3=k3_v; k4=k4_v; lambda3=lam3_v; lambda4=lam4_v; beta2=beta_v; gamma2=gamma_v;
    SC = SC0; SH = SH0; SA = SA0; IC = IC0; IA = IA0; RC = RC0; RA = RA0;
    % Compute the Matrix F/V
    FV = F/V;
    FV_val = subs(FV);
    FV_val = double(FV_val);
    % Compute the largest magnitude eigenvalue
    E_0 = eigs(FV_val, 1, 'largestabs');
    % Compute also the others epidemic numbers
    R_0 = beta_v/(gamma_v+delta_v);
    B_1 = k3_v/lam3_v;
    B_2 = k4_v/lam4_v;
end

