% clc;
% clear;
close all;
fig = 0;
%% Set the scenario; first number scenario, second what F/V matrix compute
caso =[1,3];
%% Run the code to se R_0, B1, B2 and the three resulting E_0 values
[lambda_1,lambda_2,k_1,k_2,time,title_fig, beta,gamma,delta_v,rho_v,epsilon_v,psi_v,phi_n]=scenario(caso(1));
% [fig, E_0_matrix_1] = heat_map_E0_1(fig, caso(2),beta, gamma,psi_v, rho_v,delta_v,phi_n);


%% Plot the heat map
lambda_1 = linspace(1/2,1/50,10);
lambda_2 = linspace(1/2,1/50,10);
epsilon_v = linspace(0,1,10);

titles = "E0_heatmap_1_lam1_B1_epsilon.pdf";
[xx, yy,zz] = ndgrid(lambda_1,lambda_2,epsilon_v);

    aga = E_0_matrix_1;
    FF = aga; %FF = transpose(FF);
    
    hold on
    for i = 2:2:length(lambda_1)
        fig = fig+1;
        figure(fig)
        surf(xx(:,:,i),yy(:,:,i),zz(:,:,i),FF(:,:,i),'EdgeColor','none') ;
        view(45,45)
        colorbar;
        colormap default;
        xlabel('\lambda_1') 
        ylabel('\lambda_2')
        zlabel('\epsilon ')
    end
    % fontsize(20,"points")
    % set(gcf, 'PaperUnits', 'centimeters');
    % set(gcf, 'PaperPosition', [0 0 25 16]);
    % set(gcf, 'PaperSize', [26 16]); % dimension on x axis and y axis resp.
    %  print(gcf,titles,'-dpdf')


%% Function Disease Free Equilibrium
function FDE = equilibri_FDE(psi_v,rho_v,epsilon_v, k_1, k_2, lambda_1, lambda_2, beta, gamma, delta_v,phi_n_val)
    syms psi2 SC SA RA RC SH IC IA rho epsilon k1 k2 k3 k4 k5 k6 lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 beta2 gamma2 SC0 SH0 SA0 RC0 RA0 delta phi_n   
    psi2 =psi_v;  rho=rho_v; epsilon=epsilon_v; 
    k1 = k_1; k2 = k_2; k3=k_1; k4=k_2; k5 = k_1; k6 = k_2;
    lambda1 = lambda_1; lambda2= lambda_2; lambda3=lambda_1; lambda4=lambda_2; lambda5 = lambda_1; lambda6 =lambda_2;
    beta2=beta; gamma2=gamma; delta = delta_v; phi_n = phi_n_val;
    % DFE Conditions
    IA = 0; IC = 0;
    % Equations to solve
    eq1 = - beta2*SH*(epsilon*IC+IA) - k1*psi2*SH*(SC+IC+RC) -k2*SH*(SA+IA+RA) +lambda1*SC + lambda2*SA + delta*(1-phi_n)*RC == 0;
    eq2 = k1*psi2*SH*(SC+IC+RC) - lambda1*SC -beta2*rho*SC*(epsilon*IC+IA) + delta*phi_n*RC==0;
    eq3 = k2*SH*(SA+IA+RA) - lambda2*SA - beta2*SA*(epsilon*IC+IA) + delta*RA==0;
    eq4 = beta2*rho*SC*(epsilon*IC+IA) + beta2*SH*(epsilon*IC+IA) - k4 *IC*(SA+IA+RA) + k3*psi2*IA*(SC+IC+RC) - lambda3*IC- gamma2*IC + lambda4*IA==0;
    eq5 = +beta2*SA*(epsilon*IC+IA) + k4*IC*(SA+IA+RA)+lambda3*IC - k3*psi2*IA*(SC+IC+RC)-gamma2*IA - lambda4*IA==0;
    eq6 = gamma2 * IC - k6*RC*(SA+IA+RA) + psi2*k5*RA*(SC+IC+RC) - delta*RC + lambda6*RA - lambda5 * RC==0;
    eq7 = + gamma2*IA - psi2*k5*RA*(SC+IC+RC)+ k6*RC*(SA+IA+RA) -delta*RA - lambda6*RA + lambda5 * RC==0;
    eq8 = SC + SA + SH + IA + IC + RC + RA ==1;
    % Substitute and solve the system
    eq1 = subs(eq1); eq2 = subs(eq2); eq3 = subs(eq3); 
    eq4 = subs(eq4); eq5 = subs(eq5); eq6 = subs(eq6); 
    eq7 = subs(eq7); eq8= subs(eq8); 
    sol = solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7, eq8],[SC,SA,SH,RC,RA]);
    % Creat a matrix with the DFE solutions
    FDE = [double(sol.SC), double(sol.SH), double(sol.SA), double(sol.RC),  double(sol.RA)];
        
end
  
%% Function "Scenari"
function [lambda_1,lambda_2,k1,k2,time,title_fig, beta,gamma,delta,rho,epsilon,psi,phi_n] = scenario(caso)
    % Parameters equal in all situations
    beta = 0.40;
    gamma = 1/9;
    delta = 1/90;
    psi = 1;
    rho = 0.65;
    epsilon = 0.15;
    phi_n = 0.5;
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
            lambda_2 = 1/20; %fatigue to mantain A behaviour
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
        lambda_2 = 1/5; %fatigue to mantain A behaviour
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
%% Function TNG Matrix
function [E_0, R_0, B_1, B_2] = calcolo_E_0_Van(caso, beta_v, gamma_v,psi_v, rho_v, epsilon_v, k3_v, k4_v,lam3_v, lam4_v, SH0, SC0,SA0, IC0,IA0,delta_v,RA0,RC0 )
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

function [fig, E_0] = heat_map_E0_1(fig, caso_2,beta, gamma,psi_v, rho_v,delta_v,phi_n)
    
    cnt = 0; progress = 0;
    lambda_1 = linspace(1/2,1/50,10);
    lambda_2 = linspace(1/2,1/50,10);
    epsilon_v = linspace(0,1,10);
    totale_it = (length(lambda_1))^3/100;
    E_0 = zeros(length(lambda_1), length(lambda_2), length(epsilon_v));
    for i = 1: length(lambda_1)
        for j = 1: length(lambda_2)
            for k = 1: length(epsilon_v)
                % Define parameters
                B1 = 5; k_1 = B1*lambda_1(i);
                B2 = 5; k_2 = B2*lambda_2(j);
                % Set the DFE
                SC_0 = 0; SA_0 = 0; SH_0 = 1; IC_0 =0; IA_0 = 0; RC_0 = 0; RA_0 = 0;
                % Compute E_0
                [E_0(i,j,k), R_0, B_1, B_2] = calcolo_E_0_Van(caso_2,beta, gamma,psi_v, rho_v, epsilon_v(k), k_1, k_2,lambda_1(i), lambda_2(j), SH_0, SC_0,SA_0, IC_0,IA_0,delta_v,RA_0,RC_0);
                cnt = cnt+1;
                if mod(cnt, totale_it) == 0
                    progress = progress +1
                end
            end
        end
    end
    % % Plot the heat map
    % fig = fig+1;
    % figure(fig)
    % titles = "E0_heatmap_1_lam1_B1_epsilon.pdf";
    % [xx, yy,zz] = ndgrid(lambda_1,B1,epsilon_v);
    % 
    % aga = E_0;
    % FF = aga; %FF = transpose(FF);
    % 
    % hold on
    % for i = 1:length(lambda_1)
    %     surf(xx(:,:,i),yy(:,:,i),zz(:,:,i),FF(:,:,i),'EdgeColor','none') ;
    % end
    % view(45,45)
    % colorbar;
    % colormap default;
    % xlabel('\lambda_1') 
    % ylabel('B1')
    % zlabel('\epsilon ')
    % fontsize(20,"points")
    % set(gcf, 'PaperUnits', 'centimeters');
    % set(gcf, 'PaperPosition', [0 0 25 16]);
    % set(gcf, 'PaperSize', [26 16]); % dimension on x axis and y axis resp.
    %  print(gcf,titles,'-dpdf')

end
%% Function Simulation and plots
function fig = epi_behaviour(fig,beta,gamma,delta,rho,psi,k1,k2,k3,k4,k5,k6,lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,epsilon,phi_n,C_zero,A_zero,H_zero, time, E_0,title_fig)
    a1 = 1/2; a2=1/2; p1 =1; q11 =1; % Heuns method
    dt=0.01; %un centesimo di secondo per dt è ottimo con runge kutta
    
    s_c = C_zero;    % susceptible compliant
    s_a = A_zero;    % susceptible against
    i_c  = 10/60e6;    % infected compliant and careless
    i_a = 10/60e6;    % infected against
    r_c  = 0;          % recovered compliant and careless 
    r_a = 0;          % recovered against
    s_h = H_zero; % susceptible careless
    t = 0; % initialize the time counter
    cnt=0;
    R_0 = beta/(gamma+delta);
    B1 = k1/lambda1;
    B2= k2/lambda2;
    %Array creation and inititialization
    time_epi=[]; time_epi(1) =0; 
    x_sh   =[]; x_sh(1)    = s_h; s_h3 = s_h;
    x_sc   =[]; x_sc(1)    = s_c; s_c3 = s_c;
    x_sa   =[]; x_sa(1)  = s_a;   s_a3 = s_a;
    y_c=[];   y_c(1) = i_c;       i_c3 = i_c;
    y_a=[];  y_a(1)= i_a;         i_a3 = i_a;
    z_rc = []; z_rc(1) = r_c;     r_c3 = r_c;
    z_ra=[];  z_ra(1)= r_a;       r_a3 = r_a;
    while t < time
        %every 100 millisecond I save the result
        if mod(cnt,100) == 0 && cnt ~=0 
            time_epi = cat(2,time_epi,t);
            x_sh = cat(2, x_sh, s_h3);
            x_sc = cat(2, x_sc, s_c3);
            x_sa = cat(2, x_sa, s_a3); 
            y_c= cat(2, y_c, i_c3);
            y_a= cat(2, y_a, i_a3);
            z_rc= cat(2, z_rc, r_c3);
            z_ra= cat(2, z_ra, r_a3);
        end
        
        % step 1
        im = epsilon*i_c3+i_a3;  %group that participate in the infection
        % system equation definition
        ksh1 = - beta*s_h3*im - k1*psi*s_h3*(s_c3+i_c3+r_c3) -k2*s_h3*(s_a3+i_a3+r_a3) +lambda1*s_c3 + lambda2*s_a3 + delta*(1-phi_n)*r_c3;
        ksc1 = k1*psi*s_h3*(s_c3+i_c3+r_c3) - lambda1*s_c3 -beta*rho*s_c3*im + delta*phi_n*r_c3;
        ksa1 = k2*s_h3*(s_a3+i_a3+r_a3) - lambda2*s_a3 - beta*s_a3*im + delta*r_a3;
        kic1 = beta*rho*s_c3*im + beta*s_h3*im - k4 *i_c3*(s_a3+i_a3+r_a3) + k3*psi*i_a3*(s_c3+i_c3+r_c3) - lambda3*i_c3- gamma*i_c3 + lambda4*i_a3;
        kia1 = +beta*s_a3*im + k4*i_c3*(s_a3+i_a3+r_a3)+lambda3*i_c3 - k3*psi*i_a3*(s_c3+i_c3+r_c3)-gamma*i_a3 - lambda4*i_a3;
        krc1 = gamma * i_c3 - k6*r_c3*(s_a3+i_a3+r_a3) + psi*k5*r_a3*(s_c3+i_c3+r_c3) - delta*r_c3 + lambda6*r_a3 - lambda5 * r_c3;
        kra1 = + gamma*i_a3 - psi*k5*r_a3*(s_c3+i_c3+r_c3)+ k6*r_c3*(s_a3+i_a3+r_a3) -delta*r_a3 - lambda6*r_a3 + lambda5 * r_c3;
        % step 2
        t2 = t+p1*dt;
        s_h2 = s_h3 + q11*ksh1*dt;
        s_c2 = s_c3 + q11*ksc1*dt;
        s_a2 = s_a3 + q11*ksa1*dt;
        i_c2 = i_c3 + q11*kic1*dt;
        i_a2 = i_a3 + q11*kia1*dt;
        r_c2 = r_c3 + q11*krc1*dt;
            if r_c2 <= 0
                r_c2 = 0;
            end
        r_a2 = r_a + q11*kra1*dt;
            if r_a2 <= 0
                r_a2 = 0;
            end
        im = epsilon*i_c2 + i_a2;
        %
        ksh2 = - beta*s_h2*im - k1*psi*s_h2*(s_c2+i_c2+r_c2) -k2*s_h2*(s_a2+i_a2+r_a2) +lambda1*s_c2 + lambda2*s_a2 + delta*(1-phi_n)*r_c2;
        ksc2 = k1*psi*s_h2*(s_c2+i_c2+r_c2) -lambda1*s_c2 -beta*rho*s_c2*im + delta*phi_n * r_c2 ;
        ksa2 = k2*s_h2*(s_a2+i_a2+r_a2) - lambda2*s_a2 - beta*s_a2*im + delta*r_a2;
        kic2 = beta*rho*s_c2*im + beta*s_h2*im - k4 *i_c2*(s_a2+i_a2+r_a2) + k3*psi*i_a2*(s_c2+i_c2+r_c2) - lambda3*i_c2- gamma*i_c2 + lambda4*i_a2;
        kia2 = +beta*s_a2*im + k4*i_c2*(s_a2+i_a2+r_a2)+lambda3*i_c2 - k3*psi*i_a2*(s_c2+i_c2+r_c2)-gamma*i_a2 - lambda4*i_a2;
        krc2 = gamma * i_c2 - k6*r_c2*(s_a2+i_a2+r_a2) + psi*k5*r_a2*(s_c2+i_c2+r_c2) - delta*r_c2 - lambda5*r_c2 + lambda6*r_a2;
        kra2 = + gamma*i_a2 - psi*k5*r_a2*(s_c2+i_c2+r_c2)+ k6*r_c2*(s_a2+i_a2+r_a2) - delta*r_a2 - lambda6 * r_a2 + lambda5 * r_c2;
        % update
        s_h3 = s_h3 + (a1*ksh1+a2*ksh2)*dt;
        s_c3 = s_c3 + (a1*ksc1+a2*ksc2)*dt;
        s_a3 = s_a3 + (a1*ksa1+a2*ksa2)*dt;
        i_c3 = i_c3 + (a1*kic1+a2*kic2)*dt;
        i_a3 = i_a3 + (a1*kia1+a2*kia2)*dt;
        r_c3 = r_c3 + (a1*krc1+a2*krc2)*dt;
        r_a3 = r_a3 + (a1*kra1+a2*kra2)*dt;
        if r_c3 <= 0
            r_c3 = 0;
        end
        if r_a3 <= 0
            r_a3 = 0;
        end
        t = t + dt;
        cnt = cnt + 1;
    end  
fig = fig+1;
figure(fig)
tiledlayout(3,1);
% Tile 1
nexttile
hold on
box on
plot(time_epi,x_sh, 'linewidth',2, 'Color',[0.6350 0.0780 0.1840] )
plot(time_epi,x_sc, 'linewidth',2,'Color',[0 0.4470 0.7410] )
plot(time_epi,x_sa, 'linewidth',2,'Color',[0.9290 0.6940 0.1250] )
plot(time_epi, x_sh+x_sc+x_sa, '--','linewidth',2, 'Color',[0.4940 0.1840 0.5560] )
% plot(time_epi, xs_h+ xs_c+ xs_a+yi_c+ yi_a+zr_c+zr_a, 'linewidth',2)
legend('SH', 'SC', 'SA','Total S', Orientation='horizontal', Location='southoutside')
title("Susceptibles compartments")
xlabel("t[days]");
ylabel("SC[t], SH[t], SA[t]");
% txt = {['R_0=' num2str(R_0)],['B_1=' num2str(B1)], ['B_2=' num2str(B2)], ['\beta=' num2str(round(beta,2))],['\gamma=' num2str(round(gamma,2))],['k_1=' num2str(round(k1,2))],['k_2=' num2str(round(k2,2))],['\lambda_1=' num2str(round(lambda1,2))],['\lambda_2=' num2str(round(lambda2,2))]};    
% dim = [.93 .85 .1 .1];
% annotation('textbox',dim, ...
%     'String',txt,'EdgeColor','none')
% hold off
% fontsize(20,"points")
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 22 17]);
% set(gcf, 'PaperSize', [24 17]); % dimension on x axis and y axis resp.
% % print(gcf,'-dpdf', "susceptible_"+title_fig)
% Tile 2
nexttile
hold on
box on
plot(time_epi,y_c, 'linewidth',2, 'Color',[0 0.4470 0.7410])
plot(time_epi,y_a, 'linewidth',2, 'Color',[0.9290 0.6940 0.1250] )
plot(time_epi, y_c+y_a, '--','linewidth',2, 'Color',[0.4940 0.1840 0.5560] )
legend('IC', 'IA','Total I', Orientation='horizontal', Location='southoutside')
title("Infected compartments")
xlabel("t[days]");
ylabel("IC[t], IA[t]");
% txt = {['R_0=' num2str(R_0)],['B_1=' num2str(B1)], ['B_2=' num2str(B2)], ['\beta=' num2str(round(beta,2))],['\gamma=' num2str(round(gamma,2))],['k_1=' num2str(round(k1,2))],['k_2=' num2str(round(k2,2))],['\lambda_1=' num2str(round(lambda1,2))],['\lambda_2=' num2str(round(lambda2,2))]};    
% dim = [.93 .85 .1 .1];
% annotation('textbox',dim, ...
%     'String',txt,'EdgeColor','none')
% hold off
% fontsize(20,"points")
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 22 17]);
% set(gcf, 'PaperSize', [24 17]); % dimension on x axis and y axis resp.
% % print(gcf,'-dpdf', "infected_"+title_fig)

% Tile 1
nexttile
hold on
box on
plot(time_epi,z_rc, 'linewidth',2, 'Color',[0 0.4470 0.7410] )
plot(time_epi,z_ra, 'linewidth',2, 'Color',[0.9290 0.6940 0.1250] )
plot(time_epi,z_rc+z_ra , '--','linewidth',2, 'Color',[0.4940 0.1840 0.5560] )
legend('RC', 'RA','Total R', Orientation='horizontal', Location='southoutside')
title("Recovered compartments, E_0")
xlabel("t[days]");
ylabel("RC[t], RA[t]");
txt = {['E_0=' num2str(E_0,3)],['R_0=' num2str(R_0,3)],['B_1=' num2str(B1)], ['B_2=' num2str(B2)],...
    ['\beta=' num2str(round(beta,2))],['\gamma=' num2str(round(gamma,2))],['\delta=' num2str(round(delta,2))],['k_1=' num2str(round(k1,2))],['k_2=' num2str(round(k2,2))],...
    ['\lambda_1=' num2str(round(lambda1,2))],['\lambda_2=' num2str(round(lambda2,2))],...
    ['S_{C0}=' num2str(C_zero,3)], ['S_{H0}=' num2str(H_zero,3)], ['S_{A0}=' num2str(A_zero,3)]};    
dim = [.95 .85 .1 .1];
annotation('textbox',dim, ...
    'String',txt,'EdgeColor','none')
hold off
fontsize(26,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 40]);
set(gcf, 'PaperSize', [28 40]); % dimension on x axis and y axis resp.
 print(gcf,'-dpdf', "E_0_model_"+fig+title_fig)
end
