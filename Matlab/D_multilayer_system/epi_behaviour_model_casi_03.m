clc;
close all;
%% Epidemiological Behavioural model

caso = 4;

[lambda_1,lambda_2,k1,k2, B1, B2, C_zero,A_zero, time,title_fig] = scenario(caso);

% Multi layer network implementation
beta = 0.40;  % infectivity rate
gamma = 1/9; % rate from I to R, 9 days
delta = 1/90; % rate from R to S, 90 days 

k3 = k1; k5 = k1; 
k4 = k2; k6 = k2;
lambda_3 = lambda_1;  lambda_5 = lambda_1;
lambda_4 = lambda_2;  lambda_6 = lambda_2;
R_i = [beta/gamma, k1/lambda_1, k2/lambda_2] % the basic reproduction number of disease and each opinione SEPARATELY!
epsilon = 0.15; % il 15 percento della popolazione compliant infetta partecipa all'infezione 
omega = 0.005; % non viene usato

[s_h, s_c, s_a, i_c, i_a, r_c,r_a, time_epi] = epi_behaviour(beta,gamma,delta,k1,k2,k3,k4,k5,k6,lambda_1,lambda_2,lambda_3,lambda_4,lambda_5,lambda_6,epsilon,omega,C_zero,A_zero,time);

figure(1)
hold on
box on
plot(time_epi,s_h, 'linewidth',1.3 )
plot(time_epi,s_c, 'linewidth',1.3 )
plot(time_epi,s_a, 'linewidth',1.3 )
legend('SC', 'SH', 'SA', Orientation='horizontal', Location='southoutside')
title('Susceptibles compartments dynamic')
xlabel("t[days]");
ylabel("SC[t], SH[t], SA[t]");
txt = {['R_0=' num2str(R_i(1))],['B_1=' num2str(B1)], ['B_2=' num2str(B2)], ['\beta=' num2str(round(beta,2))],['\gamma=' num2str(round(gamma,2))],['k_1=' num2str(round(k1,2))],['k_2=' num2str(round(k2,2))],['\lambda_1=' num2str(round(lambda_1,2))],['\lambda_2=' num2str(round(lambda_2,2))]};    
dim = [.93 .85 .1 .1];
annotation('textbox',dim, ...
    'String',txt,'EdgeColor','none')
hold off
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 22 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', "susceptible_"+title_fig)
%
figure(2)
hold on
box on
plot(time_epi,i_c, 'linewidth',1.3, 'Color',[0.8500 0.3250 0.0980])
plot(time_epi,i_a, 'linewidth',1.3, 'Color',[0.9290 0.6940 0.1250] )
legend('IC', 'IA', Orientation='horizontal', Location='southoutside')
title('Infected compartments dynamic')
xlabel("t[days]");
ylabel("IC[t], IA[t]");
txt = {['R_0=' num2str(R_i(1))],['B_1=' num2str(B1)], ['B_2=' num2str(B2)], ['\beta=' num2str(round(beta,2))],['\gamma=' num2str(round(gamma,2))],['k_1=' num2str(round(k1,2))],['k_2=' num2str(round(k2,2))],['\lambda_1=' num2str(round(lambda_1,2))],['\lambda_2=' num2str(round(lambda_2,2))]};    
dim = [.93 .85 .1 .1];
annotation('textbox',dim, ...
    'String',txt,'EdgeColor','none')
hold off
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 22 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', "infected_"+title_fig)
%
figure(3)
hold on
box on
plot(time_epi,r_c, 'linewidth',1.3, 'Color',[0.8500 0.3250 0.0980] )
plot(time_epi,r_a, 'linewidth',1.3, 'Color',[0.9290 0.6940 0.1250]  )
legend('RC', 'RA', Orientation='horizontal', Location='southoutside')
title('Recovered compartments dynamic')
xlabel("t[days]");
ylabel("RC[t], RA[t]");
txt = {['R_0=' num2str(R_i(1))],['B_1=' num2str(B1)], ['B_2=' num2str(B2)], ['\beta=' num2str(round(beta,2))],['\gamma=' num2str(round(gamma,2))],['k_1=' num2str(round(k1,2))],['k_2=' num2str(round(k2,2))],['\lambda_1=' num2str(round(lambda_1,2))],['\lambda_2=' num2str(round(lambda_2,2))]};    
dim = [.93 .85 .1 .1];
annotation('textbox',dim, ...
    'String',txt,'EdgeColor','none')
hold off
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 22 15]);
set(gcf, 'PaperSize', [24 15]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', "recovered_"+title_fig)



%% Function section 
function [x_sh, x_sc, x_sa, y_c, y_a, z_rc,z_ra, time_epi] = epi_behaviour(beta,gamma,delta,k1,k2,k3,k4,k5,k6,lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,epsilon,omega,C_zero,A_zero, time)
    a1 = 1/2; a2=1/2; p1 =1; q11 =1; % Heuns method
    dt=0.01; %un centesimo di secondo per dt è ottimo con runge kutta
    % Other model parameters 
    rho = 0.65; % grado di protezione asociato all'essere compliant
    gamma_le = gamma;
    mu = 0.2; % mortality? !!!!
    s_c = C_zero/60e6;    % susceptible compliant
    s_a = A_zero/60e6;    % susceptible against
    i_c  = 10/60e6;    % infected compliant and careless
    i_a = 10/60e6;    % infected against
    r_c  = 0;          % recovered compliant and careless 
    r_a = 0;          % recovered against
    s_h = 1- s_c - s_a - i_c -i_a; % susceptible careless
    t = 0; % initialize the time counter
    cnt=0;
    %Array creation and inititialization
    time_epi=[]; time_epi(1) =0; 
    x_sh   =[]; x_sh(1)    = s_h;
    x_sc   =[]; x_sc(1)    = s_c;
    x_sa   =[]; x_sa(1)  = s_a;
    y_c=[];   y_c(1) = i_c;
    y_a=[];  y_a(1)= i_a;
    z_rc = []; z_rc(1) = r_c;
    z_ra=[];  z_ra(1)= r_a; 
    while t < time
        %every 100 millisecond I save the result
        if mod(cnt,100) == 0 && cnt ~=0 
            time_epi = cat(2,time_epi,t);
            x_sh = cat(2, x_sh, s_h);
            x_sc = cat(2, x_sc, s_c);
            x_sa = cat(2, x_sa, s_a); 
            y_c= cat(2, y_c, i_c);
            y_a= cat(2, y_a, i_a);
            z_rc= cat(2, z_rc, r_c);
            z_ra= cat(2, z_ra, r_a);
        end
        psi = 1;  % Awareness / exposure / mean field forcing 
        phi_n = psi; %deve diventare awareness normalizzata, cioè tra 0 e 1.
        if phi_n > 1 && phi_n < 0
            print("Error")
        end
        % step 1
        im = epsilon*i_c+i_a;  %group that participate in the infection
        % system equation definition
        ksh1 = - beta*s_h*im - k1*psi*s_h*(s_c+i_c+r_c) -k2*s_h*(s_a+i_a+r_a) +lambda1*s_c + lambda2*s_a + delta*(1-phi_n)*r_c;
        ksc1 = k1*psi*s_h*(s_c+i_c+r_c) - lambda1*s_c -beta*rho*s_c*im + delta*phi_n*r_c;
        ksa1 = k2*s_h*(s_a+i_a+r_a) - lambda2*s_a - beta*s_a*im + delta*r_a;
        kic1 = beta*rho*s_c*im + beta*s_h*im - k4 *i_c*(s_a+i_a+r_a) + k3*psi*i_a*(s_c+i_c+r_c) - lambda3*i_c- gamma*i_c - lambda4*i_a;
        kia1 = +beta*s_a*im + k4*i_c*(s_a+i_a+r_a)+lambda3*i_c - k3*psi*i_a*(s_c+i_c+r_c)-gamma*i_a - lambda4*i_a;
        krc1 = gamma * i_c - k6*r_c*(s_a+i_a+r_a) + psi*k5*r_a*(s_c+i_c+r_c) - delta*r_c + lambda6*r_a - lambda5 * r_c;
        kra1 = + gamma*i_a - psi*k5*r_a*(s_c+i_c+r_c)+ k6*r_c*(s_a+i_a+r_a) -delta*r_a - lambda6*r_a + lambda5 * r_c;
        % step 2
        t2 = t+p1*dt;
        s_h2 = s_h + q11*ksh1*dt;
        s_c2 = s_c + q11*ksc1*dt;
        s_a2 = s_a + q11*ksa1*dt;
        i_c2 = i_c   + q11*kic1*dt;
        i_a2 = i_a + q11*kia1*dt;
        r_c2 = r_c + q11*krc1*dt;
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
        s_h = s_h + (a1*ksh1+a2*ksh2)*dt;
        s_c = s_c + (a1*ksc1+a2*ksc2)*dt;
        s_a = s_a + (a1*ksa1+a2*ksa2)*dt;
        i_c = i_c + (a1*kic1+a2*kic2)*dt;
        i_a = i_a + (a1*kia1+a2*kia2)*dt;
        r_c = r_c + (a1*krc1+a2*krc2)*dt;
        r_a = r_a + (a1*kra1+a2*kra2)*dt;
        if r_c <= 0
            r_c = 0;
        end
        if r_a <= 0
            r_a = 0;
        end
        t = t + dt;
        cnt = cnt + 1;
    end    
end

%% Function to change the initial conditions
function [lambda_1,lambda_2,k1,k2, B1, B2, C_zero,A_zero, time,title_fig] = scenario(caso)
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
    time = 10000;
    title_fig = 'epi_behav_sim_B1_B2_less_1.pdf';
    end
    if caso == 2
            B1 = 8.5;
            B2 = 8.5;
            lambda_1 = 1/25; %fatigue to mantain C behaviour
            lambda_2 = 1/30; %fatigue to mantain A behaviour
            k1 = B1*lambda_1 ; %from H to C
            k2 = B2*lambda_2; %from H to A
            time = 600;
            % Population initial condition
            C_zero = 100;
            A_zero = 100;     
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
        time = 500;
        C_zero = 100;
        A_zero = 100;     
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
        time = 500;
        C_zero = 100;
        A_zero = 100;
        title_fig = 'epi_behav_sim_B1_mag_B2_lambda2_mag.pdf';
    end

end
