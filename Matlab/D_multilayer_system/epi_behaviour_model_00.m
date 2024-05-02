clc;
close all;
%% Epidemiological Behavioural model
% Multi layer network implementation
beta = 0.53;
gamma = 0.35;
delta = 0.45;
k1 = 0.48;
k2 = 0.243;
lambda1 = 0.143;
lambda2 = 0.143;
R_i = [beta/gamma, k1/lambda1, k2/lambda2]
epsilon = 0.66; 
omega = 0.005;

time = 1550; 
[s_ca, s_co, s_ag, i_c, i_ag, r_c,r_ag,d, time_epi] = epi_behaviour(beta,gamma,delta,k1,k2,lambda1,lambda2,epsilon,omega,time);

figure(1)
hold on
plot(time_epi,s_co, 'linewidth',1.3 )
plot(time_epi,s_ca, 'linewidth',1.3 )
plot(time_epi,s_ag, 'linewidth',1.3 )
legend('s_co', 's_{ca}', 's_{ag}')
hold off
%
figure(2)
hold on
plot(time_epi,i_c, 'linewidth',1.3 )
plot(time_epi,i_ag, 'linewidth',1.3 )
legend('i_c','i_{ag}')
%
figure(3)
hold on
plot(time_epi,r_c, 'linewidth',1.3 )
plot(time_epi,r_ag, 'linewidth',1.3 )
legend('r_c', 'r_{ag}')
hold off



%% Function section 
function [x_sca, x_sco, x_sag, y_c, y_ag, z_rc,z_rag,z_d, time_epi] = epi_behaviour(beta,gamma,delta,k1,k2,lambda1,lambda2,epsilon,omega,time)
    a1 = 1/2; a2=1/2; p1 =1; q11 =1; % Heuns method
    dt=0.01; %un centesimo di secondo per dt è ottimo con runge kutta
    
    s_co = 50/60e6;    % susceptible compliant
    s_ag = 50/60e6;    % susceptible against
    i_c  = 10/60e6;    % infected compliant and careless
    i_ag = 10/60e6;    % infected against
    r_c  = 0;          % recovered compliant and careless 
    r_ag = 0;          % recovered against
    d = 0;             % deaths
    s_ca = 1- s_co - s_ag - i_c -i_ag; % susceptible careless
    t = 0; % initialize the time counter
    cnt=0;
    %Array creation and inititialization
    time_epi=[]; time_epi(1) =0; 
    x_sca   =[]; x_sca(1)    = s_ca;
    x_sco   =[]; x_sco(1)    = s_co;
    x_sag   =[]; x_sag(1)  = s_ag;
    y_c=[];   y_c(1) = i_c;
    y_ag=[];  y_ag(1)= i_ag;
    z_rc = []; z_rc(1) = r_c;
    z_rag=[];  z_rag(1)= r_ag; 
    z_d = [];  z_d(1) = d;
    while t < time
        %every 100 millisecond I save the result
        if mod(cnt,100) == 0 && cnt ~=0 
            time_epi = cat(2,time_epi,t);
            x_sca = cat(2, x_sca, s_ca);
            x_sco = cat(2, x_sco, s_co);
            x_sag = cat(2, x_sag, s_ag); 
            y_c= cat(2, y_c, i_c);
            y_ag= cat(2, y_ag, i_ag);
            z_rc= cat(2, z_rc, r_c);
            z_rag= cat(2, z_rag, r_ag);
            z_d = cat(2,z_d, d);
        end
       
        Aw = 0.5;  % Awareness / exposure / mean field forcing 
        rho = 0.7; % grado di protezione asociato all'essere compliant
        k6 = k1; 
        k7 = k2;
        k11 = k1;
        k10 = k2;
        lambda3 = lambda2;
        gamma_le = gamma;
        delta = 0.2; % mortality? !!!!
        % step 1
        im = (1-epsilon)*i_c+i_ag;  %group that participate in the infection
        
        ksca1 = - beta*s_ca*im - k1*Aw*s_ca*(s_co+i_c+r_c) -k2*s_ca*(s_ag+i_ag+r_ag) +lambda1*s_co + lambda2*s_ag + omega*r_c;
        ksco1 = k1*Aw*s_ca*(s_co+i_c+r_c) -lambda1*s_co -beta*rho*s_co*im ;
        ksag1 = k2*s_ca*(s_ag+i_ag+r_ag) - lambda2*s_ag - beta*s_ag*im + omega*r_ag;
        kic1 = beta*rho*s_co*im + beta*s_ca*im - k6 *i_c*(s_ag+i_ag+r_ag) + k7*Aw*i_ag*(s_co+i_c+r_c) - lambda3*i_c- gamma*i_c;
        kiag1 = +beta*s_ag*im + k6*i_c*(s_ag+i_ag+r_ag)+lambda3*i_c - k7*Aw*i_ag*(s_co+i_c+r_c)-gamma_le*i_ag;
        krc1 = gamma * i_c - k11*r_c*(s_ag+i_ag+r_ag) + Aw*k10*r_ag*(s_co+i_c+r_c) - delta*r_c-omega*r_c;
        krag1 = + gamma_le*i_ag - Aw*k10*r_ag*(s_co+i_c+r_c)+ k11*r_c*(s_ag+i_ag+r_ag) -delta*r_ag - omega*r_ag;
        kd1 = delta*(r_c+r_ag);

        % step 2
        t2 = t+p1*dt;
        s_ca2 = s_ca + q11*ksca1*dt;
        s_co2 = s_co + q11*ksco1*dt;
        s_ag2 = s_ag + q11*ksag1*dt;
        i_c2 = i_c   + q11*kic1*dt;
        i_ag2 = i_ag + q11*kiag1*dt;
        r_c2 = r_c + q11*krc1*dt;
        r_ag2 = r_ag + q11*krag1*dt;
        d2= dt + q11*kd1*dt;
        im = (1-epsilon)*i_c2+i_ag2;
        %
        ksca2 = - beta*s_ca2*im - k1*Aw*s_ca2*(s_co2+i_c2+r_c2) -k2*s_ca2*(s_ag2+i_ag2+r_ag2) +lambda1*s_co2 + lambda2*s_ag2 + omega*r_c2;
        ksco2 = k1*Aw*s_ca2*(s_co2+i_c2+r_c2) -lambda1*s_co2 -beta*rho*s_co2*im ;
        ksag2 = k2*s_ca2*(s_ag2+i_ag2+r_ag2) - lambda2*s_ag2 - beta*s_ag2*im + omega*r_ag2;
        kic2 = beta*rho*s_co2*im + beta*s_ca2*im - k6 *i_c2*(s_ag2+i_ag2+r_ag2) + k7*Aw*i_ag2*(s_co2+i_c2+r_c2) - lambda3*i_c2- gamma*i_c2;
        kiag2 = +beta*s_ag2*im + k6*i_c2*(s_ag2+i_ag2+r_ag2)+lambda3*i_c2 - k7*Aw*i_ag2*(s_co2+i_c2+r_c2)-gamma_le*i_ag2;
        krc2 = gamma * i_c2 - k11*r_c2*(s_ag2+i_ag2+r_ag2) + Aw*k10*r_ag2*(s_co2+i_c2+r_c2) - delta*r_c2-omega*r_c2;
        krag2 = + gamma_le*i_ag2 - Aw*k10*r_ag2*(s_co2+i_c2+r_c2)+ k11*r_c2*(s_ag2+i_ag2+r_ag2) -delta*r_ag2 - omega*r_ag2;
        kd2 = delta*(r_c2+r_ag2);
        % update
        s_ca = s_ca + (a1*ksca1+a2*ksca2)*dt;
        s_co = s_co + (a1*ksco1+a2*ksco2)*dt;
        s_ag = s_ag + (a1*ksag1+a2*ksag2)*dt;
        i_c = i_c + (a1*kic1+a2*kic2)*dt;
        i_ag = i_ag + (a1*kiag1+a2*kiag2)*dt;
        r_c = r_c + (a1*krc1+a2*krc2)*dt;
        r_ag = r_ag + (a1*krag1+a2*krag2)*dt;
        d = d + (a1*kd1+a2*kd2)*dt;
        t = t + dt;
        cnt = cnt + 1;
    end    
end
