clc;
close all;
%% Epidemiological Behavioural model
% Multi layer network implementation
beta = 0.53;  % infectivity rate
gamma = 0.35; % rate from I to R
delta = 0.45; % rate from R to S 
k1 = 0.48;    % base opinion rate
k2 = 0.243;   % base opinion rate
k3 = k1; k5 = k1; 
k4 = k2; k6 = k2;
lambda1 = 0.143; % recovery rate from opinion
lambda2 = 0.143; % recovery rate from opinion
lambda3 = lambda1;  lambda5 = lambda1;
lambda4 = lambda2;  lambda6 = lambda2;
R_i = [beta/gamma, k1/lambda1, k2/lambda2] % the basica repdoction number of disease and each opinione SEPARATELY!
epsilon = 0.15; %il 15 percento della popolazione compliant infetta partecipa all'infezione 
omega = 0.005;

time = 550; 
[s_h, s_c, s_a, i_c, i_a, r_c,r_a,d, time_epi] = epi_behaviour(beta,gamma,delta,k1,k2,k3,k4,k5,k6,lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,epsilon,omega,time);

figure(1)
hold on
plot(time_epi,s_c, 'linewidth',1.3 )
plot(time_epi,s_h, 'linewidth',1.3 )
plot(time_epi,s_a, 'linewidth',1.3 )
legend('s_c', 's_{h}', 's_{a}')
hold off
%
figure(2)
hold on
plot(time_epi,i_c, 'linewidth',1.3 )
plot(time_epi,i_a, 'linewidth',1.3 )
legend('i_c','i_{a}')
%
figure(3)
hold on
plot(time_epi,r_c, 'linewidth',1.3 )
plot(time_epi,r_a, 'linewidth',1.3 )
legend('r_c', 'r_{a}')
hold off



%% Function section 
function [x_sh, x_sc, x_sa, y_c, y_a, z_rc,z_ra,z_d, time_epi] = epi_behaviour(beta,gamma,mu,k1,k2,k3,k4,k5,k6,lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,epsilon,omega,time)
    a1 = 1/2; a2=1/2; p1 =1; q11 =1; % Heuns method
    dt=0.01; %un centesimo di secondo per dt è ottimo con runge kutta
    
    s_c = 50/60e6;    % susceptible compliant
    s_a = 50/60e6;    % susceptible against
    i_c  = 10/60e6;    % infected compliant and careless
    i_a = 10/60e6;    % infected against
    r_c  = 0;          % recovered compliant and careless 
    r_a = 0;          % recovered against
    d = 0;             % deaths
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
    z_d = [];  z_d(1) = d;
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
            z_d = cat(2,z_d, d);
        end
       
        Aw = 0.5;  % Awareness / exposure / mean field forcing 
        rho = 0.7; % grado di protezione asociato all'essere compliant
        
        
        gamma_le = gamma;
        mu = 0.2; % mortality? !!!!
        % step 1
        im = epsilon*i_c+i_a;  %group that participate in the infection
        
        ksh1 = - beta*s_h*im - k1*Aw*s_h*(s_c+i_c+r_c) -k2*s_h*(s_a+i_a+r_a) +lambda1*s_c + lambda2*s_a + omega*r_c;
        ksc1 = k1*Aw*s_h*(s_c+i_c+r_c) -lambda1*s_c -beta*rho*s_c*im ;
        ksa1 = k2*s_h*(s_a+i_a+r_a) - lambda2*s_a - beta*s_a*im + omega*r_a;
        kic1 = beta*rho*s_c*im + beta*s_h*im - k3 *i_c*(s_a+i_a+r_a) + k4*Aw*i_a*(s_c+i_c+r_c) - lambda3*i_c- gamma*i_c;
        kia1 = +beta*s_a*im + k3*i_c*(s_a+i_a+r_a)+lambda3*i_c - k4*Aw*i_a*(s_c+i_c+r_c)-gamma_le*i_a;
        krc1 = gamma * i_c - k5*r_c*(s_a+i_a+r_a) + Aw*k6*r_a*(s_c+i_c+r_c) - mu*r_c -omega*r_c;
        kra1 = + gamma_le*i_a - Aw*k6*r_a*(s_c+i_c+r_c)+ k5*r_c*(s_a+i_a+r_a) -mu*r_a - omega*r_a;
        kd1 = mu*(r_c+r_a);

        % step 2
        t2 = t+p1*dt;
        s_h2 = s_h + q11*ksh1*dt;
        s_c2 = s_c + q11*ksc1*dt;
        s_a2 = s_a + q11*ksa1*dt;
        i_c2 = i_c   + q11*kic1*dt;
        i_a2 = i_a + q11*kia1*dt;
        r_c2 = r_c + q11*krc1*dt;
        r_a2 = r_a + q11*kra1*dt;
        d2= dt + q11*kd1*dt;
        im = (1-epsilon)*i_c2+i_a2;
        %
        ksh2 = - beta*s_h2*im - k1*Aw*s_h2*(s_c2+i_c2+r_c2) -k2*s_h2*(s_a2+i_a2+r_a2) +lambda1*s_c2 + lambda2*s_a2 + omega*r_c2;
        ksc2 = k1*Aw*s_h2*(s_c2+i_c2+r_c2) -lambda1*s_c2 -beta*rho*s_c2*im ;
        ksa2 = k2*s_h2*(s_a2+i_a2+r_a2) - lambda2*s_a2 - beta*s_a2*im + omega*r_a2;
        kic2 = beta*rho*s_c2*im + beta*s_h2*im - k3 *i_c2*(s_a2+i_a2+r_a2) + k4*Aw*i_a2*(s_c2+i_c2+r_c2) - lambda3*i_c2- gamma*i_c2;
        kia2 = +beta*s_a2*im + k3*i_c2*(s_a2+i_a2+r_a2)+lambda3*i_c2 - k4*Aw*i_a2*(s_c2+i_c2+r_c2)-gamma_le*i_a2;
        krc2 = gamma * i_c2 - k5*r_c2*(s_a2+i_a2+r_a2) + Aw*k6*r_a2*(s_c2+i_c2+r_c2) - mu*r_c2-omega*r_c2;
        kra2 = + gamma_le*i_a2 - Aw*k6*r_a2*(s_c2+i_c2+r_c2)+ k5*r_c2*(s_a2+i_a2+r_a2) -mu*r_a2 - omega*r_a2;
        kd2 = mu*(r_c2+r_a2);
        % update
        s_h = s_h + (a1*ksh1+a2*ksh2)*dt;
        s_c = s_c + (a1*ksc1+a2*ksc2)*dt;
        s_a = s_a + (a1*ksa1+a2*ksa2)*dt;
        i_c = i_c + (a1*kic1+a2*kic2)*dt;
        i_a = i_a + (a1*kia1+a2*kia2)*dt;
        r_c = r_c + (a1*krc1+a2*krc2)*dt;
        r_a = r_a + (a1*kra1+a2*kra2)*dt;
        d = d + (a1*kd1+a2*kd2)*dt;
        t = t + dt;
        cnt = cnt + 1;
    end    
end
