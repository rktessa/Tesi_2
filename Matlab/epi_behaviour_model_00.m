clc;
clear all;
%% Epidemiological Behavioural model
% Multi layer network implementation







%% Function section 
function [s_co, s_ca, s_ag, i_c, i_ag, r_c,r_ag, time_epi] = epi_behaviour(beta,gamma,delta,k1,k2,lambda1,lambda2,time)
    a1 = 1/2; a2=1/2; p1 =1; q11 =1; % Heuns method
    dt=0.01; %un centesimo di secondo per dt è ottimo con runge kutta
    s_ca = 1-200/60e6; % susceptible careless
    s_co = 50/60e6;    % susceptible compliant
    s_ag = 50/60e6;    % susceptible against
    i_c  = 10/60e6;    % infected compliant and careless
    i_ag = 10/60e6;    % infected against
    r_c  = 0;          % recovered compliant and careless 
    r_ag = 0;          % recovered against
    
    t = 0; % initialize the time counter
    cnt=0;
    %Array creation and inititialization
    time_epi=[]; time_epi(1) =0; 
    x_sca   =[]; x_sca(1)    = s_ca;
    x_sco   =[]; x_sco(1)    = s_co;
    x_sag   =[]; x_sag_r(1)  = s_ag;
    y_c=[];   y_c(1) = i_c;
    y_ag=[];  y_ag(1)= i_ag;
    z_rc = []; z_rc(1) = r_c;
    z_rag=[];  z_rag(1)= r_ag; 
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
        end
        % step 1
        kx1 = - beta*x*y + delta*z;
        ky1 = beta*x*y - gamma*y;
        kz1 = + gamma*y - delta*z;
        % step 2
        t2 = t+p1*dt;
        x2 = x + q11*kx1*dt;
        y2 = y + q11*ky1*dt;
        z2 = z + q11*kz1*dt;
        %
        kx2 = - beta*x2*y2 + delta*z2;
        ky2 = beta*x2*y2 - gamma*y2;
        kz2 = gamma*y2 - delta*z2;
        % update
        x = x + (a1*kx1+a2*kx2)*dt;
        y = y + (a1*ky1+a2*ky2)*dt;
        z = z + (a1*kz1+a2*kz2)*dt;
       
        t = t + dt;
        cnt = cnt + 1;
    end    
end
