clc;
clear all;
%% Epidemiological Behavioural model
% Multi layer network implementation







%% Function section 
function [s_co, s_ca, s_ag, i_c, i_ag, r_c,r_ag, time_epi] = epi_behaviour(beta,gamma,delta,time)
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
    x_s =[]; x_s(1) = s;
    y_i =[]; y_i(1) = i;
    z_r =[]; z_r(1) = r;
    x_ca=[]; x_ca(1)= ca;
    y_co=[]; y_co(1)= co;
    z_ag=[]; z_ag(1)= ag; 
    while t < time
        %every 100 millisecond I save the result
        if mod(cnt,100) == 0 && cnt ~=0 
            time_epi = cat(2,time_epi,t);
            x_s = cat(2, x_s, s);
            y_i = cat(2, y_i, i);
            z_r = cat(2, z_r, r); 
            x_ca= cat(2, x_ca, ca);
            y_co= cat(2, y_co, co);
            z_ag= cat(2, z_ag, ag);
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
