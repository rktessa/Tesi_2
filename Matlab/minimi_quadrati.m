close all;
clear all;
clc;
%% Minimi quadrati
% Primo esempio  
dist = [50, 100, 200, 300, 500];
costo = [4300, 8200, 16000, 23500, 42000];

% Modello lineare
[a_treni, b_treni] = LS_linear(dist, costo);

costo_fit = a_treni*dist + b_treni;

i = 0;

i = i+1;
figure(i)
hold on
plot(dist,costo,'b*')
plot(dist,costo_fit,'g')
title("Tariffe ferroviare")
xlabel('Distanza [km]')
ylabel('Costo [L.]')
legend('dati', 'modello lineare')


% Secondo esempio 
anno = [0, 1, 2, 3,4];
pop = [7.24, 23.2, 62.9, 122.8, 203.2];
log_pop = log(pop);
% Modello lineare
[a_pop, b_pop] = LS_linear(anno, log_pop);

pop_fit = a_pop*anno + b_pop;

i = i+1;
figure(i)
hold on
plot(anno,log_pop,'b*')
plot(anno,pop_fit,'r')
title("Popolazione U.S.")
xlabel('Anno')
ylabel('Log Popolazione in milioni')
legend('dati', 'modello lineare')

%% Modelli lineari

%% Modelli quadratici


%% Funzioni
function [a, b] = LS_linear(x, y)
    %with x and y data points 
    m_x = mean(x); % mean
    m_y = mean(y);
    v_x = var(x);  % variance
    v_y = var(y);
    c_xy = cov(x,y); % covariance
    rms_x_2 = rms(x)^2; % rms square 
    rms_y_2 = rms(y)^2;
    m_xy = mean(x.*y);

    M = [ 1, m_x;
          m_x, rms_x_2]
    v = [m_y; m_xy];

    sol = inv(M)*v;

    b = sol(1);
    a = sol(2);
end