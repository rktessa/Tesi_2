clc;
clear all;
close all;
%% SIRS multiple simulations for study the system network
%% to Load
 load("data/multilpiSIRS.mat")
%% To save
 % filename = "multilpiSIRS.mat";
 % save(filename)

%% Readme!
% To run this code Optimization toolbox is necessary
% Solving a SIRS model with a Runge Kutta second order 

%%  SIRS simulation model 
% call the class functions container
obj = functionsContainer;
% all coefficients are modified on each simulation

% Simulation parameters
%1
gamma=1/18.95;
beta= gamma*2.38; %Così dovrei avere R0 = 2 %initial value
delta = 1/175; 
Ro = beta/gamma;
time =1500;
%2
gamma = 0.098;
beta = 0.65174;
delta = 1/68.67;
%3
beta = 0.431034;
gamma = 0.42069;
delta = 1/50;

%% For plots
i = 0;

%% Parameters I want to observe:
% 1- value at equilibrium of each compartment
% 2- Number of peaks in the infection
% 3- Time between each infection peak
% 4- Rise-time of infection 

%% 1
[taxis,xaxis,yaxis, zaxis] = obj.SIRS(beta,gamma,delta,time);
% Taking the last value of each compartments
final_S = xaxis(end); final_I = yaxis(end); final_R = zaxis(end);
tot = final_S+final_R+final_I; %check if sum equal to 1
%% 2 # of peaks
[pks,locs] = findpeaks(yaxis);
len_peak  = length(pks)
%% 3 inter-distances between peaks
dist_peaks = zeros(len_peak-1,1);
for l = 2:length(pks)
    dist_peaks(l-1) = locs(l)-locs(l-1);
end
mean_dist_peaks = mean(dist_peaks);
%% 4 rise-time of infection
[rCo,lt,ut,ll,ul] = risetime(yaxis,taxis,PercentReferenceLevels=[10 90]);
risetime2 = ut-lt %days
%% Plotting
i = i+1;
figure(i)
hold on
plot(taxis,xaxis, 'linewidth',1.5 )
plot(taxis,yaxis, 'linewidth',1.5 )
plot(taxis,zaxis, 'linewidth',1.5 )
plot(locs,pks,'*', 'Color','b')
title("SIRS model")
legend('S','I','R')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)],['delta: ' num2str(delta)]};
txt2 = ['final S: ' num2str(final_S)];
txt3 = ['final I: ' num2str(final_I)];
txt4 = ['final R: ' num2str(final_R)];
text(800,0.75, txt)
text((time-200),(final_S-0.025), txt2)
text((time-200),(final_I-0.025), txt3)
text((time-200),(final_R-0.025), txt4)


i = i+1;
figure(i)
hold on
plot(taxis,yaxis, 'linewidth',1.5 )
if  ~isempty(ut)
    plot(lt,ll, 'o', 'color', [0.4940 0.1840 0.5560])
    plot(ut,ul, 'o', 'color', [0.4940 0.1840 0.5560])
    yline([ll ul],'--',{'low level','upper level'})
end
plot(locs,pks,'*', 'Color','b')
title("SIRS model")
legend('Infected')
txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)],['delta: ' num2str(delta)]};
txt3 = ['final I: ' (final_I)];
text(800,final_I-final_I*0.5, txt)
text((time-200),(final_I-0.025), txt3)

%% Experiment with multiple values SIRS
 % [I_max, T_imax,R0, final_S, final_I, final_R, num_peak,mean_dist_peaks] = multipleSIRS();

%% Figure plot
beta_vec = linspace(0.1, 0.9,30);      d1 = length(beta_vec);
gamma_vec = linspace(1,1/15,30);       d2 = length(gamma_vec);
delta_vec = linspace(1/50, 1/360, 20); %d3 = length(delta_vec);
[xx, yy] = meshgrid(beta_vec,gamma_vec);

%% Time I max first peak
i =i+1; 
figure(i)
t = tiledlayout(1,1);
title(t,'Time where there is the first peak of infected')
for it = 1:1
    nexttile
    zz = reshape(T_imax(:,:,it),d1,d2,[]);
    zz = transpose(zz);
    for zi =1:length(zz)
        for zj = 1:length(zz)
            if zz(zi,zj) > 1499
                zz(zi,zj) = 0;
            end
        end
    end
    surf(xx,yy,zz)
    colorbar;
    colormap default;
    xlabel('\beta') 
    ylabel('\gamma') 
    txt = ['\delta = '  num2str(1/delta_vec(it))];
    title(txt)
end
%% R0
i =i+1; 
figure(i)
t = tiledlayout(1,1);
title(t,'R0')
for it = 1:1
    nexttile
    zz = reshape(R0(:,:,it),d1,d2,[]);
    zz = transpose(zz);
    surf(xx,yy,zz)
    colorbar;
    colormap default;
    xlabel('\beta') 
    ylabel('\gamma') 
    txt = ['\delta = '  num2str(1/delta_vec(it))];
    title(txt)
end
%% Max number of infected
i =i+1; 
figure(i)
t = tiledlayout(1,1);
title(t,'Max number of infected')
for it = 1:1
    nexttile
    zz = reshape(I_max(:,:,it),d1,d2,[]);
    zz = transpose(zz);
    surface(xx,yy,zz)
    colorbar;
    colormap default;
    xlabel('\beta') 
    ylabel('\gamma') 
    txt = ['\delta = '  num2str(1/delta_vec(it))];
    title(txt)
end
%% Numero di picchi nel periodo 
i =i+1; 
figure(i)
t = tiledlayout(4,5);
title(t,"Number of peaks in the time")
for it = 1:20
    nexttile
    zz = reshape(num_peak(:,:,it),d1,d2,[]);
    zz = transpose(zz);
    surface(xx,yy,zz)
    colorbar;
    colormap default;
    xlabel('\beta') 
    ylabel('\gamma') 
    txt = ['\delta = '  num2str(1/delta_vec(it))];
    title(txt)
end


%% Functions section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [taxis,xaxis,yaxis, zaxis] = SIRS(beta,gamma,delta,time)
            a1 = 1/2; a2=1/2; p1 =1; q11 =1; % Heuns method
            dt=0.01; %un centesimo di secondo per dt è ottimo con runge kutta
            x = 1-200/60e6; % susceptible
            y = 200/60e6;  % infected
            z = 0; % recovered
            t = 0;
            cnt=0;
            %Array creation and inititialization
            taxis=[]; taxis(1) =0; 
            xaxis=[]; xaxis(1) = x;
            yaxis=[]; yaxis(1) = y;
            zaxis=[]; zaxis(1) = z;
            while t < time
                if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
                    taxis = cat(2,taxis,t);
                    xaxis = cat(2,xaxis,x);
                    yaxis = cat(2,yaxis,y);
                    zaxis = cat(2,zaxis,z); 
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

        function [I_max, T_imax,R0, final_S, final_I, final_R, num_peak,mean_dist_peaks] = multipleSIRS()
            time =1500;
            % Initialize the vector of value span
            beta_vec = linspace(0.1, 0.9,30);     d1 = length(beta_vec);
            gamma_vec = linspace(1,1/15,30);       d2 = length(gamma_vec);
            delta_vec = linspace(1/50, 1/360, 20); d3 = length(delta_vec);
            % Initialize outputs
            final_S = zeros(d1,d2,d3); % all final value for the simulations
            final_I = zeros(d1,d2,d3);
            final_R = zeros(d1,d2,d3);
            num_peak = zeros(d1,d2,d3); % number of peaks in sirs
            mean_dist_peaks = zeros(d1,d2,d3);
            
            I_max = zeros(d1,d2,d3);
            T_imax = zeros(d1,d2,d3);
            R0 = zeros(d1,d2,d3);
            num = 0;
            for i = 1:length(beta_vec)
                beta = beta_vec(i);
                for j = 1:length(gamma_vec)
                    gamma = gamma_vec(j);
                    for k = 1:length(delta_vec)
                        num = num+1
                        delta = delta_vec(k);
                        %Esecuzione della trasformazione, variando beta e gamma
                        [taxis,xaxis,yaxis, zaxis] = SIRS(beta,gamma,delta,time);
                        % Taking the last value of each compartments
                        final_S(i,j,k) = xaxis(end); 
                        final_I(i,j,k) = yaxis(end); 
                        final_R(i,j,k) = zaxis(end);   
                        %% 2 # of peaks
                        [pks,locs] = findpeaks(yaxis);
                        %len = length(pks)
                        num_peak(i,j,k)  = length(pks); 
                        %% 3 inter-distances between peaks
                        if num_peak(i,j,k)>1
                            dist_peaks = zeros(num_peak(i,j,k)-1,1);
                            
                            for l = 2:length(pks)
                                dist_peaks(l-1) = locs(l)-locs(l-1);
                            end
                            mean_dist_peaks(i,j,k) = mean(dist_peaks);
                        end
                        %% 4 rise-time of infection
                        
                        %Compute max value of infected
                        [max_i, index] = max(yaxis);
                        time_i_max = taxis(index);
                        I_max(i,j,k) = max_i;
                        T_imax(i,j,k) = time_i_max;
                        R0(i,j,k) = (beta/gamma);
                    end
                end
            end
        end





