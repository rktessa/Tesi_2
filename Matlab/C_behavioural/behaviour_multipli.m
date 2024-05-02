close all;
clear all;
clc;
addpath('..\')
% addpath('..\Matlab')
%%
% filename = "multilpiCareless.mat";
% save(filename, "Careless", '-v7.3');
% %%
% filename = "multilpiCompliant.mat";
% save(filename, "Compliant", '-v7.3');
% 
% filename = "multilpiAgainst.mat";
% save(filename, "Against", '-v7.3');
%%  Behaviour simulation model
% call the class functions container
% obj = fopen('Matlab\functionsContainer.m');
 obj = functionsContainer
% Simulation parameters

k1 =  1.1/4.03 %from Ca to Co
k2 = 1/3.97 %from Ca to Ag

lambda_1 = 1/10; %fatigue to mantain Co behaviour
lambda_2 = 1/10; %fatigue to mantain Ag behaviour
time = 2000;
% For plots
i = 0;

%% Parameters I want to observe:
% 1- value at equilibrium of each compartment
% 2- Number of peaks in the infection
% 3- Time between each infection peak
% 4- Rise-time of infection 

%% 1
[taxisRK,xaxisRK,yaxisRK,zaxisRK,awareness] = obj.Behaviour_RK(k1,k2,lambda_1,lambda_2,time);
% Taking the last value of each compartments
final_Ca = xaxisRK(end); final_Co = yaxisRK(end); final_Ag = zaxisRK(end);
tot = final_Ca + final_Co + final_Ag; %check if sum equal to 1
% 
%% 2 # of peaks
[pks,locs] = findpeaks(yaxisRK);
len_peak  = length(pks)
%% 3 inter-distances between peaks
dist_peaks = zeros(len_peak-1,1);
for l = 2:length(pks)
    dist_peaks(l-1) = locs(l)-locs(l-1);
end
mean_dist_peaks = mean(dist_peaks)
%% 4 rise-time of infection
[rCo,lt,ut,ll,ul] = risetime(xaxisRK,taxisRK,PercentReferenceLevels=[10 90]);
risetime2 = ut-lt %days

%% Plotting
    i = i+1;
    figure(i)
    hold on
    plot(taxisRK,xaxisRK, 'r', 'linewidth',1.0 )
    plot(taxisRK,yaxisRK, 'b', 'linewidth',1.0 )
    plot(taxisRK,zaxisRK, 'g', 'linewidth',1.0 )
    % plot(lt,ll, 'o', 'color', [0.4940 0.1840 0.5560])
    % plot(ut,ul, 'o', 'color', [0.4940 0.1840 0.5560])
    % yline([ll ul],'--',{'low level','upper level'})
    title("Behaviour model Runge Kutta")
    legend('Careless','Compliant','Against')
    txt = {['k1: ' num2str(k1)],['k2: ' num2str(k2)],['lambda_1: ' num2str(lambda_1)],['lambda_2: ' num2str(lambda_2)]};
    txt2 = ['final Ca: ' num2str(final_Ca)];
    txt3 = ['final Co: ' num2str(final_Co)];
    txt4 = ['final Ag: ' num2str(final_Ag)];
    text(time-200,0.75,txt)
    text((time-200),(final_Ca-0.025), txt2)
    text((time-200),(final_Co+0.025), txt3)
    text((time-200),(final_Ag+0.025), txt4)
    % Compliant graph
    i = i+1;
    figure(i)
    hold on
    plot(taxisRK,yaxisRK, 'linewidth',1.5 )
    if  ~isempty(ut)
        plot(lt,ll, 'o', 'color', [0.4940 0.1840 0.5560])
        plot(ut,ul, 'o', 'color', [0.4940 0.1840 0.5560])
        yline([ll ul],'--',{'low level','upper level'})
    end
    plot(locs,pks,'*', 'Color','b')
    title("Behavioural model")
    legend('Compliant')
    txt = {['k1: ' num2str(k1)],['k2: ' num2str(k2)],['lambda_1: ' num2str(lambda_1)],['lambda_2: ' num2str(lambda_2)]};
    text(800,final_Co-final_Co*0.5, txt)
    text((time-200),(final_Co-0.025), txt3)
 
%% Varying alpha and alpha_prime parameters
% The Initial condition are the same in all the case, but alpha and
% alpha_prime coefficients vary in an interval [0.01,1] 


% [Co_max,Ag_max, T_comax, T_agmax, final_Ca, final_Co, final_Ag, Careless,Compliant, Against] = multipleBehaviour();

%% Plot the results

% alpha_pl = (0.005:step:(step*d))/N;
% alpha_prime_pl = (0.005:step:(step*d))/N;
% [xx, yy] = meshgrid(alpha_pl,alpha_prime_pl);
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,Co_max)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Max value of Compliants for alpha and alpha_{prime}')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,A_max)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Max value of Anti for alpha and alpha_{prime}')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,rCo_m)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Rise time of Compliant for alpha and alpha_{prime}')
% 
% i =i+1; 
% figure(i)
% surf(xx,yy,rA_m)
% colorbar;
% colormap default;
% xlabel('alpha_{prime}') 
% ylabel('alpha') 
% title('Rise time of Anti for alpha and alpha_{prime}')

%% Multiple simulazioni

function [Co_max,Ag_max, T_comax, T_agmax, final_Ca, final_Co, final_Ag, Careless,Compliant, Against] = multipleBehaviour()
            
            time =2000;
            % Initialize the vector of value span
            k1_vec = linspace(0.1, 0.99,30);    d1 = length(k1_vec);
            k2_vec = linspace(0.1, 0.99, 30);   d2 = length(k2_vec);
            lam1_vec = linspace(1/5, 1/30, 20); d3 = length(lam1_vec);
            lam2_vec = linspace(1/5, 1/30, 20); d4 = length(lam2_vec);
            % Initialize outputs
            final_Ca = zeros(d1,d2,d3,d4); % all final value for the simulations
            final_Co = zeros(d1,d2,d3,d4);
            final_Ag = zeros(d1,d2,d3,d4);
           
            Ag_max = zeros(d1,d2,d3, d4);
            T_agmax = zeros(d1,d2,d3, d4);
            Co_max = zeros(d1,d2,d3, d4);
            T_comax = zeros(d1,d2,d3, d4);

            % Salvo le simulazioni
            Careless = zeros(d1,d2,d3, d4, time+1);
            Compliant = zeros(d1,d2,d3, d4, time+1);
            Against = zeros(d1,d2,d3, d4, time+1);
            
            num = 0;
            for i = 1:d1
                k1 = k1_vec(i);
                for j = 1:d2
                    k2 = k2_vec(j);
                    for k = 1:d3
                        lam1 = lam1_vec(k);
                        for u = 1:d4
                            lam2 = lam2_vec(u);
                            num = num+1
                            
                            %Esecuzione della trasformazione, variando beta e gamma
                            [taxis,xaxis,yaxis,zaxis,awareness] = Behaviour_RK(k1,k2,lam1,lam2,time);
                            % Taking the last value of each compartments
                            final_Ca(i,j,k,u) = xaxis(end); 
                            final_Co(i,j,k,u) = yaxis(end); 
                            final_Ag(i,j,k,u) = zaxis(end);   
                            % Save the simualtion
                            
                            Careless(i,j,k,u,:) = xaxis;
                            Compliant(i,j,k,u,:) = yaxis;
                            Against(i,j,k,u,:) = zaxis;
                            
                            %Compute max value of Co and Ag
                            [max_i, index] = max(yaxis);
                            time_i_max = taxis(index);
                            Co_max(i,j,k,u) = max_i;
                            T_comax(i,j,k,u) = time_i_max;

                            [max_i2, index2] = max(zaxis);
                            time_i_max2 = taxis(index);
                            Ag_max(i,j,k,u) = max_i2;
                            T_agmax(i,j,k,u) = time_i_max2;
                            
                        end
                    end
                end
            end
        end


function [taxis,xaxis,yaxis,zaxis,awareness] = Behaviour_RK(k1,k2,lambda_1,lambda_2,time)
    
            a1 = 1/2; a2=1/2; p1 =1; q11 =1;
            x = 1-400/60e6; % Careless
            y = 200/60e6;  % Compliant
            z = 200/60e6; % Against
            dt = 0.01;
            t = 0;
            cnt=0;
            awareness = 0;
            %Array creation and inititialization
            taxis=[]; taxis(1) = 0; 
            xaxis=[]; xaxis(1) = x;
            yaxis=[]; yaxis(1) = y;
            zaxis=[]; zaxis(1) = z;
            while t < time
                if mod(cnt,100) == 0 && cnt ~=0 %every 100 iterations I save the result
                    taxis = cat(2,taxis,t);
                    xaxis = cat(2,xaxis,x);
                    yaxis = cat(2,yaxis,y);
                    zaxis = cat(2,zaxis,z);
                end
                % step 1
                kx1 =  -k1*x*y - k2*x*z + lambda_1*y + lambda_2*z;
                ky1 =   k1*x*y - lambda_1*y;
                kz1 =   k2*x*z - lambda_2*z;
                % step 2
                t2 = t+p1*dt;
                x2 = x + q11*kx1*dt;
                y2 = y + q11*ky1*dt;
                z2 = z + q11*kz1*dt;
                kx2 =  -k1*x2*y2 - k2*x2*z2 + lambda_1*y2 + lambda_2*z2;
                ky2 =   k1*x2*y2 - lambda_1*y2;
                kz2 =   k2*x2*z2 - lambda_2*z2;
                % update
                x = x + (a1*kx1+a2*kx2)*dt;
                y = y + (a1*ky1+a2*ky2)*dt;
                z = z + (a1*kz1+a2*kz2)*dt;
                t = t + dt;
                cnt = cnt + 1;
            end
        end




