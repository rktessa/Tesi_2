clc;
clear;
close;
fig = 0;
%% CODE TO CALCULATE THE E_0 EXPLICITLY  14/11/2024
caso =5;
[lambda_1,lambda_2,k_1,k_2,time,title_fig, beta,gamma,delta,rho,epsilon,psi, SH0,SC0,SA0,IC0,IA0]=scenario(caso);
[E_0, R_0, B_1, B_2] = calcolo_E_0(beta, gamma,psi, rho, epsilon, k_1, k_2,lambda_1, lambda_2, SH0, SC0,SA0, IC0,IA0 )
 fig = Heat_map2(fig, beta, gamma,psi, rho, epsilon, k_1, k_2, SH0, SC0,SA0, IC0,IA0)



%% Function section
% Function for E_0
function [E_0, R_0, B_1, B_2] = calcolo_E_0(beta, gamma,psi, rho, epsilon, k_1, k_2,lam_1, lam_2, SH0, SC0,SA0, IC0,IA0 )
C = SC0 + IC0;
A = SA0 + IA0;
x1 = lam_1 + A*k_2;
x2 = lam_2 + psi*C*k_1;
sigma = x1+x2+gamma;
% Basic reproduction number Epidemic
R_0 = beta/gamma;
B_1 = k_1/lam_1;
B_2 = k_2 / lam_2;
% E_0 = R_0*(SA0*((gamma+x1+epsilon*x2)/sigma) + (SH0+ rho*SC0)*((x1+epsilon*(x2+gamma))/sigma));

E_0= (beta*gamma)/(gamma^2 - lam_1*lam_2) + (beta*lam_1)/(gamma^2-lam_1*lam_2);



end

% Function "Scenari"
function [lambda_1,lambda_2,k1,k2,time,title_fig, beta,gamma,delta,rho,epsilon,psi, SH0,SC0,SA0,IC0,IA0] = scenario(caso)
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
    SC0 = 10e6/60e6;
    SA0 = 10e6/60e6;
    IC0 = 10/60e6;
    IA0 = 10/60e6;
    
    SH0 = 1- SC0-SA0-IC0-IA0;
    time = 3000;
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
        SC0 = 50/60e6;
        SA0 = 50/60e6;
        IC0 = 10/60e6;
        IA0 = 10/60e6;
        % IC0 = 0;
        % IA0 = 0;
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
        time = 500;
        C_zero = 100;
        A_zero = 100;
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
        time = 500;
        SC0 = 50/60e6;
        SA0 = 50/60e6;
        IC0 = 10/60e6;
        IA0 = 10/60e6;
        % IC0 = 0;
        % IA0 = 0;
        SH0 = 1- SC0-SA0-IC0-IA0;
        title_fig = 'epi_behav_sim_B2_mag_B1.pdf';
    end
end

function fig = Heat_map(fig, beta, gamma,psi, rho, epsilon, k_1, k_2, SH0, SC0,SA0, IC0,IA0)
    lam1_vec = linspace(1/2, 1/50, 100);  d3 = length(lam1_vec);
    lam2_vec = linspace(1/2, 1/50, 100);  d4 = length(lam2_vec);
    epsilon_vec = linspace(0,1,100); d5 = length(epsilon_vec);
    E0_matrix = zeros(d3,d4,d5); 
    for i = 1:d3
        for j = 1:d4
            for k = 1:d5
                [E_0, R_0, B_1, B_2] = calcolo_E_0(beta, gamma,psi, rho, epsilon_vec(k), k_1, k_2,lam1_vec(i), lam2_vec(j), SH0, SC0,SA0, IC0,IA0 );
                E0_matrix(i,j,k) = E_0;
            end
        end
    end

    fig = fig+1;
    figure(fig)
    titles = "E0_heatmap_1.pdf";
    [xx, yy,zz] = ndgrid(lam1_vec,lam2_vec,epsilon_vec);

    aga = E0_matrix;
    FF = aga; %FF = transpose(FF);
    
    hold on
    for i = 1:100
        surf(xx(:,:,i),yy(:,:,i),zz(:,:,i),FF(:,:,i),'EdgeColor','none') ;
    end
    view(235,10)
    
    colorbar;
    colormap default;
    xlabel('\lambda_1') 
    ylabel('\lambda_2')
    zlabel('\epsilon ')
    fontsize(20,"points")
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 25 16]);
    set(gcf, 'PaperSize', [26 16]); % dimension on x axis and y axis resp.
     print(gcf,titles,'-dpdf')

end


function fig = Heat_map2(fig, beta, gamma,psi, rho, epsilon, k_1, k_2, SH0, SC0,SA0, IC0,IA0)
    lam1_vec = linspace(1/2, 1/50, 100);  d3 = length(lam1_vec);
    lam2_vec = linspace(1/2, 1/50, 100);  d4 = length(lam2_vec);
    
    E0_matrix = zeros(d3,d4); 
    for i = 1:d3
        for j = 1:d4
           
                [E_0, R_0, B_1, B_2] = calcolo_E_0(beta, gamma,psi, rho, epsilon, k_1, k_2,lam1_vec(i), lam2_vec(j), SH0, SC0,SA0, IC0,IA0 );
                E0_matrix(i,j) = E_0;
            end
        end
    

    fig = fig+1;
    figure(fig)
    titles = "E0_heatmap_1.pdf";
    [xx, yy] = ndgrid(lam1_vec,lam2_vec);

    aga = E0_matrix;
    FF = aga; FF = transpose(FF);
    
    hold on
    surface(xx,yy,FF, 'edgecolor','none')  
    colorbar;
    colormap default;
    xlabel('\lambda_1') 
    ylabel('\lambda_2')
   
    fontsize(20,"points")
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 25 16]);
    set(gcf, 'PaperSize', [26 16]); % dimension on x axis and y axis resp.
     print(gcf,titles,'-dpdf')

     end
     