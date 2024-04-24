clc;
 clear all;
close all;
%% Bifurcation plot of behavioural model

% Load the data
 load("data/multilpiBehaviour.mat")
% What it is load:
% Co_max,Ag_max, T_comax, T_agmax, final_Ca, final_Co, final_Ag
% The parameter range of variation is:
% k1_vec   = k2_vec  = linspace(0.1, 0.99,30);    
% lam1_vec = lam2_vec= linspace(1/5, 1/30, 20);         
%%

k1_vec = linspace(0.1, 0.99,30);
lam1_vec = linspace(1/5, 1/30, 20);
vec_R2 = [2,5,8,12,16,19]; % six values to calculated R2

%% Plot Section
fig = 0; %initialized for figures

%% Max value function of R1 and R2
fig = fig+1;
figure(fig)
t = tiledlayout(1,2);
title(t, "Sensitivity plot")
%1,1
d1 = 10; d2 = 10;
nexttile
r1 = k1_vec./lam1_vec(d1);
r2 = k1_vec./lam1_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(Co_max(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Max Compliant "+ '\lambda_1: '+ num2str(lam1_vec(d1)) +'\lambda_2: '+ num2str(lam1_vec(d2));
title(txt3)
xlabel('R1') 
ylabel('R2') 
% 1,2
nexttile
r1 = k1_vec./lam1_vec(d1);
r2 = k1_vec./lam1_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(Ag_max(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Max Against "+ '\lambda_1: '+ num2str(lam1_vec(d1)) +'\lambda_2: '+ num2str(lam1_vec(d2));
title(txt3)
xlabel('R1') 
ylabel('R2') 

%% Equilibrium plot of R1 and R2
fig = fig+1;
figure(fig)
t = tiledlayout(1,2);
title(t, "Sensitivity plot at equilibrium")
%1,1
d1 = 10; d2 = 10;
nexttile
r1 = k1_vec./lam1_vec(d1);
r2 = k1_vec./lam1_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(final_Co(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Final Compliant "+ '\lambda_1: '+ num2str(lam1_vec(d1)) +'\lambda_2: '+ num2str(lam1_vec(d2));
title(txt3)
xlabel('R1') 
ylabel('R2') 
% 1,2
nexttile
r1 = k1_vec./lam1_vec(d1);
r2 = k1_vec./lam1_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(final_Ag(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Final against "+ '\lambda_1: '+ num2str(lam1_vec(d1)) +'\lambda_2: '+ num2str(lam1_vec(d2));
title(txt3)
xlabel('R1') 
ylabel('R2') 
% 1,3


fig = fig+1;
figure(fig)
r1 = k1_vec./lam1_vec(d1);
r2 = k1_vec./lam1_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(final_Ca(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Final Careless "+ '\lambda_1: '+ num2str(lam1_vec(d1)) +'\lambda_2: '+ num2str(lam1_vec(d2));
title(txt3)
xlabel('R1') 
ylabel('R2') 


%% Max value function of R1
fig = fig+1;
figure(fig)
t = tiledlayout(1,2);
title(t, "Bifurcation plot for R1")
%1,1
d1 = 1; d2 = 19; d3 = 15;
nexttile
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(Co_max(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "R2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
txt2 ="Max Compliant "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt2);
xlabel('R1') 
ylabel("Co ")
hold off
% 1,2
nexttile
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(Ag_max(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "R2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
hold off
txt2 ="Max Against "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt2);
xlabel('R1') 
ylabel('Ag') 


%% Equilibrium value function of R1
fig = fig+1;
figure(fig)
t = tiledlayout(1,2);
title(t, "Sensitivity plot at equilibrium for R1")
%1,1
d1 = 10; d3 = 15;
nexttile
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(final_Co(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "R2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
txt3 = "Final Compliant "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt3)
xlabel('R1') 
ylabel("Co ") 
hold off
% 1,2
nexttile
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(final_Ag(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "R2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
txt3 = "Final Against "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt3)
xlabel('R1') 
ylabel('Ag') 
hold off
% 1,3
fig = fig+1;
figure(fig)
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(final_Ca(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "R2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
txt3 = "Final Careless "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt3)
xlabel('R1') 
ylabel('Ca') 
hold off
