close all;
% Uncomment and load only the first time
% load("multipli_dati_behavior.mat")
 %load("multilpiAgainst.mat")

fig = 0;
lam1_vec = linspace(1/2, 1/40, 20); d3 = length(lam1_vec);
lam2_vec = linspace(1/2, 1/40, 20); d4 = length(lam2_vec);
B1_vec = linspace(0.1, 10, 20);    
B2_vec = linspace(0.1, 10, 20);   
k1_vec = B1_vec.*lam1_vec;    d1 = length(k1_vec);
k2_vec = B2_vec.*lam2_vec;   d2 = length(k2_vec);
vec_R2 = [2,5,8,12,16,19]; % six values to calculated R2

%% Equilibrium plot of R1 and R2 first case
d1 = 11; d2 = 13;
fig = fig+1;
figure(fig)
titles = "max_againste_heatmap_1.pdf";
[xx, yy] = meshgrid(lam1_vec,lam2_vec);
aga = reshape(Ag_max(d1,d2,:,:),20,20,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
xlabel('\lambda_1') 
ylabel('\lambda_2')
txt3 ="k_1= "+ num2str(round(k1_vec(d1),2))+ ", k_2= "+ num2str(round(k2_vec(d2),2));
title(txt3)
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 24]);
set(gcf, 'PaperSize', [24 24]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')
%% second case
d1 = 11; d2 = 18;
fig = fig+1;
figure(fig)
titles = "max_againste_heatmap_2.pdf";
[xx, yy] = meshgrid(lam1_vec,lam2_vec);
aga = reshape(Ag_max(d1,d2,:,:),20,20,[]);
 zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 ="k_1= "+ num2str(round(k1_vec(d1),2))+ ", k_2= "+ num2str(round(k2_vec(d2),2));
title(txt3)
xlabel('\lambda_1') 
ylabel('\lambda_2')
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 24]);
set(gcf, 'PaperSize', [24 24]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')

%% Against dynamic first example
fig = fig+1;
figure(fig)
titles = "against_evolution_vari_1.pdf";
%1,1
d1 = 11; d2 = 13;d3 =10;
hold on
for i = 3:17
    d4 = i;
    aga = reshape(Against(d1,d2,d3,d4,:),1,2001,[]);
    zz = aga; %zz = transpose(zz);
    txt = "\lambda_2: "+ num2str(round(lam1_vec(d4),2))+" B_2: "+num2str(round(k2_vec(d2)/lam2_vec(d4),2));
    plot(zz(1:400),'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
legend(Orientation='vertical', Location='eastoutside')
txt3 = "Against evolution "+ "B_1: "+num2str(round(k1_vec(d1)/lam1_vec(d3),2))+ " \lambda_1: "+ num2str(round(lam1_vec(d3),2)) +" k_1: "+ num2str(round(k1_vec(d1),2))+" k_2: "+ num2str(round(k2_vec(d2),2));
title(txt3)
xlabel('time') 
ylabel("Against ") 
hold off
fontsize(13,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 23 16]);
set(gcf, 'PaperSize', [24 16]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')
%% Against dynamic second example
fig = fig+1;
figure(fig)
titles = "against_evolution_vari_2.pdf";
%1,1
d1 = 11; d2 = 18;d3 =10;
hold on
for i = 3:17
    d4 = i;
    aga = reshape(Against(d1,d2,d3,d4,:),1,2001,[]);
    zz = aga; %zz = transpose(zz);
  txt = "\lambda_2: "+ num2str(round(lam1_vec(d4),2))+" B_2: "+num2str(round(k2_vec(d2)/lam2_vec(d4),2));
    plot(zz(1:400),'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
legend(Orientation='vertical', Location='eastoutside')
txt3 = "Against evolution "+ "B_1: "+num2str(round(k1_vec(d1)/lam1_vec(d3),2))+ " \lambda_1: "+ num2str(round(lam1_vec(d3),2)) +" k_1: "+ num2str(round(k1_vec(d1),2))+" k_2: "+ num2str(round(k2_vec(d2),2));
title(txt3)
xlabel('time') 
ylabel("Against ") 
hold off
fontsize(13,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 23 16]);
set(gcf, 'PaperSize', [24 16]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')