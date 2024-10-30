close all;
% Uncomment and load only the first time
% load("multipli_dati_behavior.mat")
 %load("multilpiAgainst.mat")
addpath(genpath('..\'))
% load("..\data/multilpiAgainst.mat")
% load("..\data/multilpiBehaviour.mat")
fig = 0;
k1_vec = linspace(0.1, 0.99,30);  d1 = length(k1_vec);
k2_vec = linspace(0.1, 0.99,30);  d2 = length(k2_vec);
lam1_vec = linspace(1/5, 1/30, 20);  d3 = length(lam1_vec);
lam2_vec = linspace(1/5, 1/30, 20);  d4 = length(lam2_vec);
%% Equilibrium plot of R1 and R2 first case
d1 = 11; d2 = 13;
fig = fig+1;
figure(fig)
titles = "max_againste_heatmap_1.pdf";
[xx, yy] = meshgrid(lam1_vec,lam2_vec);
aga = reshape(Ag_max(d1,d2,:,:),20,20,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
xline(lam1_vec(10), 'r', 'LineWidth',2)
colorbar;
colormap default;
xlabel('\lambda_1') 
ylabel('\lambda_2')
txt3 ="k_1= "+ num2str(round(k1_vec(d1),2))+ ", k_2= "+ num2str(round(k2_vec(d2),2));
title(txt3)
ylim([lam1_vec(20) lam1_vec(1)])
fontsize(20,"points")
legend("","section line", Orientation="horizontal",Location="southoutside")
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
xline(lam1_vec(10), 'r', 'LineWidth',2)
colorbar;
colormap default;
txt3 ="k_1= "+ num2str(round(k1_vec(d1),2))+ ", k_2= "+ num2str(round(k2_vec(d2),2));
title(txt3)
xlabel('\lambda_1') 
ylabel('\lambda_2')
ylim([lam1_vec(20) lam1_vec(1)])
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