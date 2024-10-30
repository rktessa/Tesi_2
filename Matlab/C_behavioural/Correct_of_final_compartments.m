close all;
% Uncomment and load only the first time
% load("multipli_dati_behavior.mat")
addpath(genpath('..\'))
%load("..\data/multilpiCompliant.mat")
%load("..\data/multilpiBehaviour.mat")
fig = 0;
% lam1_vec = linspace(1/2, 1/40, 30); d3 = length(lam1_vec);
% lam2_vec = linspace(1/2, 1/40, 30); d4 = length(lam2_vec);
% 
% B1_vec = linspace(0.1, 10, 20);    
% B2_vec = linspace(0.1, 10, 20);   
% k1_vec = B1_vec.*lam1_vec;    d1 = length(k1_vec);
% k2_vec = B2_vec.*lam2_vec;   d2 = length(k2_vec);

k1_vec = linspace(0.1, 0.99,30);  d1 = length(k1_vec);
k2_vec = linspace(0.1, 0.99,30);  d2 = length(k2_vec);
lam1_vec = linspace(1/5, 1/30, 20);  d3 = length(lam1_vec);
lam2_vec = linspace(1/5, 1/30, 20);  d4 = length(lam2_vec);
vec_R2 = [2,5,8,12,16,19]; % six values to calculated R2
%% Heat map final values first
% Compliant
titles= "Final_Compliant_heat_map.pdf";
%1,1
d1 = 5; d2 = 10;
fig = fig+1;
figure(fig)
r1 = k1_vec./lam1_vec(d1);
r2 = k2_vec./lam2_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(final_Co(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Final Compliant "+ '\lambda_1: '+ num2str(round(lam1_vec(d1),2)) +' \lambda_2: '+ num2str(round(lam1_vec(d2),2));
title(txt3)
xlabel('B_1') 
ylabel('B_2') 
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 24]);
set(gcf, 'PaperSize', [24 24]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')

%% Against
titles= "Final_Against_heat_map.pdf";
%1,1
d1 = 5; d2 = 10;
fig = fig+1;
figure(fig)
r1 = k1_vec./lam1_vec(d1);
r2 = k2_vec./lam2_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(final_Ag(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Final Against "+ '\lambda_1: '+ num2str(round(lam1_vec(d1),2)) +' \lambda_2: '+ num2str(round(lam1_vec(d2),2));
title(txt3)
xlabel('B_1') 
ylabel('B_2') 
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 24]);
set(gcf, 'PaperSize', [24 24]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')

 %% Heedless
titles= "Final_Heedless_heat_map.pdf";
%1,1
d1 = 5; d2 = 10;
fig = fig+1;
figure(fig)
r1 = k1_vec./lam1_vec(d1);
r2 = k2_vec./lam2_vec(d2);
[xx, yy] = meshgrid(r1,r2);
aga = reshape(final_Ca(:,:,d1,d2),30,30,[]);
zz = aga; zz = transpose(zz);
surface(xx,yy,zz, 'edgecolor','none')
colorbar;
colormap default;
txt3 = "Final Heedless "+ '\lambda_1: '+ num2str(round(lam1_vec(d1),2)) +' \lambda_2: '+ num2str(round(lam1_vec(d2),2));
title(txt3)
xlabel('B_1') 
ylabel('B_2') 
fontsize(20,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 24 24]);
set(gcf, 'PaperSize', [24 24]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')

 %% Final compliant as B1 function
fig = fig+1;
figure(fig)
titles = "final_compliant_B1.pdf";
%1,1
d1 = 10; d3 = 15;
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(final_Co(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "B_2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
legend(Orientation="vertical", Location="bestoutside")
txt3 = "Final Compliant "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt3)
xlabel('B_1') 
ylabel("Compliant") 
hold off
fontsize(13,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 16]);
set(gcf, 'PaperSize', [21 16]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')


%% Final Against as B1 function
fig = fig+1;
figure(fig)
titles = "final_against_B1.pdf";
%1,1
d1 = 10; d3 = 15;
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(final_Ag(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "B_2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
legend(Orientation="vertical", Location="bestoutside")
txt3 = "Final Against "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt3)
xlabel('B_1') 
ylabel("Against") 
hold off
fontsize(13,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 16]);
set(gcf, 'PaperSize', [21 16]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')

%% Final heedless as B1 function
fig = fig+1;
figure(fig)
titles = "final_heedless_B1.pdf";
%1,1
d1 = 10; d3 = 15;
hold on
r1 = k1_vec./lam1_vec(d1);
for i = 1:6
    d2 = vec_R2(i);
    aga = reshape(final_Ca(:,d3,d1,d2),1,30,[]);
    zz = aga; %zz = transpose(zz);
    txt =  "B_2= "+ num2str(k1_vec(d3) /lam1_vec(vec_R2(i)))+ ' \lambda_2: '+ num2str(lam1_vec(d2));
    plot(r1,zz,'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
legend(Orientation="vertical", Location="bestoutside")
txt3 = "Final Heedless "+ "\lambda_1: "+ num2str(lam1_vec(d1)) +" k_2: "+ num2str(k1_vec(d3));
title(txt3)
xlabel('B_1') 
ylabel("Heedless") 
hold off
fontsize(13,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 16]);
set(gcf, 'PaperSize', [21 16]); % dimension on x axis and y axis resp.
 print(gcf,titles,'-dpdf')













 %% OTHERS PRINTS NOT MORE IN USE


%% Compliant dynamic first example
fig = fig+1;
figure(fig)
titles = "compliant_evolution_vari_1.pdf";
%1,1
d1 = 18; d2 = 5;d3 =10;
hold on
for i = 1:20
    d4 = i;
    aga = reshape(Compliant(d1,d2,d3,d4,:),1,2001,[]);
    zz = aga; %zz = transpose(zz);
    txt = "\lambda_2: "+ num2str(lam2_vec(d4))+" B_2: "+num2str(k2_vec(d2)/lam2_vec(d4));
    plot(zz(1:1500),'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
legend(Orientation='vertical', Location='eastoutside')
txt3 = "Compliant evolution "+ "B_1: "+num2str(k1_vec(d1)/lam1_vec(d3))+ "\lambda_1: "+ num2str(lam1_vec(d3)) +" k_1: "+ num2str(k1_vec(d1))+" k_2: "+ num2str(k2_vec(d2));
title(txt3)
xlabel('time') 
ylabel("Against ") 
hold off
fontsize(13,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 23 16]);
set(gcf, 'PaperSize', [24 16]); % dimension on x axis and y axis resp.
  % print(gcf,titles,'-dpdf')
 %% Compliant dynamic second example
fig = fig+1;
figure(fig)
titles = "compliant_evolution_vari_2.pdf";
%1,1
d1 = 5; d2 = 3;d3 =5;
hold on
for i = 1:20
    d4 = i;
    aga = reshape(Compliant(d1,d2,d3,d4,:),1,2001,[]);
    zz = aga; %zz = transpose(zz);
    txt = "\lambda_2: "+ num2str(lam2_vec(d4))+" B_2: "+num2str(k2_vec(d2)/lam2_vec(d4));
    plot(zz(1:1500),'LineWidth',1.5 ,'DisplayName',txt)
end
legend show
legend(Orientation='vertical', Location='eastoutside')
txt3 = "compliant evolution "+ "B_1: "+num2str(k1_vec(d1)/lam1_vec(d3))+ "\lambda_1: "+ num2str(lam1_vec(d3)) +" k_1: "+ num2str(k1_vec(d1))+" k_2: "+ num2str(k2_vec(d2));
title(txt3)
xlabel('time') 
ylabel("Against ") 
hold off
fontsize(13,"points")
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 23 16]);
set(gcf, 'PaperSize', [24 16]); % dimension on x axis and y axis resp.
 % print(gcf,titles,'-dpdf')
