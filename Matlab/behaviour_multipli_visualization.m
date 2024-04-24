clc;
close all;
% clear all;
%% Visualization of the multiple behaviour data
 % load("data/multilpiAgainst.mat")
% load("data/multilpiBehaviour.mat")
   % load("data/multilpiCareless.mat")
 % load("data/multilpiCompliant.mat")
%% 
% Against = Compliant;
 Against = Careless;
%% Against plot
% plot Against wrt the variation of the lambda2 coefficient
i = 0;

d1 = 20;
d2 = 2001;
% length(Against(3,4,4,2,:))
aga = reshape(Against(3,4,4,2,:),1,d2,[]);
% Parameter under analysis:
k1_vec   = linspace(0.1, 0.99, 30);    
k2_vec   = linspace(0.1, 0.99, 30);  
lam1_vec = linspace(1/5, 1/30, 20);
lam2_vec = linspace(1/5, 1/30, 20);


% Variation of k1
i = i+1;
figure(i)
t = tiledlayout(1,3);
title(t,'Different simualation of behaviour model')
    %1,1
    nexttile
    hold on
    l1 = 10; l2 = 8; l3 =10; l4 = 10;
    for ii = 1:30
        aga = reshape(Against(ii,l2,l3,l4,:),1,d2,[]);
        txt =  "\k1-"+ num2str(ii)+ "= "+ num2str(k1_vec(ii)); 
        plot(aga, 'LineWidth',1.5 ,'DisplayName',txt)
    end
    hold off
%  txt2 = {"k1= "+ num2str(k1_vec(l1)) +"k2= " + num2str(k2_vec(l2))} ;
    title("  k2= " + num2str(k2_vec(l2))+ "  \lambda1= "+ num2str(lam1_vec(l3))+ "  \lambda2= "+ num2str(lam2_vec(l4)) )
    legend show

    %1,2
    nexttile
    hold on
    l1 = 10; l2 = 17; l3 =10; l4 = 10;
    for ii = 1:30
        aga = reshape(Against(ii,l2,l3,l4,:),1,d2,[]);
        txt =  "\k1-"+ num2str(ii)+ "= "+ num2str(k1_vec(ii)); 
        plot(aga, 'LineWidth',1.5 ,'DisplayName',txt)
    end
    hold off
    title("  k2= " + num2str(k2_vec(l2))+ "  \lambda1= "+ num2str(lam1_vec(l3))+ "  \lambda2= "+ num2str(lam2_vec(l4)) )
    legend show

     %1,3
    nexttile
    hold on
    l1 = 10; l2 = 25; l3 =10; l4 = 10;
    for ii = 1:30
        aga = reshape(Against(ii,l2,l3,l4,:),1,d2,[]);
        txt =  "\k1-"+ num2str(ii)+ "= "+ num2str(k1_vec(ii)); 
        plot(aga, 'LineWidth',1.5 ,'DisplayName',txt)
    end
    hold off
    title("  k2= " + num2str(k2_vec(l2))+ "  \lambda1= "+ num2str(lam1_vec(l3))+ "  \lambda2= "+ num2str(lam2_vec(l4)) )
    legend show



    %%%% varying lambdas
    % Variation of k1, k2 medium 
i = i+1;
figure(i)
t = tiledlayout(1,3);
title(t,'Different simualtion of behaviour model')
    %1,1
    nexttile
    hold on
    l1 = 10; l2 = 9; l3 =15; l4 = 5; %high lam1, low lam2
    for ii = 1:30
        aga = reshape(Against(ii,l2,l3,l4,:),1,d2,[]);
        txt =  "\k1-"+ num2str(ii)+ "= "+ num2str(k1_vec(ii)); 
        plot(aga, 'LineWidth',1.5 ,'DisplayName',txt)
    end
    hold off
%  txt2 = {"k1= "+ num2str(k1_vec(l1)) +"k2= " + num2str(k2_vec(l2))} ;
    title("  k2= " + num2str(k2_vec(l2))+ "  \lambda1= "+ num2str(lam1_vec(l3))+ "  \lambda2= "+ num2str(lam2_vec(l4)) )
    legend show

    %1,2
    nexttile
    hold on
    l1 = 10; l2 = 9; l3 =9; l4 = 11; %equal lambdas, high k2
    for ii = 1:30
        aga = reshape(Against(ii,l2,l3,l4,:),1,d2,[]);
        txt =  "\k1-"+ num2str(ii)+ "= "+ num2str(k1_vec(ii)); 
        plot(aga, 'LineWidth',1.5 ,'DisplayName',txt)
    end
    hold off
    title("  k2= " + num2str(k2_vec(l2))+ "  \lambda1= "+ num2str(lam1_vec(l3))+ "  \lambda2= "+ num2str(lam2_vec(l4)) )
    legend show

     %1,3
    nexttile
    hold on
    l1 = 10; l2 = 9; l3 =8; l4 = 18;%low lam1, high lam2
    for ii = 1:30
        aga = reshape(Against(ii,l2,l3,l4,:),1,d2,[]);
        txt =  "\k1-"+ num2str(ii)+ "= "+ num2str(k1_vec(ii)); 
        plot(aga, 'LineWidth',1.5 ,'DisplayName',txt)
    end
    hold off
    title("  k2= " + num2str(k2_vec(l2))+ "  \lambda1= "+ num2str(lam1_vec(l3))+ "  \lambda2= "+ num2str(lam2_vec(l4)) )
    legend show
