clc;
clear;
close all;
%% Figure of the epidemic over the year

%year= [165, 541, 735, 1346, 1520, 1548,1570, 1629, 1656, 1772, 1846, 1855,  1889, 1918, 1918,  1957, 1968, 1981, 2019 ];
year = [730, 800, 850, 1146, 1220, 1288, 1325, 1400, 1430, 1540, 1630, 1665, 1700, 1750, 1830, 1880, 1940, 2000 ];
%x_year = [165 541 735 1346 1520 1548 1570 1629 1656 1772 1846 1855 1889 1917 1919 1957 1968 1981 2019 ];

x_year = [730 800 850 1146 1220 1288 1325 1400 1430 1540 1630 1665 1700 1749 1750 1830 1880 1940 2000 ];
epidemic_deaths = [7.5,50,2, 50, 7, 16, 2.5, 1, 1.25, 2, 1, 15,1, 50, 2, 2, 43, 20 ];
y_deaths = [7.5 50 2 50 7 16 2.5 1 1.25 2 1 15 1 50 3 2 2 43 20];

epidemic_name = {'  Antonine Plague', '    Plague of Justinian', '    Japanese smallpox', '    Black death', '    Mexico smallpox', '    Cocoliztli epidemic',...
    '    Cocoliztli epidemic II','   Italian plague', '    Naples plague', '    Persian plague', '    Cholera', '    The third plague', '    Influenza',...
    '    Spanish Flu', '    Russian typhus','    Influenza pandemic', '    Hong Kong flu', '    HIV/AIDS', '    Covid-19' };

x_2 = [100,200,1750];
epi = [1, 1, 3]
figure(1)
sz = 25;
hold on
c = linspace(1,10,length(year));
% scatter(year,epidemic_deaths,sz,c,'filled', LineWidth=5)

bar(year,diag(epidemic_deaths),'stacked', 'BarWidth',1.3)
bar(x_2, epi,'BarWidth',0.4, FaceColor='red')

 colororder("reef")
box on
ylim([0, 85])
xlim([695, 2075])
xticks(x_year)
xticklabels({'165', '541', '735', '1346', '1520', '1548', '1570', '1629', '1656', '1772', '1846', '1855', '1889', '1918', ' ', '1957', '1968', '1981', '2019' })
xtickangle(60)
h = text(x_year,y_deaths,epidemic_name, HorizontalAlignment="left")
set(h,'Rotation',90);
fontsize(26,"points")
xlabel('Year')
ylabel('Death toll [millions]')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [-3 0 63 28]);
set(gcf, 'PaperSize', [60 28]); % dimension on x axis and y axis resp.
print(gcf,'-dpdf', ['worst_epidemic.pdf'])
%xticklabels(epidemic_name)
