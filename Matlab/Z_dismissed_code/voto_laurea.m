% La mia media
clc;
clear all;
close all;

media_attuale_ponderata = 26.33;

media_2 = (29*6 + 27*9 + 25*6 + 29*6 + 25*6 + 26*9 + 27*6 + 28*6 +27*6 +22*6 +26*9 +24*9 +29*6 +27*6 +23*6)/102

crediti_attuali = 120-3-15;

voto_tesi = [18,19,20,21,22,23,24,25,26,27,28,29,30,31];

possibili_medie = (voto_tesi*15+ media_attuale_ponderata*crediti_attuali)/117 

possibili_medie_2 = (voto_tesi*15+ media_2*crediti_attuali)/117


voto_laurea_finale = possibili_medie*3.85
voto_laurea_2 = possibili_medie_2*3.85

figure(1)
hold on
plot(voto_tesi, voto_laurea_finale, '*r', 'LineWidth',1.5)
plot(voto_tesi, voto_laurea_2, '*b', 'LineWidth',1.5)
legend("stima app", "stima con matlab")