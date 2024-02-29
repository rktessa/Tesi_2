% try to read covif italia csv file

M = readmatrix('covid_italia.csv');

figure(3)
plot(M(:,7)/60e6, LineWidth=2)
legend(' Totali positivi covid italia')