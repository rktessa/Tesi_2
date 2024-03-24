clc;
clear all;
close all;

%% 

matrix = zeros(30,20);

for i = 1:30
    for j = 1:20
        matrix(i,j) = 10*i + j;
    end
end

% matrix = transpose(matrix);
x_0 = linspace(1,30,30);
y_0 = linspace(1,20,20);
[xx, yy] = meshgrid(x_0,y_0);

figure(1)
surf(xx,yy,matrix)