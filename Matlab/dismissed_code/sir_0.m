clear all;
close all;
clc;


N = 10000;
n = 4; r=3;
beta= 0.5/N; 
gamma = 0.2;
I = 10; R=0; S = N-I;
T = 150;
x = [S;I;R];
xplot = x; 

% Euler method for systems
for t = 1:T
    dx(1) = - beta * x(1) * x(2);
    dx(2) = beta* x(1) * x(2) - gamma*x(2);
    dx(3) = gamma* x(2);
    x(1) = x(1)+  dx(1);
    x(2) = x(2) + dx(2);
    x(3) = x(3) + dx(3);
     
  % Append the result  
  xplot = [xplot x];
end

p = 0;
p =p+1; 
figure(p)
hold on
plot( 0:T, xplot, 'LineWidth', 0.75) 
legend('S', "I", "R")
title('Evolution of SIR model')

%% Varying beta and gamma and calculate the behaviour of SIR

I = 0.001; R=0; S = 1;
beta = 0; gamma = 0;
T = 150;
x2 = [S;I;R;beta;gamma];
xplot2 = x2(1); 
%range [0,2]
range_step = 0.1;
result = [];
result2 = [];
for j = 0:20
    gamma = j* range_step;
    for l = 0:20 
        beta =  l* range_step;
        for t = 1:T
            dx2(1) = - beta * x2(1) * x2(2);
            dx2(2) = beta* x2(1) * x2(2) - gamma*x2(2);
            dx2(3) = gamma* x2(2);
            x2(1) = x2(1)+  dx2(1);
            x2(2) = x2(2) + dx2(2);
            x2(3) = x2(3) + dx2(3);
            x2(4) = beta;
            x2(5) = gamma;

          % Append the result  
          xplot2 = [xplot2; x2(1)];
        end
        result = cat(2,result,xplot2); %dimensione 2 è fare le colonne
        x2 = [S;I;R;beta;gamma];
        xplot2 = x2(1);
    end
    result2 = cat(3,result2,result);
    result = [];
end

%% Plot the results

% p =p+1; 
% figure(p)
% hold on
% plot( result2(:,5,3)) 
% legend('S')
% title('Evolution of SIR model')
% hold off
% %surface(result2(

% Taking the numeber of S at the end of the epydemic
final_S =  result2(length(result2),:,:);

% Reshape in a matrix the final size of susceptiple
res = reshape(final_S,[21,21]);

%Realize and plot the surface

beta_pl = 0:0.1:2;
gamma_pl = 0:0.1:2;

[xx, yy] = meshgrid(beta_pl,gamma_pl);


p =p+1; 
figure(p)
surf(xx,yy,res)
colorbar;
colormap default;
xlabel('beta') 
ylabel('gamma') 
title('final S for beta and gamma')

%% Solve the problem with Matlab dsolve functions

%Define the equation symbolically
syms s(t) in(t) r(t) 
N = 100001;

range_step_b = 0.001;
range_step_g = 0.1;
result_ODE = []; 
result_ODE2 = [];
for jj = 0:10
    gam = jj* range_step_g;
    for ll = 0:200 
        bet =  ll* range_step_b;

        ode1 = diff(s) == -bet*s*in/N;
        ode2 = diff(in) == bet*s*in/N - gam*in;
        ode3 = diff(r) == gam*in;

        odes = [ode1; ode2; ode3];

        [V] = odeToVectorField(odes);

        M = matlabFunction(V,'vars', {'t','Y',});

        %ode45 solver input: function to solve, t_span, initial condition 
        
        [t,y] = ode45(M,[0:1:150],[100000 1 0]);
        
        
        result_ODE = cat(2,result_ODE,y(:,1)); %dimensione 2 è fare le colonne
        
        
    end
    result_ODE2 = cat(3,result_ODE2,result_ODE);
    result_ODE = [];
end

%Plot the solution using linspace to generate 100 points in the interval [0,20]
%and deval to evaluate the solution for each point.

% p =p+1; 
% figure(p)
% hold on
% plot(t, y(:,1),'LineWidth', 0.75)
% plot(t, y(:,2),'LineWidth', 0.75)
% plot(t, y(:,3),'LineWidth', 0.75)
% legend('S','I','R')




%% Realize and plot the surface
% Taking the numeber of S at the end of the epydemic
final_S_ODE =  result_ODE2(end,:,:);


% Reshape in a matrix the final size of susceptiple
res_ODE = reshape(final_S_ODE,[201,11]);


beta_pl_ODE = (0:range_step_b:200*range_step_b);
gamma_pl_ODE = (0:range_step_g:10*range_step_g);

[xx, yy] = meshgrid(gamma_pl_ODE,beta_pl_ODE);


p =p+1; 
figure(p)
surf(xx,yy,res_ODE)
colorbar;
colormap default;
xlabel('gamma') 
ylabel('beta') 
title('final S for beta and gamma')

