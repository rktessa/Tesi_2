classdef functionsContainer
   methods
        %% Runge Kutta second order solution of SIR      
        function [taxis,xaxis,yaxis,zaxis,i] = SIR(obj,a1,a2,p1,q11,N,beta,gamma,time,dt,i)
            x = 0.999*N; %susceptible
            y = N-0.999*N; %infected
            z = 0; % recovered
            t = 0;
            cnt=0;
            %Array creation and inititialization
            taxis=[]; taxis(1) =0; 
            xaxis=[]; xaxis(1) = x;
            yaxis=[]; yaxis(1) = y;
            zaxis=[]; zaxis(1) = z;
            while t < time
                if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
                    taxis = cat(2,taxis,t);
                    xaxis = cat(2,xaxis,x);
                    yaxis = cat(2,yaxis,y);
                    zaxis = cat(2,zaxis,z);
                end
                % step 1
                kx1 = - beta*x*y;
                ky1 = beta*x*y - gamma*y;
                % step 2
                t2 = t+p1*dt;
                x2 = x + q11*kx1*dt;
                y2 = y + q11*ky1*dt;
                kx2 = - beta*x2*y2;
                ky2 = beta*x2*y2 - gamma*y2;
                % update
                x = x + (a1*kx1+a2*kx2)*dt;
                y = y + (a1*ky1+a2*ky2)*dt;
                z = N - x - y;
                t = t + dt;
                cnt = cnt + 1;


            end


        %     i = i+1;
        %     figure(i)
        %     hold on
        %     plot(taxis,xaxis, 'r', 'linewidth',1.0 )
        %     plot(taxis,yaxis, 'b', 'linewidth',1.0 )
        %     plot(taxis,zaxis, 'g', 'linewidth',1.0 )
        %     title("SIR MODEL")
        %     legend('S','I','R')

        end

        %% Solve the problem with Matlab ODE functions
        function [taxis,xaxis,yaxis,zaxis,i] = SIR_ODE(obj,N,beta,gamma,time,dt,i)

            s0 = 0.999*N; % susceptible
            i0 = 1-0.999*N;  % infected
            r0 = 0; % recovered
            tspan = [0:dt:time];
            y0 = [s0,i0,r0];
            pars = [beta, gamma];

            %ode45 solver input: function to solve, t_span, initial condition 

            [t,y] = ode15s(@sir_rhs, tspan, y0, [], pars);

            cnt=0;
            p  =0;
            %Array creation and inititialization
            taxis=[]; taxis(1) =t(1); 
            xaxis=[]; xaxis(1) = y(1,1);
            yaxis=[]; yaxis(1) = y(1,2);
            zaxis=[]; zaxis(1) = y(1,3);


            while p < length(t) %Scorro tutta la soluzione
                if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result

                    taxis = cat(2,taxis,t(cnt));
                    xaxis = cat(2,xaxis,y(cnt,1));
                    yaxis = cat(2,yaxis,y(cnt,2));
                    zaxis = cat(2,zaxis,y(cnt,3));
                end
            cnt = cnt + 1;  
            p = p+1;
            end

        %      i = i+1;
        %     figure(i)
        %     hold on
        %     plot(taxis,xaxis, 'r', 'linewidth',1.0 )
        %     plot(taxis,yaxis, 'b', 'linewidth',1.0 )
        %     plot(taxis,zaxis, 'g', 'linewidth',1.0 )
        %     title("SIR MODEL")
        %     legend('S','I','R')

        end

        %function used to the ODE45 solver to compute correctly the system
        %evolution
        function f = sir_rhs(obj,t,y,pars)
        f = zeros(3,1);
        f(1) = -pars(1)*y(1)*y(2);
        f(2) = pars(1)*y(1)*y(2) - pars(2)*y(2);
        f(3) = pars(2) * y(2);
        end

        %% multiple SIR execution varying beta and gamma
        function [I_max, T_imax,R0] = multipleSIR(obj,d,step)
            I_max = zeros(d,d);
            T_imax = zeros(d,d);
            R0 = zeros(d,d);
            N = 1;
            time =70;
            dt=0.001; %un millesimo di secondo di dt
            a1 = 1/2; a2=1/2; p1 =1; q11 =1;

            for i = 1:d
                beta = i*step/N;
                for j = 1:d
                    gamma = j*step;

                    %Esecuzione della trasformazione, variando beta e gamma
                    [taxis,xaxis,yaxis,zaxis,i] = SIR(obj,a1,a2,p1,q11,N,beta,gamma,time,dt,i);

                    %Compute max value of infected
                    [max_i, index] = max(yaxis);
                    time_i_max = taxis(index);
                    I_max(i,j) = max_i;
                    T_imax(i,j) = time_i_max;
                    R0(i,j) = (beta/gamma)*N;

                end
            end
        end

        %% Function that calculate the percentage difference between two values
        function [perc_diff] = percentage_difference(obj,v1,v2)
        perc_diff = abs(v1-v2)/((v1+v2)/2) *100;
        end
     
        %% Simple SIR to produce plots for rapid experiments
        function [i, taxis, xaxis, yaxis, zaxis] = SIR_RK_figure(obj,N,beta,gamma,time,dt,i)
            a1 = 1/2; a2=1/2; p1 =1; q11 =1;
            x = 0.999*N; % susceptible normalized
            y = N-0.999*N;  % infected
            z = 0; % recovered
            t = 0;
            cnt=0;
            
            max_I = x+y-gamma/beta*(1+log(x*beta/gamma))
            %Array creation and inititialization
            taxis=[]; taxis(1) =0; 
            xaxis=[]; xaxis(1) = x;
            yaxis=[]; yaxis(1) = y;
            zaxis=[]; zaxis(1) = z;
            while t < time
                if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
                    taxis = cat(2,taxis,t);
                    xaxis = cat(2,xaxis,x);
                    yaxis = cat(2,yaxis,y);
                    zaxis = cat(2,zaxis,z);
                end
                % step 1
                kx1 = - beta*x*y;
                ky1 = beta*x*y - gamma*y;
                % step 2
                t2 = t+p1*dt;
                x2 = x + q11*kx1*dt;
                y2 = y + q11*ky1*dt;
                kx2 = - beta*x2*y2;
                ky2 = beta*x2*y2 - gamma*y2;
                % update
                x = x + (a1*kx1+a2*kx2)*dt;
                y = y + (a1*ky1+a2*ky2)*dt;
                z = N - x - y;
                t = t + dt;
                cnt = cnt + 1;


            end


            i = i+1;
            figure(i)
            hold on
            plot(taxis,xaxis, 'r', 'linewidth',1.0 )
            plot(taxis,yaxis, 'b', 'linewidth',1.0 )
            plot(taxis,zaxis, 'g', 'linewidth',1.0 )
            title("SIR MODEL")
            legend('S','I','R')
            txt = {['beta: ' num2str(beta)],['gamma: ' num2str(gamma)]};
            text(50,6000,txt)
        end
                
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Opinion functions
        %% Runge Kutta second order solution of Opinion     
        function [taxis,xaxis,yaxis,zaxis,i] = Opinion_RK(obj,a1,a2,p1,q11,N,alpha,alpha_prime,time,dt,i)
            
            x = 1-0.998*N;  % Compliant
            y = 1-0.998*N;  % Anti
            z = 0.998*N;% Careless
            t = 0;
            cnt=0;
            %Array creation and inititialization
            taxis=[]; taxis(1) =0; 
            xaxis=[]; xaxis(1) = x;
            yaxis=[]; yaxis(1) = y;
            zaxis=[]; zaxis(1) = z;
            while t < time
                if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result
                    taxis = cat(2,taxis,t);
                    xaxis = cat(2,xaxis,x);
                    yaxis = cat(2,yaxis,y);
                    zaxis = cat(2,zaxis,z);
                end
                % step 1
                kx1 =  alpha*x*z;
                ky1 = alpha_prime*y*z;
                kz1 = - alpha*x*z - alpha_prime*y*z;
                % step 2
                t2 = t+p1*dt;
                x2 = x + q11*kx1*dt;
                y2 = y + q11*ky1*dt;
                z2 = z + q11*kz1*dt;
                kx2 = alpha*x2*z2;
                ky2 = alpha_prime*y2*z2;
                kz2 = - alpha*x2*z2 - alpha_prime*y2*z2;
                % update
                x = x + (a1*kx1+a2*kx2)*dt;
                y = y + (a1*ky1+a2*ky2)*dt;
                z = z + (a1*kz1+a2*kz2)*dt;
                t = t + dt;
                cnt = cnt + 1;


            end
        end

        %% Solve the problem with Matlab ODE functions
        function [taxis,xaxis,yaxis,zaxis,i] = Opinion_ODE(obj,N,alpha,alpha_prime,time,dt,i)

            Ca = 0.998*N;   % careless
            Co = 1-0.998*N; % compliant
            A =  1-0.998*N; % anti social rules
            tspan = [0:dt:time];
            y0 = [Co,A,Ca];
            pars = [alpha, alpha_prime,N];

            %ode45 solver input: function to solve, t_span, initial condition 

            [t,y] = ode45(@opi_eqs, tspan, y0, [], pars);

            cnt=0;
            p  =0;
            %Array creation and inititialization
            taxis=[]; taxis(1) =t(1); 
            xaxis=[]; xaxis(1) = y(1,1);
            yaxis=[]; yaxis(1) = y(1,2);
            zaxis=[]; zaxis(1) = y(1,3);


            while p < length(t) %Scorro tutta la soluzione
                if mod(cnt,100) == 0 && cnt ~=0 %every 100 millisecond I save the result

                    taxis = cat(2,taxis,t(cnt));
                    xaxis = cat(2,xaxis,y(cnt,1));
                    yaxis = cat(2,yaxis,y(cnt,2));
                    zaxis = cat(2,zaxis,y(cnt,3));
                end
            cnt = cnt + 1;  
            p = p+1;
            end

        %     i = i+1;
        %     figure(i)
        %     hold on
        %     plot(taxis,xaxis, 'r', 'linewidth',1.0 )
        %     plot(taxis,yaxis, 'b', 'linewidth',1.0 )
        %     plot(taxis,zaxis, 'g', 'linewidth',1.0 )
        %     title("OPINION MODEL")
        %     legend('Co','A','Ca')
        %     txt = {['alpha: ' num2str(alpha)],['alpha_{prime}: ' num2str(alpha_prime)]};
        %     text(70,6000,txt)

        end

        % Function used to the ODE45 solver to compute correctly the system
        %evolution
        function f = opi_eqs(obj,t,y,pars)
        % 1- alpha - alpha' < 1 Condition
        f = zeros(3,1);
        f(1) = pars(1)*y(1)*y(3) ;                       % Co' = alpha *Co*Ca
        f(2) = pars(2)*y(2)*y(3) ;                       % A'  = alpha'*A*Ca
        f(3) = - pars(1)*y(1)*y(3) - pars(2)*y(2)*y(3); %+pars(3)*y(3); % Ca' = -alpha*Co*Ca -alpha'*A*Ca 
        end

        function [Co_max, A_max,max_gr_Co, rCo_m, rA_m] = multiple_Opinion(obj,d,step)

            Co_max = zeros(d,d);
            A_max = zeros(d,d);
            max_gr_Co = zeros(d,d);
            rCo_m = zeros(d,d);
            rA_m = zeros(d,d);
            N = 1;
            time =2500;
            dt=0.001; %un millesimo di secondo di dt
            a1 = 1/2; a2=1/2; p1 =1; q11 =1;

            for i = 1:d
                alpha = i*step/N;
                for j = 1:d
                    alpha_prime = j*step/N;

                    %Esecuzione della trasformazione, variando beta e gamma

                    % Usando ODE
                    %[taxis,xaxis,yaxis,zaxis,i] = Opinion_ODE(N,alpha,alpha_prime,time,dt,i);

                    % Usando Runge Kutta
                    [taxis,xaxis,yaxis,zaxis,i] = Opinion_RK(obj,a1,a2,p1,q11,N,alpha,alpha_prime,time,dt,i);
                    %Compute max value of infected
                    max_Co = max(xaxis);
                    max_A = max(yaxis);
                    Co_max(i,j) = max_Co;
                    A_max(i,j) = max_A;
                    gr_Co = gradient(xaxis); %Calcolo il valore massimo del gradiente di ogni curva
                    max_gr_Co(i,j) = max(gr_Co);
                    rCo_m(i,j) = risetime(xaxis,taxis,PercentReferenceLevels=[10 90]);
                    rA_m(i,j) = risetime(yaxis,taxis,PercentReferenceLevels=[10 90]);
            
                end
            end
        end
        
        
        %% Minimi quadrati
        function [a, b] = LS_linear(obj, x, y)
            %with x and y data points 
            m_x = mean(x);
            m_y = mean(y);
            v_x = var(x);
            v_y = var(y);
            c_xy = cov(x,y);
            rms_x_2 = rms(x)^2;
            rms_y_2 = rms(y)^2;
            m_xy = mean(x.*y);
            
            M = [ 1, m_x;
                  m_x, rms_x_2];
            v = [m_y; m_xy];
            
            sol = v\M;
            
            b = sol(1);
            a = sol(2);
        end
   end
end