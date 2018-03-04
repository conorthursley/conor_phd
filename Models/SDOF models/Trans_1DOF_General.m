%% Plots the transient response of a forced mass spring damper system 
clear all
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:15*pi;      % [s] time scale

%% Initial conditions: x(0) = 0, x'(0)=0 
initial_x    = 3;
initial_dxdt = 4;

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45( @rhs, t, [initial_x initial_dxdt],options );

%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1));
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45' 
set(gca,'fontsize',20) 
%%% Plot Phase plane
figure
plot2=plot(x(:,1),x(:,2));
set(plot1,'LineWidth',2)
xlabel('x'); ylabel('x^{.}');
title('Time Series')
grid on
legend 'ODE45' 
set(gca,'fontsize',20) 

%%% Calculate the PSD of the time series
FFTsize=1024;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x(:,1),hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');


%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x)
        mass1=1;		% [kg]
        stiff1=1;    % [N/m]
        damp=0.05;     % [Ns/m] keep as a small number to fix solver errors
        f=7.5; %*sin(2*pi*10*t);            % [N] amplitude of driving force
      
        
        dxdt_1 = x(2);
        dxdt_2 = -(damp/(mass1))*x(2) - (stiff1/mass1).*x(1)^3 + (f/mass1)*cos(t);

        dxdt=[dxdt_1; dxdt_2];
    end