%% Plots the transient response of a forced mass spring damper system 
clear all
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:14*pi;      % [s] time scale
tic
%% Initial conditions: x(0) = 0, x'(0)=0 
initial_x    = 0;
initial_dxdt = 0;

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45(@(t,y) rhs(t,y), t, [initial_x initial_dxdt],options );
toc
%% Import comparison 
M='U:\_PhD\Datathief\SDOF_nonlinear_ThompsonStewartNonlinearDyn\SDOF_phase_tilt.csv';
data=csvread(M,1,0);
M1='U:\_PhD\Datathief\SDOF_nonlinear_ThompsonStewartNonlinearDyn\SDOF_nonlinear_disp.csv';
data1=csvread(M1,1,0);
%% Plot the results
figure
plot1=plot(t,x(:,1)); %,data1(:,1),data1(:,2)); %,data(:,1),data(:,2),'g');
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45' 'DataThief' 
set(gca,'fontsize',20) 

figure
plot2=plot(x(:,1),x(:,2),data(:,1),data(:,2));
set(plot2,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Phase potrait')
grid on
legend 'ODE45' 'DataThief' 
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
        stiff1=0;    % [N/m]
        stiff2=1;
        damp=0.05;     % [Ns/m] keep as a small number to fix solver errors
        f=7.5; %*sin(2*pi*10*t);            % [N] amplitude of driving force
        w=1;
        
        dxdt_1 = x(2);
        dxdt_2 = -(damp/(mass1))*x(2) - (stiff1/mass1)*x(1) -(stiff2/mass1).*x(1)^3 + (f/mass1)*cos(w*t);

        dxdt=[dxdt_1; dxdt_2];
    end