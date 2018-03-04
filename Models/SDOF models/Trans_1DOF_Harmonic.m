%% Plots the transient response of a mass spring damper system excited by a harmonic signal
% Carl Howard 27/1/2018
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:3;      % [s] time scale

%% Initial conditions: x(0) = 0, x'(0)=0
initial_x    = 0;
initial_dxdt = 0.2;

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45( @rhs, t, [initial_x initial_dxdt],options );
%% Import comparison 
M='U:\_PhD\Datathief\Figure3-InmanEngvib\figure3.csv';
data=csvread(M,1,0);
%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),'b*',data(:,1),data(:,2),'r');
xlabel('t'); ylabel('x');
set(plot1,'LineWidth',2)
title('Time Series')
legend 'ODE45' 'DataThief'
grid on
set(gca,'fontsize',20) 

%%% Calculate the PSD of the time series
FFTsize=1024;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x(:,1),hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');

%%% TF estimate 
% tfestimate(input signal,x(:,1))

%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x)
        mass1=10;		% [kg]
        stiff1=1000;    % [N/m]
        damp=0.00000001;     % [Ns/m] keep as a small number to fix solver errors
        f=23;            % [N] amplitude of driving force

        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - (stiff1/mass1)*x(1) + (f/mass1)*cos(2*sqrt(stiff1/mass1)*t);

        dxdt=[dxdt_1; dxdt_2];
    end
