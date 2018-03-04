%% Plots the transient response of a forced mass spring damper system 
clear all
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:10;      % [s] time scale
tspan=t;
%% Initial conditions: x(0) = 0, x'(0)=0 
initial_x    = 0;
initial_dxdt = 0;
%% Whitenoise signal generation
whitenoise = -1 + (2).*rand(length(t),1);
figure
plot(t, whitenoise)
%%%% Use lowpass signal filter
% global WN
WN=doFilter(whitenoise);
plot(t,WN)
tpsi=t;

fs_white=1000;
FFTsize=1024;
[P_LP50Hz_force,Freq_LP50Hz_force]=pwelch(WN,hanning(FFTsize),[],FFTsize,fs_white);
figure
plot(Freq_LP50Hz_force,10*log10(abs(P_LP50Hz_force)))

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode23s(@(t,x) rhs(t,x,tpsi,WN),tspan,[initial_x initial_dxdt],options );

%% Plot the results
% Plot the time series
figure
plot(t,x(:,1));
xlabel('t'); ylabel('x');
title('Time Series')

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
    function dxdt=rhs(t,x,tpsi,WN)
%         global WN
        mass1=0.1;		% [kg]
        stiff1=1000;    % [N/m]
        damp=0.0001;     % [Ns/m] keep as a small number to fix solver errors
        signal=interp1(tpsi,WN,t);
        f=1*signal; %*sin(2*pi*10*t);            % [N] amplitude of driving force
        
        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - (stiff1/mass1)*x(1) + (f/mass1);

        dxdt=[dxdt_1; dxdt_2];
    end