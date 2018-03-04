%% Plots the transient response of a forced 2DOF mass spring damper system 

%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:10;      % [s] time scale

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0;
initial_dxdt = 0;
initial_y    = 0;
initial_dydt = 0;
z=[initial_x initial_dxdt initial_y initial_dydt];
%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45(@rhs, t, z, options);

%% Plot the results
% Plot the time series
figure
plot(t,x(:,1),t,x(:,3));
xlabel('t'); ylabel('x');
title('Time Series')
grid

%%% Calculate the PSD of the time series
FFTsize=1024;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x,hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');


%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x)
        mass1=0.1;		% [kg]
        mass2=0.05;
        stiff1=1000;    % [N/m]
        stiff2=1500;
        damp=0.000002;     % [Ns/m] keep as a small number to fix solver errors
        f=1;
     
        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - ((stiff1+stiff2)/mass1)*x(1) +(stiff2/mass1)*x(3) + (f/mass1)*sin(2*pi*15*t);
        dydt_1= x(4);
        dydt_2= (stiff2/mass2)*x(1) - (stiff2/mass2)*x(3);

        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2];
    end