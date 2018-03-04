%% Plots the transient response of a free mass spring damper system
% Carl Howard 27/1/2018


%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:12;      % [s] time scale

%% Initial conditions: x(0) = 5mm, x'(0)=2root2 mm/s
initial_x    = 5e-4;
initial_dxdt = 2*sqrt(2)*10^-3;
x0=initial_x;
v0=initial_dxdt;
wn=2;

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45( @rhs, t, [initial_x initial_dxdt],options );


%% Import comparison 
M='U:\_PhD\Datathief\Figure6-Inman EngVib\figure6.csv';
data=csvread(M,1,0);


%% Plot the results
% theoretical result
xt=(sqrt(wn^2*x0^2+v0^2)/wn)*sin(wn*t+atan((wn*x0)/v0));
% Plot the time series
figure
plot1=plot(t,x(:,1),'b*',t,xt,'r',data(:,1),data(:,2),'g');
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45' 'Theoretical' 'DataThief' 
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
        mass1=0.1;		% [kg]
        stiff1=1000;    % [N/m]
        damp=0.0001;     % [Ns/m] keep as a small number to fix solver errors
        f=0; %*sin(2*pi*10*t);            % [N] amplitude of driving force

        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - (4)*x(1) + (f/mass1);

        dxdt=[dxdt_1; dxdt_2];
    end
