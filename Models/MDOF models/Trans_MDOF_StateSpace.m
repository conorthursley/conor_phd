%% Plots the transient response of a forced 2DOF mass spring damper system 

%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:20;      % [s] time scale

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0;
initial_dxdt = 0;
initial_y    = 0;
initial_dydt = 0;
z=[initial_x initial_dxdt initial_y initial_dydt];
%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45(@rhs, t, z, options);
%% Import comparison 
% M='U:\_PhD\Datathief\MDOF_freeResponse-InmanEngVib\figure12_mode1.csv';
% data=csvread(M,1,0);
% N='U:\_PhD\Datathief\MDOF_freeResponse-InmanEngVib\figure12_mode2.csv';
% data1=csvread(N,1,0);
%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),'r',t,x(:,2),'b--');
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' 
set(gca,'fontsize',20) 
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
        mass1=4;		% [kg]
        mass2=9;
        M=[mass1 0;0 mass2];  %mass matrix
        
        stiff1=25;    % [N/m]
        stiff2=5;
        K=[stiff1+stiff2 -stiff2;-stiff2 stiff2]; %stiffness matrix
        
        F=[0.23500;2.97922]; %force vector, B
        w=2.75655; % natural frequency 
        % State space formulation
        A1=[zeros(2) eye(2); -inv(M)*K zeros(2)];
        f=inv(M)*F;
        dxdt=A1*x+[0;0;f]*sin(w*t);
     
    end