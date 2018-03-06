%% Plots the transient response of a forced 2DOF mass spring damper system 

%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:10;      % [s] time scale

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 1;
initial_dxdt = -2;
initial_y    = -1;
initial_dydt = 2;
initial_z    = 0;
initial_dzdt = -1;
z=[initial_x initial_dxdt initial_y initial_dydt initial_z initial_dzdt];
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
plot1=plot(t,x(:,1),t,x(:,3),t,x(:,4));
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' 'Datathief mode1' 'Datathief mode2'
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
        mass1=1;		% [kg]
        mass2=1;
        mass3=1;
        stiff1=1;    % [N/m]
        stiff2=1;
        damp1=0;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0;
        f1=0;
        f2=0;
        f3=0;
     
        dxdt_1 = x(2);
        dxdt_2 = -(damp1/mass1)*x(2) + (damp1/mass1)*x(4) - ((stiff1)/mass1)*x(1) +(stiff1/mass1)*x(3) + (f1/mass1);
        dydt_1= x(4);
        dydt_2=(damp1/mass2)*x(2)-((damp1+damp2)/mass2)*x(4)+(damp2/mass2)*x(6) + ((stiff1)/mass2)*x(1) -((stiff1+stiff2)/mass2)*x(3) + ((stiff1)/mass2)*x(5)+(f2/mass2);
        dzdt_1 = x(6);
        dzdt_2 = (damp2/mass3)*x(4) - (damp2/mass3)*x(6) + ((stiff2)/mass3)*x(3) -(stiff2/mass3)*x(5) + (f3/mass3);

        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2; dzdt_1; dzdt_2];
    end