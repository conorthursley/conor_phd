%% Plots the transient response of a forced 2DOF mass spring damper system 

%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:20;      % [s] time scale

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 1e-3;
initial_dxdt = 0;
initial_y    = 0;
initial_dydt = 0;
x=[initial_x initial_dxdt];
y=[initial_y initial_dydt];
z=[x y];
% initial disp and velocity
xt0=[initial_x initial_y]';
dxdt0=[initial_dxdt initial_dydt]';
%% Modal Decoupling 
k1=24;
k2=3;
m1=9;
m2=1;
%---------------------------
M=[m1 0;0 m2];
K=[k1+k2 -k2;-k2 k2];
%---------------------------
% spectral matrix
K_dash=M^-0.5*K*M^-0.5;
[modes,nat]=eig(K_dash);
W=diag(nat);
w1=sqrt(W(1));
w2=sqrt(W(2));
% supplementary matrices
P=modes;
S=M^-0.5*P;
S_dash=(P'*M^0.5)';
% initial disp and velo -> uncoupled
rxt=S_dash*xt0;
rdxdt=S_dash*dxdt0;

r1=((sqrt(w1^2*rxt(1)^2+rdxdt(1)^2))/w1)*sin(w1*t+atan((w1*rxt(1))/rdxdt(1)));
r2=((sqrt(w2^2*rxt(2)^2+rdxdt(2)^2))/w2)*sin(w2*t+atan((w2*rxt(2))/rdxdt(2)));
r=[r1;r2];

% final solution in physical coordinates
X=S*r;

figure
plot(t,X)
% 
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
plot1=plot(t,x(:,1),t,x(:,3));
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
        mass1=9;		% [kg]
        mass2=1;
        stiff1=24;    % [N/m]
        stiff2=3;
        damp1=0;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0;
        
     
        dxdt_1 = x(2);
        dxdt_2 = ((stiff2)/mass1)*x(3) -((stiff1+stiff2)/mass1)*x(1);
        dydt_1= x(4);
        dydt_2= ((stiff2)/mass2)*x(1) -((stiff2)/mass2)*x(3);
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2];
    end