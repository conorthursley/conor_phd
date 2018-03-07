%% Plots the transient response of a forced 2DOF mass spring damper system 

%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:300;      % [s] time scale

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0;
initial_dxdt = 0;
initial_y    = 0;
initial_dydt = 0;
z=[initial_x initial_y initial_dxdt initial_dydt];
%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt,'RelTol',1e-4);
[t,x]=ode45(@rhs, t, z, options);
%% Import comparison 
M='U:\_PhD\Datathief\2DOF_impulse_InmanEngVib\figure10_3_x1.csv';
data=csvread(M,1,0);
N='U:\_PhD\Datathief\2DOF_impulse_InmanEngVib\figure10_3_x2.csv';
data1=csvread(N,1,0);
%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),'r*',t,x(:,2),'b+',data(:,1),data(:,2),'g',data1(:,1),data1(:,2),'c');
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

%% FFT
figure
dt=abs(mean(diff(t)));  %average time step done 
Fs=1/(dt);
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(t); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(x(:,1),m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft1); 
f=Fs*(0:(n/2))/n;
freq1=fr(1:floor(m/2));
P1=fourier(1:floor(m/2));
plot(freq1,P1)
% plot(flop,abs(y),'LineWidth',2)
title('FFT')
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',20)

%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x)
        mass1=2;		% [kg]
        mass2=1;
        stiff1=2;    % [N/m]
        stiff2=1;
        damp1=0.3-0.05;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.05;
        w=0; % driving frequency
        % matrix allocation
        M=[mass1 0; 0 mass2];   % mass matrix
        C=[damp1+damp2 -damp2;-damp2 damp2];  % Damping matrix
        K=[stiff1+stiff2 -stiff2;-stiff2 stiff2]; %stiffness matrix
        B=[0;1]; % forcing input vector
        t1=1; t2=1.1; % pulse duration
        
        
        
        % state space formulation
        A1=[zeros(2) eye(2); -inv(M)*K -inv(M)*C];
        f=inv(M)*B;
        
        % solution 
        dxdt=A1*x+[0;0;f]*(stepfun(t,t1)-stepfun(t,t2));
    end