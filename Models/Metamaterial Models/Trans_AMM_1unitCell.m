%% Plots the transient response of a forced 2DOF mass spring damper system 
clear all
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:60;      % [s] time scale

mass1=0.1;		% [kg]
mass2=mass1*0.5;
stiff1=1000;    % [N/m]
stiff2=1.5*stiff1;
w2=sqrt(stiff2/mass2)/(2*pi);

%% Initial conditions: x(0) = 0, x'(0)=0 ,y(0)=0, y'(0)=0
initial_x    = 0e-3;
initial_dxdt = 0;
initial_y    = 0e-3;
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
plot1=plot(t,x(:,1),t,x(:,3));
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
grid on
legend 'ODE45 mode1' 'ODE45 mode2' 
set(gca,'fontsize',20) 
%% Calculate the PSD of the time series
FFTsize=1024;
[PSD_theory_f10Hz,F_theory_f10Hz]=pwelch(x,hanning(FFTsize),[],FFTsize,fs);
figure
p3=plot(F_theory_f10Hz,10*log10(abs(PSD_theory_f10Hz)));
xlabel('Frequency (Hz)');
ylabel('Displacement (dB re 1m)');
title('PSD of Displacement of Mass');
%% Normalised Displacement to mass 1

% max value of m1 to normalise
max_m1=max(x(:,1));

m1_x=x(:,1);
m2_x=x(:,3);
% normalise both 
m1_nx=m1_x./max_m1;
m2_nx=m2_x./max_m1;

% figure
figure
plot1=plot(t,(m1_nx),'b',t,(m2_nx),'r');
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('Normalised displacemtn, x');
title('Time Series')
grid on
legend 'ODE45 mass1' 'ODE45 mass2' 
set(gca,'fontsize',20) 
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

%% Amplitude and resonance 
% Xamp=zeros(length(omega),2)';
% Wn=zeros(length(omega),1)';
% for jj=1:length(omega)
%     Wn(jj)=[omega(jj)]; %omega(jj)];
%     [amp]=forced_vibration(K,M,f,Wn(jj));
%     
%     Xamp(jj)=amp;
% end
% figure
% plot(omega,(Xamp))

%% Transmission
% U=Xn/X1
U=abs(x(:,3))./abs(x(:,1));
figure
% plot(t,20*log10(U))
periodogram(x(:,3))

%% TF estimate
% takes in input and output signal
% bandpass filter input signal
wind = kaiser(length(x(:,1)),35);
w=5;
signal=1*sin(2*pi*w*t);
% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(signal,x(:,1),[],[],[],1/(dt));
% plot
figure
% hold on
graph=plot(frequencies,(20*log10(abs(txy))),'r');
% axis([0 12 -Inf Inf])
grid on
% title('Transfer function of metamaterial configurations','FontSize',14)
xlabel('Frequency, \omega','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',14)
% axis([0 10 -Inf 0])
legend({'linear AMM','noise insulation panel','nonlinear case 1'},'FontSize',14)
set(graph,'LineWidth',1.5);
%% semilog txy
figure
% hold on
graph1=semilogx(frequencies,abs(txy),'b');
% axis([0 3 0 Inf])
grid on
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',20)
legend({'linear AMM numerical result','nonlinear case 1 AMM numerical result','phononic crystal numerical result','nonlinear case 2 AMM numerical result'})
set(graph1,'LineWidth',2);
%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
function dxdt=rhs(t,x)
        mass1=0.1;		% [kg]
        mass2=mass1*0.5;
        stiff1=1000;    % [N/m]
        stiff2=1.5*stiff1;
        damp1=0.002;     % [Ns/m] keep as a small number to fix solver errors
        damp2=0.002;
        f=1; %*(stepfun(t,0)-stepfun(t,0.01));
        w=5; % Hz, forcing frequency 
     
        %---------------------------------------
        % first unit cell
        % first mass
        dxdt_1 = x(2);
        dxdt_2 = -((2*damp1+damp2)/mass1)*x(2)- ((2*stiff1+stiff2)/mass1)*x(1) +(stiff2/mass1)*x(3)+(damp2/mass1)*x(4)...
          +(f/mass1)*sin(2*pi*w*t);
        % second mass
        dydt_1= x(4);
        dydt_2= -(stiff2/mass2)*x(3) - (damp2/mass2)*x(4) + (stiff2/mass2)*x(1) + (damp2/mass2)*x(2);
        %---------------------------------------
                
        % final solution 
        dxdt=[dxdt_1; dxdt_2; dydt_1; dydt_2];
end
%% forced vibration 
function X = forced_vibration(K,M,f,omega)
% Function to calculate steady state amplitude of
% a forced linear system.
% K is nxn the stiffness matrix
% M is the nxn mass matrix
% f is the n dimensional force vector
% omega is the forcing frequency, in radians/sec.
% The function computes a vector X, giving the amplitude of
% each degree of freedom
h=eig(K,M);
s=h/(2*pi);
r=(omega./s(2));
% X = (K-M*omega^2)\f;
X = 1/sqrt((1-(r).^2).^2+(2*0.000002*r).^2);


end