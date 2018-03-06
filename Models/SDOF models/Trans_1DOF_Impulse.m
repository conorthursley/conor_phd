%% Plots the transient response of a forced mass spring damper system 
clear all
%% simulation parameters
fs=1000;        % [Hz] sampling frequency
dt=1/fs;        % [s] delta t
t=0:dt:10;      % [s] time scale
omega=logspace(-1,3,length(t));
%% Initial conditions: x(0) = 0, x'(0)=0 
initial_x    = 1e-3;
initial_dxdt = -1e-3;

%% Solve the model
options=odeset('InitialStep',dt,'MaxStep',dt);
[t,x]=ode45( @rhs, t, [initial_x initial_dxdt],options );
%% Import comparison 
M='U:\_PhD\Datathief\ForcedResponse-InmanEngVib\figure4.csv';
data=csvread(M,1,0);

%% Plot the results
% Plot the time series
figure
plot1=plot(t,x(:,1),data(:,1),data(:,2),'g');
set(plot1,'LineWidth',2)
xlabel('t'); ylabel('x');
title('Time Series')
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
Xamp=zeros(length(omega),1)';
Wn=zeros(length(omega),1)';
for jj=1:length(omega)
    Wn(jj)=omega(jj); %omega(jj)];
    [amp]=forced_vibration(4,1,1,Wn(jj));
    Xamp(jj)=amp;
end
figure
plot(omega/(2*pi),abs(Xamp))
xlabel('frequency \omega'); ylabel('Amplitude x');
title('Amplitude v Frequency ')
% legend 'ODE45' 'DataThief'
grid on
set(gca,'fontsize',20) 
%% Mass-Spring-Damper system
% The equations for the mass spring damper system have to be defined
% separately so that the ODE45 solver can call it.
    function dxdt=rhs(t,x)
        mass1=1;		% [kg]
        stiff1=4;    % [N/m]
        damp=2;     % [Ns/m] keep as a small number to fix solver errors
        if abs(t-0)<1e-3  %impulse forcing 
            f=1; %*sin(2*pi*10*t);            % [N] amplitude of driving force
        elseif abs(t-4)<1e-3
            f=-0.5;
        else                
            f=0;
        end
        
        dxdt_1 = x(2);
        dxdt_2 = -(damp/mass1)*x(2) - (stiff1/mass1)*x(1) + (f/mass1);

        dxdt=[dxdt_1; dxdt_2];
    end
%% Force Vibration Amplitude function
function X = forced_vibration(K,M,f,omega)
% Function to calculate steady state amplitude of
% a forced linear system.
% K is nxn the stiffness matrix
% M is the nxn mass matrix
% f is the n dimensional force vector
% omega is the forcing frequency, in radians/sec.
% The function computes a vector X, giving the amplitude of
% each degree of freedom
%
X = (K-M*omega^2)\f;

end