%% Simple 1D mass-in-mass system 
% dependicies - sys.m
% 11/09/2017 - Conor MacDonald
% GitHub comment to test
 
clear all
close all
tic
% Time specification
tspan = [0 1];

% Initial conditions
% velocity is initially 0, v1=0, v2=0
% displacement for m1 is 1e-3 from the fixed wall (L spacing)
% displacement for m2 is 1.5E-3 from the wall (5e-4 away from m1, within)
% y=[u1;v1;u2;v2]
y=[1e-3;0;1e-3;0];
%opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
% System simulation
% [t, y] = ode45(@sys, t, y1);
[t, result] = ode45(@sys, tspan, y);
toc

% u1=result(:,1); u2=result(:,2);
% 
% plot(t, u1,t,u2)


% Plots the input signal
% plot(t,H)
% xlabel('Time');
% ylabel('Input Signal');
% Plots the state
% plot(t, result)
% xlabel('Time');
% ylabel('State');

%% Plots displacement and velocity of displacement of 2 masses
u1=result(:,1);
u2=result(:,3);

figure
subplot(2,1,1);
plot(t,u1)
% Format plot
xlabel('time'); % Insert the x-axis label
ylabel('displacement'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('u1')
grid on
subplot(2,1,2);
plot(t,u2)
% Format plot
xlabel('time'); % Insert the x-axis label
ylabel('displacement'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('u2')
grid on

%% Ratio of displacement between masses
U1=abs(u1);
U2=abs(u2);
Ux=U2/U1;

figure
plot(t,Ux)

%% Pwelch function
% figure
% [pxx,freq] = pwelch(u2,500,300,500,max(t));
% plot(freq,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('Magnitude (dB)')
% title('Pwelch Function')
% grid on

%% FFT
%Single sided amplitude spectrum of U1(t)
% 
figure
subplot(2,1,1);
dt=mean(diff(t));  %average time step done in the ode45 computation
Fs=1/dt;
n=length(t);  %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft=fft(u1,m)/n; % DFT of signal
fr = (0:m-1)*(Fs/m)/10;
fourier = abs(dft);     
plot(fr(1:floor(m/2)),fourier(1:floor(m/2)))
title('Single-Sided Amplitude Spectrum of U1(t)')
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
%Single sided amplitude spectrum of U2(t)
% 
subplot(2,1,2);
dft2=fft(u2,m)/n; % DFT of signal
fourier2 = abs(dft2);     
plot(fr(1:floor(m/2)),fourier2(1:floor(m/2)))
title('Single-Sided Amplitude Spectrum of U2(t)')
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% Phase plane
%
% Plots trajectory 
figure
% U1
subplot(2,1,1);
plot(result(:,1),result(:,2))
xlabel('U_1'); % Insert the x-axis label
ylabel('dU_1/dt'); % Inserts the y-axis label
title('Phase plane') % Inserts the title in the plot
grid on
% U2
subplot(2,1,2);
plot(result(:,3),result(:,4))
xlabel('U_2'); % Insert the x-axis label
ylabel('dU_2/dt'); % Inserts the y-axis label
title('Phase plane') % Inserts the title in the plot
grid on

%% Frequency Response Function plot
% FRF = fft(output)./fft(input);
% input from H = sin(2*pi*1*t);
% output from result(:,1);
% plot(t,FRF);

%% Impulse Reaction Plot
% displacement over time [0-0.1]s

toc
