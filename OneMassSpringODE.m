%% Simple 1D mass-spring system with dynamic system outlined by the funciton in the ode45 
% 11/09/2017 - Conor MacDonald
%---------------------------------------------------
% clear all
% close all
tic
% Time specification
tspan = [0 100];
%---------------------------------------------------
% Initial conditions
% velocity is initially 0, v1=0, v2=0
% displacement for m1 is 1e-3 from the fixed wall (L spacing)
% displacement for m2 is 1E-3 from the wall (5e-4 away from m1, within)
% y=[u1;v1;u2;v2]
% y=[1e-3;0;1e-3;0];
y=[0 0];
%---------------------------------------------------
% Parameters
k1=1e3; %N/m
m1=0.1; %kg
k2=10;
%---------------------------------------------------
% harmonic input frequency 
% expressed in Hz and then converted to rad/s in the function
input =[5 1]; %Hz


opts = odeset('RelTol',1e-5,'AbsTol',1e-7, 'OutputFcn',@odeplot); %, 'Mass', mass, 'Events', @events);
%% System simulation
% [t, y] = ode45(@sys, t, y1);
[t,result] = ode15s(@(t,y)DuffingVal(t,y,k1,m1,k2), tspan, y,opts); %,te,ye,ie
toc


%% Plots displacement and velocity of displacement of 2 masses
u1=result(:,1);
u2=result(:,2);
%---------------------------------------------------
figure
plot(t,u1)
% Format plot
xlabel('time'); % Insert the x-axis label
ylabel('displacement'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('u1')
grid on
%---------------------------------------------------


%% Ratio of displacement between masses
% U1=abs(u1);
% U2=abs(u2);
% Ux=U2/U1;
% 
% figure
% plot(t,Ux)
% 
%% Pwelch function
% figure
% [pxx,freq] = pwelch(u1,500,300,500,max(t));
% plot(freq,10*log10(pxx))
% xlabel('Frequency (Hz)')
% ylabel('Magnitude (dB)')
% title('Pwelch Function')
% grid on

%% FFT
%Single sided amplitude spectrum of U1(t)
% 
figure
dt=mean(diff(t));  %average time step done in the ode45 computation
Fs=10000;
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
%---------------------------------------------------

%% Phase plane
%
% Plots trajectory 
figure
% U1
plot(result(:,1),result(:,2))
xlabel('U_1'); % Insert the x-axis label
ylabel('dU_1/dt'); % Inserts the y-axis label
title('Phase plane') % Inserts the title in the plot
grid on
%---------------------------------------------------

%% Frequency Response Function plot
% input= sin(2*pi*5*t);
% output=result(:,1);
% 
% FRF = fft(output)./fft(input);
% 
% plot(t,FRF);

%% LogLog Plot
figure
loglog(fr(1:floor(m/2)),fourier(1:floor(m/2)))
y1=get(gca,'ylim');
grid on
title('Frequency response magnitudes for a NL one unit-cell AMM','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('magnitude','FontSize',14)
% legend({'mass_1','mass_2'},'FontSize',14)
toc





