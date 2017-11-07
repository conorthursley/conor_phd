%% Simple 1D mass-in-mass system with dynamic system outlined by the funciton in the ode45 
% 11/09/2017 - Conor MacDonald
%---------------------------------------------------
clear all
close all
tic
% Time specification
tspan = [0 0.1];
%---------------------------------------------------
% Initial conditions
% velocity is initially 0, v1=0, v2=0
% displacement for m1 is 1e-3 from the fixed wall (L spacing)
% displacement for m2 is 1E-3 from the wall (5e-4 away from m1, within)
% y=[u1;v1;u2;v2]
% y=[1e-3;0;1e-3;0];
y=zeros(8,1);
%---------------------------------------------------
% Parameters
k1=1000; %N/m
m1=1; %kg
k2=0.1*k1;
m2=500*m1;
%---------------------------------------------------
% harmonic input frequency 
% expressed in Hz and then converted to rad/s in the function
input =2; %Hz


opts = odeset('RelTol',1e-5,'AbsTol',1e-7, 'OutputFcn',@odeplot, 'Events', @events); %, 'Mass', mass);
%% System simulation
% [t, y] = ode45(@sys, t, y1);
[t,result,te,ye,ie] = ode45(@(t,y)TwoCell(t,y,input,k1,m1,k2,m2), tspan, y,opts);
toc


%% Plots displacement and velocity of displacement of 2 masses
u1=result(:,1);
u2=result(:,3);
%---------------------------------------------------
figure
ax1=subplot(2,1,1);
plot(t,u1)
% Format plot
xlabel('time'); % Insert the x-axis label
ylabel('displacement'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('u1')
grid on
%---------------------------------------------------
ax2=subplot(2,1,2);
plot(t,u2)
% Format plot
xlabel('time'); % Insert the x-axis label
ylabel('displacement'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('u2')
grid on
%---------------------------------------------------
linkaxes([ax1,ax2],'y')

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
ax3=subplot(2,1,1);
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
% axis([0 100 0 Inf])
%---------------------------------------------------
%Single sided amplitude spectrum of U2(t)
%---------------------------------------------------
ax4=subplot(2,1,2);
dft2=fft(u2,m)/n; % DFT of signal
fourier2 = abs(dft2);     
plot(fr(1:floor(m/2)),fourier2(1:floor(m/2)))
title('Single-Sided Amplitude Spectrum of U2(t)')
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
% axis([0 100 0 Inf])
%---------------------------------------------------
linkaxes([ax3,ax4],'y')


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
%---------------------------------------------------
% U2
%---------------------------------------------------
subplot(2,1,2);
plot(result(:,3),result(:,4))
xlabel('U_2'); % Insert the x-axis label
ylabel('dU_2/dt'); % Inserts the y-axis label
title('Phase plane') % Inserts the title in the plot
grid on

%% Frequency Response Function plot
% input= sin(2*pi*5*t);
% output=result(:,1);
% 
% FRF = fft(output)./fft(input);
% 
% plot(t,FRF);

%% Transfer function estimate
% Using the TFESTIMATE function to compare input and output signals
% txy = tfestimate(x,y) finds a transfer function estimate, txy, given an input signal, x, and an output signal, y.
% input signal is our sine wave (or harmonic input) across the time length
x=0.1*(sin(input*2*pi*t));
Fsample=500;
[txu1,fu1]=tfestimate(x,u1,1024,[],[],Fsample);
[txu2,fu2]=tfestimate(x,u2,1024,[],[],Fsample);
%---------------------------------------------------
% subplot creation
figure
subplot(2,1,1);
plot(fu1,20*log10(abs(txu1)))
title('Transfer function of U1(t)')
grid on
xlabel('f (Hz)')
ylabel('Mag (dB)')
%---------------------------------------------------
subplot(2,1,2);
plot(fu2,mag2db(abs(txu2)))
title('Transfer function of U2(t)')
grid on
xlabel('f (Hz)')
ylabel('Mag (dB)')


toc
%% Event Function
% Using the ODE events function to trigger when the smaller mass reaches
% the bounds of the larer mass. As the smaller mass is inside the larger
% mass
%---------------------------------------------------

function [value,isterminal,direction] = events(~,y)
      value = [double(y(3)<=(y(1)-(5e-4))); double(y(3)>=y(1)+(5e-4))]; %need to use double to convert logical to numerical
           % detect when the bounds gets crossed
      isterminal = [0;0]; % halt integration, reverse direction 
      % 1 if the integration is to terminate when the ith event occurs. Otherwise, it is 0.
      direction = [-1;1]; % approaching the event from any which way
      %0 if all zeros are to be located (the default). A value of +1 locates only zeros where the event function is increasing, ...
      % ... and -1 locates only zeros where the event function is decreasing
end

%% Mass function
% function M = mass
% 
% M=[1 0 0 0; 0 0.5 0 0; 0 0 1 0; 0 0 0 0.5]; %mass matrix
% 
% end



