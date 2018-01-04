%% Experimental code for ODE event solver to find a solution to the secular cubic function

% Simple 1D mass-in-mass system with dynamic system outlined by the funciton in the ode45 
% 11/03/2017 - Conor MacDonald MMDDYYYY
%---------------------------------------------------
clear all
close all
%---------------------------------------------------
tic
%---------------------------------------------------
% Time specification
tstart = 0;
tfinal = 100;
tout = tstart;
%---------------------------------------------------


%---------------------------------------------------
% Initial conditions
% velocity is initially 0, v1=0, v2=0
% displacement for m1 is 1e-3 from the fixed wall (L spacing)
% displacement for m2 is 1E-3 from the wall (5e-4 away from m1, within)
% y=[u1;v1;u2;v2]
y=[0;0;0;0];
% y=zeros(8,1);
refine = 4;
yout = y.';
teout = [];
yeout = [];
ieout = [];
%---------------------------------------------------
% Parameters
k1=1000; %N/m
m1=1; %kg
k2=0.1/k1;
m2=0.3*m1;
%---------------------------------------------------
% harmonic input frequency 
% expressed in Hz and then converted to rad/s in the function

bandpassFile='C:\ANSYS\Temp\Validation\DuffingValDec17\HigherAmp.csv';
bandpass1=csvread(bandpassFile);
bandpass=bandpass1(:,2);
freq=bandpass; %Hz
Amp=bandpass1(:,1);
input=[freq Amp]; 
%---------------------------------------------------
% ode options  - see 'odeset'
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,  'Events', @events); %'OutputFcn',@odeplot,
%% System simulation

% solve continuosly from tstart to tend at each terminal event
% during a for loop from 0 to tfinish
for i=1:10000
    %---------------------------------------------------
    %previous simuluations showed iterations of the event was around 400
    %could probably write a code to keep iterating between the event (dont
    %stop) and event (stop) with a flag system.
    %---------------------------------------------------
    [t,result,te,ye,ie] = ode45(@(t,y)DuffingOsc(t,y,input,k1,m1,k2,m2), [tstart tfinal], y,opts);
%     if ~ishold
%         hold on
%     end %check that the graph is on hold
    %---------------------------------------------------
    %accumulate the output from the ODE solver once the event has triggered
    nt = length(t); %new time
    tout = [tout; t(2:nt)];
    yout = [yout; result(2:nt,:)];
    teout = [teout; te];          % Events at tstart are never reported.
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    %---------------------------------------------------
    % Set the new initial conditions,
    y=[result(nt,1);0;result(nt,3);0]; %;result(nt,5);0;result(nt,7);0]; % displacement should be the same, but once the second mass has hit the 'bounds'
    % then the velocity would have reached zero, hence setting velocity to
    % zero 
    %---------------------------------------------------
    % setup the new 'initial step' to be put back into ODE solver once this
    % for loop has ended
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    options = odeset(opts,'InitialStep',t(nt)-t(nt-refine),...
        'MaxStep',t(nt)-t(1)); %, 'Events', @event);
    %reset the time start for the next ode solve
    tstart = t(nt);
    if tstart==tfinal
        break
    end
    
end

%---------------------------------------------------
toc
%---------------------------------------------------

%% Plots displacement and velocity of displacement of 2 masses
% extract displacement from the results 
%---------------------------------------------------
u1=yout(:,1);
u2=yout(:,3);
eT=tout;
%---------------------------------------------------
figure
ax1=subplot(2,1,1);
plot(eT,u1)
% Format plot
xlabel('time'); % Insert the x-axis label
ylabel('displacement'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('u1')
grid on
%---------------------------------------------------
ax2=subplot(2,1,2);
plot(eT,u2)
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
dt=mean(diff(eT));  %average time step done in the ode45 computation
Fs=1/dt;
n=length(eT);  %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(u1,m)/n; % DFT of signal
fr = (0:m-1)*(Fs/m)/10;
fourier = abs(dft1);     
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
plot(yout(:,1),yout(:,2))
xlabel('U_1'); % Insert the x-axis label
ylabel('dU_1/dt'); % Inserts the y-axis label
title('Phase plane') % Inserts the title in the plot
grid on
%---------------------------------------------------
% U2
%---------------------------------------------------
subplot(2,1,2);
plot(yout(:,3),yout(:,4))
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
% x=input(1)*sin(2*pi*input(2)*t);
% Fsample=1000;
% [txu1,fu1]=tfestimate(x,u1,[],[],[],Fsample);
% [txu2,fu2]=tfestimate(x,u2,[],[],[],Fsample);
% %---------------------------------------------------
% % subplot creation
% figure
% ax5=subplot(2,1,1);
% plot(fu1,20*log10(abs(txu1)))
% title('Transfer function of U1(t)')
% grid on
% xlabel('f (Hz)')
% ylabel('Mag (dB)')
% %---------------------------------------------------
% ax6=subplot(2,1,2);
% plot(fu2,mag2db(abs(txu2)))
% title('Transfer function of U2(t)')
% grid on
% xlabel('f (Hz)')
% ylabel('Mag (dB)')
% %---------------------------------------------------
% linkaxes([ax5,ax6],'y')

%% LogLog plots
y1=fft(u1);
y2=fft(u2);
dt=mean(diff(eT));  %average time step done 
Fs=1/dt;
n=length(eT); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(u1,m)/n; % DFT of signal
fr = (0:m-1)*(Fs/m)/10;
fourier1 = abs(dft1); 
freq=fr(1:floor(m/2));
P3=fourier1(1:floor(m/2));
figure
ax5=subplot(2,1,1);
loglog(freq,abs(P3))
title('Transfer function of U1(t)')
grid on
xlabel('f (Hz)')
ylabel('Mag (dB)')
%---------------------------------------------------
m=pow2(nextpow2(n));  %transform length
dft2=fft(u2,m)/n; % DFT of signal
fr = (0:m-1)*(Fs/m)/10;
fourier1 = abs(dft2); 
freq=fr(1:floor(m/2));
P4=fourier1(1:floor(m/2));
ax6=subplot(2,1,2);
loglog(freq,abs(P4))
title('Transfer function of U2(t)')
grid on
xlabel('f (Hz)')
ylabel('Mag (dB)')
%---------------------------------------------------
linkaxes([ax5,ax6],'y')
toc
%% Event Function
% Using the ODE events function to trigger when the smaller mass reaches
% the bounds of the larer mass. As the smaller mass is inside the larger
% mass
%---------------------------------------------------

function [value,isterminal,direction] = events(~,y)
      value = [double((y(3)-(y(1))<=-(5e-4))); double((y(3)-y(1))>=+(5e-4))]; %double((y(7)-(y(5))<=-(5e-4))); double((y(7)-y(5))>=+(5e-4));double((y(1)-y(5))>=+(1.5e-3))]; %need to use double to convert logical to numerical
           % detect when the bounds gets crossed
      isterminal = [0;0]; % halt integration, reverse direction 
      % 1 if the integration is to terminate when the ith event occurs. Otherwise, it is 0.
      direction = [0;0]; % approaching the event from any which way
      %0 if all zeros are to be located (the default). A value of +1 locates only zeros where the event function is increasing, ...
      % ... and -1 locates only zeros where the event function is decreasing
end
