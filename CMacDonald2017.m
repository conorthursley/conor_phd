%% Simple 1D mass-in-mass system with dynamic system outlined by the funciton in the ode45 
% 11/09/2017 - Conor MacDonald
%---------------------------------------------------
clear all
tic
% Time specification
tspan = [0 10];
%---------------------------------------------------
% Initial conditions
% velocity is initially 0, v1=0, v2=0
% displacement for m1 is 1e-3 from the fixed wall (L spacing)
% displacement for m2 is 1E-3 from the wall (5e-4 away from m1, within)
% y=[u1;v1;u2;v2]
% <<<<<<< HEAD
y=[1e-3;0;1e-3;0];
opts = odeset('RelTol',1e-5,'AbsTol',1e-7);
% System simulation
% [t, y] = ode45(@sys, t, y1);
[t, result] = ode15s(@sys1, tspan, y, opts);
toc
% =======
% y=[1e-3;0;1e-3;0];
y=[0 0 0 0];
%---------------------------------------------------
% Parameters
k1=1e3; %N/m
m1=0.1; %kg
k2=1500;
m2=0.5*m1;
w1=sqrt(k1/m1)/(2*pi);
w2=sqrt(k2/m2)/(2*pi);
%---------------------------------------------------
% harmonic input frequency 
% expressed in Hz and then converted to rad/s in the function
input =[0.1 40]; %amp and freq, Hz, 
% LM 1 = 10.3753 LM 2 = 19.5283
% >>>>>>> master


opts = odeset('RelTol',1e-10,'AbsTol',1e-10, 'OutputFcn',@odeplot); %, 'Mass', mass, 'Events', @events);
%% System simulation
% [t, y] = ode45(@sys, t, y1);
[t,result] = ode45(@(t,y)Linear(t,y,input,k1,m1,k2,m2), tspan, y,opts); %,te,ye,ie
toc


%% Plots displacement and velocity of displacement of 2 masses
u1=result(:,1);
u2=result(:,2);
f=input(2)*(2*pi);
H=input(1)*(sin(f*t));
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
linkaxes([ax1,ax2],'x')
%% Combined Displacement
plot(t,u1,t,u2)
xlabel('time'); % Insert the x-axis label
ylabel('displacement'); % Inserts the y-axis label
title('mass-in-mass 1D system') % Inserts the title in the plot
legend('u1','u2')
grid on

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
dft=fft(H,m)/n; % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft);     
plot(fr(1:floor(m/2)),fourier(1:floor(m/2)))
title('Single-Sided Amplitude Spectrum of input signal')
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
% axis([0 100 0 Inf])
%---------------------------------------------------
%Single sided amplitude spectrum of U2(t)
%---------------------------------------------------
ax4=subplot(2,1,2);
dft2=fft(u1,m)/n; % DFT of signal
fourier2 = abs(dft2);     
plot(fr(1:floor(m/2)),fourier2(1:floor(m/2)))
title('Single-Sided Amplitude Spectrum of U1(t)')
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
% axis([0 100 0 Inf])
%---------------------------------------------------
linkaxes([ax3,ax4],'x')


%% Phase plane
%
% Plots trajectory 
figure
% U1
subplot(2,2,1);
plot(result(:,1),result(:,2))
xlabel('U_1'); % Insert the x-axis label
ylabel('dU_1/dt'); % Inserts the y-axis label
title('Phase plane of m1') % Inserts the title in the plot
grid on
%---------------------------------------------------
% U2
%---------------------------------------------------
subplot(2,2,3);
plot(result(:,3),result(:,4))
xlabel('U_2'); % Insert the x-axis label
ylabel('dU_2/dt'); % Inserts the y-axis label
title('Phase planeof m2') % Inserts the title in the plot
grid on
subplot(2,2,[2,4]);
plot(result(:,1),result(:,3))
grid on
title('Invariant manifold','FontSize',20)
xlabel('displacement of mass1','FontSize',20)
ylabel('displacement of mass2','FontSize',20)
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
x=input(1)*(sin(input(2)*2*pi*t));
[txy,frequencies]=tfestimate(x,u1,[],[],[],1/dt);

% plot
figure
% hold on
graph=plot(frequencies/w1,20*log10(abs(txy)),'b');
% axis([0 12 -Inf Inf])
grid on
title('Transfer function of metamaterial configurations','FontSize',14)
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',14)
% axis([0 10 -Inf 0])
% legend({'linear AMM','noise insulation panel','nonlinear case 1','nonlinear case 2','nonlinear case 3','nonlinear case 4'},'FontSize',14)
legend({'1 unit cell','2 unit cells','5 unit cells','10 unit cells'},'FontSize',14)
set(graph,'LineWidth',1.5);
alldatacursors = findall(gcf,'type','hggroup');
set(alldatacursors,'FontSize',20)
%%
figure
% hold on
graph1=semilogx(frequencies/w1,abs(txy),'b');
% axis([0 3 0 Inf])
grid on
% title('TF Estimate 5unitAMM NLH chain','FontSize',14)
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',20)
legend({'linear AMM numerical result','nonlinear case 1 AMM numerical result','phononic crystal numerical result','nonlinear case 2 AMM numerical result'})
set(graph1,'LineWidth',2);
%%
figure
% hold on
[pxx,f] = periodogram(u1,[],[],1/dt);
plot(f/w1,10*log10(pxx),'b')
grid on
title('Periodogram PSD of 10 unit cell linear system','FontSize',20)
xlabel('Frequency, \omega','FontSize',20)
ylabel('Magnitude, dB','FontSize',20)

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



