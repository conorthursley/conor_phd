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
%% *******************Parameters**********************
k1=1000; %N/m
m1=0.1; %kg
k2=320;
k2NL=0;
m2=0.5*m1;
w1=sqrt(k1/m1)/(2*pi);
w2=sqrt(k2/m2)/(2*pi);
%---------------------------------------------------
%% *******************Input Signal**********************

bandpassFile='U:\_PhD\APDL\Validation\DuffingValDec17\HigherAmp.csv';
bandpass1=csvread(bandpassFile);
Values=bandpass1(:,2);
Points=10*bandpass1(:,1);
%---------------------------------------------------
% using an interpolated function for the ode solver to read in the ode
% function
input=griddedInterpolant(Points,Values,'linear');
% plot(Points,Values)
%---------------------------------------------------
%% *******************Spring Function**********************
% nonliner spring function found in nonlinearCurve.m
% need to call the function with the parameters declared earlier
F=nonlinearCurve(k2,k2NL);
nlPoints=F(:,1);
nlValues=F(:,2);
NLspring=griddedInterpolant(nlPoints,nlValues,'spline');

% --------check--------
% plot(nlPoints,nlValues,'b-o')
% hold on
% plot(nlPoints*0.5,NLspring(nlPoints*0.5),'r-*')
%---------------------------------------------------
% ode options  - see 'odeset'
opts = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Events', @events); %'OutputFcn',@odeplot);
toc
%% System simulation

% solve continuosly from tstart to tend at each terminal event
% during a for loop from 0 to tfinish
for i=1:10000
    %---------------------------------------------------
    %previous simuluations showed iterations of the event was around 400
    %could probably write a code to keep iterating between the event (dont
    %stop) and event (stop) with a flag system.
    %---------------------------------------------------
    [t,result,te,ye,ie] = ode23(@(t,y)Linear(t,y,input,k1,m1,k2,m2), [tstart tfinal], y,opts);
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
vel1=yout(:,2);
u2=yout(:,3);
vel2=yout(:,4);
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
%% *******************Acceleration**********************
%-------------u1---------------
time=tout;
velocity=vel1;
nn = length(time); %Assume velocity vector is same length
ta = [time(3),time',time(nn-2)];
ve = [velocity(3),velocity',velocity(nn-2)];
t1 = ta(1:nn); t2 = ta(2:nn+1); t3 = ta(3:nn+2);
v1 = ve(1:nn); v2 = ve(2:nn+1); v3 = ve(3:nn+2);
t21 = t2-t1; t32 = t3-t2; t31 = t3-t1;
v21 = v2-v1; v32 = v3-v2;
ac_u1 = (v21./t21.*t32+v32./t32.*t21)./t31; % Approx. acceleration values
%----------------u2---------------
velocity=vel2;
ve = [velocity(3),velocity',velocity(nn-2)];
v1 = ve(1:nn); v2 = ve(2:nn+1); v3 = ve(3:nn+2);
v21 = v2-v1; v32 = v3-v2;
ac_u2 = (v21./t21.*t32+v32./t32.*t21)./t31; % Approx. acceleration values
%*******************Poincare Section**********************
nnnn=size(time);
%set the index of poincare points to 1
%-----------------u1------------
np_u1=1;
for i=1:nnnn(1)
        % detect the cros-section of the trajectory with the plane y1-y2
        if (ac_u1(i)>=(2*pi)*np_u1)
            % store detected cross-section point y1,y2 to ps1,ps2
        ps_u1(np_u1,1)=u1(i);
        ps_u1(np_u1,2)=vel1(i);
        % increase the index of poincare point
                np_u1=np_u1+1;
        end
end
%-----------------u2------------
np_u2=1;
for i=1:nnnn(1)
        % detect the cros-section of the trajectory with the plane y1-y2
        if (ac_u2(i)>=(2*pi)*np_u2)
            % store detected cross-section point y1,y2 to ps1,ps2
        ps_u2(np_u2,1)=u1(i);
        ps_u2(np_u2,2)=vel2(i);
        % increase the index of poincare point
                np_u2=np_u2+1;
        end
end
%% plot poincare section 
figure
subplot(2,1,1);
% plot(u1,v1,'b--')
hold on
for i=1:np_u1-1
    plot(ps_u1(i,1),ps_u1(i,2),'r*')
    % use pause to folow the plot of the poincare section
%     pause(0.5);
end
grid on
subplot(2,1,2);
hold on
for i=1:np_u2-1
    plot(ps_u2(i,1),ps_u2(i,2),'r*')
    % use pause to folow the plot of the poincare section
    %pause(0.5);
end
grid on
%plot(ps(:,1),ps(:,2),'r+');
%% Displacement Time responses
figure
plot(eT,(u1),'b',eT,(u2),'r')
grid on
title('Time response magnitudes for a free-free AMM system with 7 units','FontSize',14)
xlabel('Time, s','FontSize',14)
ylabel('Magnitude, u','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',24)
set(gca,'fontsize',20)
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
dt=mean(diff(Points));  %average time step done in the ode45 computation
Fs=1/dt;
n=length(eT);  %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(Values,m)/n; % DFT of signal
fr = (0:m-1)*(Fs/m);
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
linkaxes([ax3,ax4],'x')


%% Phase plane

%
% Plots trajectory 
figure
% U1
subplot(2,2,1);
plot(u1,yout(:,2))
xlabel('U_1'); % Insert the x-axis label
ylabel('dU_1/dt'); % Inserts the y-axis label
title('Phase plane of m1') % Inserts the title in the plot
grid on
%---------------------------------------------------
% U2
%---------------------------------------------------
subplot(2,2,3);
plot(u2,yout(:,4))
xlabel('U_2'); % Insert the x-axis label
ylabel('dU_2/dt'); % Inserts the y-axis label
title('Phase plane of m2') % Inserts the title in the plot
grid on

subplot(2,2,[2,4]);
plot(u1,u2)
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
fr = (0:m-1)*(Fs/m);
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
fr = (0:m-1)*(Fs/m);
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
%% Transfer function estimate
% Using the TFESTIMATE function to compare input and output signals
% txy = tfestimate(x,y) finds a transfer function estimate, txy, given an input signal, x, and an output signal, y.
% input signal is our sine wave (or harmonic input) across the time length
T=linspace(0,tfinal,length(tout));
x=input(T);

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
set(graph,'LineWidth',1);
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
legend({'linear AMM MATLAB result','nonlinear case 1 AMM numerical result','phononic crystal numerical result','nonlinear case 2 AMM numerical result'})
set(graph1,'LineWidth',2);
alldatacursors = findall(gcf,'type','hggroup');
set(alldatacursors,'FontSize',20)
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
      value = [double((y(3)-(y(1))<=-(5e-1))); double((y(3)-y(1))>=+(5e-1))]; %double((y(7)-(y(5))<=-(5e-4))); double((y(7)-y(5))>=+(5e-4));double((y(1)-y(5))>=+(1.5e-3))]; %need to use double to convert logical to numerical
           % detect when the bounds gets crossed
      isterminal = [1;1]; % halt integration, reverse direction 
      % 1 if the integration is to terminate when the ith event occurs. Otherwise, it is 0.
      direction = [-1;1]; % approaching the event from any which way
      %0 if all zeros are to be located (the default). A value of +1 locates only zeros where the event function is increasing, ...
      % ... and -1 locates only zeros where the event function is decreasing
end
