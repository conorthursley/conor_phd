%% ODE computation of single mass spring system - BASIC
clear all

tspan = [0 100];
y0=zeros(4,1);
% broadband=wgn(1000000,1,0);


[t,y]=ode23s(@fun, tspan, y0);

a=1;

u1=y(a:end,1);
v1=y(a:end,3);
% u2=y(:,5);
% v2=y(:,7);
% u3=y(:,9);
% v3=y(:,11);
%% Plots
% Time,displacement
figure
% plot(t,y(:,1),'r',t,y(:,5),'b:')
plot(t(a:end),u1,t(a:end),v1) %,t,u2,t,v2,t,u3,t,v3)
legend({'U1','V1','U2','V2','U3','V3'},'FontSize',14)


%FFt
%% FFT of input (displacement)
figure
ax1=subplot(2,1,1);
dt=abs(mean(diff(t)));  %average time step done 
Fs=1/(dt);
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(t); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft1=fft(u1,m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft1); 
f=Fs*(0:(n/2))/n;
freq1=fr(1:floor(m/2));
P1=fourier(1:floor(m/2));
plot(freq1,P1)
% plot(flop,abs(y),'LineWidth',2)
title('FFT of input signal')
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',20)
% FFT of output (displacement)

ax2=subplot(2,1,2);
dt=mean(diff(t));  %average time step done 
Fs=1/dt;
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(t); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft=fft(v1,m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft); 
f=Fs*(0:(n/2))/n;
freq=fr(1:floor(m/2));
P=fourier(1:floor(m/2));
plot(freq,P)
% plot(flop,abs(y),'LineWidth',2)
title('FFT of output amp')
grid on
xlabel('Normalised frequency, \Omega (Hz)')
ylabel('|P1(f)|')
linkaxes([ax1,ax2],'x')
set(gca,'fontsize',20)


%% TF estimate
% takes in input and output signal
% bandpass filter input signal
wind = kaiser(length(u1),105);
% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(u1,v1,[],[],[],1/(dt));

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
%periodogram
figure
wind1 = kaiser(length(u1),105);
[pxx,f] = periodogram(v1,[],[],1/dt);
plot(f,10*log10(pxx),'b')
grid on
% title('Periodogram PSD of 10 unit cell linear system','FontSize',20)
xlabel('Frequency, \omega','FontSize',20)
ylabel('Magnitude, dB','FontSize',20)
set(gca,'FontSize',24)
%% Functions
function dy=fun(t,y)
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);
% x3=y(5);
% x4=y(6);
% y3=y(7);
% y4=y(8);
% x5=y(9);
% x6=y(10);
% y5=y(11);
% y6=y(12);
% parameters
m1=0.1;
m2=0.5*m1;
% m2=0.05;
k1=1000;
k2=1.5*k1;
F=0.1*cos(5*2*pi*t);
zeta1=0.00002;
c1=zeta1*2*sqrt(k1*m1);  
c2=zeta1*2*sqrt(k2*m2);  
% dy/dx
dx1=x2;
dx2=(F+k2*y1-2*k1*x1-x1*k2-2*c1*x2-c2*x2)/m1;
dy1=y2;
dy2=(x1*k2-y1*k2+c2*y2-c2*x2)/m2;
% dx3=x4;
% dx4=(-2*k*x3+k*x1+k*y3-c*2*x4+c*y4)/m; %+k*x5
% dy3=y4;
% dy4=(x3*k-y3*k+c*y4-c*x4)/m;
% dx5=x6;
% dx6=(-2*k*x5+k*x3+k*y3-2*c*x6+c*y6+c*x4)/m;
% dy5=y6;
% dy6=(x5*k-y5*k-y6*c+x6*c)/m;
dy=[dx1;dx2;dy1;dy2];
% dy=[dx1;dx2;dy1;dy2;dx3;dx4;dy3;dy4;dx5;dx6;dy5;dy6];
% dy=[dx1;dx2;dy1;dy2;dx3;dx4;dy3;dy4];
end


    