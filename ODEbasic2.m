%% ODE computation of single mass spring system - BASIC
clear all

tspan = [0 100];
y0=zeros(12,1);
broadband=wgn(1000000,1,0);


[t,y]=ode45(@fun, tspan, y0);

u1=y(:,1);
v1=y(:,3);
u2=y(:,5);
v2=y(:,7);
u3=y(:,9);
v3=y(:,11);
%% Plots
% Time,displacement
figure
% plot(t,y(:,1),'r',t,y(:,5),'b:')
plot(t,u1,t,v1,t,u2,t,v2,t,u3,t,v3)
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
dft=fft(u2,m); % DFT of signal
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
wind = kaiser(length(u2),105);
% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(u1,v2,[],[],[],1/(dt));

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
periodogram(v2)
%% Functions
function dy=fun(t,y)
x1=y(1);
x2=y(2);
y1=y(3);
y2=y(4);
x3=y(5);
x4=y(6);
y3=y(7);
y4=y(8);
x5=y(9);
x6=y(10);
y5=y(11);
y6=y(12);
% parameters
m=0.1;
% m2=0.05;
k=1000;
F=0.1*cos(5*2*pi*t);
zeta=0.002;
c=zeta*2*sqrt(k*m);  
% dy/dx
dx1=x2;
dx2=(F+k*y1-2*k*x1+x3*k-2*c*x2+c*y2)/m;
dy1=y2;
dy2=x1*k/m-y1*k/m+c*y2/m-c*x2/m;
dx3=x4;
dx4=(-2*k*x3+k*x1+k*y3-c*2*x4+c*x2+c*x6)/m; %+k*x5
dy3=y4;
dy4=(x3*k-y3*k+c*y4-c*x4)/m;
dx5=x6;
dx6=(-2*k*x5+k*x3+k*y3-c*x6+c*y6)/m;
dy5=y6;
dy6=(x5*k-y5*k-y6*c+x6*c)/m;
dy=[dx1;dx2;dy1;dy2;dx3;dx4;dy3;dy4;dx5;dx6;dy5;dy6];
% dy=[dx1;dx2;dy1;dy2;dx3;dx4;dy3;dy4];
end


    