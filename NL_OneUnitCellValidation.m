    clear all
%% Parameters
m1=0.1; 
m2=0.3;
k1=1000;
k2=1.0659e4;
w1=sqrt((k1)/m1)/(2*pi);
w2=sqrt((k2)/m2)/(2*pi);

% effective mass
theta=m2/m1;
% w2=sqrt(k2/m2);

%% File read from APDL simulation (numerical)

file = 'C:\ANSYS\Temp\Validation\DuffingValDec17\DuffOneUnitTrans62.csv';
M=csvread(file,1,0); %start reading from row 1, column 1

ansys_time = M((1:length(M)),1); % time
t=ansys_time;
ansys_amp_1 = M((1:length(M)),2);
ansys_amp_2 = M((1:length(M)),3);

bandpassFile='C:\ANSYS\Temp\Validation\DuffingValDec17\deleteme.csv';
bandpass1=csvread(bandpassFile);
bandpass=bandpass1(:,2);
% bandpass=10*(sin(100*2*pi*t)+sin(200*2*pi*t)+sin(300*2*pi*t)+sin(400*2*pi*t));
%% FFT of input (displacement)
figure
ax1=subplot(2,1,1);
dt=abs(mean(diff(bandpass(:,1))));  %average time step done 
Fs=1/dt;
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(bandpass); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft=fft(bandpass,m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft); 
f=Fs*(0:(n/2))/n;
freq=fr(1:floor(m/2));
P=fourier(1:floor(m/2));
plot(freq,P)
% plot(flop,abs(y),'LineWidth',2)
title('FFT of input signal (simple)')
grid on
xlabel('Frequency,  (Hz)')
ylabel('|P1(f)|')
set(gca,'fontsize',14)
% FFT of output (displacement)
ax2=subplot(2,1,2);
dt=mean(diff(ansys_time));  %average time step done 
Fs=1/dt;
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(ansys_time); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft=fft(ansys_amp_2,m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft); 
f=Fs*(0:(n/2))/n;
freq=fr(1:floor(m/2));
P=fourier(1:floor(m/2));
plot(freq,P)
% plot(flop,abs(y),'LineWidth',2)
title('FFT of output amp of a free-free AMM system with 7 units')
grid on
xlabel('Normalised frequency, \Omega (Hz)')
ylabel('|P1(f)|')
linkaxes([ax1,ax2],'x')
set(gca,'fontsize',14)

%% Displacement Time responses
figure
plot(ansys_time,(ansys_amp_1),ansys_time,(ansys_amp_2),'r','LineWidth',0.005)
grid on
title('Time response magnitudes for a free-free AMM system with 7 units','FontSize',14)
xlabel('Time, s','FontSize',14)
ylabel('Magnitude, u','FontSize',14)
legend({'mass_1','mass_2'},'FontSize',14)
set(gca,'fontsize',14)
%% Loglog plot of frequency response
figure
loglog(freq/w1,abs(P))
y1=get(gca,'ylim');
grid on
title('Frequency response magnitudes for a 10 unit NLH spring mass-spring system','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('magnitude','FontSize',14)
% legend({'mass_1','mass_2'},'FontSize',14)
%% TF estimate
% takes in input and output signal
% bandpass filter input signal

% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(ansys_amp_2,ansys_amp_1,[],[],[],1/dt);

% plot
figure
plot(frequencies,20*log10(abs(txy)))
% axis([0 12 -Inf Inf])
grid on
title('TF Estimate 5unitAMM NLH chain','FontSize',14)
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',14)
figure
loglog(frequencies,abs(txy))
% axis([0 12 0 10])
grid on
title('TF Estimate 5unitAMM NLH chain','FontSize',14)
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
set(gca,'fontsize',14)
% dt=mean(diff(bandpass(:,1)));  %average time step done 
% Fs=1/dt;
% % y = fft(ansys_amp_1);  
% % flop = (0:length(y)-1)*Fs/length(y);
% n=length(bandpass(:,1)); %length of signal = number of samples
% m=pow2(nextpow2(n));  %transform length
% dft=fft(bandpass(:,2),m); % DFT of signal
% fr = (0:m-1)*(Fs/m);
% fourier = abs(dft); 
% f=Fs*(0:(n/2))/n;
% freq=fr(1:floor(m/2));
% P1=fourier(1:floor(m/2));
% figure
% plot(freq,P,'g',freq,P1,'b')
%% M_effective
% effective mass
w=linspace(0,w2*2*pi*2.5,80000);
wT=sqrt(k2/m2);
A=w;
B=wT;
C=A./B;

meff=m1+(m2*wT^2)./(wT^2-A.^2);

% figure
% plot(meff/m1,C)
% grid
% axis([-20 20 0 3])
% Dispersion relation (theoretical prediction)
% Yao 2008 
qL=2*asin(sqrt((meff)/(4*k1).*A.^2));

% figure
% plot(A/(2*pi),real(qL)/pi)
% grid
% axis([0 12 0 1])

% Transmittance of cells
% Pls automate somehow
B=A/(2*pi);
T7=(k1)./(k1*(1)-meff.*A.^2);
T6=(k1)./(k1*(2-T7)-meff.*A.^2);
T5=(k1)./(k1*(2-T6)-meff.*A.^2);
T4=(k1)./(k1*(2-T5)-meff.*A.^2);
T3=(k1)./(k1*(2-T4)-meff.*A.^2);
T2=(k1)./(k1*(2-T3)-meff.*A.^2);
T1=(k1)./(k1*(2-T2)-meff.*A.^2);
% figure
Tall=T7+T6+T5+T4+T3+T2+T1;
% semilogy(B,abs(Tall))
% grid

%% Graph 3 plots together (dispersion relation, theoretical transmittance,
% numerical transmittance)
SP=7.671; %upper limit of bandgap - work out how to calculate - band gap upper and lower limit
VP=5.752; %lower limit
figure
%------------SubPlot1----------------%
ay1=subplot(3,1,1);
semilogy(frequencies,abs(txy))
grid
line([SP SP],ylim,'Color',[1 0 0])
line([VP VP],ylim,'Color',[1 0 0])
% axis([0.1 12 0.01 10])
title('Dispersion Relation and Transmittance of free-free AMM system with 7 units (validation of Yao 2008)','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)
legend({'Numerical Result (APDL)'},'FontSize',14)
set(gca,'fontsize',14)
%------------SubPlot2----------------%
ay2=subplot(3,1,2);
semilogy(B,abs(Tall))
grid
line([SP SP],ylim,'Color',[1 0 0])
line([VP VP],ylim,'Color',[1 0 0])
% axis([0 12 -10 100])
ylabel('Magnitude, dB','FontSize',14)
legend({'Theoretical Result'},'FontSize',14)
set(gca,'fontsize',14)
%------------SubPlot3----------------%
ay3=subplot(3,1,3);
plot(A/(2*pi),real(qL)/pi)
% axis([0 12 0 1])
grid
line([SP SP],ylim,'Color',[1 0 0])
line([VP VP],ylim,'Color',[1 0 0])
xlabel('Frequency, Hz','FontSize',14)
ylabel('Real(qa/\pi)' ,'FontSize',14)
legend({'Dispersion Relation'},'FontSize',14)
set(gca,'fontsize',14)
% linkaxes([ay1,ay2,ay3],'x')  
    