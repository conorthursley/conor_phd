    clear all
%% Parameters
m1=0.1; 
m2=0.3;
k1=1000;
k2=10;
w0=2*pi*sqrt(k2/m2);

% effective mass
theta=m2/m1;
w2=sqrt(k2/m2);

%% File read from APDL simulation (numerical)

file = 'C:\ANSYS\Temp\Validation\DuffingValDec17\DuffOneUnitTrans25.csv';
M=csvread(file,1,0); %start reading from row 1, column 1

ansys_time = M((1:length(M)),1); % time
t=ansys_time;
ansys_amp_1 = M((1:length(M)),2);
ansys_amp_2 = M((1:length(M)),3);

bandpassFile='C:\ANSYS\Temp\Validation\DuffingValDec17\lowpassFilterAmp10Length10.csv';
bandpass=csvread(bandpassFile);
%% FFT of input (displacement)
figure
ax1=subplot(2,1,1);
dt=mean(diff(ansys_time));  %average time step done 
Fs=1/dt;
% y = fft(ansys_amp_1);  
% flop = (0:length(y)-1)*Fs/length(y);
n=length(ansys_time); %length of signal = number of samples
m=pow2(nextpow2(n));  %transform length
dft=fft(bandpass(:,2),m); % DFT of signal
fr = (0:m-1)*(Fs/m);
fourier = abs(dft); 
f=Fs*(0:(n/2))/n;
freq=fr(1:floor(m/2));
P=fourier(1:floor(m/2));
plot(freq/(w2),P)
% plot(flop,abs(y),'LineWidth',2)
title('FFT of input amp (simple)')
grid on
xlabel('Normalised frequency, \Omega (Hz)')
ylabel('|P1(f)|')
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
plot(freq/(w2),P)
% plot(flop,abs(y),'LineWidth',2)
title('FFT of output amp (simple)')
grid on
xlabel('Normalised frequency, \Omega (Hz)')
ylabel('|P1(f)|')
linkaxes([ax1,ax2],'x')

%% Displacement Time responses
figure
plot(ansys_time,(ansys_amp_1),ansys_time,(ansys_amp_2),'r','LineWidth',0.005)
grid on
title('Time response magnitudes for a 10 unit linear spring mass-spring system','FontSize',14)
xlabel('Time, s','FontSize',14)
ylabel('Magnitude, u','FontSize',14)
legend({'mass_1'},'FontSize',14)
%% Loglog plot of frequency response
figure
loglog(freq,abs(P))
y1=get(gca,'ylim');
grid on
title('Frequency response magnitudes for a 10 unit linear spring mass-spring system','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('magnitude','FontSize',14)
% legend({'mass_1','mass_2'},'FontSize',14)
%% TF estimate
% takes in input and output signal
% bandpass filter input signal

% [txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_1,[],[],[],500);
[txy,frequencies]=tfestimate(bandpass(:,2),ansys_amp_2,[],[],[],1/(dt*1*t(end)));

% plot
figure
plot(frequencies/(w2),20*log10(abs(txy)))
grid on
title('TF Estimate','FontSize',14)
xlabel('Normalised frequency, \omega/\omega_0','FontSize',14)
ylabel('Magnitude, dB','FontSize',14)


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
%% amplitude vs freq
% need to get the average amplitude in each section of the frequency signal
% % each frequency is active for 20 seconds
% tend=20;
% ti=1;
% for i=1:15
%     [c, index] = min(abs(ansys_time-tend));
%     SVG=mean(abs(ansys_amp_1(ti:index)));
%     EG(i,1)=SVG;
%     tend=tend+20;
%     ti=index;
%     
% end
% figure
% gee=13:0.25:16.5;
% gg=flipud(EG);
% plot(gee,200*EG)
% grid on

    
    
    