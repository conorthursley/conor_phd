    clear all

%% File read from APDL simulation (numerical)

file = 'C:\ANSYS\Temp\Validation\NLOneUnitCell8.csv';
M=csvread(file,1,0); %start reading from row 2, column 1

ansys_time = M((1:1004),1); % time
ansys_amp_1 = M((1:1004),2);
ansys_amp_2 = M((1:1004),3);
%% FFT of time to frequency
figure
dt=mean(diff(ansys_time));  %average time step done 
Fs=1/dt;
y = fft(ansys_amp_1);  
f = (0:length(y)-1)*Fs/length(y);
% n=length(ansys_time); %length of signal = number of samples
% Y=fft(ansys_amp_1);
% P2 = abs(Y/n);
% P1 = P2(1:n/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% m=pow2(nextpow2(n));  %transform length
% dft=fft(ansys_amp_1,m)/n; % DFT of signal
% fr = (0:m-1)*(Fs/m)/10;
% fourier = abs(dft); 
% f=Fs*(0:(n/2))/n;
% freq=fr(1:floor(m/2));
% P=fourier(1:floor(m/2));
% plot(freq,P)
plot(f,abs(y),'LineWidth',2)
title('FFT of amp (simple)')
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
%% Displacement Time responses
figure
plot(ansys_time, ansys_amp_1,'r','LineWidth',0.05)
grid on
title('Time response magnitudes for a NL one unit-cell AMM','FontSize',14)
xlabel('Time, s','FontSize',14)
ylabel('Magnitude, u','FontSize',14)
% legend({'mass_1'},'FontSize',14)
%% Loglog plot of frequency response
figure
loglog(f,(abs(y)))
y1=get(gca,'ylim');
grid on
title('Frequency response magnitudes for a NL one unit-cell AMM','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('magnitude','FontSize',14)
% legend({'mass_1','mass_2'},'FontSize',14)

