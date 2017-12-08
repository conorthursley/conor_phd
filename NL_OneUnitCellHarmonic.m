    clear all

%% File read from APDL simulation (numerical)

file = 'C:\ANSYS\Temp\Validation\DuffingValDec17\DuffTenUnitHarmonic.csv';
M=csvread(file,2,0); %start reading from row 2, column 1    

ansys_freq = M((1:length(M)),1); % excitation frequency
ansys_amp_1 = M((1:length(M)),2);
% ansys_amp_2 = M((1:10000),4);

%% Displacement freq responses mass 1
figure
plot(ansys_freq,abs(ansys_amp_1),'LineWidth',0.05)
grid on
title('Frequency response magnitudes for a Duffing Oscillator','FontSize',14)
xlabel('Freq, Hz','FontSize',14)
ylabel('Magnitude, u','FontSize',14)
% legend({'mass_1'},'FontSize',14)
%% Displacement freq responses mass2
% figure
% plot(ansys_freq,abs(ansys_amp_1),'LineWidth',0.05)
% grid on
% title('frequency response magnitudes for a NL one unit-cell AMM, secondary mass','FontSize',14)
% xlabel('Freq, Hz','FontSize',14)
% ylabel('Magnitude, u','FontSize',14)
%% Loglog plot of frequency response
figure
loglog(ansys_freq,(abs(ansys_amp_1)))
y1=get(gca,'ylim');
grid on
title('LogLog frequency response magnitudes for a Duffing Oscillator','FontSize',14)
xlabel('Frequency, Hz','FontSize',14)
ylabel('magnitude','FontSize',14)
% legend({'mass_1','mass_2'},'FontSize',14)

